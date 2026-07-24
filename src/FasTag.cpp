// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#ifdef FASTAG_HAVE_MZPEAK
#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <functional>
#endif

#include "FasTagger.h"
#include "FastaFilter.h"
#include "SpectrumSampler.h"
#include "TaxIndex.h"
#include "TaxStats.h"

#include <fstream>
#include <map>
#include <limits>
#include <memory>
#include <algorithm>
#include <iterator>
#include <OpenMS/SYSTEM/File.h>

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <random>
#include <set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

using namespace OpenMS;

namespace
{
  // Deliberately OUTSIDE the FASTAG_HAVE_MZPEAK guard below. -species
  // works on any build, so a mzPeak-less build must still compile this; it
  // lived inside the guard briefly and broke every stock-OpenMS build while
  // CI (which always builds the patched OpenMS) stayed green.
  /// Directory holding the bundled taxonomy, or "" when there is none.
  ///
  /// Resolution order:
  ///   1. $FASTAG_TAXONOMY_DIR   -- a custom or larger reference set
  ///   2. <executable dir>/../share/FasTag/taxonomy   -- the release layout
  ///
  /// Deliberately anchored to the EXECUTABLE, not OPENMS_DATA_PATH: this is
  /// FasTag's own data, and a user pointing OPENMS_DATA_PATH at a system OpenMS
  /// installation must not silently lose the taxonomy that shipped beside the
  /// binary they are running.
  /// k of a taxonomy index, read from its 16-byte header alone.
  ///
  /// Exists so the tag-length precondition can be checked BEFORE tagging.
  /// Loading the index to learn k costs seconds and ~2 GB, which is exactly the
  /// work we want to refuse to waste.
  /// Layout (TaxIndex::save): "FTXI" | uint32 version | uint32 k | uint64 n.
  int peekTaxdbK_(const String& path)
  {
    std::ifstream f(path.c_str(), std::ios::binary);
    if (!f) return -1;
    char magic[4];
    if (!f.read(magic, 4) || std::memcmp(magic, "FTXI", 4) != 0) return -1;
    std::uint32_t ver = 0, k = 0;
    if (!f.read(reinterpret_cast<char*>(&ver), 4) || ver != 1) return -1;
    if (!f.read(reinterpret_cast<char*>(&k), 4)) return -1;
    return static_cast<int>(k);
  }

  String taxonomyDir_()
  {
    const char* env = std::getenv("FASTAG_TAXONOMY_DIR");
    if (env != nullptr && *env != '\0' && File::isDirectory(String(env))) return String(env);

    // Two layouts, because the project ships both:
    //   ../share/FasTag/taxonomy   a normal `cmake --install` tree (bin/ + share/)
    //   ./share-FasTag-taxonomy    the release tarball, which is flat: the
    //                              executable sits at the root next to lib/ and
    //                              share-OpenMS/, so there is no ../share to find.
    // Checking only the first would leave -species broken in every release while
    // working perfectly from a local install -- the worst way to get this wrong.
    // A directory only counts if it actually HOLDS the taxonomy. Accepting the
    // first that merely exists meant an empty ../share/FasTag/taxonomy masked a
    // complete share-FasTag-taxonomy beside it -- exactly the layout a release
    // has if the index was never installed.
    const String exe = File::getExecutablePath();
    String first_existing;
    for (const String& cand : {exe + "../share/FasTag/taxonomy",
                               exe + "share-FasTag-taxonomy",
                               exe + "taxonomy"})
    {
      if (!File::isDirectory(cand)) continue;
      if (first_existing.empty()) first_existing = cand;
      if (File::exists(cand + "/nodes.dmp") && File::exists(cand + "/names.dmp")) return cand;
    }
    // Fall back to one that exists but is incomplete, so the error names a real
    // path the user can inspect rather than reporting "nothing found anywhere".
    return first_existing;
  }
}

//-------------------------------------------------------------------------
/**
  @page TOPP_FasTag FasTag

  @brief Infers partial sequence tags from MS/MS spectra.

  Reimplements DirecTag (Tabb et al., J Proteome Res 2008, 7:3838) with three
  differences: the intensity null is computed by dynamic programming instead of
  by enumerating every C(n,k) subset of peak ranks, which is what makes tag
  lengths above four reachable; a seed tag can be extended, so the configured
  length is a minimum; and reported tags can be restricted to sequences supplied
  as FASTA, with the carrying spectra written out as mzML.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FasTag.cli
*/
//-------------------------------------------------------------------------

#ifdef FASTAG_HAVE_MZPEAK
namespace
{
  /// Buffers pushed spectra into bounded chunks and hands each to a callback.
  ///
  /// MzPeakFile is push-based -- transform() calls consumeSpectrum() once per
  /// spectrum -- while the tagging loop wants a batch to parallelise over. This
  /// adapter is the whole of the mzPeak support; everything downstream is the
  /// existing mzML code path.
  ///
  /// The chunk is bounded by PEAKS, not by spectrum count. A fixed count is the
  /// obvious choice and the wrong one: spectra range from ~100 peaks to over
  /// 130,000 across the benchmark corpus, so "2048 spectra" is somewhere between
  /// 200 K and 270 M peaks. Bounding peaks keeps the buffer flat whatever the
  /// data looks like, which is the property the streaming reader bought and this
  /// must not hand back.
  ///
  /// MS1 and precursor-less spectra are dropped HERE rather than in the tagging
  /// callback, so they never occupy the buffer. On a DDA file that is most of
  /// the input.
  class ChunkingConsumer : public Interfaces::IMSDataConsumer
  {
  public:
    using Flush = std::function<void(std::vector<MSSpectrum>&)>;

    /// @param keep_fraction  <1.0 subsamples: each MS2 spectrum is kept with this
    ///   probability (seeded, so a run is reproducible). 1.0 keeps everything.
    ///   Only a fraction is offered here -- the push interface has no index, so an
    ///   exact count would need reservoir sampling; the mzML path handles counts.
    ChunkingConsumer(Flush flush, size_t peak_budget, double keep_fraction = 1.0,
                     uint32_t seed = 1)
      : flush_(std::move(flush)), budget_(peak_budget),
        keep_(keep_fraction), rng_(seed) {}

    /// Progress hooks. @p tick is called once per spectrum OFFERED (before any
    /// filtering), and @p on_size once with the run's spectrum count, so the two
    /// share a denominator and a determinate percentage is possible. Optional:
    /// without them the consumer behaves exactly as before.
    void setProgressHooks(std::function<void()> tick, std::function<void(Size)> on_size)
    {
      tick_ = std::move(tick);
      on_size_ = std::move(on_size);
    }

    void consumeSpectrum(SpectrumType& s) override
    {
      // Counted before the filters below: progress tracks how far through the
      // INPUT we are, which is what setExpectedSize() counts too. Counting only
      // kept MS2 would stall the bar short of 100% on any file with MS1 scans.
      if (tick_) tick_();
      if (s.getMSLevel() != 2 || s.empty() || s.getPrecursors().empty()) return;
      // Subsample AFTER the MS2 filter so the fraction is of taggable spectra, and
      // BEFORE buffering so dropped spectra never occupy memory. Called serially
      // by transform(), so one RNG is safe.
      if (keep_ < 1.0 && unit_(rng_) >= keep_) return;
      peaks_ += s.size();
      buf_.push_back(s);
      if (peaks_ >= budget_) { flush_(buf_); peaks_ = 0; }
    }

    /// FasTag tags spectra; chromatograms are not input to it.
    void consumeChromatogram(ChromatogramType&) override {}

    /// Still NOT used for presizing -- that would make correctness depend on a
    /// call the interface only recommends ("expected to be called"), and the
    /// flush callback grows its own storage instead. It is forwarded to the
    /// progress hook only, where being advisory is harmless: a wrong or missing
    /// count costs an approximate percentage, never a wrong result.
    void setExpectedSize(Size n_spectra, Size) override
    {
      if (on_size_) on_size_(n_spectra);
    }

    void setExperimentalSettings(const ExperimentalSettings& e) override { settings_ = e; }
    const ExperimentalSettings& settings() const { return settings_; }

    /// Must be called after transform(): the final partial chunk is otherwise
    /// never flushed and its spectra vanish without a word.
    void finish() { if (!buf_.empty()) { flush_(buf_); peaks_ = 0; } }

  private:
    Flush flush_;
    std::function<void()> tick_;
    std::function<void(Size)> on_size_;
    size_t budget_;
    double keep_ = 1.0;
    std::mt19937 rng_;
    std::uniform_real_distribution<double> unit_{0.0, 1.0};
    size_t peaks_ = 0;
    std::vector<MSSpectrum> buf_;
    ExperimentalSettings settings_;
  };
}
#endif

class TOPPFasTag : public TOPPBase
{
public:
  // official = false: FasTag lives outside the OpenMS tree and so is not in
  // ToolHandler's list. Passing true makes the TOPPBase constructor throw
  // InvalidValue on startup.
  TOPPFasTag()
    : TOPPBase("FasTag", "Infers partial sequence tags from MS/MS spectra.", false,
               {{"Tabb DL, Ma ZQ, Martin DB, Ham AJ, Chambers MC",
                 "DirecTag: accurate sequence tags from peptide MS/MS through statistical scoring",
                 "J Proteome Res 2008; 7(9): 3838-46", "10.1021/pr800154p"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input spectra");
    // mzpeak is offered only when the OpenMS this was built against provides
    // MzPeakFile, so the accepted formats differ between builds. That is
    // deliberate -- advertising a format the binary cannot read would turn a
    // clear "unsupported" into a confusing parse failure -- and the configure
    // step prints which one you have.
#ifdef FASTAG_HAVE_MZPEAK
    setValidFormats_("in", ListUtils::create<String>("mzML,mzpeak"));
#else
    setValidFormats_("in", ListUtils::create<String>("mzML"));
#endif
    registerOutputFile_("out", "<file>", "", "Tag list (tab-separated)");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerInputFile_("fasta", "<file>", "", "Report only tags occurring in these sequences", false);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));
    registerOutputFile_("out_spectra", "<file>", "",
                        "Write spectra carrying a reported tag here (mzML, or "
                        "mzpeak if this build has it). Note this path holds one "
                        "slot per input spectrum, so unlike the default it needs "
                        "memory proportional to the FILE, not to the thread count", false);
#ifdef FASTAG_HAVE_MZPEAK
    setValidFormats_("out_spectra", ListUtils::create<String>("mzML,mzpeak"));
#else
    setValidFormats_("out_spectra", ListUtils::create<String>("mzML"));
#endif

    registerIntOption_("tag_length", "<n>", 3, "Seed tag length in residues", false);
    setMinInt_("tag_length", 1);
    registerIntOption_("extension", "<n>", 0,
                       "Maximum residues appended per terminus; 0 disables extension", false);
    setMinInt_("extension", 0);
    // tag_length + 2*extension is computed as int and indexes the null tables.
    // Unbounded, --extension 1073741824 overflows and aborts with std::length_error.
    setMaxInt_("extension", (FasTag::MAX_FILTER_LEN - 1) / 2);
    registerIntOption_("gaps", "<n>", 0,
                       "Allow a tag to cross one missing peak, spelling the two "
                       "residues either side of it from their summed mass; "
                       "0 disables", false);
    setMinInt_("gaps", 0);
    registerFlag_("deisotope",
                  "Collapse isotope clusters to their monoisotopic peak and move "
                  "multiply-charged fragments onto the singly-charged scale before "
                  "peak selection", false);
    // One gap only. Each additional gap multiplies the branching and asserts
    // another unobserved split, and a two-gap tag would be mostly inference.
    setMaxInt_("gaps", 1);

    registerDoubleOption_("fragment_tolerance", "<value>", 20.0, "Fragment mass tolerance", false);
    // Without a floor a negative or zero tolerance matches nothing and the run
    // reports zero tags with no indication why.
    setMinFloat_("fragment_tolerance", 1e-9);
    registerStringOption_("fragment_tolerance_unit", "<unit>", "ppm", "Tolerance unit", false);
    setValidStrings_("fragment_tolerance_unit", ListUtils::create<String>("ppm,Da"));

    registerIntOption_("max_peaks", "<n>", 400,
                       "Peaks retained per spectrum; 0 uses the internal ceiling "
                       "of 1024, not unlimited -- the scoring tables are built to "
                       "that bound. The ceiling for 'peaks_per_window'", false);
    setMinInt_("max_peaks", 0);
    registerIntOption_("peaks_per_window", "<n>", 10,
                       "Keep this many peaks per 100 Da window instead of the "
                       "strongest 'max_peaks' overall; 0 disables. Scales the "
                       "effective peak budget with the spectrum's m/z range, so "
                       "dense spectra are not starved and sparse ones not padded", false);
    setMinInt_("peaks_per_window", 0);
    registerIntOption_("max_tags", "<n>", 50, "Tags reported per spectrum; 0 = unlimited", false);
    setMinInt_("max_tags", 0);

    // Random subsampling of input spectra. For a taxon call (or a quick preview)
    // a run does not need every spectrum -- a uniformly random subset gives the
    // same lowest-common-ancestor at a fraction of the cost. Selection is seeded
    // (-subsample_seed) so a run is reproducible. 0 / 0.0 disables. If both are
    // set, the absolute count wins.
    registerIntOption_("subsample_spectra", "<n>", 0,
                       "Tag only this many randomly chosen input spectra; 0 = all", false);
    setMinInt_("subsample_spectra", 0);
    registerDoubleOption_("subsample_fraction", "<f>", 0.0,
                          "Tag only this random fraction (0..1] of input spectra; 0 = all", false);
    setMinFloat_("subsample_fraction", 0.0);
    setMaxFloat_("subsample_fraction", 1.0);
    registerIntOption_("subsample_seed", "<n>", 1, "Seed for -subsample_* selection", false);
    setMinInt_("subsample_seed", 0);

    // Machine-readable progress for a GUI/pipeline driving FasTag. Off by
    // default so the CLI stays quiet; when set, emit periodic lines to stderr:
    //   FASTAG_PROGRESS done=<n> total=<n>
    // total is the spectrum count for indexed input, or 0 when it is not known
    // ahead of time (the push-based mzPeak path streams and cannot count first).
    // Emitted from a single thread per update, so lines never interleave.
    registerFlag_("progress", "Emit 'FASTAG_PROGRESS done=<n> total=<n>' lines to stderr for a GUI progress bar");

    // Taxonomic / species detection: throw the tags at a prebuilt tag->taxon
    // index (buildtaxdb) and infer the taxa present by lowest-common-ancestor
    // over the tag hits. A reduced reference set gives honest genus/family
    // resolution; short peptide tags cannot resolve species. Pair with
    // -subsample_* for a fast call on a large run.
    // The switch. Off by default -- loading an index costs seconds and ~2 GB,
    // which no plain tagging run should pay for. With it on, the three inputs
    // below default to the taxonomy bundled beside the binary, so `-species` on
    // its own is a complete request.
    registerFlag_("species", "Infer which taxa are present from the tags. Uses the bundled "
                             "taxonomy unless -taxdb/-taxonomy_* are given, and writes "
                             "<out>.species.tsv unless -species_out is given");

    registerInputFile_("taxdb", "<file>", "", "Tag->taxon index (built by buildtaxdb). Default: the bundled tax_k7.taxdb", false);
    registerInputFile_("taxonomy_nodes", "<file>", "", "NCBI taxonomy nodes.dmp. Default: the bundled nodes.dmp", false);
    registerInputFile_("taxonomy_names", "<file>", "", "NCBI taxonomy names.dmp. Default: the bundled names.dmp", false);
    registerOutputFile_("species_out", "<file>", "", "Ranked taxa TSV. Default: <out> with a .species.tsv suffix", false);
    setValidFormats_("species_out", ListUtils::create<String>("tsv"));
    registerIntOption_("species_min_len", "<n>", 0,
                       "Ignore tags shorter than this for taxonomy; 0 uses the index k", false);
    setMinInt_("species_min_len", 0);
    registerStringOption_("species_rank", "<rank>", "genus",
                          "Report taxa at this NCBI rank (genus, family, ...)", false);
    registerDoubleOption_("max_evalue", "<value>", 20.0, "E-value cutoff; 0 disables", false);
    setMinFloat_("max_evalue", 0.0);
    registerDoubleOption_("gap_penalty", "<value>", 100.0,
                          "Rank gapped tags as if their E-value were this many times "
                          "worse. Affects ORDER only, never which tags are reported. "
                          "1 disables. Gapped tags are otherwise heavily over-ranked: "
                          "they take ~95% of top-1 slots while being ~3.6x less likely "
                          "to be correct", false);
    setMinFloat_("gap_penalty", 1.0);

    registerDoubleOption_("isobaric_tolerance", "<value>", 0.04,
                          "When matching against 'fasta', treat a residue as interchangeable "
                          "with an isobaric residue pair within this tolerance; 0 requires "
                          "exact strings", false);
    registerIntOption_("min_filter_length", "<n>", 0,
                       "Ignore tags shorter than this when matching against 'fasta'; "
                       "0 derives a floor from the database size", false);
    registerStringOption_("orientation", "<mode>", "both",
                          "Match a tag as written, or also reversed. Tags are stored N->C under a "
                          "y-ion assumption, so b-derived tags read reversed", false);
    setValidStrings_("orientation", ListUtils::create<String>("both,forward"));

    // Modifications, resolved through OpenMS ModificationsDB so any UniMod name
    // works, not a fixed list. Fixed mods shift a residue's mass and are not
    // annotated (the shift is implicit, as carbamidomethyl-C conventionally is);
    // variable mods add a modified alternative that is written inline in the tag
    // as X[Name] (e.g. S[Phospho]). Residue-specific mods only -- a terminal mod
    // shifts the whole precursor and so is already carried in the flanking
    // masses, with nothing per-residue to annotate.
    registerStringList_("fixed_modifications", "<mods>",
                        ListUtils::create<String>("Carbamidomethyl (C)"),
                        "Fixed modifications by OpenMS/UniMod name, e.g. "
                        "'Carbamidomethyl (C)' 'TMT6plex (K)'. Shift residue masses; "
                        "not annotated in the tag", false);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""),
                        "Variable modifications by OpenMS/UniMod name, e.g. "
                        "'Phospho (S)' 'Phospho (T)' 'Phospho (Y)' 'Formyl (K)'. "
                        "Add a modified alternative, written inline as X[Name]", false);
  }

  /// Resolve OpenMS/UniMod modification names to FasTag::ModSpec via
  /// ModificationsDB. Terminal mods are skipped with a note -- they are absorbed
  /// into the reported flanking masses, not the internal residue alphabet.
  void resolveMods_(const StringList& names, bool variable, std::vector<FasTag::ModSpec>& out)
  {
    auto* db = ModificationsDB::getInstance();
    for (const String& nm : names)
    {
      if (nm.empty()) continue;
      const ResidueModification* mod = nullptr;
      try { mod = db->getModification(nm); }
      catch (Exception::BaseException&)
      {
        OPENMS_LOG_ERROR << "Unknown modification '" << nm << "'; ignoring." << std::endl;
        continue;
      }
      if (!mod) continue;
      if (mod->getTermSpecificity() != ResidueModification::ANYWHERE)
      {
        OPENMS_LOG_INFO << "Modification '" << nm << "' is terminal; it is absorbed "
                           "into the reported flanking masses and not annotated "
                           "per-residue." << std::endl;
        continue;
      }
      const char origin = mod->getOrigin();
      if (origin < 'A' || origin > 'Z')
      {
        OPENMS_LOG_WARN << "Modification '" << nm << "' has no single-residue origin; "
                           "ignoring." << std::endl;
        continue;
      }
      FasTag::ModSpec ms;
      ms.residue = (origin == 'I') ? 'L' : origin;  // I is folded to L in the alphabet
      ms.delta = mod->getDiffMonoMass();
      ms.name = mod->getId();
      ms.variable = variable;
      out.push_back(ms);
    }
  }

  ExitCodes main_(int, const char**) override
  {
    // Species precondition, checked BEFORE the run rather than after.
    //
    // The index is keyed on k-mers, so a tag shorter than k can never be looked
    // up. With the default tag_length of 3 against the bundled k=7 index that is
    // EVERY tag: the run used to succeed, spend seconds and ~2 GB loading the
    // index, and write a header-only report -- which reads as "nothing found"
    // instead of "nothing could be found". A warning after the fact was not
    // enough; refuse up front and say exactly what to change.
    if (getFlag_("species"))
    {
      String probe_taxdb = getStringOption_("taxdb");
      if (probe_taxdb.empty())
      {
        const String d = taxonomyDir_();
        if (!d.empty()) probe_taxdb = d + "/tax_k7.taxdb";
      }
      const int kk = probe_taxdb.empty() ? -1 : peekTaxdbK_(probe_taxdb);
      if (kk > 0)
      {
        const int reach = getIntOption_("tag_length") + 2 * getIntOption_("extension");
        if (reach < kk)
        {
          OPENMS_LOG_ERROR << "-species: tag_length (" << getIntOption_("tag_length")
                           << ") + 2*extension (" << getIntOption_("extension") << ") = "
                           << reach << " cannot reach the index k = " << kk
                           << ", so every tag would be too short to look up and the "
                           << "report would be empty. Set -tag_length " << kk
                           << " (or raise -extension)." << std::endl;
          return ILLEGAL_PARAMETERS;
        }
      }
    }

    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const String fasta = getStringOption_("fasta");
    const String out_spectra = getStringOption_("out_spectra");

    FasTag::Param p;
    p.tag_length = getIntOption_("tag_length");
    p.max_extension = getIntOption_("extension");
    p.max_gaps = getIntOption_("gaps");
    p.deisotope = getFlag_("deisotope");
    p.peaks_per_window = getIntOption_("peaks_per_window");
    p.frag_tol = getDoubleOption_("fragment_tolerance");
    p.tol_ppm = getStringOption_("fragment_tolerance_unit") == "ppm";
    p.max_peak_count = static_cast<size_t>(getIntOption_("max_peaks"));
    p.max_tag_count = getIntOption_("max_tags");
    p.max_evalue = getDoubleOption_("max_evalue");
    p.gap_penalty = getDoubleOption_("gap_penalty");
    resolveMods_(getStringList_("fixed_modifications"), false, p.mods);
    resolveMods_(getStringList_("variable_modifications"), true, p.mods);
    if (!p.mods.empty())
    {
      String summary;
      for (const auto& m : p.mods)
        summary += (summary.empty() ? "" : ", ") + m.name + " (" + String(m.residue) + ", "
                 + (m.variable ? "variable" : "fixed") + ")";
      OPENMS_LOG_INFO << "Modifications: " << summary << std::endl;
    }

    // Unconditional: the realised length bounds the null tables whether or not a
    // FASTA filter is in use.
    const int max_len = p.tag_length + 2 * p.max_extension;
    if (max_len > FasTag::MAX_FILTER_LEN)
    {
      OPENMS_LOG_ERROR << "Realised tag length can reach " << max_len << ", beyond the filter's "
                       << FasTag::MAX_FILTER_LEN << "-residue encoding." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    FasTag::FastaFilter filt(getStringOption_("orientation") == "both");
    const bool filtering = !fasta.empty();
    if (filtering)
    {
      std::vector<FASTAFile::FASTAEntry> entries;
      FASTAFile().load(fasta, entries);
      std::string err;
      if (!filt.load(entries, &err))
      {
        OPENMS_LOG_ERROR << "FASTA: " << err << std::endl;
        return INPUT_FILE_EMPTY;
      }
      const int floor_ = getIntOption_("min_filter_length");
      if (floor_ > 0) filt.setMinLen(floor_);
      const double iso = getDoubleOption_("isobaric_tolerance");
      if (iso > 0) filt.deriveCollapses(iso);
      filt.build(p.tag_length, max_len);
      OPENMS_LOG_INFO << "Filter index: " << filt.indexedKeys() << " keys" << std::endl;

      OPENMS_LOG_INFO << "Filter: " << filt.sequenceCount() << " sequences, "
                      << filt.residueCount() << " residues; minimum tag length " << filt.minLen()
                      << (filt.minLenAuto() ? " (derived)" : " (set)")
                      << "; " << filt.collapseRules() << " isobaric rules" << std::endl;
      if (!filt.minLenAuto() && filt.minLen() < filt.autoMinLen())
      {
        OPENMS_LOG_WARN << "min_filter_length " << filt.minLen() << " is below the derived floor "
                        << filt.autoMinLen() << "; about " << 100 * filt.chanceRate(filt.minLen())
                        << "% of random tags that length occur in this database by chance."
                        << std::endl;
      }
    }

    // Stream from disk rather than loading the run into memory.
    //
    // loadExperiment() reads the whole file into a PeakMap first: measured at
    // 7.1 GB peak RSS for a 5.3 GB mzML, with most of the wall time spent in that
    // single-threaded load. OnDiscMSExperiment reads one spectrum at a time from
    // an indexed mzML, so memory becomes O(threads) instead of O(file).
    //
    // It is documented as NOT thread-safe -- it holds an open stream -- so each
    // worker gets its own copy, exactly as the class comment prescribes. Falls
    // back to a full load when the file carries no index, since random access is
    // then impossible.
    // Skip the up-front metadata pass; each worker reads its own spectrum's
    // metadata from the block it already parses.
    //
    // openFile() otherwise calls loadMetaData_(), a fully serial parse of every
    // spectrum's metadata before any work starts -- ~55 s of a 60 s run on S23,
    // and the reason wall time stopped improving past 16 threads however many
    // cores it was given.
    //
    // REQUIRES a patched OpenMS. Stock MzMLSpectrumDecoder::domParseSpectrum()
    // fills the binary data and native ID only, so with skipMetaData every
    // spectrum has MS level 0 and no precursor and the run silently returns
    // nothing. The patch harvests MS level, RT and precursors from the DOM that
    // function already builds; verified identical to the canonical parser across
    // all 632,677 spectra of S23. Guarded below rather than assumed.
    //
    // -out_spectra still needs getMetaData() for the run-level settings, so that
    // path keeps the full load.
    // mzPeak is read through its own push-based path, so none of the mzML
    // reader setup below applies to it.
    // FileTypes::MZPEAK is referenced only under the guard.
    //
    // An OpenMS can ship MzPeakFile.h while its FileTypes enum has no MZPEAK
    // member -- bioconda's build is exactly that -- so detecting the header is
    // not sufficient to know the enum exists. CI caught this: the reference sat
    // outside the #ifdef and failed to compile on every platform, which is the
    // one thing the compile-time gate was supposed to prevent.
#ifdef FASTAG_HAVE_MZPEAK
    const bool mzpeak_in = FileHandler::getType(in) == FileTypes::MZPEAK;
    // Which writer -out_spectra gets, decided from its own extension rather
    // than the input's: every one of the four in/out combinations is allowed,
    // so mzML->mzpeak and mzpeak->mzML both work as a side effect of tagging.
    const bool mzpeak_out = !out_spectra.empty()
                            && FileHandler::getType(out_spectra) == FileTypes::MZPEAK;
#else
    const bool mzpeak_in = false;
    const bool mzpeak_out = false;
#endif

    // A pointer, not a value: the fallback below must replace this with a
    // genuinely fresh reader, and OnDiscMSExperiment's operator= is private
    // (copy-construction only), so an in-place reassignment isn't available.
    auto ondisc = std::make_unique<OnDiscMSExperiment>();
    PeakMap exp;
    const bool streaming = !mzpeak_in && ondisc->openFile(in, out_spectra.empty());
    if (!streaming && !mzpeak_in)
    {
      OPENMS_LOG_WARN << "'" << in << "' has no usable index; reading it entirely "
                         "into memory. Run FileConverter to write an indexed mzML "
                         "if the file is large." << std::endl;
      FileHandler().loadExperiment(in, exp, {FileTypes::MZML});
    }
    // Refuse to run blind on an unpatched OpenMS.
    //
    // Without the patch every spectrum comes back at MS level 1 (the
    // default-constructed value, never overwritten) with no precursor, so the
    // MS2 test rejects all of them and the tool reports a clean run with an
    // empty output file. That is the worst possible failure -- indistinguishable
    // from a file that genuinely holds no MS2 -- so probe for it and take the
    // slow path instead of producing a confident nothing.
    //
    // The probe checks for level 2 WITH a precursor -- what tag_one() actually
    // requires -- not merely "level > 0". A bare ">0" is satisfied by the very
    // default this guard exists to catch (level 1 is > 0), so it passed on
    // spectrum 0 of every file regardless of what the reader actually reported,
    // and the fallback never engaged. Confirmed against a real 617 MB Thermo
    // file built with stock bioconda OpenMS: every spectrum read back at level 1
    // through the fast path, 0 of 53,521 MS2 spectra found, no warning printed.
    if (streaming && out_spectra.empty() && !mzpeak_in)
    {
      bool have_meta = false;
      const Size probe = std::min<Size>(ondisc->getNrSpectra(), 64);
      for (Size i = 0; i < probe && !have_meta; ++i)
      {
        const MSSpectrum s = ondisc->getSpectrum(i);
        if (s.getMSLevel() == 2 && !s.getPrecursors().empty()) have_meta = true;
      }
      if (!have_meta && probe > 0)
      {
        OPENMS_LOG_WARN << "This OpenMS does not report spectrum metadata from the "
                           "per-spectrum reader, so the fast path is unavailable; "
                           "loading metadata up front instead. Expect roughly a "
                           "minute of extra single-threaded startup on a large file."
                        << std::endl;
        // A fresh instance, not a second openFile() on this one: calling
        // openFile() again on the SAME OnDiscMSExperiment crashed (EXC_BAD_ACCESS
        // inside MSSpectrum's copy constructor, reached via getSpectrum()) on the
        // same real file this whole guard exists for -- this fallback had never
        // actually run before, since the old ">0" probe never triggered it.
        ondisc = std::make_unique<OnDiscMSExperiment>();
        ondisc->openFile(in);
      }
    }
    const size_t n_total = mzpeak_in ? 0 : (streaming ? ondisc->getNrSpectra() : exp.size());

    // Subsampling selection. For the indexed (mzML) paths the mask is exact -- it
    // is built over the spectrum indices up front. The count over ALL input
    // spectra, not MS2 only: the fast reader does not know the MS level up front,
    // and a fraction is proportional regardless. For the push-based mzPeak path,
    // where there is no index, only a FRACTION is supported, applied as a seeded
    // per-spectrum Bernoulli keep inside the consumer (below).
    const size_t subsample_n = static_cast<size_t>(getIntOption_("subsample_spectra"));
    const double subsample_frac = getDoubleOption_("subsample_fraction");
    const uint32_t subsample_seed = static_cast<uint32_t>(getIntOption_("subsample_seed"));
    const bool subsampling = subsample_n > 0 || subsample_frac > 0.0;
    if (mzpeak_in && subsample_n > 0)
    {
      OPENMS_LOG_ERROR << "-subsample_spectra (absolute count) needs an index and is "
                          "not supported for mzPeak input; use -subsample_fraction."
                       << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    FasTag::SampleMask sample_mask;
    if (subsampling && !mzpeak_in)
    {
      sample_mask = subsample_n > 0 ? FasTag::sampleByCount(n_total, subsample_n, subsample_seed)
                                    : FasTag::sampleByFraction(n_total, subsample_frac, subsample_seed);
      size_t sel = 0; for (char c : sample_mask) sel += c ? 1 : 0;
      OPENMS_LOG_INFO << "Subsampling: tagging " << sel << " of " << n_total
                      << " input spectra (seed " << subsample_seed << ")" << std::endl;
    }

    // Warn when the fragment tolerance looks far too tight for the data.
    //
    // A high-resolution tolerance on low-resolution data is silent, and looks
    // exactly like a bad file. Measured on an ion-trap MS2 run: 3,007 tags at
    // 20 ppm against 824,959 at 0.3 Da -- a factor of 274, with nothing in the
    // output to suggest the setting was at fault rather than the spectra.
    //
    // The analyser is not read from run-level metadata, because the fast reader
    // path deliberately does not load it. Infer resolution from the peaks
    // instead: an ion trap cannot place two peaks closer than roughly 0.2-0.3 Th,
    // while an Orbitrap or TOF routinely does. If the closest pair anywhere in
    // the sample is still many times the tolerance, no real fragment can be
    // matched at that tolerance either.
    {
      double tightest = mzpeak_in ? 0.0 : std::numeric_limits<double>::max();
      double at_mz = 0;
      const Size want = std::min<Size>(n_total, 200);
      const Size step = std::max<Size>(1, n_total / std::max<Size>(want, 1));
      Size seen = 0;
      for (Size i = 0; i < n_total && seen < want; i += step)
      {
        MSSpectrum s = streaming ? ondisc->getSpectrum(i) : exp[i];
        if (s.getMSLevel() != 2 || s.size() < 8) continue;
        ++seen;
        s.sortByPosition();
        for (Size k = 1; k < s.size(); ++k)
        {
          const double d = s[k].getMZ() - s[k - 1].getMZ();
          if (d > 0 && d < tightest) { tightest = d; at_mz = s[k].getMZ(); }
        }
      }
      if (at_mz > 0)
      {
        const double tol = p.tol_ppm ? at_mz * p.frag_tol * 1e-6 : p.frag_tol;
        if (tol > 0 && tightest > 20.0 * tol)
        {
          OPENMS_LOG_WARN
            << "No two peaks anywhere in this file are closer than " << tightest
            << " Th, yet the fragment tolerance is " << tol << " Th near m/z "
            << at_mz << ". That is what low-resolution (ion trap) MS2 read with a "
                        "high-resolution tolerance looks like, and it yields almost "
                        "no tags for a reason the output cannot show. If the MS2 is "
                        "ion trap, try -fragment_tolerance 0.3 "
                        "-fragment_tolerance_unit Da." << std::endl;
        }
      }
    }

    const FasTag::Tables tables(p);
    std::ofstream tsv(out.c_str());
    tsv << "spectrum\ttag\tlength\tcharge\tnterm_mass\tcterm_mass\textended\tgapped\tevalue\tmin_conf\tmean_conf\tfasta_hit\n";

    PeakMap kept;
    if (streaming)
    {
      if (auto meta = ondisc->getMetaData()) kept.getExperimentalSettings() = *meta;
    }
    else kept.getExperimentalSettings() = exp.getExperimentalSettings();

    size_t n_ms2 = 0, n_tags = 0, n_reported = 0;
    std::map<int, std::pair<size_t, size_t>> by_len;   // length -> (seen, matched)

    // Parallel over spectra, which is the only safe level: each spectrum is
    // independent, `tables` and `filt` are immutable once built, and per-spectrum
    // work is ~60 us -- far too little to amortise a fork/join inside the tagger.
    //
    // Results are collected into a per-spectrum vector and written afterwards in
    // input order, so the output is identical whatever -threads is set to.
    const SignedSize n_spec = static_cast<SignedSize>(n_total);
    std::vector<std::string> rows(n_total);
    std::vector<char> keep(n_total, 0);
    std::vector<MSSpectrum> kept_spec(out_spectra.empty() ? 0 : n_total);
    std::vector<std::map<int, std::pair<size_t, size_t>>> per_thread_len(
        static_cast<size_t>(std::max(1, omp_get_max_threads())));
    std::vector<size_t> per_thread_ms2(per_thread_len.size(), 0),
                        per_thread_tags(per_thread_len.size(), 0),
                        per_thread_rep(per_thread_len.size(), 0);

    // The per-spectrum work, shared by every input path.
    //
    // Takes its thread id rather than calling omp_get_thread_num() internally,
    // and returns its rows rather than writing them: the caller owns where a
    // result lands. That is what lets a second input path reuse this without the
    // two drifting apart, which is the failure mode that matters here -- two
    // readers that mostly agree are worse than one.
    //
    // Callers must not resize rows/keep/kept_spec while this runs.
    auto tag_one = [&](const MSSpectrum& spec, size_t tid) -> std::string
    {
      std::string buf;
      if (spec.getMSLevel() != 2 || spec.empty() || spec.getPrecursors().empty()) return buf;
      ++per_thread_ms2[tid];

      const auto& prec = spec.getPrecursors().front();
      const auto tags = FasTag::tagSpectrum(spec, prec.getMZ(), prec.getCharge(), p, tables);
      per_thread_tags[tid] += tags.size();

      for (const auto& t : tags)
      {
        const char* hit = "-";
        if (filtering)
        {
          ++per_thread_len[tid][static_cast<int>(t.n_res)].first;
          const auto h = filt.match(FasTag::baseSequence(t.seq));
          if (h == FasTag::FastaFilter::Hit::None) continue;
          ++per_thread_len[tid][static_cast<int>(t.n_res)].second;
          // A reverse-only match identifies the ion series: the tag was read off
          // the b series, so its flanking masses carry a one-water offset.
          hit = (h == FasTag::FastaFilter::Hit::Forward) ? "fwd" : "rev";
        }
        ++per_thread_rep[tid];
        char line[512];
        // Flanking masses at 4 decimals (0.1 mDa), not %g.
        //
        // %g gives 6 significant digits, so a 1234.56789 Da flank prints as
        // 1234.57 and a 2500 Da flank keeps a single decimal -- coarser than the
        // tolerance the tag was found with, and these are precisely the values a
        // downstream search uses as a precursor constraint. E-values keep %g,
        // where relative precision is what matters.
        const int m = std::snprintf(line, sizeof line,
                                    "%s\t%s\t%zu\t%d\t%.4f\t%.4f\t%d\t%d\t%g\t%.3f\t%.3f\t%s\n",
                                    spec.getNativeID().c_str(), t.seq.c_str(), t.n_res,
                                    t.charge, t.nterm_mass, t.cterm_mass,
                                    t.extended ? 1 : 0, t.gapped ? 1 : 0, t.evalue,
                                    t.min_conf, t.mean_conf, hit);
        // snprintf returns the length it WOULD have written, which a long native
        // ID can push past the buffer. Appending that many bytes reads past
        // `line`. Clamp, and re-format on the heap rather than emit a truncated
        // row -- half a row is a corrupt TSV, not a cosmetic problem.
        if (m < 0) continue;
        if (static_cast<size_t>(m) < sizeof line)
        {
          buf.append(line, static_cast<size_t>(m));
        }
        else
        {
          std::vector<char> big(static_cast<size_t>(m) + 1);
          const int m2 = std::snprintf(big.data(), big.size(),
                                       "%s\t%s\t%zu\t%d\t%.4f\t%.4f\t%d\t%d\t%g\t%s\n",
                                       spec.getNativeID().c_str(), t.seq.c_str(), t.n_res,
                                       t.charge, t.nterm_mass, t.cterm_mass,
                                       t.extended ? 1 : 0, t.gapped ? 1 : 0, t.evalue, hit);
          if (m2 > 0) buf.append(big.data(), static_cast<size_t>(m2));
        }
      }
      return buf;
    };

    // Record one spectrum's result at its global index. Serial or parallel: the
    // index is the caller's, so output order never depends on scheduling.
    auto record = [&](size_t idx, std::string&& buf, const MSSpectrum& spec)
    {
      if (buf.empty()) return;
      rows[idx].swap(buf);
      keep[idx] = 1;
      if (!out_spectra.empty()) kept_spec[idx] = spec;
    };

    // Tag a buffered chunk in parallel and append its rows in order.
    //
    // Used by the mzPeak path, which is push-based: MzPeakFile hands over one
    // spectrum at a time, so there is nothing to index into and no random access
    // to parallelise over. Buffer, then run the same per-spectrum work over the
    // buffer, appending at base+j so order follows input regardless of
    // scheduling -- the guarantee the mzML path gets from indexing by i.
    //
    // The vectors grow HERE, outside the parallel region. Growing them inside
    // would be a data race, and presizing from setExpectedSize would make
    // correctness depend on a call the interface only recommends.
    // Progress reporting (opt-in via -progress), shared across both input paths.
    // One atomic counter; the single thread that observes each step boundary
    // emits one line under a critical section, so lines never interleave. total
    // is the spectrum count for indexed input, or 0 (unknown) for the streaming
    // mzPeak path, which cannot be counted before it is read.
    const bool emit_progress = getFlag_("progress");
    std::atomic<long long> progress_done{0};
    // Not const, and 0 until known: the mzPeak path learns its spectrum count
    // from setExpectedSize() partway in (see the consumer's progress hooks), so
    // that run starts indeterminate and becomes a real percentage.
    std::atomic<long long> progress_total{mzpeak_in ? 0 : static_cast<long long>(n_spec)};
    std::atomic<int> progress_pct{-1};
    std::atomic<long long> progress_ms{0};
    const auto progress_t0 = std::chrono::steady_clock::now();
    auto elapsed_ms = [&progress_t0]() {
      return static_cast<long long>(std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - progress_t0).count());
    };

    // At most ONE line per whole percent, and at least one every 5 s.
    //
    // "At most one per percent", not "one for every percent": the emitting
    // thread reports the counter's CURRENT value, which other threads have
    // advanced past, so a busy run legitimately jumps 5% -> 8%. That is right
    // for a progress bar and wrong to describe as per-percent.
    //
    // Percent alone goes quiet for minutes on a slow file (one 5 GB run spends
    // ~a minute in metadata before the first spectrum), which reads as a hung
    // GUI; a pure time interval spams a fast run. Together the line rate is
    // bounded by ~101 + elapsed/5s regardless of file size or speed.
    auto tick = [&]()
    {
      if (!emit_progress) return;
      const long long d = ++progress_done;
      const long long tot = progress_total.load(std::memory_order_relaxed);
      const int pct = tot > 0 ? static_cast<int>((d * 100) / tot) : -1;
      // 100% is NOT emitted here. The loop finishing is not the tool finishing --
      // species classification and -out_spectra still run -- so the completion
      // line is reserved for the end of main_(). Without this the last percent
      // and the completion line print the same text twice.
      bool want = (d == 1) || (pct < 100 && pct > progress_pct.load(std::memory_order_relaxed));
      // The clock is consulted only every 64th spectrum: the time rule exists for
      // SLOW runs, where 64 spectra is a rounding error, and this keeps the hot
      // path free of a clock read per spectrum.
      // The clock is read on EVERY tick, not one in 64. Sampling made the
      // advertised 5 s guarantee false exactly where it matters: at one
      // spectrum per second the gap became ~63 s, and a single slow spectrum
      // produced no heartbeat at all. steady_clock::now() is tens of
      // nanoseconds against microseconds-to-seconds of tagging per spectrum.
      if (!want && elapsed_ms() - progress_ms.load(std::memory_order_relaxed) >= 5000)
      {
        want = true;
      }
      if (!want) return;
#pragma omp critical(fastag_progress)
      {
        // Re-check under the lock: several threads can decide to emit at once,
        // and without this they each print a line for the same percent.
        const long long dd = progress_done.load();
        const long long tt = progress_total.load();
        const int pp = tt > 0 ? static_cast<int>((dd * 100) / tt) : -1;
        const long long now = elapsed_ms();
        if (dd == 1 || (pp < 100 && pp > progress_pct.load())
            || (pp < 100 && now - progress_ms.load() >= 5000))
        {
          progress_pct.store(pp);
          progress_ms.store(now);
          std::cerr << "FASTAG_PROGRESS done=" << dd << " total=" << tt << std::endl;
        }
      }
    };

    auto flush_chunk = [&](std::vector<MSSpectrum>& buf)
    {
      if (buf.empty()) return;
      const size_t base = rows.size();
      rows.resize(base + buf.size());
      keep.resize(base + buf.size(), 0);
      if (!out_spectra.empty()) kept_spec.resize(base + buf.size());

      const SignedSize n = static_cast<SignedSize>(buf.size());
#pragma omp parallel for schedule(dynamic, 16)
      for (SignedSize j = 0; j < n; ++j)
      {
        const size_t tid = static_cast<size_t>(omp_get_thread_num());
        record(base + static_cast<size_t>(j), tag_one(buf[j], tid), buf[j]);
      }
      buf.clear();
    };

    if (mzpeak_in)
    {
#ifdef FASTAG_HAVE_MZPEAK
      // 1 M peaks per chunk, ~16 MB of Peak1D, independent of how many spectra
      // that turns out to be.
      //
      // WARNING -- mzPeak input is NOT memory-bounded, and the cause is not
      // here. MzPeakFile::transform() materialises the run rather than streaming
      // it, despite the consumer interface existing to avoid exactly that.
      // Measured: a 2.11 GB .mzpeak peaks at 23.7 GB resident, ~11x the file,
      // against 169 MB for the same data as mzML. Varying this budget over a 16x
      // range moved peak memory by under 10% (1027 / 945 / 986 MB on a 93 MB
      // file), which is what shows the buffer is not the term that matters.
      //
      // Output is correct and byte-identical to the mzML path. Fine for files
      // small relative to RAM; prefer mzML for large runs until the upstream
      // reader streams. Tracked in doc/BACKLOG-mzpeak.md.
      ChunkingConsumer consumer(flush_chunk, 1000u * 1000u,
                                subsample_frac > 0.0 ? subsample_frac : 1.0, subsample_seed);
      // Progress for the streaming path: the count arrives from the reader via
      // setExpectedSize(), so the run reports "n spectra" until it does and a
      // real percentage afterwards.
      consumer.setProgressHooks(tick, [&progress_total](Size n) {
        progress_total.store(static_cast<long long>(n));
      });
      MzPeakFile().transform(in, &consumer);
      consumer.finish();   // the final partial chunk, otherwise silently dropped
      // Run-level settings for -out_spectra come from the consumer here; the
      // mzML paths take them from getMetaData()/the loaded map above, neither
      // of which ran. Without this the written file has no run-level metadata
      // at all.
      if (!out_spectra.empty()) kept.getExperimentalSettings() = consumer.settings();
#endif
    }
    else
    {
#pragma omp parallel
    {
      // One reader per thread: OnDiscMSExperiment keeps an open stream and is
      // documented as not thread-safe.
      // Copy-constructed, not assigned: OnDiscMSExperiment's operator= is private.
      std::unique_ptr<OnDiscMSExperiment> reader;
      if (streaming) reader = std::make_unique<OnDiscMSExperiment>(*ondisc);

#pragma omp for schedule(dynamic, 64)
      for (SignedSize i = 0; i < n_spec; ++i)
      {
        if (!sample_mask.empty() && !sample_mask[static_cast<size_t>(i)]) { tick(); continue; }
        const MSSpectrum loaded = streaming ? reader->getSpectrum(static_cast<Size>(i))
                                            : MSSpectrum();
        const MSSpectrum& spec = streaming ? loaded : exp[static_cast<Size>(i)];
        const size_t tid = static_cast<size_t>(omp_get_thread_num());
        record(static_cast<size_t>(i), tag_one(spec, tid), spec);
        tick();
      }
    }
    }

    // Land the bar on 100%.
    //
    // The denominator can be an upper bound the run never reaches: for mzPeak it
    // comes from setExpectedSize(), which counts every spectrum the file's
    // metadata describes, while the reader delivers only those with point data
    // (42,092 of 53,521 on a real Lumos run -- the bar would stop at 79%).
    // Emitting done==total once at the end costs one line and avoids a GUI that
    // sits at four-fifths on a finished run.
    // NOTE: this is deliberately NOT the completion line. It reports that the
    // TAGGING LOOP finished; species classification and -out_spectra can still
    // run for many seconds after it, and an earlier version emitted 100% here
    // and then kept working -- or emitted 100% and then returned an error. The
    // real completion line is at the end of main_().
    if (false)
    {
      // Skipped when the loop already reported 100%, which the indexed path
      // does via its d == total case -- otherwise every mzML run ends with the
      // same line twice.
      //
      // total is the LARGER of what was announced and what was delivered, never
      // the delivered count alone. The mzPeak reader announces every spectrum in
      // the file but delivers only those with point data (42,092 of 53,521), and
      // rewriting total downwards made it non-monotonic -- a consumer could not
      // tell "the reader skipped some" from "the total was always smaller".
      const long long dd = progress_done.load();
      const long long tt = std::max(dd, progress_total.load());
      std::cerr << "FASTAG_PROGRESS done=" << tt << " total=" << tt << std::endl;
    }

    for (size_t t = 0; t < per_thread_len.size(); ++t)
    {
      n_ms2 += per_thread_ms2[t];
      n_tags += per_thread_tags[t];
      n_reported += per_thread_rep[t];
      for (const auto& kv : per_thread_len[t])
      {
        by_len[kv.first].first += kv.second.first;
        by_len[kv.first].second += kv.second.second;
      }
    }
    for (size_t i = 0; i < rows.size(); ++i)
    {
      if (!keep[i]) continue;
      tsv << rows[i];
      if (!out_spectra.empty()) kept.addSpectrum(kept_spec[i]);
    }
    tsv.close();

    // Taxonomic / species detection from the tags.
    //
    // For every spectrum, the SET of taxa its tags support (a tag supports a
    // taxon when every one of its k-mers is in that taxon -- intersection, so a
    // longer tag is more specific). Each taxon is counted ONCE per spectrum, not
    // per tag, so correlated tags from one spectrum cannot manufacture evidence.
    // Roll the per-leaf counts up the taxonomy, compare each node against the
    // breadth it has in the index, and rank. The q-value is a ranking aid, not a
    // calibrated FDR: the per-k-mer background is a proxy and the chimeric-null
    // problem (see F4 in doc/BACKLOG.md) applies here too.
    String taxdb = getStringOption_("taxdb");
    String species_out = getStringOption_("species_out");
    String nodes = getStringOption_("taxonomy_nodes");
    String names = getStringOption_("taxonomy_names");
    // -species is the switch. Passing -taxdb and -species_out explicitly still
    // works without it, which is how this was driven before the flag existed.
    const bool want_species = getFlag_("species") || (!taxdb.empty() && !species_out.empty());

    if (want_species)
    {
      const String tdir = taxonomyDir_();
      // An index and a taxonomy are ONE coherent set. Defaulting them
      // independently let `-taxdb custom.taxdb` silently pair a custom index
      // with the bundled 117-node dumps: rollUp() drops every taxid the dumps
      // do not know, so calls go missing and the background is biased -- with
      // no error. Either all three come from the bundle, or the user supplies
      // the taxonomy that goes with their index.
      const bool custom_index = !taxdb.empty();
      if (!custom_index && !tdir.empty()) taxdb = tdir + "/tax_k7.taxdb";
      if (!custom_index)
      {
        if (nodes.empty() && !tdir.empty()) nodes = tdir + "/nodes.dmp";
        if (names.empty() && !tdir.empty()) names = tdir + "/names.dmp";
      }
      else if (nodes.empty() || names.empty())
      {
        OPENMS_LOG_ERROR << "-taxdb was given explicitly, so -taxonomy_nodes and "
                         << "-taxonomy_names must be too: the bundled taxonomy "
                         << "describes only the bundled index, and pairing it with "
                         << "another index silently drops every taxon it does not "
                         << "know." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      if (species_out.empty())
      {
        // Strip an extension only if the dot is in the FILE NAME. `-out
        // /tmp/run.v1/tags` has its last dot in the directory, and blindly
        // cutting there wrote /tmp/run.species.tsv -- a different directory.
        const size_t slash = out.find_last_of("/\\");
        const size_t dot = out.rfind('.');
        const bool ext = dot != std::string::npos && (slash == std::string::npos || dot > slash);
        species_out = (ext ? out.substr(0, dot) : out) + ".species.tsv";
      }

      // Say which file is missing and where it was looked for. "Failed to load"
      // on an empty path is useless when the whole point of the defaults is that
      // the user never typed one.
      for (const auto& need : {std::make_pair("taxdb", taxdb),
                               std::make_pair("taxonomy_nodes", nodes),
                               std::make_pair("taxonomy_names", names)})
      {
        if (need.second.empty() || !File::exists(need.second))
        {
          OPENMS_LOG_ERROR << "-species: no " << need.first << ". ";
          if (tdir.empty())
          {
            OPENMS_LOG_ERROR << "No bundled taxonomy was found next to the executable "
                             << "(expected <bin>/../share/FasTag/taxonomy). Pass -" << need.first
                             << " explicitly, or set FASTAG_TAXONOMY_DIR." << std::endl;
          }
          else
          {
            OPENMS_LOG_ERROR << "Looked in '" << tdir << "'. Pass -" << need.first
                             << " explicitly, or set FASTAG_TAXONOMY_DIR." << std::endl;
          }
          return ILLEGAL_PARAMETERS;
        }
      }
      FasTag::TaxIndex idx;
      std::string ierr;
      if (!idx.load(taxdb, &ierr))
      {
        OPENMS_LOG_ERROR << "Failed to load -taxdb '" << taxdb << "': " << ierr << std::endl;
        return INPUT_FILE_CORRUPT;
      }
      FasTag::Taxonomy tax;
      std::string terr;
      if (!tax.load(nodes, names, &terr))
      {
        OPENMS_LOG_ERROR << "Failed to load taxonomy: " << terr << std::endl;
        return INPUT_FILE_CORRUPT;
      }
      const int kk = idx.k();
      int min_len = getIntOption_("species_min_len");
      if (min_len < kk) min_len = kk;

      // The silent-empty trap: the index is keyed on k-mers, so a tag shorter
      // than k can never be looked up. With the default tag_length of 3 against
      // a k=7 index that is EVERY tag -- the run then succeeds, reports "0
      // spectra contributed taxon evidence" and writes a header-only file, which
      // reads as "nothing found" rather than "nothing could have been found".
      const int reach = getIntOption_("tag_length") + 2 * getIntOption_("extension");
      if (reach < kk)
      {
        OPENMS_LOG_WARN << "-species: no tag can reach the index k (this run was "
                        << "driven by -taxdb/-species_out rather than -species, which "
                        << "checks up front). tag_length ("
                        << getIntOption_("tag_length") << ") + 2*extension ("
                        << getIntOption_("extension") << ") = " << reach << " < k = " << kk
                        << ", so every tag is too short to look up and the report WILL be empty. "
                        << "Raise -tag_length to " << kk << " (or use -extension to reach it)."
                        << std::endl;
      }

      // Per-spectrum taxon support -> per-leaf spectrum counts.
      std::map<uint32_t, uint64_t> hits;
      size_t n_units = 0;
      // Group the written rows by spectrum: field 0 = spectrum id, field 1 = tag.
      std::map<std::string, std::vector<std::string>> by_spec;
      for (size_t i = 0; i < rows.size(); ++i)
      {
        if (!keep[i]) continue;
        const std::string& r = rows[i];
        size_t start = 0;
        while (start < r.size())
        {
          size_t nl = r.find('\n', start);
          if (nl == std::string::npos) nl = r.size();
          size_t t1 = r.find('\t', start);
          if (t1 != std::string::npos && t1 < nl)
          {
            size_t t2 = r.find('\t', t1 + 1);
            if (t2 != std::string::npos && t2 <= nl)
              by_spec[r.substr(start, t1 - start)].push_back(r.substr(t1 + 1, t2 - t1 - 1));
          }
          start = nl + 1;
        }
      }
      for (const auto& sp : by_spec)
      {
        std::set<uint32_t> taxa;  // taxa supported anywhere in this spectrum
        for (const std::string& raw : sp.second)
        {
          const std::string seq = FasTag::baseSequence(raw);
          if (static_cast<int>(seq.size()) < min_len) continue;
          // Intersection of the tag's k-mer taxon sets: the taxon must carry the
          // whole tag, not just one window.
          std::vector<uint32_t> acc;
          bool first = true;
          const int L = static_cast<int>(seq.size());
          for (int i = 0; i + kk <= L; ++i)
          {
            const std::vector<uint32_t>& t = idx.lookup(FasTag::TaxIndex::fold(seq.substr(static_cast<size_t>(i), static_cast<size_t>(kk))));
            if (first) { acc = t; first = false; }
            else
            {
              std::vector<uint32_t> tmp;
              std::set_intersection(acc.begin(), acc.end(), t.begin(), t.end(), std::back_inserter(tmp));
              acc.swap(tmp);
            }
            if (acc.empty()) break;
          }
          for (uint32_t tx : acc) taxa.insert(tx);
        }
        if (taxa.empty()) continue;
        ++n_units;
        for (uint32_t tx : taxa) ++hits[tx];
      }

      // Roll observed hits AND the index breadth up the tree, so each node's
      // background is the breadth of its whole subtree.
      std::vector<std::pair<uint32_t, uint64_t>> hit_pairs(hits.begin(), hits.end());
      std::vector<std::pair<uint32_t, uint64_t>> post_pairs;
      for (uint32_t tx : idx.taxa()) post_pairs.emplace_back(tx, idx.postings(tx));
      const auto rolled = tax.rollUp(hit_pairs);
      const auto rolled_post = tax.rollUp(post_pairs);
      std::map<uint32_t, double> bg;
      const double nk = static_cast<double>(std::max<uint64_t>(1, idx.nKmers()));
      for (const auto& np : rolled_post) bg[np.taxid] = static_cast<double>(np.subtree_hits) / nk;
      std::vector<std::pair<uint32_t, double>> bg_pairs(bg.begin(), bg.end());

      const double qthr = 1.0;  // rank everything; the report shows q, caller thresholds
      auto calls = tax.call(rolled, bg_pairs, static_cast<uint64_t>(n_units), qthr);

      const String want_rank = getStringOption_("species_rank");
      std::ofstream so(species_out.c_str());
      so << "rank\ttaxid\tname\tobserved\texpected\tenrichment\tlog_pvalue\tqvalue\n";
      size_t shown = 0;
      for (const auto& c : calls)
      {
        if (c.rank != want_rank) continue;
        const double enr = c.expected > 0 ? c.observed / c.expected : 999.0;
        char line[512];
        std::snprintf(line, sizeof line, "%s\t%u\t%s\t%llu\t%.1f\t%.1fx\t%g\t%g\n",
                      c.rank.c_str(), c.taxid, c.name.c_str(),
                      static_cast<unsigned long long>(c.observed), c.expected, enr,
                      c.log_pvalue, c.qvalue);
        so << line;
        if (++shown >= 50) break;
      }
      so.close();
      OPENMS_LOG_INFO << "Species: " << n_units << " spectra contributed taxon evidence; "
                      << "top " << want_rank << " calls -> " << species_out << std::endl;
    }

    OPENMS_LOG_INFO << n_ms2 << " MS2 spectra, " << n_tags << " tags";
    if (filtering) OPENMS_LOG_INFO << " -> " << n_reported << " after filtering";
    OPENMS_LOG_INFO << std::endl;

    // Report matches beside the number expected by chance, so a count is never
    // mistaken for signal: below the derived floor the filter passes almost
    // everything.
    if (filtering && !by_len.empty())
    {
      char line[160];
      std::snprintf(line, sizeof line, "%5s %12s %12s %12s %10s",
                    "len", "tags", "matched", "by chance", "enrichment");
      OPENMS_LOG_INFO << line << std::endl;
      for (const auto& kv : by_len)
      {
        const size_t seen = kv.second.first, hit = kv.second.second;
        if (seen == 0) continue;
        const double expect = filt.chanceRate(kv.first) * static_cast<double>(seen);
        // Only a length that actually matched something can be uninformative;
        // zero matches is simply no evidence either way.
        const char* note = "";
        if (hit > 0 && static_cast<double>(hit) <= 2.0 * expect) note = "  <- at chance";
        char enr[24] = "-";
        if (hit > 0) std::snprintf(enr, sizeof enr, "%.1fx",
                                   expect > 0 ? hit / expect : 999.0);
        std::snprintf(line, sizeof line, "%5d %12zu %12zu %12zu %10s%s",
                      kv.first, seen, hit, static_cast<size_t>(expect + 0.5), enr, note);
        OPENMS_LOG_INFO << line << std::endl;
      }
    }

    if (!out_spectra.empty())
    {
      if (kept.empty())
      {
        OPENMS_LOG_ERROR << "No spectrum carried a reported tag; not writing " << out_spectra
                         << std::endl;
        return UNEXPECTED_RESULT;
      }
      addDataProcessing_(kept, getProcessingInfo_(DataProcessing::FILTERING));
      if (mzpeak_out)
      {
#ifdef FASTAG_HAVE_MZPEAK
        // Not routed through FileHandler: its storeExperiment() has no mzPeak
        // branch, so asking it for one silently writes something else.
        MzPeakFile().store(out_spectra, kept);
#endif
      }
      else
      {
        FileHandler().storeExperiment(out_spectra, kept, {FileTypes::MZML});
      }
      OPENMS_LOG_INFO << "Wrote " << kept.size() << " spectra to " << out_spectra << std::endl;
    }
    // Completion, emitted once everything the run promised has actually been
    // written: tags, species report and -out_spectra. total is the LARGER of the
    // announced and delivered counts, so it never shrinks (the mzPeak reader
    // announces every spectrum in the file but delivers only those with point
    // data, 42,092 of 53,521).
    if (emit_progress)
    {
      const long long dd = progress_done.load();
      const long long tt = std::max(dd, progress_total.load());
      std::cerr << "FASTAG_PROGRESS done=" << tt << " total=" << tt << std::endl;
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFasTag tool;
  return tool.main(argc, argv);
}
