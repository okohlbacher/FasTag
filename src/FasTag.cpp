// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include "FasTagger.h"
#include "FastaFilter.h"

#include <fstream>
#include <map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

using namespace OpenMS;

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
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Tag list (tab-separated)");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerInputFile_("fasta", "<file>", "", "Report only tags occurring in these sequences", false);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));
    registerOutputFile_("out_spectra", "<file>", "",
                        "Write spectra carrying a reported tag to this mzML", false);
    setValidFormats_("out_spectra", ListUtils::create<String>("mzML"));

    registerIntOption_("tag_length", "<n>", 3, "Seed tag length in residues", false);
    setMinInt_("tag_length", 1);
    registerIntOption_("extension", "<n>", 0,
                       "Maximum residues appended per terminus; 0 disables extension", false);
    setMinInt_("extension", 0);

    registerDoubleOption_("fragment_tolerance", "<value>", 20.0, "Fragment mass tolerance", false);
    registerStringOption_("fragment_tolerance_unit", "<unit>", "ppm", "Tolerance unit", false);
    setValidStrings_("fragment_tolerance_unit", ListUtils::create<String>("ppm,Da"));

    registerIntOption_("max_peaks", "<n>", 100, "Peaks retained per spectrum", false);
    registerIntOption_("max_tags", "<n>", 50, "Tags reported per spectrum; 0 = unlimited", false);
    registerDoubleOption_("max_evalue", "<value>", 20.0, "E-value cutoff; 0 disables", false);

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
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const String fasta = getStringOption_("fasta");
    const String out_spectra = getStringOption_("out_spectra");

    FasTag::Param p;
    p.tag_length = getIntOption_("tag_length");
    p.max_extension = getIntOption_("extension");
    p.frag_tol = getDoubleOption_("fragment_tolerance");
    p.tol_ppm = getStringOption_("fragment_tolerance_unit") == "ppm";
    p.max_peak_count = static_cast<size_t>(getIntOption_("max_peaks"));
    p.max_tag_count = getIntOption_("max_tags");
    p.max_evalue = getDoubleOption_("max_evalue");

    const int max_len = p.tag_length + 2 * p.max_extension;
    if (!fasta.empty() && max_len > FasTag::MAX_FILTER_LEN)
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

    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML});

    const FasTag::Tables tables(p);
    std::ofstream tsv(out.c_str());
    tsv << "spectrum\ttag\tlength\tcharge\tnterm_mass\tcterm_mass\textended\tevalue\tfasta_hit\n";

    PeakMap kept;
    kept.getExperimentalSettings() = exp.getExperimentalSettings();

    size_t n_ms2 = 0, n_tags = 0, n_reported = 0;
    std::map<int, std::pair<size_t, size_t>> by_len;   // length -> (seen, matched)

    // Parallel over spectra, which is the only safe level: each spectrum is
    // independent, `tables` and `filt` are immutable once built, and per-spectrum
    // work is ~60 us -- far too little to amortise a fork/join inside the tagger.
    //
    // Results are collected into a per-spectrum vector and written afterwards in
    // input order, so the output is identical whatever -threads is set to.
    const SignedSize n_spec = static_cast<SignedSize>(exp.size());
    std::vector<std::string> rows(exp.size());
    std::vector<char> keep(exp.size(), 0);
    std::vector<std::map<int, std::pair<size_t, size_t>>> per_thread_len(
        static_cast<size_t>(std::max(1, omp_get_max_threads())));
    std::vector<size_t> per_thread_ms2(per_thread_len.size(), 0),
                        per_thread_tags(per_thread_len.size(), 0),
                        per_thread_rep(per_thread_len.size(), 0);

#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < n_spec; ++i)
    {
      const auto& spec = exp[static_cast<Size>(i)];
      if (spec.getMSLevel() != 2 || spec.empty() || spec.getPrecursors().empty()) continue;
      const size_t tid = static_cast<size_t>(omp_get_thread_num());
      ++per_thread_ms2[tid];

      const auto& prec = spec.getPrecursors().front();
      const auto tags = FasTag::tagSpectrum(spec, prec.getMZ(), prec.getCharge(), p, tables);
      per_thread_tags[tid] += tags.size();

      std::string buf;
      for (const auto& t : tags)
      {
        const char* hit = "-";
        if (filtering)
        {
          ++per_thread_len[tid][static_cast<int>(t.seq.size())].first;
          const auto h = filt.match(t.seq);
          if (h == FasTag::FastaFilter::Hit::None) continue;
          ++per_thread_len[tid][static_cast<int>(t.seq.size())].second;
          // A reverse-only match identifies the ion series: the tag was read off
          // the b series, so its flanking masses carry a one-water offset.
          hit = (h == FasTag::FastaFilter::Hit::Forward) ? "fwd" : "rev";
        }
        ++per_thread_rep[tid];
        char line[512];
        const int m = std::snprintf(line, sizeof line, "%s\t%s\t%zu\t%d\t%g\t%g\t%d\t%g\t%s\n",
                                    spec.getNativeID().c_str(), t.seq.c_str(), t.seq.size(),
                                    t.charge, t.nterm_mass, t.cterm_mass,
                                    t.extended ? 1 : 0, t.evalue, hit);
        if (m > 0) buf.append(line, static_cast<size_t>(m));
      }
      if (!buf.empty()) { rows[static_cast<size_t>(i)].swap(buf); keep[static_cast<size_t>(i)] = 1; }
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
      if (!out_spectra.empty()) kept.addSpectrum(exp[i]);
    }
    tsv.close();

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
      FileHandler().storeExperiment(out_spectra, kept, {FileTypes::MZML});
      OPENMS_LOG_INFO << "Wrote " << kept.size() << " spectra to " << out_spectra << std::endl;
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFasTag tool;
  return tool.main(argc, argv);
}
