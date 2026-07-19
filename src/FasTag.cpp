// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <FasTag/fastag.h>
#include <FasTag/fasta_filter.h>

#include <fstream>

using namespace OpenMS;

//-------------------------------------------------------------------------
/**
  @page TOPP_FasTag FasTag

  @brief Infers partial sequence tags from MS/MS spectra.

  FasTag reimplements the DirecTag algorithm (Tabb et al., J Proteome Res 2008,
  7:3838) with three changes that matter:

  - the intensity rank-sum null is computed by dynamic programming rather than by
    enumerating every C(n, k) subset of peak ranks, which is what makes tag
    lengths above 4 reachable at all (2.1e12 subsets at length 7);
  - a seed tag can be extended, so the configured length is a minimum;
  - reported tags can be restricted to exact (optionally isobaric) substrings of
    supplied sequences, and the matching spectra written to a filtered mzML.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FasTag.cli
*/
//-------------------------------------------------------------------------

class TOPPFasTag : public TOPPBase
{
public:
  TOPPFasTag()
    : TOPPBase("FasTag", "Infers partial sequence tags from MS/MS spectra.", true,
               {{"Tabb DL, Ma ZQ, Martin DB, Ham AJ, Chambers MC",
                 "DirecTag: accurate sequence tags from peptide MS/MS through statistical scoring",
                 "J Proteome Res 2008; 7(9): 3838-46", "10.1021/pr800154p"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input spectra");
    setValidFormats_("in", ListUtils::create<std::string>("mzML"));

    registerOutputFile_("out", "<file>", "", "Tag list (tab-separated)");
    setValidFormats_("out", ListUtils::create<std::string>("tsv"));

    registerInputFile_("fasta", "<file>", "", "Report only tags that occur in these sequences", false);
    setValidFormats_("fasta", ListUtils::create<std::string>("fasta"));

    registerOutputFile_("out_spectra", "<file>", "",
                        "Write spectra carrying a reported tag to this mzML", false);
    setValidFormats_("out_spectra", ListUtils::create<std::string>("mzML"));

    registerIntOption_("tag_length", "<n>", 3, "Seed tag length in residues", false);
    setMinInt_("tag_length", 1);
    registerIntOption_("extension", "<n>", 0,
                       "Maximum residues appended per terminus; 0 disables extension", false);
    setMinInt_("extension", 0);

    registerDoubleOption_("fragment_tolerance", "<value>", 20.0, "Fragment mass tolerance", false);
    registerStringOption_("fragment_tolerance_unit", "<unit>", "ppm", "Unit of the fragment tolerance",
                          false);
    setValidStrings_("fragment_tolerance_unit", ListUtils::create<std::string>("ppm,Da"));

    registerIntOption_("max_peaks", "<n>", 100, "Peaks retained per spectrum", false);
    registerIntOption_("max_tags", "<n>", 50, "Tags reported per spectrum", false);
    registerDoubleOption_("max_evalue", "<value>", 20.0,
                          "Discard tags above this E-value; 0 disables", false);

    registerDoubleOption_("isobaric_tolerance", "<value>", 0.04,
                          "Treat a residue as interchangeable with an isobaric residue pair "
                          "within this tolerance when matching against 'fasta'; 0 requires "
                          "exact strings", false);
    registerIntOption_("min_filter_length", "<n>", 0,
                       "Ignore tags shorter than this when matching against 'fasta'; "
                       "0 derives a floor from the database size", false);
    registerStringOption_("orientation", "<mode>", "both",
                          "Match a tag as written, or also reversed. Tags are stored N->C under a "
                          "y-ion assumption, so b-derived tags read reversed relative to the "
                          "sequence", false);
    setValidStrings_("orientation", ListUtils::create<std::string>("both,forward"));
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const String fasta = getStringOption_("fasta");
    const String out_spectra = getStringOption_("out_spectra");

    fastag::Param p;
    p.tag_length = getIntOption_("tag_length");
    p.max_extension = getIntOption_("extension");
    p.frag_tol = getDoubleOption_("fragment_tolerance");
    p.tol_ppm = getStringOption_("fragment_tolerance_unit") == "ppm";
    p.max_peak_count = static_cast<size_t>(getIntOption_("max_peaks"));
    p.max_tag_count = getIntOption_("max_tags");
    p.max_tag_score = getDoubleOption_("max_evalue");

    const int max_len = p.tag_length + 2 * p.max_extension;
    if (!fasta.empty() && max_len > fastag::MAX_FILTER_LEN)
    {
      OPENMS_LOG_ERROR << "Realised tag length can reach " << max_len
                       << ", beyond the filter's " << fastag::MAX_FILTER_LEN
                       << "-residue encoding." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // sequence filter
    //-------------------------------------------------------------
    fastag::FastaFilter filt(getStringOption_("orientation") == "both");
    const bool filtering = !fasta.empty();
    if (filtering)
    {
      // Read via OpenMS so we accept whatever FASTAFile accepts, then hand the
      // plain sequences to the filter.
      std::vector<FASTAFile::FASTAEntry> entries;
      FASTAFile().load(fasta, entries);
      if (entries.empty())
      {
        OPENMS_LOG_ERROR << "No sequences in " << fasta << std::endl;
        return INPUT_FILE_EMPTY;
      }
      String tmp = File::getTemporaryFile();
      {
        std::ofstream f(tmp.c_str());
        for (const auto& e : entries) f << ">" << e.identifier << "\n" << e.sequence << "\n";
      }
      std::string err;
      if (!filt.load(tmp, &err))
      {
        OPENMS_LOG_ERROR << err << std::endl;
        return INPUT_FILE_CORRUPT;
      }
      const int floor_ = getIntOption_("min_filter_length");
      if (floor_ > 0) filt.setMinLen(floor_);
      const double iso = getDoubleOption_("isobaric_tolerance");
      if (iso > 0) filt.deriveCollapses(iso);
      filt.finalise();
      filt.buildRange(p.tag_length, max_len);

      OPENMS_LOG_INFO << "Filter: " << filt.stats().sequences << " sequences, "
                      << filt.stats().residues << " residues; minimum tag length "
                      << filt.minLen()
                      << (filt.stats().min_len_auto ? " (derived)" : " (set)") << std::endl;
      if (!filt.stats().min_len_auto && filt.minLen() < filt.autoMinLen())
      {
        OPENMS_LOG_WARN << "min_filter_length " << filt.minLen() << " is below the derived floor "
                        << filt.autoMinLen() << "; about "
                        << 100 * filt.chanceRate(filt.minLen())
                        << "% of random tags of that length occur in this database by chance."
                        << std::endl;
      }
    }

    //-------------------------------------------------------------
    // tagging
    //-------------------------------------------------------------
    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML});

    const fastag::Tables tables(p);
    std::ofstream tsv(out.c_str());
    tsv << "spectrum\ttag\tlength\tcharge\tnterm_mass\tcterm_mass\textended\tevalue\tfasta_hit\n";

    PeakMap kept;
    kept.getExperimentalSettings() = exp.getExperimentalSettings();

    size_t n_ms2 = 0, n_tags = 0, n_reported = 0;
    for (const auto& spec : exp)
    {
      if (spec.getMSLevel() != 2 || spec.empty() || spec.getPrecursors().empty()) continue;
      ++n_ms2;

      std::vector<double> mz, intensity;
      mz.reserve(spec.size());
      intensity.reserve(spec.size());
      for (const auto& peak : spec)
      {
        mz.push_back(peak.getMZ());
        intensity.push_back(peak.getIntensity());
      }
      const auto& prec = spec.getPrecursors().front();
      const int z = prec.getCharge() > 0 ? prec.getCharge() : 2;

      const auto tags = fastag::tagSpectrum(mz, intensity, prec.getMZ(), z, p, tables);
      n_tags += tags.size();

      bool any = false;
      for (const auto& t : tags)
      {
        const char* hit = "-";
        if (filtering)
        {
          const auto h = filt.match(t.seq);
          if (h == fastag::FastaFilter::Hit::None) continue;
          // A reverse-only match identifies the ion series: the tag was read off
          // the b series, so its flanking masses carry a one-water offset.
          hit = (h == fastag::FastaFilter::Hit::Forward) ? "fwd" : "rev";
        }
        any = true;
        ++n_reported;
        tsv << spec.getNativeID() << '\t' << t.seq << '\t' << t.seq.size() << '\t'
            << t.charge << '\t' << t.nterm_mass << '\t' << t.cterm_mass << '\t'
            << (t.extended ? 1 : 0) << '\t' << t.evalue << '\t' << hit << '\n';
      }
      if (any && !out_spectra.empty()) kept.addSpectrum(spec);
    }
    tsv.close();

    OPENMS_LOG_INFO << n_ms2 << " MS2 spectra, " << n_tags << " tags";
    if (filtering) OPENMS_LOG_INFO << " -> " << n_reported << " after filtering";
    OPENMS_LOG_INFO << std::endl;

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
