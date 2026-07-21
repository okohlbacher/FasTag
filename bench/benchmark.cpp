// FasTag benchmark: synthetic spectra with exact ground truth, plus the
// calibration tool that ties them to real acquisitions.
//
//   benchmark stats <file.mzML> [tol_ppm]   measure a real file
//   benchmark synth <profile> [options]     generate, tag, score
//
// WHY SYNTHETIC. Quality was previously measured one way only: SAGE PSMs on one
// diaTracer file. That gives real spectra but costs a search engine, a proteome,
// and a 5 GB input, none of which travel with the repo -- and when those
// artefacts went missing the measurement could not be repeated at all. It also
// answers only "was a correct tag found", never "was it found FIRST", because
// PSM-based truth is per spectrum, not per tag.
//
// Synthetic spectra invert both trade-offs. Truth is exact and per tag, so
// ranking is measurable; and the whole benchmark is a seed plus a profile, so it
// runs anywhere in seconds and gates CI.
//
// WHAT IT DOES NOT CLAIM. Synthetic spectra cannot discover failure modes the
// generator does not model -- co-eluting near-isobars, unmodelled PTMs,
// detector saturation, or whatever real data does that nobody thought to write
// down. This measures RANKING BEHAVIOUR, which is exactly what it is for, and it
// does not replace validation on real data. The `stats` mode exists to keep the
// generator honest about the one structural property that governs tagging.
//
// CALIBRATION. Profiles are not invented. Each is fitted so its ladder edge
// density and peak count match a real acquisition measured by `stats` -- see
// doc/TEST-DATA.md for the four cornerstones and PROFILES below for the fits.
// Ladder edge density (the fraction of peaks with another peak exactly one
// residue mass away) is the right target because it is the structural property
// tagging consumes, and it needs no identifications to measure.
#include "FasTagger.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // The 19 residues FasTag spells with: 20 standard amino acids minus I, which
  // is isobaric with L and folded onto it everywhere in this project.
  const char* ALPHABET = "AGSPVTCLNDQKEMHFRYW";

  /// Approximate UniProt residue frequencies, I folded into L.
  ///
  /// Realistic composition matters here only because it sets how often a peptide
  /// contains repeated or near-isobaric stretches, which is what makes a tag
  /// ambiguous. Uniform composition would make the benchmark measurably easier.
  const double FREQ[19] = {
      0.0826, 0.0707, 0.0665, 0.0473, 0.0686, 0.0536, 0.0138, 0.1563, 0.0406,
      0.0546, 0.0393, 0.0582, 0.0674, 0.0241, 0.0227, 0.0386, 0.0553, 0.0292,
      0.0110};

  double residueMass(char c)
  {
    static std::map<char, double> m;
    if (m.empty())
      for (const char* p = ALPHABET; *p; ++p)
        m[*p] = AASequence::fromString(std::string(1, *p)).getMonoWeight(Residue::Internal);
    auto it = m.find(c);
    return it == m.end() ? 0.0 : it->second;
  }

  /// Is `tag` a subsequence of `pep`, in either orientation, I/L folded?
  ///
  /// Both orientations because tags are stored N->C under a y-ion assumption, so
  /// a b-derived tag reads reversed. Accepting only one orientation would score
  /// roughly half of all correct tags as wrong.
  bool reads(const std::string& pep, const std::string& tag)
  {
    auto fold = [](std::string x) {
      for (char& c : x) if (c == 'I') c = 'L';
      return x;
    };
    const std::string p = fold(pep), t = fold(tag);
    if (t.empty()) return false;
    const std::string r(t.rbegin(), t.rend());
    return p.find(t) != std::string::npos || p.find(r) != std::string::npos;
  }

  // ------------------------------------------------------------------ stats

  /// Fraction of peaks having another peak exactly one residue mass away.
  ///
  /// The structural quantity tagging actually consumes, and it needs no
  /// identifications -- which is what makes it usable as a calibration target on
  /// files whose peptides are unknown.
  double ladderEdgeDensity(const MSSpectrum& s, double tol, bool ppm)
  {
    const size_t n = s.size();
    if (n < 2) return 0.0;
    size_t with_edge = 0;
    for (size_t i = 0; i < n; ++i)
    {
      bool found = false;
      for (const char* p = ALPHABET; *p && !found; ++p)
      {
        const double target = s[i].getMZ() + residueMass(*p);
        const double t = ppm ? target * tol * 1e-6 : tol;
        const int j = s.findNearest(target, t);
        if (j >= 0 && static_cast<size_t>(j) != i) found = true;
      }
      if (found) ++with_edge;
    }
    return static_cast<double>(with_edge) / n;
  }

  int statsMode(const std::string& file, double tol, bool ppm)
  {
    // On-disc, because the largest cornerstone is 13 GB and a full load would
    // not fit. Requires an indexed mzML -- if the index is missing, say so
    // rather than silently falling back to a path that cannot work here.
    OnDiscMSExperiment ondisc;
    if (!ondisc.openFile(file))
    {
      std::fprintf(stderr, "cannot open (indexed mzML required): %s\n", file.c_str());
      return 1;
    }
    const Size n = ondisc.getNrSpectra();

    // Sample evenly across the run, never a prefix. An earlier measurement in
    // this project sampled the opening 300 spectra and understated ladder
    // density by a factor of two, because the early gradient is mostly short,
    // poorly-fragmenting peptides. Retention time is strongly confounded with
    // spectrum quality; a prefix is not a sample.
    const Size want = 2000;
    const Size stride = std::max<Size>(1, n / want);

    std::vector<size_t> npeaks;
    std::vector<double> density;
    std::map<int, size_t> charges;
    size_t ms2 = 0;

    for (Size i = 0; i < n; i += stride)
    {
      MSSpectrum s = ondisc.getSpectrum(i);
      if (s.getMSLevel() != 2 || s.empty()) continue;
      ++ms2;
      npeaks.push_back(s.size());
      s.sortByPosition();
      density.push_back(ladderEdgeDensity(s, tol, ppm));
      if (!s.getPrecursors().empty())
        ++charges[s.getPrecursors()[0].getCharge()];
      if (ms2 >= want) break;
    }

    if (npeaks.empty()) { std::fprintf(stderr, "no MS2 spectra found\n"); return 1; }

    auto pct = [](std::vector<size_t> v, double q) {
      std::sort(v.begin(), v.end());
      return v[std::min(v.size() - 1, static_cast<size_t>(q * v.size()))];
    };
    const double mean_density =
        std::accumulate(density.begin(), density.end(), 0.0) / density.size();

    std::printf("file            %s\n", file.c_str());
    std::printf("spectra total   %zu\n", static_cast<size_t>(n));
    std::printf("MS2 sampled     %zu (stride %zu)\n", ms2, static_cast<size_t>(stride));
    std::printf("peaks/MS2       p25=%zu med=%zu p75=%zu\n",
                pct(npeaks, 0.25), pct(npeaks, 0.50), pct(npeaks, 0.75));
    std::printf("ladder density  %.1f%%  (tol %.4g %s)\n",
                100.0 * mean_density, tol, ppm ? "ppm" : "Da");
    std::printf("charges        ");
    for (const auto& kv : charges)
      std::printf(" %d:%.0f%%", kv.first, 100.0 * kv.second / ms2);
    std::printf("\n");
    return 0;
  }

  // ------------------------------------------------------------------ synth

  /// One synthetic acquisition class.
  ///
  /// p_ion and n_noise are the two fitted knobs; everything else is fixed by
  /// chemistry or by the instrument's stated tolerance. Fits reproduce the
  /// ladder density and peak count that `stats` measures on the corresponding
  /// real file (doc/TEST-DATA.md), so the generator is anchored rather than
  /// invented -- but the fit is to TWO summary statistics, not to the spectra,
  /// and matching them does not make a synthetic spectrum real.
  struct Profile
  {
    const char* name;
    double p_ion;        ///< probability an individual b/y ion is observed
    double p_internal;   ///< probability of each candidate internal fragment
    double n_chimera;    ///< mean co-fragmenting peptides (Poisson)
    double n_noise;      ///< mean unstructured noise peaks (Poisson)
    double tol;
    bool   ppm;
    int    pep_len_min, pep_len_max;
  };

  // Fitted so `synthstats` lands on what `stats` measures on the real file.
  //
  // The first fit was badly wrong in an instructive way. It modelled a dense
  // Orbitrap spectrum as ~30 real ions in ~900 uniform noise peaks, which
  // reproduced the PEAK COUNT exactly and the ladder density not at all -- 33%
  // against a real 69%. Real dense spectra are not mostly noise: their peaks are
  // overwhelmingly fragment-like, from secondary series, neutral losses,
  // isotopes and co-isolated peptides. Matching one statistic and not the other
  // is what caught it, which is the argument for calibrating on two.
  const Profile PROFILES[] = {
      // Orbitrap Astral, label-free DDA. Real: med 1337 peaks, 69.4% ladder.
      {"astral",    0.88, 0.55, 13.0, 330.0, 20.0, true,  8, 25},
      // timsTOF ddaPASEF. Real: med 140 peaks, 38.5% ladder.
      {"ddapasef",  0.62, 0.14,  2.0,  45.0, 20.0, true,  7, 22},
      // diaTracer pseudo-MS2. Constructed spectra with a hard 500-peak cap; the
      // real file is no longer local, so this profile is anchored to the peak
      // count and ladder density recorded in doc/TEST-DATA.md, not re-measured.
      {"diatracer", 0.75, 0.34, 6.0, 110.0, 20.0, true,  7, 22},
      // Orbitrap Eclipse, ion-trap MS2. Real: med 267 peaks, 79.2% ladder at
      // 0.3 Da -- but at that tolerance a large share of those edges are
      // coincidence, so this profile's ladder figure is the weakest anchor here.
      {"iontrap",   0.72, 0.22, 2.0,  70.0,  0.3, false, 8, 24},
  };

  const Profile* findProfile(const std::string& n)
  {
    for (const Profile& p : PROFILES) if (n == p.name) return &p;
    return nullptr;
  }

  std::string randomPeptide(std::mt19937& rng, int lo, int hi)
  {
    std::discrete_distribution<int> res(FREQ, FREQ + 19);
    std::uniform_int_distribution<int> len(lo, hi);
    const int L = len(rng);
    std::string p;
    p.reserve(static_cast<size_t>(L));
    for (int i = 0; i < L - 1; ++i) p.push_back(ALPHABET[res(rng)]);
    // Tryptic C-terminus. Not cosmetic: K and R fix the y1 ion's mass, and a
    // fixed terminus is what real tryptic data gives the tagger.
    p.push_back((rng() & 1) ? 'K' : 'R');
    return p;
  }

  const double CO = 27.99491, H2O = 18.010565, NH3 = 17.026549, NEUTRON = 1.003355;

  /// Build one spectrum from a known peptide under a profile.
  ///
  /// Models b/y plus the secondary structure that dominates a real dense
  /// spectrum: a-ions, water and ammonia losses, +1 isotopes, doubly-charged
  /// fragments on 3+ precursors, and co-isolated peptides. Those are not
  /// decoration -- they are most of the peaks, and leaving them out understates
  /// ladder density by half (see PROFILES).
  MSSpectrum synth(const std::string& pep, const Profile& pf, std::mt19937& rng,
                   double& prec_mz, int& prec_charge)
  {
    const AASequence aa = AASequence::fromString(pep);
    const double M = aa.getMonoWeight();
    // Charge 2 or 3, as tryptic peptides overwhelmingly are.
    prec_charge = (std::uniform_real_distribution<double>(0, 1)(rng) < 0.62) ? 2 : 3;
    prec_mz = (M + prec_charge * Constants::PROTON_MASS_U) / prec_charge;

    MSSpectrum s;
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::lognormal_distribution<double> ion_int(0.0, 0.9), noise_int(-1.6, 0.9);
    std::normal_distribution<double> err(0.0, 0.35);   // in units of the tolerance
    const double hi_mz = prec_mz * prec_charge;

    auto tolAt = [&](double mz) { return pf.ppm ? mz * pf.tol * 1e-6 : pf.tol; };
    auto place = [&](double mz, double intensity) {
      if (mz < 100.0 || mz > hi_mz) return;
      s.emplace_back(mz + err(rng) * tolAt(mz), intensity);
    };

    /// Add one peptide's fragments at a given intensity scale.
    ///
    /// @param scale  1.0 for the precursor of interest, lower for co-isolated
    ///               peptides -- which is what lets a chimera be present without
    ///               automatically outranking the real answer.
    auto addPeptide = [&](const AASequence& q, int z, double scale) {
      for (size_t i = 1; i < q.size(); ++i)
      {
        // Each ion independently observed or not. Independent dropout is the
        // point: it creates the single-ion holes gapped tags exist to bridge,
        // at a rate the profile controls.
        for (int which = 0; which < 2; ++which)
        {
          if (u01(rng) >= pf.p_ion) continue;
          const bool is_b = which == 0;
          const double m = is_b ? q.getPrefix(i).getMonoWeight(Residue::BIon, 1)
                                : q.getSuffix(i).getMonoWeight(Residue::YIon, 1);
          const double I = ion_int(rng) * scale;
          place(m, I);
          if (u01(rng) < 0.55) place(m + NEUTRON, I * 0.35);
          if (is_b && u01(rng) < 0.30) place(m - CO, I * 0.30);      // a-ion
          if (u01(rng) < 0.22) place(m - H2O, I * 0.28);
          if (u01(rng) < 0.15) place(m - NH3, I * 0.25);
          // Doubly-charged fragments appear once the precursor is 3+.
          if (z >= 3 && u01(rng) < 0.30)
            place((m + Constants::PROTON_MASS_U) / 2.0, I * 0.45);
        }
      }

      // Internal fragments (b-type, from a second backbone cleavage).
      //
      // Needed to reach real peak counts without inventing implausible numbers
      // of co-isolated peptides: b/y plus losses accounts for roughly a third of
      // a dense Orbitrap spectrum's peaks, and internals are most of the rest.
      // They also carry genuine ladder structure -- an internal series steps by
      // residue masses like any other -- which is why adding them raises both
      // calibration targets at once, as the real data says it should.
      for (size_t i = 1; i + 2 < q.size(); ++i)
        for (size_t L = 2; L <= 6 && i + L < q.size(); ++L)
        {
          if (u01(rng) >= pf.p_internal) continue;
          const double m = q.getSubsequence(i, L).getMonoWeight(Residue::BIon, 1);
          place(m, ion_int(rng) * scale * 0.22);
        }
    };

    addPeptide(aa, prec_charge, 1.0);

    // Co-isolated peptides. Real chimericity was measured at ~32% of spectra
    // ungated; a dense Astral window holds several. Without them the benchmark
    // would never exercise a confident tag reading the WRONG peptide, which is
    // the failure mode a prefilter most needs to be honest about.
    std::poisson_distribution<int> n_chim(pf.n_chimera);
    const int nc = n_chim(rng);
    for (int c = 0; c < nc; ++c)
    {
      const AASequence oa = AASequence::fromString(
          randomPeptide(rng, pf.pep_len_min, pf.pep_len_max));
      // Co-isolated species are weaker than the selected precursor, but not
      // uniformly so -- occasionally one dominates, which is exactly when a
      // tagger confidently reads the wrong peptide.
      addPeptide(oa, 2, 0.15 + 0.5 * u01(rng));
    }

    std::poisson_distribution<int> n_nz(pf.n_noise);
    std::uniform_real_distribution<double> nz(100.0, std::max(200.0, hi_mz));
    const int nn = n_nz(rng);
    for (int i = 0; i < nn; ++i) s.emplace_back(nz(rng), noise_int(rng));

    s.sortByPosition();
    return s;
  }

  struct Score
  {
    size_t spectra = 0, tagged = 0;
    size_t rank1_ok = 0, top5_ok = 0, any_ok = 0;
    size_t rank1_gapped = 0, rank1_gapped_ok = 0;
    size_t rank1_contig = 0, rank1_contig_ok = 0;
    size_t tags = 0, tags_ok = 0;
  };

  void report(const char* label, const Score& c)
  {
    auto pc = [](size_t a, size_t b) { return b ? 100.0 * a / b : 0.0; };
    std::printf("%-10s %7zu %8.1f%% %8.1f%% %8.1f%% %8.1f%% %9.1f%% %9.1f%% %7.1f\n",
                label, c.spectra,
                pc(c.rank1_ok, c.spectra), pc(c.top5_ok, c.spectra),
                pc(c.any_ok, c.spectra), pc(c.tags_ok, c.tags),
                pc(c.rank1_contig_ok, c.rank1_contig),
                pc(c.rank1_gapped_ok, c.rank1_gapped),
                c.spectra ? static_cast<double>(c.tags) / c.spectra : 0.0);
    std::printf("%-10s rank-1 slots held by gapped tags: %.1f%%\n", "",
                pc(c.rank1_gapped, c.rank1_gapped + c.rank1_contig));
  }
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::fprintf(stderr,
        "usage: benchmark stats <file.mzML> [tol] [ppm|da]\n"
        "       benchmark synth <profile> [n] [tag_len] [gaps] [ext] [peaks] [gap_penalty]\n"
        "profiles:");
    for (const Profile& p : PROFILES) std::fprintf(stderr, " %s", p.name);
    std::fprintf(stderr, "\n");
    return 2;
  }

  if (std::strcmp(argv[1], "stats") == 0)
  {
    if (argc < 3) { std::fprintf(stderr, "stats needs a file\n"); return 2; }
    const double tol = argc > 3 ? std::atof(argv[3]) : 20.0;
    const bool ppm = argc > 4 ? std::strcmp(argv[4], "da") != 0 : true;
    return statsMode(argv[2], tol, ppm);
  }

  const bool synthstats = std::strcmp(argv[1], "synthstats") == 0;
  if (!synthstats && std::strcmp(argv[1], "synth") != 0)
  { std::fprintf(stderr, "unknown mode\n"); return 2; }
  if (argc < 3) { std::fprintf(stderr, "needs a profile\n"); return 2; }

  const Profile* pf = findProfile(argv[2]);
  if (!pf) { std::fprintf(stderr, "unknown profile %s\n", argv[2]); return 2; }

  // Same statistics as `stats`, on generated spectra. This is what closes the
  // calibration loop: a profile is only anchored if its numbers land on the
  // real file's, and without this mode "calibrated" would be an assertion.
  if (synthstats)
  {
    std::mt19937 rng(std::hash<std::string>{}(pf->name) & 0xffffffffu);
    std::vector<size_t> npeaks;
    std::vector<double> density;
    const int n = argc > 3 ? std::atoi(argv[3]) : 300;
    for (int i = 0; i < n; ++i)
    {
      double mz = 0; int z = 2;
      const MSSpectrum s = synth(randomPeptide(rng, pf->pep_len_min, pf->pep_len_max),
                                 *pf, rng, mz, z);
      npeaks.push_back(s.size());
      density.push_back(ladderEdgeDensity(s, pf->tol, pf->ppm));
    }
    auto pct = [](std::vector<size_t> v, double q) {
      std::sort(v.begin(), v.end());
      return v[std::min(v.size() - 1, static_cast<size_t>(q * v.size()))];
    };
    std::printf("profile         %s (synthetic, n=%d)\n", pf->name, n);
    std::printf("peaks/MS2       p25=%zu med=%zu p75=%zu\n",
                pct(npeaks, 0.25), pct(npeaks, 0.50), pct(npeaks, 0.75));
    std::printf("ladder density  %.1f%%  (tol %.4g %s)\n",
                100.0 * std::accumulate(density.begin(), density.end(), 0.0) / density.size(),
                pf->tol, pf->ppm ? "ppm" : "Da");
    return 0;
  }

  const int n_spec   = argc > 3 ? std::atoi(argv[3]) : 2000;
  const int tag_len  = argc > 4 ? std::atoi(argv[4]) : 4;
  const int gaps     = argc > 5 ? std::atoi(argv[5]) : 0;
  const int ext      = argc > 6 ? std::atoi(argv[6]) : 0;
  const int peaks    = argc > 7 ? std::atoi(argv[7]) : 100;
  const double gpen  = argc > 8 ? std::atof(argv[8]) : -1.0;
  const double maxev = argc > 9 ? std::atof(argv[9]) : 20.0;
  const int max_tags = argc > 10 ? std::atoi(argv[10]) : 50;

  FasTag::Param p;
  p.tag_length = tag_len;
  p.max_gaps = gaps;
  p.max_extension = ext;
  p.max_peak_count = static_cast<size_t>(peaks);
  p.frag_tol = pf->tol;
  p.tol_ppm = pf->ppm;
  p.complement_tol = pf->ppm ? 0.02 : pf->tol;
  p.deisotope = true;
  if (gpen >= 0) p.gap_penalty = gpen;   // negative keeps the shipped default
  p.max_evalue = maxev;
  p.max_tag_count = max_tags;
  const FasTag::Tables tables(p);

  Score sc;
  // One RNG for the whole run, seeded from the profile name, so a given
  // (profile, n) is the same corpus on every machine and every run.
  std::mt19937 rng(std::hash<std::string>{}(pf->name) & 0xffffffffu);

  for (int i = 0; i < n_spec; ++i)
  {
    const std::string pep = randomPeptide(rng, pf->pep_len_min, pf->pep_len_max);
    double prec_mz = 0; int z = 2;
    const MSSpectrum s = synth(pep, *pf, rng, prec_mz, z);

    const std::vector<FasTag::Tag> tags = FasTag::tagSpectrum(s, prec_mz, z, p, tables);
    ++sc.spectra;
    if (tags.empty()) continue;
    ++sc.tagged;

    // Per-spectrum flags, set inside the loop and counted once after it. An
    // earlier version incremented the top5/any counters inside the loop under a
    // "counter < spectra" guard, which double-counts a spectrum with two correct
    // tags and silently caps at 100%.
    bool hit_top5 = false, hit_any = false;
    for (size_t t = 0; t < tags.size(); ++t)
    {
      const bool ok = reads(pep, tags[t].seq);
      ++sc.tags;
      if (ok) ++sc.tags_ok;
      if (t == 0)
      {
        if (ok) ++sc.rank1_ok;
        if (tags[t].gapped) { ++sc.rank1_gapped; if (ok) ++sc.rank1_gapped_ok; }
        else                { ++sc.rank1_contig; if (ok) ++sc.rank1_contig_ok; }
      }
      if (ok && t < 5) hit_top5 = true;
      if (ok) hit_any = true;
    }
    if (hit_top5) ++sc.top5_ok;
    if (hit_any) ++sc.any_ok;
  }

  std::printf("profile=%s n=%d TL=%d gaps=%d ext=%d peaks=%d gap_penalty=%.3g\n",
              pf->name, n_spec, tag_len, gaps, ext, peaks, gpen);
  std::printf("%-10s %7s %9s %9s %9s %9s %10s %10s %7s\n",
              "", "spectra", "rank1", "top5", "any", "tagprec",
              "r1-contig", "r1-gapped", "tags/sp");
  report("result", sc);

  // Exit non-zero on an obviously broken run, so ctest can gate on this.
  //
  // Deliberately a floor, not a threshold on accuracy: pinning a quality number
  // would make every scoring change a CI failure and train people to edit the
  // expectation. What this catches is the tool falling over or going silent on a
  // whole acquisition class -- which is exactly what happened when `-deisotope`
  // threw on the iontrap profile.
  if (sc.tagged == 0)
  {
    std::fprintf(stderr, "FAIL: profile %s produced no tags on any of %d spectra\n",
                 pf->name, n_spec);
    return 1;
  }
  return 0;
}
