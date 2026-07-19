// End-to-end: build a spectrum from a known peptide and require FasTag to read
// that peptide back. Residue masses come from OpenMS, so the test uses OpenMS to
// build the spectrum too -- if the two disagreed, this would catch it.
#include "FasTagger.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cstdio>
#include <random>
#include <string>

using namespace OpenMS;

static int failures = 0;
#define CHECK(c, ...) do { if (!(c)) { ++failures; std::printf("FAIL %d: ", __LINE__); \
  std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

/// Ideal singly-charged b and y ions, plus optional noise.
static MSSpectrum synth(const std::string& pep, int noise, unsigned seed, double& prec_mz)
{
  const AASequence aa = AASequence::fromString(pep);
  const double M = aa.getMonoWeight();
  prec_mz = (M + 2 * Constants::PROTON_MASS_U) / 2;

  MSSpectrum s;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> hi(0.6, 1.0), lo(0.02, 0.3),
      nz(120.0, std::max(200.0, M - 40.0));

  for (size_t i = 1; i < aa.size(); ++i)
  {
    s.emplace_back(aa.getPrefix(i).getMonoWeight(Residue::BIon, 1), hi(rng));
    s.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1), hi(rng));
  }
  for (int i = 0; i < noise; ++i) s.emplace_back(nz(rng), lo(rng));
  s.sortByPosition();
  return s;
}

/// Is `tag` a subsequence of `pep`? Either orientation, I/L folded -- tags are
/// stored N->C under a y-ion assumption, so b-derived tags read reversed.
static bool reads(const std::string& pep, std::string tag)
{
  auto fold = [](std::string x) { for (char& c : x) if (c == 'I') c = 'L'; return x; };
  const std::string p = fold(pep), t = fold(tag);
  const std::string r(t.rbegin(), t.rend());
  return p.find(t) != std::string::npos || p.find(r) != std::string::npos;
}

int main()
{
  const std::vector<std::string> peps = {
      "SAMPLER", "PEPTIDEK", "VGAHAGEYGAEALER", "DFTPAELR", "YLGYLEQLLR"};

  // 1. clean spectra: the top-ranked tag must read the peptide
  {
    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    const FasTag::Tables tab(p);
    int ok = 0, n = 0;
    for (const auto& pep : peps)
    {
      double pmz = 0;
      const MSSpectrum s = synth(pep, 0, 42, pmz);
      const auto tags = FasTag::tagSpectrum(s, pmz, 2, p, tab);
      CHECK(!tags.empty(), "no tags for %s", pep.c_str());
      if (tags.empty()) continue;
      ++n;
      if (reads(pep, tags[0].seq)) ++ok;
      else std::printf("   note: top tag for %s was '%s'\n", pep.c_str(), tags[0].seq.c_str());
    }
    CHECK(ok == n, "top tag correct for %d/%d", ok, n);
    std::printf("1. top-ranked tag reads the peptide for %d/%d clean spectra\n", ok, n);
  }

  // 2. noise: precision of the reported set must stay high
  {
    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02; p.max_tag_count = 10;
    const FasTag::Tables tab(p);
    int good = 0, all = 0;
    for (const auto& pep : peps)
    {
      double pmz = 0;
      const MSSpectrum s = synth(pep, 60, 7, pmz);
      for (const auto& t : FasTag::tagSpectrum(s, pmz, 2, p, tab))
      { ++all; if (reads(pep, t.seq)) ++good; }
    }
    const double prec = all ? double(good) / all : 0;
    CHECK(prec > 0.5, "precision %.2f too low with noise", prec);
    std::printf("2. with 60 noise peaks: %d/%d tags correct (%.0f%%)\n", good, all, 100 * prec);
  }

  // 3. flanking masses place the tag on the peptide
  //    nterm + tag + cterm == peptide residue mass holds by construction, so it
  //    proves nothing; the real check is that nterm lands on an actual prefix.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    const AASequence aa = AASequence::fromString(pep);
    std::vector<double> prefix{0.0};
    for (size_t i = 1; i <= aa.size(); ++i) prefix.push_back(aa.getPrefix(i).getMonoWeight(Residue::Internal));

    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    const FasTag::Tables tab(p);
    double pmz = 0;
    const MSSpectrum s = synth(pep, 0, 1, pmz);
    const auto tags = FasTag::tagSpectrum(s, pmz, 2, p, tab);
    int placed = 0;
    for (const auto& t : tags)
      for (double a : prefix)
        if (std::fabs(a - t.nterm_mass) < 0.03) { ++placed; break; }
    // Only y-derived tags can land on a prefix. Flanking masses assume the
    // low-m/z peak is a y ion, so b-derived tags carry a one-water offset --
    // measured at ~57% of output, which caps this fraction near 43%.
    CHECK(placed > 0, "no tag placed on a real prefix mass");
    std::printf("3. flanking masses: %d/%zu tags land on a real peptide prefix\n"
                "   (the rest are b-derived, offset by one water by construction)\n",
                placed, tags.size());
  }

  // 4. extension lengthens tags without wrecking them
  {
    const std::string pep = "VGAHAGEYGAEALER";
    FasTag::Param base;
    base.frag_tol = 0.02; base.tol_ppm = false; base.complement_tol = 0.02;
    FasTag::Param ext = base; ext.max_extension = 4;
    const FasTag::Tables tb(base), te(ext);
    double pmz = 0;
    const MSSpectrum s = synth(pep, 0, 1, pmz);
    size_t lb = 0, le = 0;
    for (const auto& t : FasTag::tagSpectrum(s, pmz, 2, base, tb)) lb = std::max(lb, t.seq.size());
    const auto grown = FasTag::tagSpectrum(s, pmz, 2, ext, te);
    for (const auto& t : grown) le = std::max(le, t.seq.size());
    CHECK(lb == 3, "seed length should be exactly 3, got %zu", lb);
    CHECK(le > lb, "extension produced no longer tag (%zu vs %zu)", le, lb);
    int tr = 0, tot = 0;
    for (const auto& t : grown) if (t.extended) { ++tot; if (reads(pep, t.seq)) ++tr; }
    CHECK(tot == 0 || tr > 0, "no extended tag was correct");
    std::printf("4. extension: max length %zu -> %zu; %d/%d extended tags read the peptide\n",
                lb, le, tr, tot);
  }

  // 5. determinism -- the tie-breaks exist for this
  {
    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02; p.max_extension = 4;
    const FasTag::Tables tab(p);
    double pmz = 0;
    const MSSpectrum s = synth("VGAHAGEYGAEALER", 40, 3, pmz);
    const auto a = FasTag::tagSpectrum(s, pmz, 2, p, tab);
    const auto b = FasTag::tagSpectrum(s, pmz, 2, p, tab);
    bool same = a.size() == b.size();
    for (size_t i = 0; same && i < a.size(); ++i)
      same = a[i].seq == b[i].seq && a[i].evalue == b[i].evalue;
    CHECK(same, "output is not reproducible across runs");
    std::printf("5. repeated runs produce identical output\n");
  }

  // 6. degenerate input must not crash or fabricate
  {
    FasTag::Param p;
    const FasTag::Tables tab(p);
    MSSpectrum empty;
    CHECK(FasTag::tagSpectrum(empty, 500, 2, p, tab).empty(), "empty spectrum");
    MSSpectrum tiny;
    tiny.emplace_back(100.0, 1.0);
    tiny.emplace_back(200.0, 1.0);
    CHECK(FasTag::tagSpectrum(tiny, 500, 2, p, tab).empty(), "too few peaks");
    double pmz = 0;
    const MSSpectrum s = synth("SAMPLER", 0, 3, pmz);
    FasTag::tagSpectrum(s, pmz, 0, p, tab);      // unknown charge
    FasTag::tagSpectrum(s, 0.0, 2, p, tab);      // no precursor
    std::printf("6. empty, too-few-peaks, unknown charge and missing precursor all survive\n");
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
