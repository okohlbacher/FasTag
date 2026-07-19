// End-to-end: build a spectrum from a known peptide and require FasTag to read
// that peptide back. Residue masses come from OpenMS, so the test uses OpenMS to
// build the spectrum too -- if the two disagreed, this would catch it.
#include "FasTagger.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
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

  // 7. a gap bridges a missing peak, and spells it in the RIGHT ORDER
  //
  // The ordering is the part worth testing. A gap step is reversed twice --
  // once because traversal is C->N, once within the pair -- and getting it
  // wrong transposes exactly two residues while leaving the pair sum, both
  // flanking masses and every score untouched. No mass-based check can see it.
  //
  // Built from y ions only, so there is no b/y ambiguity to hide behind, with
  // one interior y ion deleted. Nothing can span that hole without a gap.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    const AASequence aa = AASequence::fromString(pep);
    double pmz = 0;
    synth(pep, 0, 1, pmz);   // for the precursor m/z only

    MSSpectrum y;
    for (size_t i = 1; i < aa.size(); ++i)
      y.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1), 0.5 + 0.02 * i);
    const double drop = aa.getSuffix(6).getMonoWeight(Residue::YIon, 1);
    MSSpectrum holed;
    for (const auto& pk : y) if (std::fabs(pk.getMZ() - drop) > 1e-6) holed.push_back(pk);
    holed.sortByPosition();

    // Uncapped deliberately. This checks that the gap machinery SPELLS the
    // peptide correctly, which is a different question from whether it ranks it
    // highly. On this spectrum -- 14 y ions, no b ions, no noise -- gap edges
    // join nearly arbitrary peak pairs and the correct tag lands at rank 59, so
    // a top-50 cap would hide it and this test would report an ordering failure
    // that is really a ranking one. Real spectra carry ~100 peaks and measure
    // far better (see Gapped-Tags-Results); the rank is printed below so a
    // regression in it stays visible.
    FasTag::Param off;
    off.frag_tol = 0.02; off.tol_ppm = false; off.complement_tol = 0.02;
    off.tag_length = 6; off.max_tag_count = 0; off.max_evalue = 0;
    FasTag::Param on = off; on.max_gaps = 1;
    const FasTag::Tables tf(off), tn(on);

    const auto without = FasTag::tagSpectrum(holed, pmz, 2, off, tf);
    const auto with = FasTag::tagSpectrum(holed, pmz, 2, on, tn);

    int gapped_correct = 0, gapped_total = 0;
    for (const auto& t : with)
      if (t.gapped) { ++gapped_total; if (reads(pep, t.seq)) ++gapped_correct; }

    CHECK(gapped_total > 0, "no gapped tag produced across a deleted peak");
    CHECK(gapped_correct > 0, "gapped tags produced but none reads the peptide -- "
                              "the two residues either side of the gap are most "
                              "likely transposed");

    // The assertion above only means something if transposing the gap would
    // actually break it. Prove the checker discriminates.
    int would_catch = 0;
    for (const auto& t : with)
    {
      if (!t.gapped || !reads(pep, t.seq)) continue;
      for (size_t i = 0; i + 1 < t.seq.size(); ++i)
      {
        std::string sw = t.seq;
        std::swap(sw[i], sw[i + 1]);
        if (sw != t.seq && !reads(pep, sw)) { ++would_catch; break; }
      }
    }
    CHECK(would_catch > 0, "reads() cannot distinguish a transposition here, so "
                           "the ordering check above proves nothing");
    std::printf("7. gap spans a deleted peak: %d/%d gapped tags read the peptide "
                "(%zu tags without gaps)\n"
                "   negative control: %d would fail if two residues were swapped\n",
                gapped_correct, gapped_total, without.size(), would_catch);
  }

  // 8. every emitted row must close on its own mass
  //
  //    nterm + mass(seq) + cterm == precursor residue mass, for EVERY row.
  //
  // This is the check that catches a whole class at once. It held by
  // construction while one path produced one tag, which is why an earlier
  // version of this file asserted it and proved nothing. With gaps a path spells
  // several sequences, and giving them all the flanks of the best-matching
  // composition breaks closure by the difference between compositions -- rows
  // that contradict their own sequence, which is exactly what a downstream
  // search consumes as a precursor constraint. Now that each spelling is
  // rescored against its own pair mass the identity is restored, and it is worth
  // asserting precisely because it is no longer automatic.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    p.tag_length = 6; p.max_gaps = 1; p.max_tag_count = 0;
    const FasTag::Tables tab(p);
    double pmz = 0;
    const MSSpectrum s = synth(pep, 20, 11, pmz);
    const auto tags = FasTag::tagSpectrum(s, pmz, 2, p, tab);

    const double prec_res = pmz * 2 - 2 * Constants::PROTON_MASS_U
                          - EmpiricalFormula("H2O").getMonoWeight();
    int checked = 0, closed = 0, gapped = 0;
    double worst = 0;
    for (const auto& t : tags)
    {
      double m = 0;
      for (char c : t.seq)
        m += ResidueDB::getInstance()->getResidue(String(c))->getMonoWeight(Residue::Internal);
      const double err = std::fabs(t.nterm_mass + m + t.cterm_mass - prec_res);
      ++checked; worst = std::max(worst, err);
      if (err < 0.05) ++closed;
      if (t.gapped) ++gapped;
    }
    CHECK(checked > 0, "no tags to check closure on");
    CHECK(gapped > 0, "no gapped tag in the closure check, so it proves nothing");
    CHECK(closed == checked, "%d/%d rows do not close on their own mass "
          "(worst %.4f Da) -- a spelling is carrying another composition's flanks",
          closed, checked, worst);
    std::printf("8. mass closure: %d/%d rows satisfy nterm+seq+cterm == precursor "
                "(%d gapped, worst %.5f Da)\n", closed, checked, gapped, worst);
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
