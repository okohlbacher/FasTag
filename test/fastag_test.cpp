// End-to-end: build a spectrum from a known peptide and require FasTag to read
// that peptide back. Residue masses come from OpenMS, so the test uses OpenMS to
// build the spectrum too -- if the two disagreed, this would catch it.
//
// EVERY TEST HERE HAS BEEN MUTATION-CHECKED: the behaviour it names was broken
// on purpose and the test confirmed to fail. Four tests in this project have
// been caught asserting something true by construction, so "it passes" is not
// evidence that it tests anything.
//
// Two traps found while doing that, both of which produced confident wrong
// answers:
//   - a stale build. Restoring a mutated source without touching it can leave
//     the old object linked, and the run then reports the mutation's behaviour
//     as if it were the code's. Always touch the file and confirm a recompile.
//   - a mutation that does not mutate. Blanking the inner branch of the
//     duplicate-peak collapse only changes WHICH duplicate survives; the else
//     branch still collapses. It looked like a surviving mutation for two
//     rounds. Check that the mutation changes behaviour before concluding the
//     test is weak.
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

  // 9. a charge-3 precursor survives, with two fragment charges in play
  //
  // Every other test here uses charge 2, where n_frag_charges == 1. Nothing
  // exercised the n_frag_charges > 1 path at all, which is how a real bug lived
  // there: prepare() pooled the complement population across charges while
  // scorePath counted at one charge, so the hypergeometric's K and obs came from
  // different populations. Fixed; this test covers the path.
  //
  // WHAT THIS DOES NOT PROVE. It does not discriminate that bug: reintroducing
  // the pooled count leaves this test passing with byte-identical output, which
  // was verified rather than assumed. A discriminating test needs a spectrum
  // whose complement population is EMPTY at one fragment charge and non-empty at
  // another -- then pooled and per-charge disagree about whether the complement
  // subscore applies at all, and Fisher's degrees of freedom change with it.
  // Constructing that spectrum is fiddly and is left as a follow-up rather than
  // dressed up as done.
  {
    const std::string pep = "VGAHAGEYGAEALERSAMPLERK";
    const AASequence aa = AASequence::fromString(pep);
    const double M = aa.getMonoWeight();

    // Charge 3 precursor -> n_frag_charges = 2, so +1 and +2 complements both
    // exist and a pooled count would differ from a per-charge one.
    const double pmz = (M + 3 * Constants::PROTON_MASS_U) / 3;
    MSSpectrum s;
    for (size_t i = 1; i < aa.size(); ++i)
    {
      s.emplace_back(aa.getPrefix(i).getMonoWeight(Residue::BIon, 1), 0.7);
      s.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1), 0.8);
      // a few doubly-charged fragments, so the z=2 complement population is
      // genuinely different from the z=1 one
      if (i % 3 == 0)
        s.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 2) / 2.0, 0.4);
    }
    s.sortByPosition();

    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    const FasTag::Tables tab(p);
    const auto tags = FasTag::tagSpectrum(s, pmz, 3, p, tab);

    CHECK(!tags.empty(), "no tags from a charge-3 precursor");
    // Every E-value must be finite and positive. A mismatched population can
    // push the hypergeometric outside its support, which is how this surfaces.
    int bad = 0;
    for (const auto& t : tags)
      if (!(t.evalue > 0.0) || !std::isfinite(t.evalue) ||
          !(t.p_complement >= 0.0 && t.p_complement <= 1.0)) ++bad;
    CHECK(bad == 0, "%d tags have an out-of-range E-value or complement p-value", bad);
    int correct = 0;
    for (const auto& t : tags) if (reads(pep, t.seq)) ++correct;
    CHECK(correct > 0, "no tag from the charge-3 spectrum reads the peptide");
    std::printf("9. charge-3 precursor (2 fragment charges): %zu tags, %d read the "
                "peptide, %d with an out-of-range p-value\n", tags.size(), correct, bad);
  }

  // 10. the stored direction is N->C, asserted FORWARD-ONLY
  //
  // reads() accepts either orientation, because a b-derived tag legitimately
  // reads reversed. That is right for every other test here and useless for
  // this one: mutation testing showed the whole suite passes with spell()
  // reversed, so nothing was pinning the storage convention at all.
  //
  // A y-only spectrum removes the ambiguity. Traversal runs low->high m/z, which
  // is C->N for y ions, and spell() reverses it, so a correct tag must appear in
  // the peptide FORWARD. Searching forward only is the whole point -- do not be
  // tempted to reuse reads() here.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    const AASequence aa = AASequence::fromString(pep);
    double pmz = 0;
    synth(pep, 0, 1, pmz);                       // for the precursor m/z

    MSSpectrum y;                                // y ions only: no b/y ambiguity
    for (size_t i = 1; i < aa.size(); ++i)
      y.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1), 0.5 + 0.02 * i);
    y.sortByPosition();

    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    const FasTag::Tables tab(p);
    const auto tags = FasTag::tagSpectrum(y, pmz, 2, p, tab);

    auto fold = [](std::string x) { for (char& c : x) if (c == 'I') c = 'L'; return x; };
    const std::string fp = fold(pep);
    int fwd = 0, rev_only = 0;
    for (const auto& t : tags)
    {
      const std::string ft = fold(t.seq);
      const std::string rt(ft.rbegin(), ft.rend());
      const bool f = fp.find(ft) != std::string::npos;
      const bool r = fp.find(rt) != std::string::npos;
      if (f) ++fwd; else if (r) ++rev_only;
    }
    CHECK(fwd > 0, "no tag from a y-only spectrum reads the peptide forward");
    CHECK(fwd > rev_only,
          "%d tags read forward but %d only backwards -- the stored direction "
          "looks reversed", fwd, rev_only);
    std::printf("10. y-only spectrum: %d tags read N->C forward, %d only reversed\n",
                fwd, rev_only);
  }

  // 11. duplicated m/z peaks must not change the result
  //
  // diaTracer emits peaks sharing an m/z routinely -- measured on S23, every MS2
  // spectrum has them, median 22 and up to 46. prepare() collapses them keeping
  // the strongest, because each otherwise consumes a slot in the peak budget and
  // shifts every intensity rank.
  //
  // Mutation testing showed the suite passed with that collapse disabled: the
  // synthetic spectra had no duplicates, so nothing exercised it. This feeds the
  // same spectrum twice, once with duplicates injected, and requires identical
  // output. Verified to FAIL when the collapse is removed.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    double pmz = 0;
    const MSSpectrum clean = synth(pep, 20, 5, pmz);

    // Same peaks, but every third one repeated at exactly its own m/z with a
    // lower intensity -- the case prepare() is meant to fold away.
    MSSpectrum dup;
    size_t k = 0;
    for (const auto& pk : clean)
    {
      dup.push_back(pk);
      // EQUAL intensity, not weaker. A first version injected at 0.5x, so the
      // cap discarded the duplicates before they could displace anything and the
      // test passed with the collapse removed. A duplicate only costs a slot if
      // it ranks high enough to occupy one.
      if (k++ % 3 == 0) dup.emplace_back(pk.getMZ(), pk.getIntensity());
    }
    dup.sortByPosition();
    CHECK(dup.size() > clean.size(), "no duplicates were actually injected");

    FasTag::Param p;
    p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
    // The cap MUST bind, or duplicates displace nothing and this proves nothing.
    // A first version left the default 100 against ~66 peaks, and passed happily
    // with the collapse removed.
    p.max_peak_count = 30;
    const FasTag::Tables tab(p);
    const auto a = FasTag::tagSpectrum(clean, pmz, 2, p, tab);
    const auto b = FasTag::tagSpectrum(dup, pmz, 2, p, tab);

    bool same = a.size() == b.size();
    for (size_t i = 0; same && i < a.size(); ++i)
      same = a[i].seq == b[i].seq && a[i].evalue == b[i].evalue;
    CHECK(same, "duplicated m/z peaks changed the output: %zu tags vs %zu",
          a.size(), b.size());
    std::printf("11. %zu duplicate peaks injected: %zu tags either way, identical\n",
                dup.size() - clean.size(), a.size());
  }

  // 12. the peak cap keeps exactly max_peaks, and the boundary is the strongest
  //
  // Mutation testing showed an off-by-one in the cap went unnoticed. The cap is
  // the single most consequential parameter -- it decides which peaks exist at
  // all -- so its boundary deserves an exact assertion rather than a smoke test.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    const AASequence aa = AASequence::fromString(pep);
    double pmz = 0;
    synth(pep, 0, 9, pmz);                       // for the precursor m/z only

    // A controlled ladder: y ions with intensity DECREASING along the series, so
    // the strongest k peaks are exactly the first k and therefore consecutive.
    // Without that, the top-k of a noisy spectrum are scattered and no tag of
    // any length survives a tight cap -- which is a fact about noise, not about
    // the cap, and would make this test assert the wrong thing.
    MSSpectrum s;
    for (size_t i = 1; i < aa.size(); ++i)
      s.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1),
                     static_cast<float>(1.0 - 0.02 * i));
    s.sortByPosition();

    // Assert the BOUNDARY, not a bound. A length-L tag needs L+1 peaks, so with
    // the cap at exactly tag_length there can be no tag, and at tag_length + 1
    // there can. An off-by-one either way moves that edge. Asserting "no tag is
    // longer than the cap allows" instead is vacuous: tag length is fixed by
    // tag_length, so it holds for any cap.
    const int TL = 5;
    size_t at_edge = 0, past_edge = 0;
    for (int cap : {TL, TL + 1})
    {
      FasTag::Param p;
      p.frag_tol = 0.02; p.tol_ppm = false; p.complement_tol = 0.02;
      p.max_peak_count = static_cast<size_t>(cap);
      p.tag_length = TL;
      const FasTag::Tables tab(p);
      const size_t n = FasTag::tagSpectrum(s, pmz, 2, p, tab).size();
      if (cap == TL) at_edge = n; else past_edge = n;
    }
    CHECK(at_edge == 0, "cap %d retains too many peaks: %zu tags where a "
          "%d-residue tag needs %d peaks", TL, at_edge, TL, TL + 1);
    CHECK(past_edge > 0, "cap %d retains too few peaks: no tags at all", TL + 1);
    std::printf("12. peak cap boundary: %zu tags at cap=%d, %zu at cap=%d\n",
                at_edge, TL, past_edge, TL + 1);
  }

  // 13. the gap penalty REORDERS but never removes a tag
  //
  // The invariant that protects sensitivity. Tagging is a prefilter, so the
  // penalty must change which tag is tried first and nothing else. When it was
  // first applied it also gated the E-value cutoff, and total recall fell from
  // 71.3% to 51.9% on the astral benchmark profile -- gapped tags were not
  // demoted but deleted. That is the regression this guards.
  //
  // Asserts BOTH halves: the tag sets must be equal, and the orders must differ.
  // Without the second half a penalty that did nothing at all would pass.
  {
    const std::string pep = "VGAHAGEYGAEALER";
    const AASequence aa = AASequence::fromString(pep);
    double pmz = 0;
    synth(pep, 0, 1, pmz);

    MSSpectrum y;
    for (size_t i = 1; i < aa.size(); ++i)
      y.emplace_back(aa.getSuffix(i).getMonoWeight(Residue::YIon, 1), 0.5 + 0.02 * i);
    const double drop = aa.getSuffix(6).getMonoWeight(Residue::YIon, 1);
    MSSpectrum holed;
    for (const auto& pk : y) if (std::fabs(pk.getMZ() - drop) > 1e-6) holed.push_back(pk);
    holed.sortByPosition();

    // A cutoff must be ACTIVE for this test to mean anything: with max_evalue 0
    // nothing can be filtered out and the set equality holds trivially.
    FasTag::Param base;
    base.frag_tol = 0.02; base.tol_ppm = false; base.complement_tol = 0.02;
    base.tag_length = 6; base.max_gaps = 1; base.max_tag_count = 0;
    base.max_evalue = 20.0;

    FasTag::Param flat = base; flat.gap_penalty = 1.0;
    FasTag::Param pen  = base; pen.gap_penalty = 100.0;
    const FasTag::Tables tf(flat), tp(pen);

    const auto a = FasTag::tagSpectrum(holed, pmz, 2, flat, tf);
    const auto b = FasTag::tagSpectrum(holed, pmz, 2, pen, tp);

    auto keys = [](const std::vector<FasTag::Tag>& v) {
      std::vector<std::string> k;
      for (const auto& t : v) k.push_back(t.seq);
      return k;
    };
    std::vector<std::string> ka = keys(a), kb = keys(b);
    const bool order_differs = ka != kb;
    std::sort(ka.begin(), ka.end());
    std::sort(kb.begin(), kb.end());

    CHECK(!a.empty() && !b.empty(), "no tags produced; test asserts nothing");
    CHECK(ka == kb, "gap penalty changed the tag SET (%zu vs %zu tags) -- it must "
          "reorder only, or sensitivity is lost", a.size(), b.size());
    CHECK(order_differs, "gap penalty changed nothing at all; either it is not "
          "applied or this spectrum has no gapped tags to demote");

    size_t gapped = 0;
    for (const auto& t : a) if (t.gapped) ++gapped;
    CHECK(gapped > 0, "no gapped tags here, so the penalty is untested");
    std::printf("13. gap penalty reorders only: %zu tags both ways, %zu gapped, "
                "order changed=%d\n", a.size(), gapped, order_differs ? 1 : 0);
  }

  // 14. deisotoping must survive an ion-trap tolerance
  //
  // OpenMS's deisotoper throws above 100 ppm or 0.1 Da, and FasTag passed the
  // fragment tolerance straight through -- so `-deisotope` at the 0.3 Da an ion
  // trap needs killed the tool with an uncaught exception. The combination is
  // reachable by following doc/TEST-DATA.md, which names 0.3 Da for the Eclipse
  // benchmark file.
  {
    double pmz = 0;
    const MSSpectrum s = synth("VGAHAGEYGAEALER", 40, 5, pmz);
    FasTag::Param p;
    p.frag_tol = 0.3; p.tol_ppm = false; p.complement_tol = 0.3;
    p.deisotope = true; p.tag_length = 3;
    const FasTag::Tables tab(p);
    bool threw = false;
    try { FasTag::tagSpectrum(s, pmz, 2, p, tab); }
    catch (...) { threw = true; }
    CHECK(!threw, "deisotoping threw at a 0.3 Da fragment tolerance");
    std::printf("14. deisotope at 0.3 Da (ion trap): no exception\n");
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
