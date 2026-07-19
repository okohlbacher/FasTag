// Self-check for fastag.h. Build & run:
//   c++ -std=c++17 -O2 -o fastag_test fastag_test.cpp && ./fastag_test
//
// The load-bearing test is end-to-end: synthesise a b/y spectrum from a known
// peptide and require FasTag to recover that peptide's real subsequences.
#include <FasTag/fastag.h>

#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <string>

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++failures; std::printf("  FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static double residueMass(char c)
{
  for (const auto& r : fastag::residues()) if (r.code == c) return r.mass;
  if (c == 'I') return 113.08406;
  return 0.0;
}

static double peptideMass(const std::string& pep)
{
  double m = fastag::WATER;
  for (char c : pep) m += residueMass(c);
  return m;
}

/// Ideal b/y ion series, singly charged, with optional noise peaks.
static void synth(const std::string& pep, std::vector<double>& mz,
                  std::vector<double>& intensity, int noise, uint32_t seed)
{
  const double M = peptideMass(pep);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> noise_mz(150.0, M - 50.0);
  std::uniform_real_distribution<double> lo(0.05, 0.35);

  double b = 0.0;
  for (size_t i = 0; i + 1 < pep.size(); ++i)
  {
    b += residueMass(pep[i]);
    mz.push_back(b + fastag::PROTON);
    intensity.push_back(1.0);
  }
  double y = fastag::WATER;
  for (size_t i = pep.size(); i-- > 1;)
  {
    y += residueMass(pep[i]);
    mz.push_back(y + fastag::PROTON);
    intensity.push_back(1.0);
  }
  for (int i = 0; i < noise; ++i) { mz.push_back(noise_mz(rng)); intensity.push_back(lo(rng)); }

  // sort by m/z, as a centroided spectrum would be
  std::vector<size_t> idx(mz.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b2) { return mz[a] < mz[b2]; });
  std::vector<double> m2, i2;
  for (size_t i : idx) { m2.push_back(mz[i]); i2.push_back(intensity[i]); }
  mz = m2; intensity = i2;
}

/// Is `tag` a real subsequence of `pep` (in either orientation, I/L-insensitive)?
static bool isTrueTag(const std::string& pep, const std::string& tag)
{
  auto norm = [](std::string s) { for (char& c : s) if (c == 'I') c = 'L'; return s; };
  const std::string p = norm(pep), t = norm(tag);
  std::string rt(t.rbegin(), t.rend());
  return p.find(t) != std::string::npos || p.find(rt) != std::string::npos;
}

int main()
{
  // ---------------------------------------------------------------- unit bits
  {
    // chi2 survival, even DOF, against known values.
    CHECK(std::fabs(fastag::chi2SfEven(0.0, 2) - 1.0) < 1e-12, "chi2 sf at 0");
    // DOF=2: P(X>x) = exp(-x/2)
    CHECK(std::fabs(fastag::chi2SfEven(4.0, 1) - std::exp(-2.0)) < 1e-12, "chi2 dof2");
    // DOF=4: P(X>x) = exp(-x/2)(1+x/2)
    CHECK(std::fabs(fastag::chi2SfEven(4.0, 2) - std::exp(-2.0) * 3.0) < 1e-12, "chi2 dof4");
    CHECK(fastag::chi2SfEven(1e5, 3) >= 0.0, "chi2 huge x must not go negative");
    // monotone decreasing in x
    double prev = 1.0;
    for (double x = 0.5; x < 60; x += 0.5)
    { const double v = fastag::chi2SfEven(x, 3); CHECK(v <= prev + 1e-15, "chi2 monotone at %.1f", x); prev = v; }
    std::printf("1. chi-square survival (even DOF) correct and monotone\n");
  }
  {
    // Hypergeometric tail: P(X>=k) must be 1 at k=0, and sum to sane values.
    CHECK(std::fabs(fastag::hyperSf(100, 20, 4, 0) - 1.0) < 1e-12, "hyper k=0");
    CHECK(fastag::hyperSf(100, 20, 4, 5) == 0.0, "hyper k>n must be 0");
    // The paper's own worked example: 100 peaks, 20 complemented, 3 of 4 by chance.
    const double p3of4 = fastag::hyperSf(100, 20, 4, 3);
    CHECK(p3of4 > 0.001 && p3of4 < 0.05, "paper example 3-of-4 = %.4f, expected ~0.023", p3of4);
    std::printf("2. hypergeometric tail: paper's 3-of-4 example = %.4f (paper says 2.3%%)\n", p3of4);
    // monotone decreasing in k
    double prev = 1.0;
    for (int k = 0; k <= 4; ++k)
    { const double v = fastag::hyperSf(100, 20, 4, k); CHECK(v <= prev + 1e-15, "hyper monotone k=%d", k); prev = v; }
  }
  {
    const auto v = fastag::expandIL("PLK");
    CHECK(v.size() == 2, "expandIL PLK -> %zu", v.size());
    CHECK(fastag::expandIL("LLL").size() == 8, "expandIL LLL -> 8");
    CHECK(fastag::expandIL("PEK").size() == 1, "expandIL no L -> 1");
    std::printf("3. I/L expansion at output only (LLL -> 8 variants)\n");
  }

  // ------------------------------------------------- end-to-end: clean spectra
  const std::vector<std::string> peps = {
    "SAMPLER", "PEPTIDEK", "VGAHAGEYGAEALER", "DFTPAELR", "YLGYLEQLLR"};

  {
    int total = 0, with_true_top = 0;
    for (const auto& pep : peps)
    {
      std::vector<double> mz, in;
      synth(pep, mz, in, 0, 42);
      const double M = peptideMass(pep);
      fastag::Param p;
      p.frag_tol = 0.02; p.complement_tol = 0.02;   // clean synthetic data
      fastag::Tables tab(p);
      const auto tags = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, p, tab);
      CHECK(!tags.empty(), "no tags for %s", pep.c_str());
      if (tags.empty()) continue;
      ++total;
      bool ok = false;
      for (const auto& v : fastag::expandIL(tags[0].seq)) if (isTrueTag(pep, v)) ok = true;
      if (ok) ++with_true_top;
      else std::printf("  note: top tag for %s was '%s' (E=%.3g)\n",
                       pep.c_str(), tags[0].seq.c_str(), tags[0].evalue);
    }
    CHECK(with_true_top == total, "top tag correct for %d/%d peptides", with_true_top, total);
    std::printf("4. top-ranked tag is a true subsequence for %d/%d clean spectra\n",
                with_true_top, total);
  }

  // --------------------------------------------- end-to-end: noisy, precision
  {
    int true_tags = 0, all_tags = 0;
    for (const auto& pep : peps)
    {
      std::vector<double> mz, in;
      synth(pep, mz, in, 60, 7);          // 60 noise peaks
      const double M = peptideMass(pep);
      fastag::Param p;
      p.frag_tol = 0.02; p.complement_tol = 0.02; p.max_tag_count = 10;
      fastag::Tables tab(p);
      const auto tags = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, p, tab);
      for (const auto& t : tags)
      {
        ++all_tags;
        for (const auto& v : fastag::expandIL(t.seq)) if (isTrueTag(pep, v)) { ++true_tags; break; }
      }
    }
    const double prec = all_tags ? static_cast<double>(true_tags) / all_tags : 0.0;
    std::printf("5. with 60 noise peaks: %d/%d top-10 tags correct (precision %.0f%%)\n",
                true_tags, all_tags, 100 * prec);
    CHECK(prec > 0.5, "precision %.2f too low on noisy spectra", prec);
  }

  // ---------------------------------------------------- flanking-mass accuracy
  {
    // The identity nterm + tag + cterm == precursor residue mass holds BY
    // CONSTRUCTION for every tag, so it validates nothing. Real placement means
    // the N-terminal flank equals an actual peptide prefix mass -- which is what
    // downstream database search matches on.
    const std::string pep = "VGAHAGEYGAEALER";
    std::vector<double> mz, in;
    synth(pep, mz, in, 0, 1);
    const double M = peptideMass(pep);
    std::vector<double> prefix{0.0};
    { double m = 0; for (char c : pep) { m += residueMass(c); prefix.push_back(m); } }

    fastag::Param p;
    p.frag_tol = 0.02; p.complement_tol = 0.02;
    fastag::Tables tab(p);
    const auto tags = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, p, tab);
    CHECK(!tags.empty(), "no tags for flank test");

    int identity_ok = 0, placed = 0; double worst_placed = 1e9;
    for (const auto& t : tags)
    {
      double tagmass = 0;
      for (char c : t.seq) tagmass += residueMass(c);
      if (std::fabs(t.nterm_mass + tagmass + t.cterm_mass - (M - fastag::WATER)) < 0.03) ++identity_ok;
      // y-derived: nterm flank lands on a real prefix, and nterm+tag on another
      for (double a : prefix)
      {
        if (std::fabs(a - t.nterm_mass) > 0.03) continue;
        for (double b : prefix)
          if (std::fabs(b - (a + tagmass)) < 0.03)
          { ++placed; worst_placed = std::min(worst_placed, std::fabs(a - t.nterm_mass)); goto next; }
      }
      next:;
    }
    CHECK(identity_ok == static_cast<int>(tags.size()),
          "conservation identity should hold for ALL tags, got %d/%zu", identity_ok, tags.size());
    CHECK(placed > 0, "no tag placed on a real prefix mass");
    std::printf("6. flanking masses: conservation identity holds for %d/%zu (by construction);\n"
                "   %d/%zu placed on a real peptide prefix (y-derived, exact to %.4f Da).\n"
                "   The rest are b-derived: flanks carry a one-water offset -- see fastag.h.\n",
                identity_ok, tags.size(), placed, tags.size(),
                worst_placed > 1e8 ? 0.0 : worst_placed);
  }

  // --------------------------------------------------------------- extension
  {
    const std::string pep = "VGAHAGEYGAEALER";
    std::vector<double> mz, in;
    synth(pep, mz, in, 0, 1);
    const double M = peptideMass(pep);

    fastag::Param base;
    base.frag_tol = 0.02; base.complement_tol = 0.02; base.tag_length = 3;
    fastag::Tables t1(base);
    const auto seeds = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, base, t1);

    fastag::Param ext = base;
    ext.max_extension = 4;
    fastag::Tables t2(ext);
    const auto grown = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, ext, t2);

    CHECK(!seeds.empty() && !grown.empty(), "extension test produced no tags");
    size_t max_seed = 0, max_grown = 0;
    for (const auto& t : seeds) max_seed = std::max(max_seed, t.seq.size());
    for (const auto& t : grown) max_grown = std::max(max_grown, t.seq.size());
    CHECK(max_seed == 3, "seed length should be exactly 3, got %zu", max_seed);
    CHECK(max_grown > max_seed, "extension produced no longer tag (%zu vs %zu)", max_grown, max_seed);

    // Extended tags must still be real subsequences on a clean spectrum.
    // Note a "wrong" tag is often an isobaric alternative reading, not an error:
    // GA and Q are both 128.0586, GG and N are both 114.0429, so the graph
    // legitimately contains a Q edge spanning the same peaks as a G->A path.
    int tr = 0, tot = 0;
    for (const auto& t : grown)
    {
      if (!t.extended) continue;
      ++tot;
      for (const auto& v : fastag::expandIL(t.seq)) if (isTrueTag(pep, v)) { ++tr; break; }
    }
    std::printf("7. extension: seed len %zu -> max len %zu; %d/%d extended tags are exact\n"
                "   subsequences (most of the remainder are isobaric readings, e.g. Q for GA)\n",
                max_seed, max_grown, tr, tot);
    CHECK(tot == 0 || tr > 0, "no extended tag was correct");

    // Isobaric ambiguity is real and must be demonstrated, not assumed.
    CHECK(std::fabs((residueMass('G') + residueMass('A')) - residueMass('Q')) < 1e-4,
          "GA vs Q should be isobaric");
    CHECK(std::fabs((residueMass('G') + residueMass('G')) - residueMass('N')) < 1e-4,
          "GG vs N should be isobaric");

    // Determinism: same input, same output, twice.
    fastag::Tables t3(ext);
    const auto again = fastag::tagSpectrum(mz, in, (M + 2 * fastag::PROTON) / 2, 2, ext, t3);
    bool same = again.size() == grown.size();
    for (size_t i = 0; same && i < again.size(); ++i)
      same = (again[i].seq == grown[i].seq && again[i].evalue == grown[i].evalue);
    CHECK(same, "output is not deterministic across runs");
    std::printf("8. output is deterministic across runs\n");
  }

  // ------------------------------------------------------------- edge cases
  {
    fastag::Param p;
    fastag::Tables tab(p);
    std::vector<double> e;
    CHECK(fastag::tagSpectrum(e, e, 500, 2, p, tab).empty(), "empty spectrum");
    std::vector<double> tiny_mz{100, 200}, tiny_in{1, 1};
    CHECK(fastag::tagSpectrum(tiny_mz, tiny_in, 500, 2, p, tab).empty(), "too few peaks");
    // NaN / non-positive intensities must be dropped, not crash
    std::vector<double> bad_mz{100, std::nan(""), 200, 300}, bad_in{1, 1, -1, 0};
    fastag::tagSpectrum(bad_mz, bad_in, 500, 2, p, tab);
    // unknown precursor charge
    std::vector<double> mz, in;
    synth("SAMPLER", mz, in, 0, 3);
    fastag::tagSpectrum(mz, in, 400, 0, p, tab);
    fastag::tagSpectrum(mz, in, 400, 1, p, tab);
    std::printf("9. edge cases survive: empty, too-few-peaks, NaN/zero intensity, charge 0/1\n");
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
