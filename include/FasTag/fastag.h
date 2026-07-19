// FasTag -- DirecTag-style peptide sequence tagging.
//
// Standalone and dependency-free so the algorithm can be tested without an
// OpenMS build; the TOPP tool is a thin adapter over this.
//
// Reimplemented from the published algorithm (Tabb et al., J Proteome Res 2008,
// 7:3838) and from reading the reference implementation for semantics -- not
// transcribed. freicore is Apache-2.0, OpenMS is BSD-3-Clause.
//
// Deviations from DirecTag, all deliberate and reviewed:
//  - rank-sum null by DP, not C(n,k) enumeration            (see ranksum.h)
//  - complement null computed exactly, not via an MVH table
//  - Fisher DOF = 2 * (number of active subscores)          (compat mode available)
//  - optional tag extension: seed length becomes a minimum
#pragma once

#include <FasTag/ranksum.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <atomic>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace fastag
{
  constexpr double PROTON = 1.00727646677;
  constexpr double WATER = 18.0105646863;

  /// The 19 distinct residue masses. I and L are isobaric: we traverse on the
  /// canonical 'L' and expand to 'I' at output, as DirecTag's TagExploder does.
  /// Expanding during traversal (as OpenMS::Tagger does) costs 2^k tags.
  struct Residue { char code; double mass; };
  inline const std::vector<Residue>& residues()
  {
    static const std::vector<Residue> r = {
      {'G',  57.02146}, {'A',  71.03711}, {'S',  87.03203}, {'P',  97.05276},
      {'V',  99.06841}, {'T', 101.04768}, {'C', 103.00919}, {'L', 113.08406},
      {'N', 114.04293}, {'D', 115.02694}, {'Q', 128.05858}, {'K', 128.09496},
      {'E', 129.04259}, {'M', 131.04049}, {'H', 137.05891}, {'F', 147.06841},
      {'R', 156.10111}, {'Y', 163.06333}, {'W', 186.07931}};
    return r;
  }

  struct Param
  {
    int    tag_length      = 3;      ///< seed length in residues
    double frag_tol        = 0.5;    ///< fragment m/z tolerance
    bool   tol_ppm         = false;
    double complement_tol  = 0.5;
    double precursor_tol   = 1.5;
    size_t max_peak_count  = 100;
    double tic_cutoff      = 1.0;    ///< keep peaks up to this fraction of TIC
    int    max_tag_count   = 50;     ///< output cap, applied after scoring
    double max_tag_score   = 20.0;   ///< E-value cutoff; <=0 disables

    /// Extension: seed length becomes a minimum. Default OFF -- at a 0.5 Da
    /// fragment tolerance ~67% of peaks have a spurious residue-gap neighbour,
    /// so a greedy extension adds ~2 wrong residues per seed, and a wrong
    /// extension destroys an identification the seed would have made.
    /// Only enable at high resolution.
    int    max_extension   = 0;      ///< residues added per terminus; 0 = off

    /// Reproduce DirecTag's rank-sum off-by-one (observed sum initialised to 1
    /// while the null is enumerated from 0).
    bool   compat_ranksum  = false;
    /// Reproduce DirecTag's DOF = floor(sum of weights)*2 instead of 2k.
    bool   compat_fisher_dof = false;
    double w_intensity = 1.0, w_mzfidelity = 1.0, w_complement = 1.0;

    uint32_t seed = 20080717u;       ///< the m/z-fidelity null is Monte Carlo
    int    mzfidelity_samples = 10000;
  };

  struct Tag
  {
    std::string seq;          ///< N->C orientation
    double nterm_mass = 0;    ///< residues missing from the N-terminal side
    double cterm_mass = 0;    ///< residues missing from the C-terminal side
    int    charge = 1;        ///< fragment charge this tag was read at
    double p_intensity = 1, p_mzfidelity = 1, p_complement = 1;
    double chisq = 0;
    double pvalue = 1;
    double evalue = 1;        ///< pvalue * seeds enumerated
    int    length = 0;        ///< realised length in residues (seed + extension)
    bool   extended = false;
    int    ext_nterm = 0;     ///< residues added walking toward higher m/z
    int    ext_cterm = 0;     ///< residues added walking toward lower m/z
    int    ext_total = 0;
    double low_mz = 0, high_mz = 0;
  };

  // ---------------------------------------------------------------- internals

  /// Survival function of chi-square with even DOF = 2m, in closed form:
  ///   P(X > x) = exp(-x/2) * sum_{i<m} (x/2)^i / i!
  /// Fisher's method always yields even DOF, so no incomplete-gamma is needed.
  inline double chi2SfEven(double x, int m)
  {
    if (m <= 0) return 1.0;
    if (x <= 0.0) return 1.0;
    const double h = x / 2.0;
    // Sum in log space to survive large x without underflowing to 0 prematurely.
    double term = 1.0, sum = 1.0;
    for (int i = 1; i < m; ++i) { term *= h / i; sum += term; }
    const double lg = -h + std::log(sum);
    return lg < -700.0 ? 0.0 : std::exp(lg);
  }

  inline double lnChoose(int n, int k)
  {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
  }

  /// P(X >= k) for X ~ Hypergeometric(N population, K successes, n draws).
  /// DirecTag builds an MVH table for this; the exact tail is cheaper and clearer.
  inline double hyperSf(int N, int K, int n, int k)
  {
    if (k <= 0) return 1.0;
    if (K <= 0 || n <= 0) return (k <= 0) ? 1.0 : 0.0;
    const double ln_den = lnChoose(N, n);
    double p = 0.0;
    const int hi = std::min(K, n);
    for (int i = k; i <= hi; ++i)
    {
      const double ln_num = lnChoose(K, i) + lnChoose(N - K, n - i);
      if (std::isfinite(ln_num)) p += std::exp(ln_num - ln_den);
    }
    return std::min(1.0, std::max(0.0, p));
  }

  /// Empirical null for the m/z-fidelity SSE, built by simulation exactly as
  /// DirecTag does. Deterministically seeded -- the reference is not, so two of
  /// its runs disagree with each other.
  class MzFidelityNull
  {
  public:
    MzFidelityNull() = default;
    MzFidelityNull(int npeaks, double tol, uint32_t seed, int samples)
    {
      std::mt19937 rng(seed + static_cast<uint32_t>(npeaks) * 2654435761u);
      std::uniform_real_distribution<double> err(-tol, tol);
      sse_.reserve(samples);
      std::vector<double> e(npeaks);
      for (int s = 0; s < samples; ++s)
      {
        double mean = 0.0;
        for (int i = 0; i < npeaks; ++i) { e[i] = err(rng); mean += e[i]; }
        mean /= npeaks;
        double sse = 0.0;
        for (int i = 0; i < npeaks; ++i) { const double d = e[i] - mean; sse += d * d; }
        sse_.push_back(sse);
      }
      std::sort(sse_.begin(), sse_.end());
    }
    /// P(SSE <= observed) -- small is good, as for the other subscores.
    double pvalue(double sse) const
    {
      if (sse_.empty()) return 1.0;
      const size_t r = static_cast<size_t>(
          std::lower_bound(sse_.begin(), sse_.end(), sse) - sse_.begin());
      return std::max(1.0 / sse_.size(), static_cast<double>(r) / sse_.size());
    }
  private:
    std::vector<double> sse_;
  };

  /// Per-run caches, built EAGERLY in the constructor and immutable thereafter.
  ///
  /// Eager construction is what makes multi-threading safe: a lazily-filled cache
  /// would be written from every worker thread. Everything a run can need is
  /// known from Param up front -- tag_peaks spans [tag_length+1,
  /// tag_length+1+2*max_extension] because extension changes the realised length,
  /// and peak counts span [0, max_peak_count]. Both lookups are then const and
  /// lock-free.
  ///
  /// Cost at defaults (seed 3, extension 6, 100 peaks): a few MB, tens of ms.
  /// It grows as O(k_max^2 * n_max^2), so a large max_peak_count is expensive --
  /// see the scaling note in ranksum.h.
  class Tables
  {
  public:
    explicit Tables(const Param& p)
      : k_min_(std::max(2, p.tag_length + 1)),
        k_max_(std::max(2, p.tag_length + 1 + 2 * std::max(0, p.max_extension))),
        n_max_(p.max_peak_count > 0 ? p.max_peak_count : 1024)
    {
      RanksumTables rs(p.compat_ranksum);
      cdf_.resize(static_cast<size_t>(k_max_ - k_min_ + 1));
      for (int k = k_min_; k <= k_max_; ++k)
      {
        auto& per_n = cdf_[static_cast<size_t>(k - k_min_)];
        per_n.resize(n_max_ + 1);
        for (size_t n = static_cast<size_t>(k); n <= n_max_; ++n) per_n[n] = rs.get(k, static_cast<int>(n));
      }
      mz_.reserve(static_cast<size_t>(k_max_ - k_min_ + 1));
      for (int k = k_min_; k <= k_max_; ++k)
        mz_.emplace_back(k, p.frag_tol, p.seed, p.mzfidelity_samples);
    }

    /// p-value for an observed rank sum. Conservative (1.0) outside the built range.
    double intensityP(int tag_peaks, int npeaks, int ranksum) const
    {
      if (tag_peaks < k_min_ || tag_peaks > k_max_) return 1.0;
      if (npeaks < 0 || static_cast<size_t>(npeaks) > n_max_) return 1.0;
      const std::vector<float>& c = cdf_[static_cast<size_t>(tag_peaks - k_min_)][npeaks];
      if (c.empty()) return 1.0;
      if (ranksum < 0) return 0.0;
      if (static_cast<size_t>(ranksum) >= c.size()) return 1.0;
      return c[ranksum];
    }

    const MzFidelityNull& mzNull(int tag_peaks) const
    {
      const int k = std::min(std::max(tag_peaks, k_min_), k_max_);
      return mz_[static_cast<size_t>(k - k_min_)];
    }

  private:
    int k_min_, k_max_;
    size_t n_max_;
    std::vector<std::vector<std::vector<float>>> cdf_;   ///< [k][npeaks] -> CDF
    std::vector<MzFidelityNull> mz_;                     ///< [k]
  };

  /// Preprocessed spectrum: peaks sorted by m/z, with intensity ranks and
  /// per-charge complement flags.
  struct Spectrum
  {
    std::vector<double> mz;
    std::vector<double> intensity;
    std::vector<int>    rank;        ///< 1 = most intense
    std::vector<uint8_t> has_compl;  ///< bit z-1 set if a complement exists at charge z
    double precursor_mz = 0;
    int    precursor_charge = 0;
    double precursor_mass = 0;       ///< neutral
    int    n_frag_charges = 1;
    int    n_with_compl = 0;         ///< peaks with a complement at any charge
  };

  inline double tolAt(const Param& p, double mz)
  {
    return p.tol_ppm ? mz * p.frag_tol * 1e-6 : p.frag_tol;
  }

  /// Index of the peak nearest to @p target within @p tol, or -1.
  inline int findNear(const std::vector<double>& mz, double target, double tol)
  {
    if (mz.empty()) return -1;
    auto it = std::lower_bound(mz.begin(), mz.end(), target);
    int best = -1; double bestd = tol;
    for (int k = -1; k <= 0; ++k)
    {
      auto j = it + k;
      if (j < mz.begin() || j >= mz.end()) continue;
      const double d = std::fabs(*j - target);
      if (d <= bestd) { bestd = d; best = static_cast<int>(j - mz.begin()); }
    }
    return best;
  }

  /// Filter, rank, and find complements. Mirrors DirecTag's FilterPeaks +
  /// Preprocess, minus the precursor-mass grid search (off by default there too,
  /// and ~51 complement searches per spectrum).
  inline Spectrum preprocess(const std::vector<double>& raw_mz,
                             const std::vector<double>& raw_int,
                             double precursor_mz, int precursor_charge,
                             const Param& p)
  {
    Spectrum s;
    s.precursor_mz = precursor_mz;
    s.precursor_charge = precursor_charge;
    s.precursor_mass = precursor_mz * precursor_charge - precursor_charge * PROTON;
    s.n_frag_charges = std::max(1, precursor_charge - 1);

    // Order peaks by index so we can filter on intensity then restore m/z order.
    std::vector<size_t> idx;
    const double max_mass = s.precursor_mass + PROTON + p.precursor_tol;
    for (size_t i = 0; i < raw_mz.size(); ++i)
    {
      if (!std::isfinite(raw_mz[i]) || !std::isfinite(raw_int[i])) continue;
      if (raw_int[i] <= 0.0) continue;
      if (raw_mz[i] > max_mass) continue;   // nothing above the precursor
      idx.push_back(i);
    }
    // Tie-break order is load-bearing and must match DirecTag exactly.
    // freicore/PeakSpectrum.h:507 FilterByPeakCount builds a multimap keyed by
    // intensity, inserting in ascending m/z (peakPreData is an m/z-keyed map, and
    // multimap::insert appends at the upper bound of an equal range), then walks
    // it with a reverse_iterator. The effective key is therefore
    // (intensity DESC, m/z DESC). A plain unstable sort on intensity alone picks
    // a different peak set whenever intensities tie -- and with integer detector
    // counts they tie often, which is where real-data parity was being lost.
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      if (raw_int[a] != raw_int[b]) return raw_int[a] > raw_int[b];
      return raw_mz[a] > raw_mz[b];
    });

    // TIC cutoff (freicore/PeakSpectrum.h:463). DirecTag accumulates relative
    // intensity descending until the cutoff is reached, then keeps every peak at
    // or above that threshold intensity -- including ties past the count, via
    // lower_bound. A no-op at the default cutoff of 1.0.
    if (p.tic_cutoff < 1.0)
    {
      double total = 0.0;
      for (size_t i : idx) total += raw_int[i];
      if (total > 0.0)
      {
        double run = 0.0; size_t keep = 0;
        for (; keep < idx.size() && run < p.tic_cutoff * total; ++keep) run += raw_int[idx[keep]];
        if (keep > 0 && keep < idx.size())
        {
          const double thresh = raw_int[idx[keep - 1]];
          while (keep < idx.size() && raw_int[idx[keep]] == thresh) ++keep;  // keep ties
        }
        idx.resize(keep);
      }
    }
    if (p.max_peak_count > 0 && idx.size() > p.max_peak_count) idx.resize(p.max_peak_count);

    // rank by intensity (1-based) before restoring m/z order
    std::vector<int> rank_of(idx.size());
    for (size_t r = 0; r < idx.size(); ++r) rank_of[r] = static_cast<int>(r) + 1;
    std::vector<std::pair<size_t, int>> keep;
    keep.reserve(idx.size());
    for (size_t r = 0; r < idx.size(); ++r) keep.emplace_back(idx[r], rank_of[r]);
    std::sort(keep.begin(), keep.end(),
              [&](const auto& a, const auto& b) { return raw_mz[a.first] < raw_mz[b.first]; });

    s.mz.reserve(keep.size()); s.intensity.reserve(keep.size()); s.rank.reserve(keep.size());
    for (const auto& kv : keep)
    {
      s.mz.push_back(raw_mz[kv.first]);
      s.intensity.push_back(raw_int[kv.first]);
      s.rank.push_back(kv.second);
    }
    s.has_compl.assign(s.mz.size(), 0);

    // Complements: b + y = M + 2*proton, so complement m/z at charge z is
    // M/z + 2*proton - mz.
    for (size_t i = 0; i < s.mz.size(); ++i)
    {
      for (int z = 1; z <= s.n_frag_charges; ++z)
      {
        const double cmz = s.precursor_mass / z + 2.0 * PROTON - s.mz[i];
        if (findNear(s.mz, cmz, p.complement_tol) >= 0)
          s.has_compl[i] |= static_cast<uint8_t>(1u << (z - 1));
      }
    }
    for (size_t i = 0; i < s.mz.size(); ++i) if (s.has_compl[i]) ++s.n_with_compl;
    return s;
  }

  /// Gap graph in CSR form, per fragment charge. Edges are peak indices, so the
  /// whole structure is a few tens of KB and never invalidated during traversal.
  struct Graph
  {
    // forward: edges out of peak i, toward higher m/z
    std::vector<uint32_t> off;   ///< size npeaks+1
    std::vector<uint32_t> dst;
    std::vector<uint8_t>  res;   ///< index into residues()
    // reverse: edges into peak i, so C-terminal extension is also O(degree).
    // Without this, walking backwards costs a full O(n * E) scan per step.
    std::vector<uint32_t> roff;
    std::vector<uint32_t> rsrc;
    std::vector<uint8_t>  rres;
  };

  inline Graph buildGraph(const Spectrum& s, const Param& p, int charge)
  {
    const size_t n = s.mz.size();
    Graph g;
    g.off.assign(n + 1, 0);
    g.roff.assign(n + 1, 0);
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> adj(n), radj(n);
    const auto& R = residues();
    for (size_t i = 0; i < n; ++i)
    {
      for (uint8_t r = 0; r < R.size(); ++r)
      {
        const double target = s.mz[i] + R[r].mass / charge;
        const int j = findNear(s.mz, target, tolAt(p, target));
        if (j > static_cast<int>(i))
        {
          adj[i].emplace_back(static_cast<uint32_t>(j), r);
          radj[j].emplace_back(static_cast<uint32_t>(i), r);
        }
      }
    }
    uint32_t total = 0, rtotal = 0;
    for (size_t i = 0; i < n; ++i)
    {
      g.off[i] = total;   total  += adj[i].size();
      g.roff[i] = rtotal; rtotal += radj[i].size();
    }
    g.off[n] = total; g.roff[n] = rtotal;
    g.dst.reserve(total); g.res.reserve(total);
    g.rsrc.reserve(rtotal); g.rres.reserve(rtotal);
    for (size_t i = 0; i < n; ++i)
      for (const auto& e : adj[i]) { g.dst.push_back(e.first); g.res.push_back(e.second); }
    for (size_t i = 0; i < n; ++i)
      for (const auto& e : radj[i]) { g.rsrc.push_back(e.first); g.rres.push_back(e.second); }
    return g;
  }

  /// Score one path (peak indices + residue indices) and fill a Tag.
  inline Tag scorePath(const Spectrum& s, const Param& p, const Tables& tab,
                       const std::vector<uint32_t>& peaks,
                       const std::vector<uint8_t>& res, int charge)
  {
    const auto& R = residues();
    const int tag_peaks = static_cast<int>(peaks.size());
    Tag t;
    t.charge = charge;
    t.length = tag_peaks - 1;

    // Back-project every peak to an estimate of the first peak's position, then
    // take the mean. SSE about that mean is the m/z-fidelity statistic.
    std::vector<double> est(tag_peaks);
    double cum = 0.0;
    est[0] = s.mz[peaks[0]];
    for (int i = 0; i < tag_peaks - 1; ++i)
    {
      cum += R[res[i]].mass / charge;
      est[i + 1] = s.mz[peaks[i + 1]] - cum;
    }
    const double avg = std::accumulate(est.begin(), est.end(), 0.0) / tag_peaks;
    double sse = 0.0;
    for (double e : est) { const double d = e - avg; sse += d * d; }
    t.p_mzfidelity = tab.mzNull(tag_peaks).pvalue(sse);

    // Intensity rank sum. DirecTag seeds this at 1 while enumerating its null
    // from 0, shifting every lookup by a bin; reproduced only in compat mode.
    int ranksum = p.compat_ranksum ? 1 : 0;
    int n_compl = 0;
    for (int i = 0; i < tag_peaks; ++i)
    {
      ranksum += s.rank[peaks[i]];
      if (s.has_compl[peaks[i]] & (1u << (charge - 1))) ++n_compl;
    }
    t.p_intensity = tab.intensityP(tag_peaks, static_cast<int>(s.mz.size()), ranksum);

    // Complementarity: probability of seeing at least this many complemented
    // peaks among tag_peaks drawn from the spectrum. Dropped entirely when the
    // spectrum has no complement pairs at all, as in DirecTag.
    const bool use_compl = s.n_with_compl > 0;
    t.p_complement = use_compl
        ? hyperSf(static_cast<int>(s.mz.size()), s.n_with_compl, tag_peaks, n_compl)
        : 1.0;

    // Fisher. Clamp before log: an underflowed subscore would give -inf.
    auto safe = [](double v) { return std::log(std::max(v, 1e-300)); };
    const double wi = p.w_intensity, wm = p.w_mzfidelity;
    const double wc = use_compl ? p.w_complement : 0.0;
    t.chisq = -2.0 * (safe(t.p_intensity) * wi + safe(t.p_mzfidelity) * wm
                      + (use_compl ? safe(t.p_complement) * wc : 0.0));
    const int k_active = 2 + (use_compl ? 1 : 0);
    const int m = p.compat_fisher_dof
        ? std::max(1, static_cast<int>(std::floor(wi + wm + wc)))  // DirecTag's DOF bug
        : k_active;
    t.pvalue = chi2SfEven(t.chisq, m);

    // Flanking masses from the FITTED positions, not the raw ones. These feed
    // downstream database search as precursor-mass constraints and are the only
    // consumer of the fit.
    //
    // ION-TYPE CONVENTION (inherited from DirecTag, and load-bearing downstream):
    // the low-m/z peak is assumed to be a y ion, so cterm_mass is the residue
    // mass below it and nterm_mass is the residue mass above the high peak.
    // Peaks are not labelled by ion type, so a tag read off the b series gets
    // flanks that are offset by exactly one water: its cterm_mass equals
    // (N-terminal prefix mass - WATER). Measured on synthetic b/y spectra,
    // ~43% of tags are y-derived with exact flanks and the rest are b-derived
    // with that offset.
    //
    // Note the identity nterm + tag + cterm == precursor_mass - WATER holds for
    // BOTH cases by construction, so it is *not* a validity check -- placement
    // must be tested against prefix masses instead.
    // Downstream search must try both interpretations; that is standard for
    // tag-based search and is why TagRecon filters on flank consistency.
    const double fitted_first = avg;
    const double fitted_last  = avg + cum;
    t.cterm_mass = std::max(0.0, fitted_first * charge - PROTON * charge - WATER);
    t.nterm_mass = std::max(0.0, s.precursor_mass - (fitted_last * charge - PROTON * charge));

    // Traversal runs low->high m/z, which is C->N; store N->C.
    std::string seq;
    for (int i = tag_peaks - 2; i >= 0; --i) seq.push_back(R[res[i]].code);
    t.seq = seq;
    t.low_mz = s.mz[peaks.front()];
    t.high_mz = s.mz[peaks.back()];
    return t;
  }

  /// Extend a path as far as the graph allows, greedily by smallest |m/z error|.
  /// Ties are broken by intensity rank then peak index so output is deterministic.
  /// Returns the number of residues appended.
  /// Greedily extend a path by one terminus, as far as the graph and the budget
  /// allow. Returns the number of residues appended.
  ///
  /// Selection is by smallest |m/z error|, tie-broken by better (lower) intensity
  /// rank, then by peak index. All three are needed: ties at equal error are
  /// common, and without a total order the output is nondeterministic.
  ///
  /// Extension stays within the seed's fragment charge -- gaps are charge-scaled,
  /// and DirecTag forbids mixed-charge paths by design.
  inline int extendPath(const Spectrum& s, const Param& p, const Graph& g,
                        std::vector<uint32_t>& peaks, std::vector<uint8_t>& res,
                        int charge, bool forward)
  {
    const auto& R = residues();
    int added = 0;
    std::vector<bool> used(s.mz.size(), false);
    for (uint32_t i : peaks) used[i] = true;   // never revisit a peak

    while (added < p.max_extension)
    {
      const uint32_t cur = forward ? peaks.back() : peaks.front();
      const uint32_t e0 = forward ? g.off[cur]  : g.roff[cur];
      const uint32_t e1 = forward ? g.off[cur + 1] : g.roff[cur + 1];

      int best = -1; uint8_t best_r = 0; double best_err = 0;
      for (uint32_t e = e0; e < e1; ++e)
      {
        const uint32_t cand = forward ? g.dst[e] : g.rsrc[e];
        const uint8_t  rr   = forward ? g.res[e] : g.rres[e];
        if (used[cand]) continue;
        // Expected position of the *later* peak in the pair, either way round.
        const double err = forward
            ? std::fabs(s.mz[cand] - (s.mz[cur] + R[rr].mass / charge))
            : std::fabs(s.mz[cur]  - (s.mz[cand] + R[rr].mass / charge));
        const bool better =
            best < 0 ||
            err < best_err - 1e-12 ||
            (std::fabs(err - best_err) <= 1e-12 &&
             (s.rank[cand] < s.rank[best] ||
              (s.rank[cand] == s.rank[best] && cand < static_cast<uint32_t>(best))));
        if (better) { best = static_cast<int>(cand); best_r = rr; best_err = err; }
      }
      if (best < 0) break;

      used[best] = true;
      if (forward) { peaks.push_back(static_cast<uint32_t>(best)); res.push_back(best_r); }
      else { peaks.insert(peaks.begin(), static_cast<uint32_t>(best));
             res.insert(res.begin(), best_r); }
      ++added;
    }
    return added;
  }

  /// Tag one spectrum. Returns tags sorted by E-value, best first.
  inline std::vector<Tag> tagSpectrum(const std::vector<double>& raw_mz,
                                      const std::vector<double>& raw_int,
                                      double precursor_mz, int precursor_charge,
                                      const Param& p, const Tables& tab)
  {
    std::vector<Tag> out;
    if (p.tag_length < 1) return out;
    const Spectrum s = preprocess(raw_mz, raw_int, precursor_mz, precursor_charge, p);
    if (s.mz.size() < static_cast<size_t>(p.tag_length) + 1) return out;

    std::vector<Tag> scored;
    size_t n_seeds = 0;

    for (int charge = 1; charge <= s.n_frag_charges; ++charge)
    {
      const Graph g = buildGraph(s, p, charge);
      std::vector<uint32_t> peaks;
      std::vector<uint8_t>  res;

      // Explicit-stack DFS to exactly tag_length edges.
      struct Frame { uint32_t node; uint32_t edge; };
      for (uint32_t start = 0; start < s.mz.size(); ++start)
      {
        std::vector<Frame> st;
        peaks.clear(); res.clear();
        peaks.push_back(start);
        st.push_back({start, g.off[start]});
        while (!st.empty())
        {
          Frame& f = st.back();
          if (static_cast<int>(res.size()) == p.tag_length)
          {
            ++n_seeds;
            std::vector<uint32_t> pk = peaks;
            std::vector<uint8_t>  rs = res;
            int a = 0, b = 0;
            if (p.max_extension > 0)
            {
              // Budgets are independent per terminus, not shared.
              a = extendPath(s, p, g, pk, rs, charge, true);
              b = extendPath(s, p, g, pk, rs, charge, false);
            }
            // Rescore at the realised length: the seed's statistic does not
            // describe an extended tag, and the per-length nulls are cheap.
            Tag t = scorePath(s, p, tab, pk, rs, charge);
            // Traversal runs low->high m/z = C->N, so a forward (higher-m/z)
            // extension grows the N-terminal end of the stored N->C sequence.
            t.ext_nterm = a;
            t.ext_cterm = b;
            t.ext_total = a + b;
            t.extended = (a + b) > 0;
            scored.push_back(std::move(t));
            res.pop_back(); peaks.pop_back(); st.pop_back();
            continue;
          }
          if (f.edge >= g.off[f.node + 1]) { st.pop_back(); if (!res.empty()) { res.pop_back(); peaks.pop_back(); } continue; }
          const uint32_t e = f.edge++;
          const uint32_t nxt = g.dst[e];
          peaks.push_back(nxt);
          res.push_back(g.res[e]);
          st.push_back({nxt, g.off[nxt]});
        }
      }
    }

    // E-value: expected number of seeds this good by chance. Valid under
    // arbitrary dependence between tags by linearity of expectation.
    for (Tag& t : scored) t.evalue = t.pvalue * static_cast<double>(std::max<size_t>(n_seeds, 1));

    std::sort(scored.begin(), scored.end(), [](const Tag& a, const Tag& b) {
      if (a.evalue != b.evalue) return a.evalue < b.evalue;
      if (a.seq != b.seq) return a.seq < b.seq;
      return a.low_mz < b.low_mz;
    });

    for (Tag& t : scored)
    {
      if (p.max_tag_score > 0 && t.evalue > p.max_tag_score) break;
      out.push_back(std::move(t));
      if (p.max_tag_count > 0 && static_cast<int>(out.size()) >= p.max_tag_count) break;
    }
    return out;
  }

  /// One spectrum's worth of input for the batch API.
  struct SpectrumInput
  {
    std::vector<double> mz, intensity;
    double precursor_mz = 0;
    int    charge = 0;
    std::string id;
  };

  /// Tag a batch of spectra across @p nthreads worker threads.
  ///
  /// Parallelism is over SPECTRA and nowhere else. Each spectrum is independent,
  /// `Tables` is immutable after construction, and every worker owns its own
  /// scratch, so there is no shared mutable state and no locking on the hot path.
  ///
  /// Deliberately not parallelised inside tagSpectrum: the per-spectrum work is
  /// ~60 us, far too small to amortise a fork/join, and nesting it under a caller
  /// that already parallelises over spectra is the oversubscription trap flagged
  /// in review.
  ///
  /// Output is written to a pre-sized vector indexed by input position, so the
  /// result is **bit-identical to the single-threaded path regardless of thread
  /// count** -- no append races, no ordering dependence.
  ///
  /// @param nthreads 0 selects hardware_concurrency().
  inline std::vector<std::vector<Tag>> tagSpectra(const std::vector<SpectrumInput>& in,
                                                  const Param& p, const Tables& tab,
                                                  unsigned nthreads = 0)
  {
    std::vector<std::vector<Tag>> out(in.size());
    if (in.empty()) return out;

    if (nthreads == 0) nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    nthreads = std::min<unsigned>(nthreads, static_cast<unsigned>(in.size()));

    if (nthreads == 1)
    {
      for (size_t i = 0; i < in.size(); ++i)
        out[i] = tagSpectrum(in[i].mz, in[i].intensity, in[i].precursor_mz,
                             in[i].charge ? in[i].charge : 2, p, tab);
      return out;
    }

    // Dynamic scheduling via one atomic counter. Spectra differ several-fold in
    // cost (peak count, charge, graph density), so a static split leaves threads
    // idle at the tail.
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; ++t)
    {
      workers.emplace_back([&]() {
        for (;;)
        {
          const size_t i = next.fetch_add(1, std::memory_order_relaxed);
          if (i >= in.size()) return;
          out[i] = tagSpectrum(in[i].mz, in[i].intensity, in[i].precursor_mz,
                               in[i].charge ? in[i].charge : 2, p, tab);
        }
      });
    }
    for (auto& w : workers) w.join();
    return out;
  }

  /// Expand 'L' to every I/L variant. DirecTag does this at output; doing it
  /// during traversal (as OpenMS::Tagger does) costs 2^k tags.
  inline std::vector<std::string> expandIL(const std::string& tag)
  {
    std::vector<std::string> out{tag};
    for (size_t i = 0; i < tag.size(); ++i)
    {
      if (tag[i] != 'L') continue;
      const size_t n = out.size();
      for (size_t j = 0; j < n; ++j) { std::string v = out[j]; v[i] = 'I'; out.push_back(v); }
    }
    return out;
  }
}
