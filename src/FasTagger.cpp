// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "FasTagger.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <tuple>
#include <random>

using namespace OpenMS;

namespace FasTag
{
  /// Peak ceiling used when max_peak_count is 0. Tables and prepare() must agree:
  /// intensityP() returns 1.0 above the table extent, silently disabling the
  /// intensity subscore rather than erroring.
  static constexpr size_t PEAK_CEILING = 1024;

  namespace
  {
    /// Residue masses from OpenMS, not a hardcoded table.
    ///
    /// "Natural19WithoutI" is exactly the canonical set the graph needs: I and L
    /// are isobaric, so a spectrum cannot distinguish them and enumerating both
    /// would only double the output. OpenMS::Tagger uses the same set.
    ///
    /// Measured against the table this replaced: worst difference 0.005 mDa
    /// (cysteine), i.e. 4000x smaller than the 20 mDa tolerance at m/z 1000.
    struct ResidueTable
    {
      std::vector<char>   code;
      std::vector<double> mass;

      ResidueTable()
      {
        for (const Residue* r : ResidueDB::getInstance()->getResidues("Natural19WithoutI"))
        {
          code.push_back(r->getOneLetterCode()[0]);
          mass.push_back(r->getMonoWeight(Residue::Internal));
        }
        // Ascending mass keeps the graph edges ordered, which makes the DFS emit
        // tags deterministically without a later sort.
        std::vector<size_t> idx(mass.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return mass[a] < mass[b]; });
        std::vector<char> c2; std::vector<double> m2;
        for (size_t i : idx) { c2.push_back(code[i]); m2.push_back(mass[i]); }
        code.swap(c2); mass.swap(m2);
      }
    };

    const ResidueTable& residues()
    {
      static const ResidueTable t;
      return t;
    }

    const double WATER = EmpiricalFormula("H2O").getMonoWeight();
    const double PROTON = Constants::PROTON_MASS_U;

    inline double tolAt(const Param& p, double mz)
    {
      return p.tol_ppm ? mz * p.frag_tol * 1e-6 : p.frag_tol;
    }

    /// Combine independent-ish p-values by Fisher's method.
    ///
    /// The three subscores are computed from the same peaks and are therefore not
    /// independent, so the chi-square tail is optimistic. That is inherited from
    /// DirecTag and is why the result is reported as an E-value for ranking rather
    /// than as a calibrated probability.
    double fisher(double p1, double p2, double p3, int k)
    {
      auto safe = [](double v) { return std::log(std::max(v, 1e-300)); };
      const double x = -2.0 * (safe(p1) + safe(p2) + (k == 3 ? safe(p3) : 0.0));
      if (x <= 0 || k <= 0) return 1.0;
      return boost::math::cdf(boost::math::complement(
          boost::math::chi_squared(2.0 * k), x));
    }
  }

  // ---------------------------------------------------------------- Tables

  Tables::Tables(const Param& p)
    : k_min_(std::max(2, p.tag_length + 1)),
      k_max_(std::max(2, p.tag_length + 1 + 2 * std::max(0, p.max_extension))),
      n_max_(p.max_peak_count > 0 ? p.max_peak_count : PEAK_CEILING)
  {
    ranksum_.resize(static_cast<size_t>(k_max_ - k_min_ + 1));
    for (int k = k_min_; k <= k_max_; ++k)
    {
      auto& per_n = ranksum_[static_cast<size_t>(k - k_min_)];
      per_n.resize(n_max_ + 1);
      for (size_t n = static_cast<size_t>(k); n <= n_max_; ++n)
        per_n[n] = ranksumCdf(k, static_cast<int>(n));
    }

    // The m/z-fidelity null has no closed form: it is the distribution of the
    // sum of squared deviations of k uniform errors about their own mean.
    // Sampled, with an explicit seed -- DirecTag leaves this unseeded, so two of
    // its runs disagree with each other.
    //
    // Sampled in UNITS OF THE TOLERANCE, on [-1, 1], and the observed SSE is
    // divided by tol^2 before lookup. Sampling on [-frag_tol, frag_tol] instead
    // is wrong whenever the tolerance is in ppm: frag_tol is then 20, so the null
    // spans +-20 Da while observed deviations are ~1e-4 Da. Every tag lands below
    // the whole null, mzFidelityP returns its floor, and one of the three Fisher
    // subscores becomes a constant -- contributing nothing. Working in tolerance
    // units also makes one null valid at every m/z.
    mz_null_.resize(static_cast<size_t>(k_max_ - k_min_ + 1));
    for (int k = k_min_; k <= k_max_; ++k)
    {
      std::mt19937 rng(p.seed + static_cast<unsigned>(k) * 2654435761u);
      std::uniform_real_distribution<double> err(-1.0, 1.0);
      auto& v = mz_null_[static_cast<size_t>(k - k_min_)];
      v.reserve(static_cast<size_t>(p.mzfidelity_samples));
      std::vector<double> e(static_cast<size_t>(k));
      for (int s = 0; s < p.mzfidelity_samples; ++s)
      {
        double mean = 0;
        for (int i = 0; i < k; ++i) { e[i] = err(rng); mean += e[i]; }
        mean /= k;
        double sse = 0;
        for (int i = 0; i < k; ++i) { const double d = e[i] - mean; sse += d * d; }
        v.push_back(sse);
      }
      std::sort(v.begin(), v.end());
    }
  }

  double Tables::intensityP(int tag_peaks, int npeaks, int ranksum) const
  {
    if (tag_peaks < k_min_ || tag_peaks > k_max_) return 1.0;
    if (npeaks < 0 || static_cast<size_t>(npeaks) > n_max_) return 1.0;
    const auto& c = ranksum_[static_cast<size_t>(tag_peaks - k_min_)][static_cast<size_t>(npeaks)];
    if (c.empty()) return 1.0;
    if (ranksum < 0) return 0.0;
    if (static_cast<size_t>(ranksum) >= c.size()) return 1.0;
    return c[static_cast<size_t>(ranksum)];
  }

  double Tables::mzFidelityP(int tag_peaks, double sse) const
  {
    const int k = std::min(std::max(tag_peaks, k_min_), k_max_);
    const auto& v = mz_null_[static_cast<size_t>(k - k_min_)];
    if (v.empty()) return 1.0;
    const size_t r = static_cast<size_t>(std::lower_bound(v.begin(), v.end(), sse) - v.begin());
    return std::max(1.0 / v.size(), static_cast<double>(r) / v.size());
  }

  // ------------------------------------------------------------- internals

  namespace
  {
    /// Preprocessed spectrum: peaks sorted by m/z with intensity ranks and
    /// per-charge complement flags.
    struct Prepared
    {
      MSSpectrum spec;                 ///< filtered, m/z sorted
      std::vector<int> rank;           ///< 1 = most intense
      std::vector<uint8_t> has_compl;  ///< bit z-1 set if a complement exists at charge z
      double precursor_mass = 0;
      int n_frag_charges = 1;
      int n_with_compl = 0;
    };

    Prepared prepare(const MSSpectrum& in, double precursor_mz, int charge, const Param& p)
    {
      Prepared s;
      s.precursor_mass = precursor_mz * charge - charge * PROTON;
      // +2 fragments are only sought from precursors of charge >= 3, as the paper
      // specifies.
      s.n_frag_charges = std::max(1, charge - 1);

      // Drop everything above the precursor, then keep the n most intense.
      //
      // NOT OpenMS::NLargest, deliberately. It does sortByIntensity() -- an
      // unstable sort with no secondary key -- and keeps the first n. Peaks
      // carry integer detector counts, so ties at the cut are common: measured
      // on real timsTOF data, 23% of spectra have at least one peak tied at the
      // 100-peak boundary. Which of them survives would then depend on the
      // standard library's sort, and peak selection provably drives tag output.
      // The total order below (intensity desc, m/z desc) makes it reproducible.
      MSSpectrum work;
      work.reserve(in.size());
      const double max_mz = s.precursor_mass + PROTON + p.precursor_tol;
      for (const auto& pk : in)
      {
        if (!std::isfinite(pk.getMZ()) || !std::isfinite(pk.getIntensity())) continue;
        if (pk.getIntensity() <= 0 || pk.getMZ() > max_mz) continue;
        work.push_back(pk);
      }

      // Collapse peaks sharing an m/z, keeping the strongest.
      //
      // diaTracer emits them routinely -- measured on S23, every MS2 spectrum has
      // them, median 22 and up to 46 per spectrum, some with equal and some with
      // differing intensity. They are indistinguishable to the graph (findNearest
      // returns one of them) but each still consumes a slot in the peak budget and
      // shifts every intensity rank, and they let the same tag be enumerated more
      // than once. Collapsing is a property of the input, so it belongs here
      // rather than in a deduplication pass over the output.
      if (!work.empty())
      {
        work.sortByPosition();
        std::vector<Peak1D> uniq;
        uniq.reserve(work.size());
        for (const auto& pk : work)
        {
          if (!uniq.empty() && pk.getMZ() == uniq.back().getMZ())
          {
            if (pk.getIntensity() > uniq.back().getIntensity()) uniq.back() = pk;
          }
          else uniq.push_back(pk);
        }
        work.clear(false);
        for (const auto& pk : uniq) work.push_back(pk);
      }
      std::vector<size_t> order(work.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (work[a].getIntensity() != work[b].getIntensity())
          return work[a].getIntensity() > work[b].getIntensity();
        return work[a].getMZ() > work[b].getMZ();
      });
      // max_peak_count == 0 means "use the safety ceiling", NOT "unlimited".
      //
      // The rank-sum tables are built for peak counts up to n_max_, and
      // intensityP() returns 1.0 above that -- silently disabling one of the three
      // subscores rather than erroring. So an genuinely uncapped run would quietly
      // lose a third of the scoring on any spectrum with more peaks than the
      // tables cover. Both places must agree on the same ceiling.
      const size_t cap = p.max_peak_count > 0 ? p.max_peak_count : PEAK_CEILING;
      if (order.size() > cap) order.resize(cap);

      // Keep the rank alongside the peak, then restore m/z order: the graph and
      // findNearest both need a sorted spectrum.
      std::vector<std::pair<Peak1D, int>> kept;
      kept.reserve(order.size());
      for (size_t r = 0; r < order.size(); ++r)
        kept.emplace_back(work[order[r]], static_cast<int>(r) + 1);
      std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) {
        return a.first.getMZ() < b.first.getMZ();
      });

      s.spec.reserve(kept.size());
      s.rank.reserve(kept.size());
      for (const auto& kv : kept) { s.spec.push_back(kv.first); s.rank.push_back(kv.second); }
      const size_t n = s.spec.size();

      // Complements: b + y = M + 2*proton, so a fragment at m/z m has its partner
      // at M/z + 2*proton - m.
      s.has_compl.assign(n, 0);
      for (size_t i = 0; i < n; ++i)
        for (int z = 1; z <= s.n_frag_charges; ++z)
        {
          const double cmz = s.precursor_mass / z + 2.0 * PROTON - s.spec[i].getMZ();
          if (s.spec.findNearest(cmz, p.complement_tol) >= 0)
            s.has_compl[i] |= static_cast<uint8_t>(1u << (z - 1));
        }
      for (size_t i = 0; i < n; ++i) if (s.has_compl[i]) ++s.n_with_compl;
      return s;
    }

    /// Residue-gap graph in compressed-sparse-row form, per fragment charge.
    struct Graph
    {
      std::vector<uint32_t> off, dst;   ///< forward edges, toward higher m/z
      std::vector<uint8_t>  res;
      std::vector<uint32_t> roff, rsrc; ///< reverse edges, so C-terminal
      std::vector<uint8_t>  rres;       ///< extension is O(degree) too
    };

    Graph buildGraph(const Prepared& s, const Param& p, int charge)
    {
      const size_t n = s.spec.size();
      const auto& R = residues();
      Graph g;
      g.off.assign(n + 1, 0);
      g.roff.assign(n + 1, 0);
      std::vector<std::vector<std::pair<uint32_t, uint8_t>>> adj(n), radj(n);
      for (size_t i = 0; i < n; ++i)
        for (uint8_t r = 0; r < R.mass.size(); ++r)
        {
          const double target = s.spec[i].getMZ() + R.mass[r] / charge;
          const int j = s.spec.findNearest(target, tolAt(p, target));
          if (j > static_cast<int>(i))
          {
            adj[i].emplace_back(static_cast<uint32_t>(j), r);
            radj[static_cast<size_t>(j)].emplace_back(static_cast<uint32_t>(i), r);
          }
        }
      uint32_t t = 0, rt = 0;
      for (size_t i = 0; i < n; ++i)
      { g.off[i] = t; t += adj[i].size(); g.roff[i] = rt; rt += radj[i].size(); }
      g.off[n] = t; g.roff[n] = rt;
      g.dst.reserve(t); g.res.reserve(t); g.rsrc.reserve(rt); g.rres.reserve(rt);
      for (size_t i = 0; i < n; ++i)
        for (const auto& e : adj[i]) { g.dst.push_back(e.first); g.res.push_back(e.second); }
      for (size_t i = 0; i < n; ++i)
        for (const auto& e : radj[i]) { g.rsrc.push_back(e.first); g.rres.push_back(e.second); }
      return g;
    }

    Tag scorePath(const Prepared& s, const Param& p, const Tables& tab,
                  const std::vector<uint32_t>& peaks, const std::vector<uint8_t>& res,
                  int charge)
    {
      const auto& R = residues();
      const int k = static_cast<int>(peaks.size());
      Tag t;
      t.charge = charge;

      // Back-project every peak onto an estimate of the first peak's position.
      // The spread of those estimates is the m/z-fidelity statistic.
      std::vector<double> est(static_cast<size_t>(k));
      double cum = 0;
      est[0] = s.spec[peaks[0]].getMZ();
      for (int i = 0; i < k - 1; ++i)
      {
        cum += R.mass[res[static_cast<size_t>(i)]] / charge;
        est[static_cast<size_t>(i) + 1] = s.spec[peaks[static_cast<size_t>(i) + 1]].getMZ() - cum;
      }
      const double avg = std::accumulate(est.begin(), est.end(), 0.0) / k;
      double sse = 0;
      for (double e : est) { const double d = e - avg; sse += d * d; }
      // Express the deviation in tolerance units, so the null is scale-free and
      // one table serves every m/z and both tolerance modes.
      const double tol = tolAt(p, avg);
      t.p_mzfidelity = tab.mzFidelityP(k, tol > 0 ? sse / (tol * tol) : sse);

      int ranksum = 0, n_compl = 0;
      for (int i = 0; i < k; ++i)
      {
        ranksum += s.rank[peaks[static_cast<size_t>(i)]];
        if (s.has_compl[peaks[static_cast<size_t>(i)]] & (1u << (charge - 1))) ++n_compl;
      }
      t.p_intensity = tab.intensityP(k, static_cast<int>(s.spec.size()), ranksum);

      // Probability of seeing at least this many complemented peaks among k drawn
      // from the spectrum. Dropped when the spectrum has no complements at all.
      const bool use_compl = s.n_with_compl > 0;
      if (use_compl && n_compl > 0)
      {
        const unsigned N = static_cast<unsigned>(s.spec.size());
        const unsigned K = static_cast<unsigned>(s.n_with_compl);
        const unsigned n = static_cast<unsigned>(k);
        // The support is [max(0, n+K-N), min(K, n)], not [0, n]: when most peaks
        // are complemented, a tag CANNOT have few. Boost throws on an argument
        // outside it, so resolve the ends analytically instead of clamping into
        // a distribution that does not cover them.
        const unsigned lo = (n + K > N) ? n + K - N : 0u;
        const unsigned hi = std::min(K, n);
        const unsigned obs = static_cast<unsigned>(n_compl);
        if (obs <= lo)      t.p_complement = 1.0;   // cannot do worse; no evidence
        else if (obs > hi)  t.p_complement = 0.0;   // impossible under the null
        else
        {
          boost::math::hypergeometric_distribution<double> h(K, n, N);
          t.p_complement = boost::math::cdf(boost::math::complement(h, obs - 1u));
        }
      }

      t.evalue = fisher(t.p_intensity, t.p_mzfidelity, t.p_complement, use_compl ? 3 : 2);

      // Flanking masses from the FITTED positions, not the raw ones: they are the
      // only consumer of the fit, and they are what a downstream search matches on.
      // The low peak is assumed to be a y ion, so a b-derived tag reads reversed
      // and its flanks carry a one-water offset.
      t.cterm_mass = std::max(0.0, avg * charge - PROTON * charge - WATER);
      t.nterm_mass = std::max(0.0, s.precursor_mass - ((avg + cum) * charge - PROTON * charge));

      // Traversal runs low->high m/z, which is C->N; store N->C.
      t.seq.reserve(static_cast<size_t>(k - 1));
      for (int i = k - 2; i >= 0; --i) t.seq.push_back(R.code[res[static_cast<size_t>(i)]]);
      t.low_mz = s.spec[peaks.front()].getMZ();
      return t;
    }

    /// Greedily extend one terminus. Selection is by smallest |m/z error|, tied
    /// on better intensity rank, then peak index -- all three are needed, since
    /// ties at equal error are common and a partial order leaves the output
    /// dependent on traversal accidents.
    int extendPath(const Prepared& s, const Param& p, const Graph& g,
                   std::vector<uint32_t>& peaks, std::vector<uint8_t>& res,
                   int charge, bool forward)
    {
      const auto& R = residues();
      int added = 0;
      std::vector<bool> used(s.spec.size(), false);
      for (uint32_t i : peaks) used[i] = true;

      while (added < p.max_extension)
      {
        const uint32_t cur = forward ? peaks.back() : peaks.front();
        const uint32_t e0 = forward ? g.off[cur] : g.roff[cur];
        const uint32_t e1 = forward ? g.off[cur + 1] : g.roff[cur + 1];
        int best = -1; uint8_t best_r = 0; double best_err = 0;
        for (uint32_t e = e0; e < e1; ++e)
        {
          const uint32_t cand = forward ? g.dst[e] : g.rsrc[e];
          const uint8_t rr = forward ? g.res[e] : g.rres[e];
          if (used[cand]) continue;
          const double err = forward
              ? std::fabs(s.spec[cand].getMZ() - (s.spec[cur].getMZ() + R.mass[rr] / charge))
              : std::fabs(s.spec[cur].getMZ() - (s.spec[cand].getMZ() + R.mass[rr] / charge));
          const bool better = best < 0 || err < best_err - 1e-12 ||
              (std::fabs(err - best_err) <= 1e-12 &&
               (s.rank[cand] < s.rank[static_cast<size_t>(best)] ||
                (s.rank[cand] == s.rank[static_cast<size_t>(best)] &&
                 cand < static_cast<uint32_t>(best))));
          if (better) { best = static_cast<int>(cand); best_r = rr; best_err = err; }
        }
        if (best < 0) break;
        used[static_cast<size_t>(best)] = true;
        if (forward) { peaks.push_back(static_cast<uint32_t>(best)); res.push_back(best_r); }
        else { peaks.insert(peaks.begin(), static_cast<uint32_t>(best));
               res.insert(res.begin(), best_r); }
        ++added;
      }
      return added;
    }
  }

  // ------------------------------------------------------------------ API

  std::vector<Tag> tagSpectrum(const MSSpectrum& in, double precursor_mz, int charge,
                               const Param& p, const Tables& tables)
  {
    std::vector<Tag> out;
    if (p.tag_length < 1 || in.empty() || precursor_mz <= 0) return out;
    if (charge <= 0) charge = 2;

    const Prepared s = prepare(in, precursor_mz, charge, p);
    if (s.spec.size() < static_cast<size_t>(p.tag_length) + 1) return out;

    std::vector<Tag> scored;
    size_t n_seeds = 0;

    for (int z = 1; z <= s.n_frag_charges; ++z)
    {
      const Graph g = buildGraph(s, p, z);
      std::vector<uint32_t> peaks;
      std::vector<uint8_t> res;
      struct Frame { uint32_t node, edge; };

      for (uint32_t start = 0; start < s.spec.size(); ++start)
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
            std::vector<uint8_t> rs = res;
            int grew = 0;
            if (p.max_extension > 0)
              grew = extendPath(s, p, g, pk, rs, z, true)
                   + extendPath(s, p, g, pk, rs, z, false);
            // Rescore at the realised length: a seed's statistic does not
            // describe an extended tag, and the per-length nulls are already built.
            Tag t = scorePath(s, p, tables, pk, rs, z);
            t.extended = grew > 0;
            scored.push_back(std::move(t));
            res.pop_back(); peaks.pop_back(); st.pop_back();
            continue;
          }
          if (f.edge >= g.off[f.node + 1])
          {
            st.pop_back();
            if (!res.empty()) { res.pop_back(); peaks.pop_back(); }
            continue;
          }
          const uint32_t e = f.edge++;
          const uint32_t nxt = g.dst[e];
          peaks.push_back(nxt);
          res.push_back(g.res[e]);
          st.push_back({nxt, g.off[nxt]});
        }
      }
    }

    // E-value: expected number of tags this good by chance. Valid under arbitrary
    // dependence between tags, by linearity of expectation.
    for (Tag& t : scored) t.evalue *= static_cast<double>(std::max<size_t>(n_seeds, 1));

    // stable_sort, not sort: the comparator is not total (it ignores charge and
    // the flanks), and 99.4% of duplicate groups are identical in every emitted
    // field today -- but a future field would make the survivor choice
    // implementation-defined. Stability costs nothing at these sizes.
    std::stable_sort(scored.begin(), scored.end(), [](const Tag& a, const Tag& b) {
      if (a.evalue != b.evalue) return a.evalue < b.evalue;
      if (a.seq != b.seq) return a.seq < b.seq;
      return a.low_mz < b.low_mz;
    });

    // Drop tags identical in every reported field, before the output cap so the
    // cap counts distinct results.
    //
    // Extension reaches the same physical tag by different N/C splits: measured
    // with extension 6, 28% of rows are duplicates and 99.4% of duplicate groups
    // agree in every emitted field. Keying on the exact reported values -- rather
    // than on rounded ones -- keeps the 0.6% that genuinely differ (flanks
    // 0.6 mDa apart with different E-values) and makes the survivor unambiguous,
    // since the duplicates are indistinguishable.
    std::set<std::tuple<std::string, int, double, double>> seen;
    for (Tag& t : scored)
    {
      if (p.max_evalue > 0 && t.evalue > p.max_evalue) break;
      if (!seen.emplace(t.seq, t.charge, t.nterm_mass, t.cterm_mass).second) continue;
      out.push_back(std::move(t));
      if (p.max_tag_count > 0 && static_cast<int>(out.size()) >= p.max_tag_count) break;
    }
    return out;
  }
}
