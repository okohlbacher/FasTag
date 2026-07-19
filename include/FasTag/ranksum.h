// FasTag -- intensity rank-sum p-value tables for DirecTag-style sequence tagging.
//
// Replaces DirecTag's CalculateIRBins_R (directagAPIConfig.h:193), which
// enumerates every C(n, k) subset of peak ranks. That is ~8.8e13 leaves for
// tag length 7 at 150 peaks -- it does not terminate. The same distribution is
// a restricted-partition count computable by DP in O(k*n*smax).
//
// Reimplemented from the algorithm description, not transcribed: freicore is
// Apache-2.0, OpenMS is BSD-3-Clause.
#pragma once

#include <cstddef>
#include <vector>

namespace fastag
{
  /// Number of distinct-rank subsets summing to each value.
  /// counts[s] = #{ S subset of {1..n} : |S| = k, sum(S) = s }, for s in [0, smax].
  ///
  /// Counts are held as double. Exact while C(n,k) <= 2^53; C(150,8) = 5.3e12 is
  /// comfortably inside that. Beyond it the ratios stay accurate but the integer
  /// counts do not -- only the ratios are ever used.
  inline std::vector<double> ranksumCounts(int k, int n)
  {
    if (k <= 0 || n <= 0 || k > n) return {};

    // Max sum picks the top k ranks: n + (n-1) + ... + (n-k+1).
    const int smax = k * n - k * (k - 1) / 2;

    // dp[j][s] = ways to choose j distinct values from those considered so far
    // summing to s. Values are introduced ascending; j and s descend so each
    // value is used at most once (0/1 knapsack over distinct parts).
    std::vector<std::vector<double>> dp(k + 1, std::vector<double>(smax + 1, 0.0));
    dp[0][0] = 1.0;

    for (int v = 1; v <= n; ++v)
    {
      for (int j = (k < v ? k : v); j >= 1; --j)
      {
        std::vector<double>& cur = dp[j];
        const std::vector<double>& prev = dp[j - 1];
        for (int s = smax; s >= v; --s)
        {
          if (prev[s - v] != 0.0) cur[s] += prev[s - v];
        }
      }
    }
    return dp[k];
  }

  /// Cumulative p-values: cdf[s] = P(ranksum <= s) under the null that the k
  /// peaks are drawn uniformly at random from the n ranked peaks.
  ///
  /// Mirrors CalculateIRBins (directagAPIConfig.h:212-221): running sum divided
  /// by the total. Stored as float -- these index a binned lookup, so double
  /// precision buys nothing and costs half the cache.
  ///
  /// @param compat_offset reproduces DirecTag's off-by-one. DirecTag enumerates
  ///   sums from 0 but initialises the observed sum to 1 (directagSpectrum.cpp:641),
  ///   shifting every lookup by one bin. Pass true for bit-compatibility with the
  ///   reference, false for the corrected statistic.
  inline std::vector<float> ranksumCdf(int k, int n, bool compat_offset = false)
  {
    const std::vector<double> counts = ranksumCounts(k, n);
    if (counts.empty()) return {};

    double total = 0.0;
    for (double c : counts) total += c;
    if (total <= 0.0) return {};

    std::vector<float> cdf(counts.size() + (compat_offset ? 1 : 0), 0.0f);
    double running = 0.0;
    const size_t base = compat_offset ? 1 : 0;
    for (size_t s = 0; s < counts.size(); ++s)
    {
      running += counts[s];
      cdf[s + base] = static_cast<float>(running / total);
    }
    if (compat_offset && !cdf.empty()) cdf[0] = 0.0f;
    return cdf;
  }

  /// Lazily-built cache of CDF tables keyed by (k, n).
  ///
  /// Built on demand rather than over the full cross-product: memory is
  /// O(k_max^2 * n_max^2) if built eagerly (2.5 MB at k<=8/n<=150, but 29 MB at
  /// n<=500 and 125 MB for lengths 3-15), whereas real spectra cluster around a
  /// few peak counts. Lookups are bounds-checked -- DirecTag reads out of bounds
  /// if MaxPeakCount is raised without regenerating its on-disk cache.
  class RanksumTables
  {
  public:
    explicit RanksumTables(bool compat_offset = false) : compat_offset_(compat_offset) {}

    const std::vector<float>& get(int k, int n)
    {
      Key key{k, n};
      for (size_t i = 0; i < keys_.size(); ++i)
      {
        if (keys_[i].k == key.k && keys_[i].n == key.n) return tables_[i];
      }
      keys_.push_back(key);
      tables_.push_back(ranksumCdf(k, n, compat_offset_));
      return tables_.back();
    }

    /// p-value for an observed rank sum; 1.0 if out of range (conservative).
    double pvalue(int k, int n, int ranksum)
    {
      const std::vector<float>& cdf = get(k, n);
      if (cdf.empty()) return 1.0;
      if (ranksum < 0) return 0.0;
      if (static_cast<size_t>(ranksum) >= cdf.size()) return 1.0;
      return cdf[ranksum];
    }

  private:
    struct Key { int k; int n; };
    bool compat_offset_;
    std::vector<Key> keys_;                  // ponytail: linear scan, a few hundred
    std::vector<std::vector<float>> tables_; // entries; swap for a hash map if it shows in a profile
  };
}
