// Distribution of a sum of k distinct ranks drawn from 1..n.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
//
// This is the one piece of mathematics FasTag keeps in-house: neither OpenMS nor
// Boost provides the exact null of a rank sum over distinct ranks. Boost offers
// the normal approximation to Mann-Whitney, which is not usable here because the
// interesting tags sit far in the tail, where that approximation is worst.
//
// DirecTag computes the same distribution by enumerating every C(n, k) subset.
// Measured on the reference binary: 0.1 s at k=4, 173 s at k=7, ~30 min at k=8,
// since the work grows ~10x per residue. The distribution is a restricted-
// partition count, so a DP gets it in O(k*n*smax) -- all tables for tag lengths
// 3-16 build in 0.15 s. That is what makes tag lengths above four usable at all.
#pragma once

#include <cstddef>
#include <vector>

namespace FasTag
{
  /// Number of distinct-rank subsets summing to each value:
  /// counts[s] = #{ S subset of {1..n} : |S| = k, sum(S) = s }.
  ///
  /// Counts are held as double, exact while C(n,k) <= 2^53. C(150,8) = 5.3e12 is
  /// comfortably inside; beyond it the ratios stay accurate, and only ratios are
  /// ever used.
  inline std::vector<double> ranksumCounts(int k, int n)
  {
    if (k <= 0 || n <= 0 || k > n) return {};
    const int smax = k * n - k * (k - 1) / 2;   // sum of the top k ranks

    // dp[j][s] = ways to choose j distinct values summing to s from those seen so
    // far. Values ascend; j and s descend, so each value is used at most once.
    std::vector<std::vector<double>> dp(static_cast<size_t>(k) + 1,
                                        std::vector<double>(static_cast<size_t>(smax) + 1, 0.0));
    dp[0][0] = 1.0;
    for (int v = 1; v <= n; ++v)
      for (int j = (k < v ? k : v); j >= 1; --j)
      {
        auto& cur = dp[static_cast<size_t>(j)];
        const auto& prev = dp[static_cast<size_t>(j) - 1];
        for (int s = smax; s >= v; --s)
          if (prev[static_cast<size_t>(s - v)] != 0.0)
            cur[static_cast<size_t>(s)] += prev[static_cast<size_t>(s - v)];
      }
    return dp[static_cast<size_t>(k)];
  }

  /// Cumulative form: cdf[s] = P(rank sum <= s) under the null that the k peaks
  /// are drawn uniformly from the n ranked peaks.
  ///
  /// Stored as float: these index a binned lookup, so double precision buys
  /// nothing and costs twice the cache.
  inline std::vector<float> ranksumCdf(int k, int n)
  {
    const std::vector<double> counts = ranksumCounts(k, n);
    if (counts.empty()) return {};
    double total = 0.0;
    for (double c : counts) total += c;
    if (total <= 0.0) return {};

    std::vector<float> cdf(counts.size(), 0.0f);
    double running = 0.0;
    for (size_t s = 0; s < counts.size(); ++s)
    {
      running += counts[s];
      cdf[s] = static_cast<float>(running / total);
    }
    return cdf;
  }
}
