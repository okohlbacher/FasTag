// Deterministic spectrum subsampling: exact count, uniformity, reproducibility.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "SpectrumSampler.h"

#include <iostream>
#include <numeric>
#include <string>

using namespace FasTag;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    if (!ok) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
  }
  size_t count(const SampleMask& m) { return std::accumulate(m.begin(), m.end(), size_t(0),
                                        [](size_t a, char b){ return a + (b ? 1 : 0); }); }
}

int main()
{
  // Exact count.
  {
    auto m = sampleByCount(1000, 100, 42);
    check(m.size() == 1000, "mask length");
    check(count(m) == 100, "exactly n_want selected");
  }
  // Edge cases.
  check(count(sampleByCount(1000, 0, 1)) == 0, "n_want 0 selects nothing");
  check(count(sampleByCount(1000, 5000, 1)) == 1000, "n_want > n_total selects all");
  check(sampleByCount(0, 10, 1).empty(), "empty input");

  // Fraction rounds to nearest and is proportional.
  check(count(sampleByFraction(1000, 0.1, 7)) == 100, "10% of 1000 == 100");
  check(count(sampleByFraction(1000, 0.0, 7)) == 0, "0% selects nothing");
  check(count(sampleByFraction(1000, 1.0, 7)) == 1000, "100% selects all");
  check(count(sampleByFraction(3, 0.01, 7)) == 1, "tiny fraction still selects at least one");

  // Reproducible: same seed -> identical mask; different seed -> (almost surely) different.
  {
    auto a = sampleByCount(10000, 500, 123);
    auto b = sampleByCount(10000, 500, 123);
    auto c = sampleByCount(10000, 500, 124);
    check(a == b, "same seed reproduces the exact subset");
    check(a != c, "different seed gives a different subset");
    check(count(c) == 500, "different seed still exact count");
  }

  // Uniformity smoke test: over many seeds each index is picked ~ n_want/n_total.
  {
    const size_t N = 200, K = 20, TRIALS = 4000;
    std::vector<size_t> hits(N, 0);
    for (uint32_t s = 0; s < TRIALS; ++s)
    {
      auto m = sampleByCount(N, K, s + 1);
      for (size_t i = 0; i < N; ++i) if (m[i]) ++hits[i];
    }
    // Expected hit rate K/N = 0.1; check every index lands within a loose band.
    double exp_rate = double(K) / double(N);
    bool uniform = true;
    for (size_t i = 0; i < N; ++i)
    {
      double r = double(hits[i]) / double(TRIALS);
      if (r < exp_rate * 0.6 || r > exp_rate * 1.4) uniform = false;
    }
    check(uniform, "selection is roughly uniform across indices");
  }

  if (failures == 0) std::cout << "sampler_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
