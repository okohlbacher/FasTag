// Self-check for ranksum.h. Build & run:
//   c++ -std=c++17 -O2 -o ranksum_test ranksum_test.cpp && ./ranksum_test
#include <FasTag/ranksum.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++failures; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// Exhaustive ground truth: enumerate every k-subset of {1..n}, the same thing
// DirecTag's CalculateIRBins_R does. Only tractable for tiny inputs -- which is
// the entire point of the DP.
static std::map<int, double> bruteForce(int k, int n)
{
  std::map<int, double> hist;
  std::vector<int> idx(k);
  for (int i = 0; i < k; ++i) idx[i] = i + 1;
  while (true)
  {
    int sum = 0;
    for (int i = 0; i < k; ++i) sum += idx[i];
    hist[sum] += 1.0;

    int i = k - 1;
    while (i >= 0 && idx[i] == n - (k - 1 - i)) --i;
    if (i < 0) break;
    ++idx[i];
    for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
  }
  return hist;
}

static double binom(int n, int k)
{
  if (k < 0 || k > n) return 0.0;
  double r = 1.0;
  for (int i = 0; i < k; ++i) r = r * (n - i) / (i + 1);
  return r;
}

int main()
{
  // 1. DP must equal exhaustive enumeration, bin for bin.
  int compared = 0;
  for (int k = 1; k <= 5; ++k)
  {
    for (int n = k; n <= 12; ++n)
    {
      const std::vector<double> dp = fastag::ranksumCounts(k, n);
      const std::map<int, double> bf = bruteForce(k, n);
      for (size_t s = 0; s < dp.size(); ++s)
      {
        auto it = bf.find(static_cast<int>(s));
        const double expect = (it == bf.end()) ? 0.0 : it->second;
        CHECK(dp[s] == expect, "k=%d n=%d s=%zu dp=%.0f brute=%.0f", k, n, s, dp[s], expect);
      }
      ++compared;
    }
  }
  std::printf("1. DP == exhaustive enumeration over %d (k,n) pairs\n", compared);

  // 2. Total count must be C(n,k).
  for (int k = 1; k <= 8; ++k)
  {
    for (int n = k; n <= 60; n += 7)
    {
      const std::vector<double> c = fastag::ranksumCounts(k, n);
      double tot = 0.0;
      for (double v : c) tot += v;
      const double want = binom(n, k);
      CHECK(std::fabs(tot - want) <= want * 1e-9, "k=%d n=%d total=%.0f want=%.0f", k, n, tot, want);
    }
  }
  std::printf("2. sum(counts) == C(n,k)\n");

  // 3. CDF is monotone non-decreasing and reaches 1.
  for (int k = 2; k <= 8; ++k)
  {
    for (int n = 20; n <= 150; n += 43)
    {
      const std::vector<float> cdf = fastag::ranksumCdf(k, n);
      CHECK(!cdf.empty(), "empty cdf k=%d n=%d", k, n);
      for (size_t s = 1; s < cdf.size(); ++s)
      {
        CHECK(cdf[s] >= cdf[s - 1], "cdf not monotone k=%d n=%d s=%zu", k, n, s);
      }
      CHECK(std::fabs(cdf.back() - 1.0f) < 1e-5f, "cdf end %.6f k=%d n=%d", cdf.back(), k, n);
    }
  }
  std::printf("3. CDF monotone, ends at 1.0\n");

  // 4. The minimum possible rank sum is k(k+1)/2; nothing below it may have mass.
  for (int k = 2; k <= 6; ++k)
  {
    const std::vector<double> c = fastag::ranksumCounts(k, 40);
    const int smin = k * (k + 1) / 2;
    for (int s = 0; s < smin; ++s) CHECK(c[s] == 0.0, "mass below min sum k=%d s=%d", k, s);
    CHECK(c[smin] == 1.0, "exactly one way to hit min sum k=%d", k);
  }
  std::printf("4. support starts at k(k+1)/2 with a single subset\n");

  // 5. compat_offset shifts the table by exactly one bin.
  {
    const std::vector<float> a = fastag::ranksumCdf(4, 50, false);
    const std::vector<float> b = fastag::ranksumCdf(4, 50, true);
    CHECK(b.size() == a.size() + 1, "compat size %zu vs %zu", b.size(), a.size());
    for (size_t s = 0; s < a.size(); ++s) CHECK(a[s] == b[s + 1], "compat shift at %zu", s);
  }
  std::printf("5. compat_offset shifts lookups by one bin\n");

  // 6. Gate 1: the full table set must build in under a second.
  {
    const auto t0 = std::chrono::steady_clock::now();
    size_t entries = 0;
    fastag::RanksumTables tables;
    for (int len = 3; len <= 7; ++len)
    {
      for (int n = 4; n <= 150; ++n) entries += tables.get(len + 1, n).size();
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("6. lengths 3-7 x peaks 4-150: %zu entries in %.3f s\n", entries, secs);
    CHECK(secs < 1.0, "Gate 1 build time %.3f s exceeds 1 s", secs);
  }

  // 7. Out-of-range lookups are conservative, never out-of-bounds reads.
  {
    fastag::RanksumTables tables;
    CHECK(tables.pvalue(4, 50, 1 << 20) == 1.0, "huge ranksum should give p=1");
    CHECK(tables.pvalue(4, 50, -1) == 0.0, "negative ranksum should give p=0");
    CHECK(tables.pvalue(9, 4, 10) == 1.0, "k>n should give p=1");
  }
  std::printf("7. out-of-range lookups clamp instead of reading out of bounds\n");

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
