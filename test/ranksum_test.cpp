// The one proof worth running in isolation: the DP equals exhaustive enumeration.
//   The DP is what replaces DirecTag's C(n,k) subset walk, so if it is wrong,
//   every intensity p-value is wrong and nothing downstream would notice.
#include <FasTag/ranksum.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

static int failures = 0;
#define CHECK(c, ...) do { if (!(c)) { ++failures; std::printf("FAIL %d: ", __LINE__); \
  std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

/// Exhaustive ground truth -- exactly what DirecTag's CalculateIRBins_R does.
/// Tractable only for tiny inputs, which is the whole point.
static std::map<int, double> brute(int k, int n)
{
  std::map<int, double> h;
  std::vector<int> idx(k);
  for (int i = 0; i < k; ++i) idx[i] = i + 1;
  for (;;)
  {
    int s = 0;
    for (int i = 0; i < k; ++i) s += idx[i];
    h[s] += 1.0;
    int i = k - 1;
    while (i >= 0 && idx[i] == n - (k - 1 - i)) --i;
    if (i < 0) break;
    ++idx[i];
    for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
  }
  return h;
}

static double binom(int n, int k)
{
  double r = 1.0;
  for (int i = 0; i < k; ++i) r = r * (n - i) / (i + 1);
  return r;
}

int main()
{
  int pairs = 0;
  for (int k = 1; k <= 5; ++k)
    for (int n = k; n <= 12; ++n)
    {
      const auto dp = FasTag::ranksumCounts(k, n);
      const auto bf = brute(k, n);
      for (size_t s = 0; s < dp.size(); ++s)
      {
        auto it = bf.find(static_cast<int>(s));
        const double want = (it == bf.end()) ? 0.0 : it->second;
        CHECK(dp[s] == want, "k=%d n=%d s=%zu dp=%.0f brute=%.0f", k, n, s, dp[s], want);
      }
      ++pairs;
    }
  std::printf("1. DP == exhaustive enumeration over %d (k,n) pairs\n", pairs);

  for (int k = 1; k <= 8; ++k)
    for (int n = k; n <= 60; n += 7)
    {
      double tot = 0;
      for (double v : FasTag::ranksumCounts(k, n)) tot += v;
      const double want = binom(n, k);
      CHECK(std::fabs(tot - want) <= want * 1e-9, "k=%d n=%d total mismatch", k, n);
    }
  std::printf("2. sum(counts) == C(n,k)\n");

  for (int k = 2; k <= 8; ++k)
    for (int n = 20; n <= 150; n += 43)
    {
      const auto c = FasTag::ranksumCdf(k, n);
      CHECK(!c.empty(), "empty cdf k=%d n=%d", k, n);
      for (size_t s = 1; s < c.size(); ++s)
        CHECK(c[s] >= c[s - 1], "cdf not monotone k=%d n=%d", k, n);
      CHECK(std::fabs(c.back() - 1.0f) < 1e-5f, "cdf ends at %.6f", c.back());
    }
  std::printf("3. CDF monotone and reaches 1\n");

  // Support starts at the sum of the k smallest ranks, reachable exactly one way.
  for (int k = 2; k <= 6; ++k)
  {
    const auto c = FasTag::ranksumCounts(k, 40);
    const int smin = k * (k + 1) / 2;
    for (int s = 0; s < smin; ++s) CHECK(c[s] == 0.0, "mass below min sum k=%d", k);
    CHECK(c[smin] == 1.0, "min sum should be reachable exactly once, k=%d", k);
  }
  std::printf("4. support starts at k(k+1)/2 with a single subset\n");

  CHECK(FasTag::ranksumCounts(0, 5).empty(), "k=0 must be empty");
  CHECK(FasTag::ranksumCounts(9, 4).empty(), "k>n must be empty");
  std::printf("5. degenerate (k, n) return empty rather than misbehaving\n");

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
