// The filter must never report a tag that is not in the sequences, and never
// drop one that is. Both directions are checked against a naive oracle.
#include "FastaFilter.h"

#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace OpenMS;

static int failures = 0;
#define CHECK(c, ...) do { if (!(c)) { ++failures; std::printf("FAIL %d: ", __LINE__); \
  std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static std::vector<FASTAFile::FASTAEntry> entries(const std::vector<std::string>& seqs)
{
  std::vector<FASTAFile::FASTAEntry> v;
  for (size_t i = 0; i < seqs.size(); ++i)
    v.emplace_back(String("p") + String(i), "", seqs[i]);
  return v;
}

static std::string fold(std::string s)
{
  for (char& c : s) if (c == 'I') c = 'L';
  return s;
}

int main()
{
  std::mt19937 rng(1234);
  const std::string alpha = "GASPVTCLNDQKEMHFRYW";
  std::uniform_int_distribution<int> pick(0, int(alpha.size()) - 1);

  std::string prot;
  for (int i = 0; i < 4000; ++i) prot.push_back(alpha[pick(rng)]);

  // 1. soundness and completeness against a naive substring oracle
  {
    FasTag::FastaFilter f(true);
    std::string err;
    CHECK(f.load(entries({prot}), &err), "load failed: %s", err.c_str());
    f.setMinLen(8);
    f.build(8, 8);
    const std::string p = fold(prot);
    int checked = 0, wrong = 0;
    for (int i = 0; i < 4000; ++i)
    {
      std::string t;
      for (int j = 0; j < 8; ++j) t.push_back(alpha[pick(rng)]);
      const std::string ft = fold(t), rev(ft.rbegin(), ft.rend());
      const bool want = p.find(ft) != std::string::npos || p.find(rev) != std::string::npos;
      const bool got = f.match(t) != FasTag::FastaFilter::Hit::None;
      if (want != got) ++wrong;
      ++checked;
    }
    // and every genuine 8-mer of the protein must be found
    int missed = 0;
    for (size_t i = 0; i + 8 <= prot.size(); i += 37)
      if (f.match(prot.substr(i, 8)) == FasTag::FastaFilter::Hit::None) ++missed;
    CHECK(wrong == 0, "%d/%d random tags misclassified", wrong, checked);
    CHECK(missed == 0, "%d genuine 8-mers not found", missed);
    std::printf("1. soundness+completeness vs naive oracle: %d random tags, 0 errors\n", checked);
  }

  // 2. orientation: a reversed tag is found only when both orientations are on
  {
    const std::string sub = prot.substr(100, 9);
    std::string rev(sub.rbegin(), sub.rend());
    FasTag::FastaFilter both(true), fwd(false);
    both.load(entries({prot})); both.setMinLen(9); both.build(9, 9);
    fwd.load(entries({prot}));  fwd.setMinLen(9);  fwd.build(9, 9);
    CHECK(both.match(sub) == FasTag::FastaFilter::Hit::Forward, "forward hit expected");
    CHECK(both.match(rev) == FasTag::FastaFilter::Hit::Reverse, "reverse hit expected");
    CHECK(fwd.match(rev) == FasTag::FastaFilter::Hit::None, "forward-only must reject a reversal");
    std::printf("2. orientation: fwd/rev distinguished; forward-only rejects reversals\n");
  }

  // 3. isobaric collapse: a tag reading N where the sequence spells GG is a
  //    mass-correct reading and must be accepted; without collapse it is not.
  {
    const std::string seq = "ARNDCQEGGHKLMFPSTWYV";   // contains GG
    FasTag::FastaFilter exact(true), iso(true);
    exact.load(entries({seq})); exact.setMinLen(4); exact.build(4, 6);
    iso.load(entries({seq}));   iso.setMinLen(4);   iso.deriveCollapses(0.04); iso.build(4, 6);
    CHECK(iso.collapseRules() > 0, "no collapse rules derived");
    // "QEGGH" read with N in place of GG is "QENH" -- 4 residues, so the
    // length floor does not interfere.
    CHECK(exact.match("QENH") == FasTag::FastaFilter::Hit::None,
          "exact matching must reject the N reading");
    CHECK(iso.match("QENH") != FasTag::FastaFilter::Hit::None,
          "collapse must accept N where the sequence spells GG");
    std::printf("3. isobaric collapse: %zu rules; N accepted where the sequence spells GG\n",
                iso.collapseRules());
  }

  // 4. the derived floor scales with database size, and short tags are the
  //    reason it exists: at length 3 the filter is the identity function.
  {
    FasTag::FastaFilter small(true), big(true);
    small.load(entries({prot.substr(0, 392)}));
    // Independently random, not prot repeated: repetition adds residues but no
    // new k-mers, so a repeated sequence would still hold only ~66% of 3-mers.
    std::string huge;
    huge.reserve(1000000);
    for (int i = 0; i < 1000000; ++i) huge.push_back(alpha[pick(rng)]);
    big.load(entries({huge}));
    CHECK(small.autoMinLen() < big.autoMinLen(), "floor must grow with database size (%d vs %d)",
          small.autoMinLen(), big.autoMinLen());
    CHECK(small.chanceRate(small.autoMinLen()) < 0.05, "derived floor should keep chance < 5%%");
    FasTag::FastaFilter noop(true);
    noop.load(entries({huge})); noop.setMinLen(3); noop.build(3, 3);
    int hits = 0;
    for (int i = 0; i < 500; ++i)
    {
      std::string t;
      for (int j = 0; j < 3; ++j) t.push_back(alpha[pick(rng)]);
      if (noop.match(t) != FasTag::FastaFilter::Hit::None) ++hits;
    }
    // Against a proteome-sized database every 3-mer is present, so the filter
    // degenerates to the identity function. That is why the floor exists.
    CHECK(hits > 490, "at length 3 the filter should accept nearly everything, got %d/500", hits);
    std::printf("4. derived floor %d (392 aa) vs %d (1M aa); at length 3 %d/500 random tags pass\n",
                small.autoMinLen(), big.autoMinLen(), hits);
  }

  // 5. degenerate input
  {
    FasTag::FastaFilter f(true);
    std::string err;
    CHECK(!f.load(entries({}), &err), "empty FASTA must be rejected");
    FasTag::FastaFilter g(true);
    CHECK(g.load(entries({"acdefghk"})), "lowercase must be accepted");
    g.setMinLen(4); g.build(4, 4);
    CHECK(g.match("ACDE") != FasTag::FastaFilter::Hit::None, "lowercase should have been folded");
    FasTag::FastaFilter h(true);
    h.load(entries({"AAXAA"})); h.setMinLen(3); h.build(3, 3);
    CHECK(h.match("AXA") == FasTag::FastaFilter::Hit::None, "ambiguity codes must not match");
    std::printf("5. degenerate input: empty rejected, lowercase folded, X not matchable\n");
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
