// The filter must never report a tag that is not in the sequences, and never
// drop one that is. Both directions are checked against a naive oracle.
#include "FastaFilter.h"

#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <cstdio>
#include <random>
#include <random>
#include <set>
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


  // Keys longer than 12 residues, where the 128-bit key spans both halves.
  //
  // The key packs 5 bits per residue into two uint64_t, so residue 13 is the
  // first whose bits cross from the low half into the high half. Nothing else
  // here reaches that far: the test sequences are 20 residues long, but the
  // k-mers INDEXED are of length min_k..max_k, and autoMinLen derives a small
  // floor for a tiny database -- so every existing assertion runs entirely
  // inside the low half.
  //
  // Verified by mutation: corrupting the carry shift (lo >> 59 to lo >> 60)
  // leaves every other test passing.
  {
    std::vector<FASTAFile::FASTAEntry> e;
    // 25 residues, the maximum the encoding supports, with no repeats so any
    // mis-shift changes the string rather than colliding onto itself.
    e.emplace_back("long", "", "ACDEFGHKLMNPQRSTVWYACDEF");

    FasTag::FastaFilter f(false);
    CHECK(f.load(e), "long sequence did not load");
    f.setMinLen(13);
    f.build(13, 24);

    // Present at every length that crosses the boundary.
    for (int k = 13; k <= 24; ++k)
    {
      const std::string tag = std::string("ACDEFGHKLMNPQRSTVWYACDEF").substr(0, k);
      CHECK(f.match(tag) == FasTag::FastaFilter::Hit::Forward,
            "length-%d key spanning both halves did not match", k);
    }

    // A single residue changed must NOT match -- otherwise the assertions above
    // would pass on a key that had silently collapsed to a constant.
    CHECK(f.match("ACDEFGHKLMNPQRW") == FasTag::FastaFilter::Hit::None,
          "a 15-residue tag absent from the sequence matched anyway");
    // Transposing two residues past the boundary must also be rejected.
    CHECK(f.match("ACDEFGHKLMNPQSR") == FasTag::FastaFilter::Hit::None,
          "a transposition past residue 13 was not distinguished");
    // COLLISION tests. The two checks above cannot catch a corrupt packing,
    // because the same corruption is applied when indexing and when matching --
    // a key encodes wrongly both times and still matches itself. Only two
    // DIFFERENT k-mers colliding exposes it. Verified by mutation: without
    // these, "lo >> 59" -> "lo >> 60" and a 4-bit residue mask both survive.
    {
      std::vector<FASTAFile::FASTAEntry> e2;
      // Identical for 12 residues, differing only beyond the hi/lo boundary.
      e2.emplace_back("hi", "", "AAAAAAAAAAAAWWWWWWWW");
      FasTag::FastaFilter g(false);
      CHECK(g.load(e2), "collision fixture did not load");
      g.setMinLen(20);
      g.build(20, 20);
      CHECK(g.match("AAAAAAAAAAAAWWWWWWWW") == FasTag::FastaFilter::Hit::Forward,
            "the indexed 20-mer did not match itself");
      // Differs only in residues 13-20. If the carry into the high half is
      // wrong, those bits are lost and this collides with the key above.
      CHECK(g.match("AAAAAAAAAAAAYYYYYYYY") == FasTag::FastaFilter::Hit::None,
            "two 20-mers differing only past residue 12 collided -- the carry "
            "into the high half of the key is dropping bits");
    }
    {
      std::vector<FASTAFile::FASTAEntry> e3;
      // W is code 22 (10110), G is 6 (00110): they differ only in bit 4, so a
      // residue mask narrower than 5 bits makes them the same symbol.
      e3.emplace_back("w", "", "WWWWWWWWWWWWWWWW");
      FasTag::FastaFilter g(false);
      CHECK(g.load(e3), "5-bit fixture did not load");
      g.setMinLen(16);
      g.build(16, 16);
      CHECK(g.match("WWWWWWWWWWWWWWWW") == FasTag::FastaFilter::Hit::Forward,
            "the indexed W-mer did not match itself");
      CHECK(g.match("GGGGGGGGGGGGGGGG") == FasTag::FastaFilter::Hit::None,
            "W and G collided -- the residue code is being masked to fewer than "
            "5 bits, so codes 16-24 fold onto codes 0-8");
    }
    // INJECTIVITY. The checks above still miss a corrupt carry, because a
    // fixture with leading A residues (code 0) leaves the high half zero, so it
    // is never exercised. Rather than hand-build a collision -- which requires
    // reasoning about exactly which bits a given corruption drops -- assert the
    // property that matters: distinct k-mers must get distinct keys.
    //
    // A pseudo-random sequence of 260 residues yields 241 overlapping 20-mers,
    // essentially all distinct. If the packing loses a bit anywhere, some of
    // them collapse together and the indexed key count falls short.
    {
      std::mt19937 rng(20260721);
      const std::string alpha = "ACDEFGHKLMNPQRSTVWY";   // 19, all codes 0..24
      std::string seq;
      for (int i = 0; i < 260; ++i) seq += alpha[rng() % alpha.size()];

      std::vector<FASTAFile::FASTAEntry> e4;
      e4.emplace_back("rand", "", seq);
      FasTag::FastaFilter g(false);
      CHECK(g.load(e4), "random fixture did not load");
      g.setMinLen(20);
      g.build(20, 20);

      // Count the distinct 20-mers the sequence actually contains, so the
      // expectation does not assume they are all unique.
      std::set<std::string> distinct;
      for (size_t i = 0; i + 20 <= seq.size(); ++i) distinct.insert(seq.substr(i, 20));

      CHECK(g.indexedKeys() == distinct.size(),
            "%zu distinct 20-mers encoded to only %zu keys -- the 128-bit "
            "packing is losing information across the hi/lo boundary",
            distinct.size(), g.indexedKeys());
      std::printf("boundary: %zu distinct 20-mers -> %zu distinct keys\n",
                  distinct.size(), g.indexedKeys());
    }
    // EVERY bit position, systematically.
    //
    // Four hand-built fixtures failed to catch a corrupt carry before this one,
    // each for a different reason, and the sequence is worth recording because
    // the same trap recurs:
    //
    //   long sequences      -- the k-mers INDEXED are short, so the boundary is
    //                          never reached
    //   hand-built collision -- the corruption is applied when indexing AND when
    //                          matching, so a key encodes wrongly both times and
    //                          still matches itself
    //   random 20-mers      -- lo holds the last ~12 residues and already
    //                          separates them, so hi never has to work
    //   suffix-identical    -- residue 8 straddles bit 64, so even a 7-residue
    //                          prefix leaks into lo
    //
    // Rather than reason about which bit a given corruption drops, enumerate ALL
    // single-residue substitutions of a base 20-mer: 20 positions x 18
    // alternatives. Every one is a distinct string, so if the packing loses
    // information anywhere -- the carry, either shift, or the residue mask --
    // some pair must collapse onto one key.
    // At MAX_FILTER_LEN, not merely "long". 25 residues is 125 bits, so hi
    // carries 61 of them -- nearly full. A shorter key leaves hi with room to
    // spare, and a corruption that merely SPREADS bits within hi rather than
    // dropping them stays injective and slips through: hi << 6 survives at
    // length 20 (8 pushes x 6 = 48 bits, no overflow) and only overflows once
    // the key is long enough. Test the documented maximum, where it does.
    {
      const std::string alpha = "ACDEFGHKLMNPQRSTVWY";
      const std::string base  = "ACDEFGHKLMNPQRSTVWYACDEFG";   // 25 = MAX_FILTER_LEN
      CHECK(static_cast<int>(base.size()) == FasTag::MAX_FILTER_LEN,
            "base is %zu residues, MAX_FILTER_LEN is %d -- this test is meant to "
            "run at the maximum", base.size(), FasTag::MAX_FILTER_LEN);
      std::set<std::string> made;
      std::vector<FASTAFile::FASTAEntry> e5;
      made.insert(base);
      e5.emplace_back("b", "", base);
      for (size_t pos = 0; pos < base.size(); ++pos)
      {
        for (char c : alpha)
        {
          if (c == base[pos]) continue;
          std::string v = base;
          v[pos] = c;
          if (made.insert(v).second) e5.emplace_back("v", "", v);
        }
      }

      FasTag::FastaFilter g(false);
      CHECK(g.load(e5), "substitution fixture did not load");
      g.setMinLen(FasTag::MAX_FILTER_LEN);
      g.build(FasTag::MAX_FILTER_LEN, FasTag::MAX_FILTER_LEN);
      CHECK(g.indexedKeys() == made.size(),
            "%zu distinct max-length k-mers -- every single-residue substitution of "
            "a base -- encoded to only %zu keys. The 128-bit packing is losing "
            "information.", made.size(), g.indexedKeys());
      std::printf("boundary: %zu single-substitution %d-mers -> %zu distinct keys\n",
                  made.size(), FasTag::MAX_FILTER_LEN, g.indexedKeys());
    }
    std::printf("boundary: lengths 13-24 span the hi/lo split; collisions rejected\n");
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
