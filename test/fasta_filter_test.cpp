// Correctness properties for the FASTA filter.
//
//   c++ -std=c++17 -O2 -o fasta_filter_test fasta_filter_test.cpp && ./fasta_filter_test
//
// Every property runs at length >= 8. At the default seed length 3 all of them
// pass trivially -- the tagger emits all 19^3 = 6859 3-mers and every one occurs
// in any real proteome, so a matcher that returns `true` would score green.
#include <FasTag/fasta_filter.h>

#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(c, ...) do { if(!(c)) { ++failures; std::printf("  FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while(0)

static const std::string AA = "GASPVTCLNDQKEMHFRYW";

static std::string randomSeq(std::mt19937& rng, int n)
{
  std::uniform_int_distribution<int> d(0, static_cast<int>(AA.size()) - 1);
  std::string s;
  for (int i = 0; i < n; ++i) s.push_back(AA[d(rng)]);
  return s;
}

static void writeFasta(const std::string& path, const std::vector<std::string>& seqs,
                       int wrap = 60)
{
  FILE* f = fopen(path.c_str(), "w");
  for (size_t i = 0; i < seqs.size(); ++i)
  {
    fprintf(f, ">seq%zu test sequence\n", i);
    for (size_t j = 0; j < seqs[i].size(); j += static_cast<size_t>(wrap))
      fprintf(f, "%s\n", seqs[i].substr(j, static_cast<size_t>(wrap)).c_str());
  }
  fclose(f);
}

/// Independent oracle: naive search over PROPERLY PARSED sequences.
/// Deliberately does not share the hash index. It necessarily shares the I/L and
/// orientation rules -- so it validates the index, not those decisions. The
/// independent check for orientation is the decoy test below.
static bool naiveMatch(const std::vector<std::string>& seqs, const std::string& tag,
                       bool both)
{
  std::string t = tag;
  for (char& c : t) if (c == 'I') c = 'L';
  std::string r(t.rbegin(), t.rend());
  for (const auto& s0 : seqs)
  {
    std::string s = s0;
    for (char& c : s) if (c == 'I') c = 'L';
    if (s.find(t) != std::string::npos) return true;
    if (both && s.find(r) != std::string::npos) return true;
  }
  return false;
}

int main()
{
  std::mt19937 rng(20260719u);
  const std::string F = "/tmp/fastag_filter_test.fasta";

  // ---------------------------------------------------------------- encoding
  {
    // The bug this replaces: 5 bits x 13 = 65 > 64, so a uint64 packing collides
    // silently at length >= 13, which extension reaches (3 + 2*6 = 15).
    CHECK(fastag::MAX_FILTER_LEN >= 15, "encoder must cover seed 3 + extension 6");
    std::set<__uint128_t> seen;
    int collisions = 0;
    for (int trial = 0; trial < 20000; ++trial)
    {
      const std::string s = randomSeq(rng, 15);
      const auto e = fastag::FastaFilter::encode(s.data(), 15);
      if (!seen.insert(e).second) ++collisions;
    }
    CHECK(collisions == 0, "%d collisions among 20000 random 15-mers", collisions);
    // and distinct 15-mers differing only in the FIRST residue must differ
    std::string a = randomSeq(rng, 15), b = a;
    b[0] = (a[0] == 'G') ? 'W' : 'G';
    CHECK(fastag::FastaFilter::encode(a.data(), 15) != fastag::FastaFilter::encode(b.data(), 15),
          "leading residue must survive the encoding");
    std::printf("1. encoding: 15-mers distinct, leading residue preserved (uint64 would collide)\n");
  }

  // --------------------------------------------------- soundness+completeness
  {
    std::vector<std::string> seqs{randomSeq(rng, 900), randomSeq(rng, 500)};
    writeFasta(F, seqs);
    fastag::FastaFilter filt(true);
    std::string err;
    CHECK(filt.load(F, &err), "load: %s", err.c_str());
    filt.setMinLen(8);
    filt.finalise();
    filt.buildRange(8, 12);

    int sound = 0, complete = 0, n = 0, pos = 0;
    for (int L = 8; L <= 12; ++L)
    {
      // real substrings (must all match), plus random tags (must match iff naive says so)
      for (int t = 0; t < 400; ++t)
      {
        const auto& s = seqs[t % 2];
        const size_t off = rng() % (s.size() - static_cast<size_t>(L));
        const std::string sub = s.substr(off, static_cast<size_t>(L));
        if (!filt.matches(sub)) ++complete;          // completeness violation
        ++pos;
        // reversed real substring must also match under orientation=both
        std::string rev(sub.rbegin(), sub.rend());
        if (!filt.matches(rev)) ++complete;
        ++pos;
      }
      for (int t = 0; t < 2000; ++t)
      {
        const std::string tag = randomSeq(rng, L);
        const bool got = filt.matches(tag);
        const bool want = naiveMatch(seqs, tag, true);
        if (got != want) ++sound;
        ++n;
      }
    }
    CHECK(sound == 0, "%d/%d disagreements with the naive oracle", sound, n);
    CHECK(complete == 0, "%d/%d true substrings missed", complete, pos);
    std::printf("2. soundness+completeness vs naive oracle: %d random tags, %d known substrings, 0 errors\n", n, pos);
  }

  // ----------------------------------------------------------- specificity
  {
    // A tag set drawn from an ABSENT protein must match at ~the chance rate,
    // not at the rate of the real one. Nothing in a "soundness" property forces
    // the filter to ever reject.
    std::vector<std::string> target{randomSeq(rng, 900)};
    writeFasta(F, target);
    fastag::FastaFilter filt(true);
    filt.load(F);
    filt.setMinLen(8); filt.finalise(); filt.buildRange(8, 10);

    int hit_real = 0, hit_absent = 0;
    const std::string absent = randomSeq(rng, 900);
    for (int t = 0; t < 1000; ++t)
    {
      const size_t o1 = rng() % (target[0].size() - 9);
      if (filt.matches(target[0].substr(o1, 9))) ++hit_real;
      const size_t o2 = rng() % (absent.size() - 9);
      if (filt.matches(absent.substr(o2, 9))) ++hit_absent;
    }
    CHECK(hit_real == 1000, "present-protein 9-mers matched %d/1000", hit_real);
    CHECK(hit_absent <= 5, "absent-protein 9-mers matched %d/1000 (expected ~0)", hit_absent);
    std::printf("3. specificity: present 9-mers %d/1000, absent-protein 9-mers %d/1000\n",
                hit_real, hit_absent);
  }

  // ------------------------------------------------------ monotonic min-len
  {
    std::vector<std::string> seqs{randomSeq(rng, 2000)};
    writeFasta(F, seqs);
    std::vector<std::string> tags;
    for (int i = 0; i < 4000; ++i) tags.push_back(randomSeq(rng, 4 + static_cast<int>(rng() % 7)));
    size_t prev = SIZE_MAX;
    bool mono = true;
    for (int floor_ = 4; floor_ <= 10; ++floor_)
    {
      fastag::FastaFilter f(true);
      f.load(F); f.setMinLen(floor_); f.finalise(); f.buildRange(4, 10);
      size_t n = 0;
      for (const auto& t : tags) if (f.matches(t)) ++n;
      if (n > prev) mono = false;
      prev = n;
    }
    CHECK(mono, "match count must fall monotonically as --min-filter-len rises");
    std::printf("4. |matches| falls monotonically as the length floor rises\n");
  }

  // --------------------------------------------------------- auto floor
  {
    // Single small protein vs a large database must get different floors.
    writeFasta(F, {randomSeq(rng, 392)});
    fastag::FastaFilter a(true); a.load(F); a.finalise();
    std::vector<std::string> big;
    for (int i = 0; i < 400; ++i) big.push_back(randomSeq(rng, 2500));   // 1M residues
    writeFasta(F, big);
    fastag::FastaFilter b(true); b.load(F); b.finalise();
    CHECK(a.minLen() < b.minLen(), "auto floor must grow with database size (%d vs %d)",
          a.minLen(), b.minLen());
    CHECK(a.chanceRate(a.minLen()) < 0.05, "auto floor must hold chance rate under 5%% (got %.3f)",
          a.chanceRate(a.minLen()));
    CHECK(b.chanceRate(b.minLen()) < 0.05, "auto floor must hold chance rate under 5%% (got %.3f)",
          b.chanceRate(b.minLen()));
    std::printf("5. auto floor: 392aa -> %d, 1M aa -> %d; chance rate < 5%% in both\n",
                a.minLen(), b.minLen());
  }

  // ------------------------------------------------------ orientation semantics
  {
    std::vector<std::string> seqs{randomSeq(rng, 1200)};
    writeFasta(F, seqs);
    fastag::FastaFilter both(true), fwd(false);
    both.load(F); both.setMinLen(9); both.finalise(); both.buildRange(9, 9);
    fwd.load(F);  fwd.setMinLen(9);  fwd.finalise();  fwd.buildRange(9, 9);

    int f_only = 0, r_only = 0, rev_missed_by_fwd = 0;
    for (int t = 0; t < 600; ++t)
    {
      const size_t off = rng() % (seqs[0].size() - 9);
      const std::string sub = seqs[0].substr(off, 9);
      std::string rev(sub.rbegin(), sub.rend());
      if (both.match(sub) == fastag::FastaFilter::Hit::Forward) ++f_only;
      if (both.match(rev) == fastag::FastaFilter::Hit::Reverse) ++r_only;
      if (!fwd.matches(rev)) ++rev_missed_by_fwd;
    }
    CHECK(f_only == 600, "forward substrings must report Hit::Forward (%d/600)", f_only);
    CHECK(r_only == 600, "reversed substrings must report Hit::Reverse (%d/600)", r_only);
    CHECK(rev_missed_by_fwd == 600, "orientation=forward must reject reversed (%d/600)", rev_missed_by_fwd);
    std::printf("6. orientation: fwd/rev correctly distinguished; forward-only rejects all 600 reversals\n");
  }

  // ------------------------------------------------------------ degenerate input
  {
    fastag::FastaFilter f(true);
    std::string err;
    // empty file
    { FILE* h = fopen(F.c_str(), "w"); fclose(h); }
    CHECK(!f.load(F, &err), "empty FASTA must fail to load");
    // header only
    { FILE* h = fopen(F.c_str(), "w"); fprintf(h, ">only a header\n"); fclose(h); }
    fastag::FastaFilter g(true);
    CHECK(!g.load(F, &err), "header-only FASTA must fail to load");
    // lowercase + CRLF + wrapping + a sequence shorter than the tag
    { FILE* h = fopen(F.c_str(), "w");
      fprintf(h, ">a\r\nggaassppvv\r\nttccllnndd\r\n>b\r\nGA\r\n"); fclose(h); }
    fastag::FastaFilter k(true);
    CHECK(k.load(F, &err), "lowercase/CRLF/wrapped FASTA must load: %s", err.c_str());
    k.setMinLen(4); k.finalise(); k.buildRange(4, 8);
    CHECK(k.matches("GGAA"), "lowercase must be upcased and wrapped lines joined");
    CHECK(k.matches("VVTT"), "line-wrap boundary must be joined, not broken");
    CHECK(!k.matches("QQQQ"), "absent 4-mer must not match");
    // ambiguity codes must not become matchable residues
    { FILE* h = fopen(F.c_str(), "w"); fprintf(h, ">x\nGGXXAAAA\n"); fclose(h); }
    fastag::FastaFilter x(true);
    x.load(F); x.setMinLen(4); x.finalise(); x.buildRange(4, 8);
    CHECK(!x.matches("GGXX"), "a tag spelling X must not match an indexed X");
    CHECK(x.matches("AAAA"), "unambiguous region must still match");
    std::printf("7. degenerate FASTA: empty/header-only rejected; lowercase, CRLF,\n"
                "   line-wrap joined; ambiguity codes not matchable\n");
  }

  // ------------------------------------------------------- TL=3 is a no-op
  {
    // Documented behaviour, asserted rather than left implicit: at length 3 a
    // filter against any reasonably sized sequence set accepts nearly everything.
    std::vector<std::string> big;
    for (int i = 0; i < 200; ++i) big.push_back(randomSeq(rng, 2500));
    writeFasta(F, big);
    fastag::FastaFilter f(true);
    f.load(F); f.setMinLen(3); f.finalise(); f.buildRange(3, 3);
    int hit = 0;
    for (int t = 0; t < 2000; ++t) if (f.matches(randomSeq(rng, 3))) ++hit;
    CHECK(hit > 1900, "at length 3 the filter should accept ~everything (got %d/2000)", hit);
    CHECK(f.autoMinLen() > 3, "auto floor must never be 3 for a database this size");
    std::printf("8. TL=3 no-op confirmed: %d/2000 random 3-mers accepted; auto floor is %d\n",
                hit, f.autoMinLen());
  }

  std::printf(failures ? "\n%d FAILURES\n" : "\nall checks passed\n", failures);
  return failures ? 1 : 0;
}
