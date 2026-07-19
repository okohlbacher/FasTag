// FasTag -- restrict reported tags to exact substrings of supplied sequences.
//
// Encoding is __uint128_t at 5 bits/residue, holding 25 residues. A uint64_t
// holds only 12 (5*13 = 65 > 64), while the realised tag length with a seed of 3
// and 6 residues of extension per terminus reaches 15 -- so a 64-bit packing
// silently collides on exactly the long tags that carry the most information.
// Lengths beyond the encoder's capacity are rejected, never truncated.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <FasTag/fastag.h>   // residues()

namespace fastag
{
  /// Residues per key at 5 bits each in a 128-bit word.
  constexpr int MAX_FILTER_LEN = 25;

  /// Sentinel for ambiguity codes; no tag can contain it.
  constexpr char AMBIG = '#';

  struct Kmer128Hash
  {
    size_t operator()(__uint128_t v) const noexcept
    {
      // splitmix64 over both halves; std::hash has no __uint128_t specialisation.
      auto mix = [](uint64_t x) {
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
      };
      return static_cast<size_t>(mix(static_cast<uint64_t>(v)) ^
                                 (mix(static_cast<uint64_t>(v >> 64)) * 3));
    }
  };

  class FastaFilter
  {
  public:
    struct Stats
    {
      size_t sequences = 0;
      size_t residues = 0;
      int    min_len = 0;          ///< effective floor (auto or explicit)
      bool   min_len_auto = true;
      bool   il_equiv = true;
      bool   both_orientations = true;
    };

    /// I/L folding is UNCONDITIONAL and there is deliberately no flag for it.
    /// The tagger traverses on canonical 'L' (I and L are isobaric) and
    /// `expandIL` is never called in the pipeline, so an emitted 'L' carries zero
    /// information about which residue it was. A "literal" match would therefore
    /// not be strict -- it would silently reject every match landing on an Ile,
    /// biasing against Ile-rich proteins.
    ///
    /// @param both_orientations also match the reversed tag. FasTag stores tags
    ///   N->C under a y-ion assumption; b-derived tags (~57% of output) appear
    ///   reversed relative to the protein, so forward-only discards most of them.
    ///   Costs a factor of 2 in chance rate = 0.26 residues of length; buys 2.33x
    ///   recall. Unambiguously worth it.
    explicit FastaFilter(bool both_orientations = true)
    {
      st_.il_equiv = true;
      st_.both_orientations = both_orientations;
    }

    /// Normalise one residue for indexing: upper-case, I folded to L, ambiguity
    /// codes replaced by a sentinel no tag can contain. Returns 0 to skip.
    static char normResidue(char c)
    {
      if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
      if (c < 'A' || c > 'Z') return 0;                 // digits, gaps, '*'
      // Ambiguity codes must not become matchable residues: an indexed 'X' would
      // match a tag spelling X, and B/Z/J each stand for two residues.
      if (c == 'X' || c == 'B' || c == 'Z' || c == 'J' || c == 'U' || c == 'O') return AMBIG;
      return c == 'I' ? 'L' : c;
    }

    /// Adopt sequences that a caller already parsed (e.g. OpenMS FASTAFile), so
    /// the tool need not round-trip through a temporary file.
    bool loadSequences(const std::vector<std::string>& seqs, std::string* err = nullptr)
    {
      for (const auto& raw : seqs)
      {
        std::string cur;
        cur.reserve(raw.size());
        for (char c : raw) { const char n = normResidue(c); if (n) cur.push_back(n); }
        if (cur.empty()) continue;
        st_.residues += cur.size();
        ++st_.sequences;
        seqs_.push_back(std::move(cur));
      }
      if (seqs_.empty()) { if (err) *err = "no usable sequence"; return false; }
      return true;
    }

    /// Parse FASTA. Headers are excluded and wrapped lines joined -- matching the
    /// raw bytes would both match header text and miss sequences split by newline.
    bool load(const std::string& path, std::string* err = nullptr)
    {
      std::ifstream f(path);
      if (!f) { if (err) *err = "cannot open FASTA: " + path; return false; }
      std::string line, cur;
      auto flush = [&]() {
        if (cur.empty()) return;
        seqs_.push_back(cur);
        st_.residues += cur.size();
        ++st_.sequences;
        cur.clear();
      };
      while (std::getline(f, line))
      {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') { flush(); continue; }
        for (char c : line) { const char n = normResidue(c); if (n) cur.push_back(n); }
      }
      flush();
      if (seqs_.empty()) { if (err) *err = "FASTA contains no sequence: " + path; return false; }
      return true;
    }

    /// Effective alphabet size, NOT 19.
    ///
    /// Amino acids are far from uniform: with I/L folded, sum(f^2) = 0.068 over
    /// the human proteome, so the effective alphabet is 1/0.068 = 14.7. Using 19
    /// understates the chance-match rate and picks a floor one residue too low
    /// (7 instead of 8 for a whole proteome).
    static constexpr double EFF_ALPHABET = 14.7;

    /// Smallest length whose expected chance-match rate is below 5%.
    ///
    /// P(random length-k tag matches) = 1 - exp(-c*M / A^k), c = 2 for both
    /// orientations, A = effective alphabet. Solving for 5%:
    /// k = ceil(log_A(20*c*M)). Gives 4 for a single 392-residue protein, 8 for
    /// the human proteome. Without a floor the filter is the IDENTITY FUNCTION at
    /// k=3 -- measured: 99.5% of spectra retained against a single 400-aa target.
    int autoMinLen() const
    {
      const double c = st_.both_orientations ? 2.0 : 1.0;
      const double m = static_cast<double>(st_.residues);
      if (m <= 0) return 1;
      const int k = static_cast<int>(
          std::ceil(std::log(20.0 * c * m) / std::log(EFF_ALPHABET)));
      return std::max(1, std::min(k, MAX_FILTER_LEN));
    }

    /// Expected fraction of random length-k tags that match. For reporting the
    /// chance rate alongside the observed one, so the two are never conflated.
    double chanceRate(int k) const
    {
      if (k < 1 || k > MAX_FILTER_LEN) return 0.0;
      const double c = st_.both_orientations ? 2.0 : 1.0;
      double space = 1.0;
      for (int i = 0; i < k; ++i) space *= EFF_ALPHABET;
      return 1.0 - std::exp(-c * static_cast<double>(st_.residues) / space);
    }

    /// Isobaric collapse rules, DERIVED from the residue table rather than
    /// hardcoded: any residue whose mass equals the sum of two residues within
    /// @p tol is a pair the tagger can read as one.
    ///
    /// Measured at tol = 0.04 Da (20 ppm on two peaks at m/z 1000):
    ///   N = GG, Q = GA/AG      (0.01 mDa -- truly isobaric, always confusable)
    ///   R = GV/VG              (11.2 mDa)
    ///   W = AD/DA/GE/EG/SV/VS  (15-21 mDa)
    ///   K = GA/AG              (36.4 mDa, only at high m/z)
    /// No residue equals a sum of three, so one level of collapse suffices.
    ///
    /// Without this, a tag reading `N` where the protein spells `GG` is a
    /// mass-correct reading of the true peptide that exact string matching
    /// rejects -- a recall loss that is sequence-correlated, not random.
    void deriveCollapses(double tol)
    {
      collapse_tol_ = tol;
      rules_.clear();
      const auto& R = residues();
      for (const auto& one : R)
      {
        const char c1 = one.code == 'I' ? 'L' : one.code;
        for (const auto& a : R)
          for (const auto& b : R)
            if (std::fabs(one.mass - (a.mass + b.mass)) <= tol)
              rules_.push_back({a.code == 'I' ? 'L' : a.code,
                                b.code == 'I' ? 'L' : b.code, c1});
      }
    }

    size_t collapseRuleCount() const { return rules_.size(); }

    /// Emit every length-k reading starting at @p i, applying collapses.
    /// Depth is k and branching is at most 2 per step, but real sequences give
    /// ~1, so the realised fan-out is small. Capped defensively.
    void emitReadings_(const std::string& seq, size_t i, int k,
                       std::string& cur, std::unordered_set<__uint128_t, Kmer128Hash>& out,
                       size_t& budget) const
    {
      if (static_cast<int>(cur.size()) == k)
      { out.insert(encode(cur.data(), k)); return; }
      if (i >= seq.size() || budget == 0) return;
      if (seq[i] == AMBIG) return;

      // literal residue
      cur.push_back(seq[i]);
      emitReadings_(seq, i + 1, k, cur, out, budget);
      cur.pop_back();

      // collapsed pair -> single residue
      if (i + 1 < seq.size() && seq[i + 1] != AMBIG)
      {
        for (const auto& r : rules_)
          if (r.a == seq[i] && r.b == seq[i + 1])
          {
            if (budget == 0) break;
            --budget;
            cur.push_back(r.one);
            emitReadings_(seq, i + 2, k, cur, out, budget);
            cur.pop_back();
          }
      }
    }

    /// Build the k-mer index for one length. Lazy: only lengths above the floor
    /// that actually occur are ever built. Indexes every isobaric *reading* of
    /// each window, not just its literal spelling.
    void build(int k)
    {
      if (k < 1 || k > MAX_FILTER_LEN) return;
      if (static_cast<size_t>(k) < sets_.size() && built_[k]) return;
      if (sets_.size() <= static_cast<size_t>(k)) { sets_.resize(k + 1); built_.resize(k + 1, 0); }
      auto& s = sets_[k];
      for (const auto& seq : seqs_)
      {
        if (static_cast<int>(seq.size()) < k) continue;
        for (size_t i = 0; i + k <= seq.size(); ++i)
        {
          if (rules_.empty())
          {
            if (memchr(seq.data() + i, AMBIG, static_cast<size_t>(k))) continue;
            s.insert(encode(seq.data() + i, k));
          }
          else
          {
            std::string cur;
            cur.reserve(static_cast<size_t>(k));
            size_t budget = 64;               // per start position
            emitReadings_(seq, i, k, cur, s, budget);
          }
        }
      }
      built_[k] = true;
    }

    void setMinLen(int n) { st_.min_len = n; st_.min_len_auto = false; }
    void finalise() { if (st_.min_len_auto) st_.min_len = autoMinLen(); }
    const Stats& stats() const { return st_; }
    int minLen() const { return st_.min_len; }

    /// Is @p tag an exact substring of any sequence, under the active rules?
    /// Tags shorter than the floor return false: they cannot carry evidence, and
    /// reporting them as matches is precisely the false-positive generator.
    enum class Hit { None = 0, Forward = 1, Reverse = 2 };

    /// Which orientation matched -- not merely whether one did.
    ///
    /// A reverse-only match is diagnostic: FasTag stores tags N->C under a y-ion
    /// assumption, so a tag that matches only reversed was read off the b series,
    /// and its reported flanking masses carry the known one-water offset. That is
    /// free information the caller should not have to rediscover.
    Hit match(const std::string& tag) const
    {
      const int k = static_cast<int>(tag.size());
      if (k < st_.min_len || k > MAX_FILTER_LEN) return Hit::None;
      if (static_cast<size_t>(k) >= sets_.size() || !built_[k]) return Hit::None;
      std::string t = tag;
      for (char& c : t) if (c == 'I') c = 'L';   // tags always carry canonical L
      const auto& s = sets_[k];
      if (s.count(encode(t.data(), k))) return Hit::Forward;
      if (st_.both_orientations)
      {
        std::string r(t.rbegin(), t.rend());
        if (s.count(encode(r.data(), k))) return Hit::Reverse;
      }
      return Hit::None;
    }

    bool matches(const std::string& tag) const { return match(tag) != Hit::None; }

    /// Which lengths need an index, given the tag lengths a run can produce.
    void buildRange(int lo, int hi)
    {
      finalise();
      for (int k = std::max(lo, st_.min_len); k <= std::min(hi, MAX_FILTER_LEN); ++k) build(k);
    }

    static __uint128_t encode(const char* s, int k)
    {
      __uint128_t v = 0;
      for (int i = 0; i < k; ++i)
      {
        char c = s[i];
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        const unsigned d = (c >= 'A' && c <= 'Z') ? static_cast<unsigned>(c - 'A') : 26u;
        v = (v << 5) | d;
      }
      return v;
    }

    const std::vector<std::string>& sequences() const { return seqs_; }

  private:
    struct Collapse { char a, b, one; };   ///< a+b reads as `one` within tolerance
    std::vector<Collapse> rules_;
    double collapse_tol_ = 0.0;
    std::vector<std::string> seqs_;
    std::vector<std::unordered_set<__uint128_t, Kmer128Hash>> sets_;
    std::vector<char> built_;   ///< not vector<bool>: bit-packing buys nothing for a few entries
    Stats st_;
  };
}
