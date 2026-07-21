// Restrict reported tags to those occurring in supplied sequences.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#pragma once

#include <OpenMS/FORMAT/FASTAFile.h>

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace FasTag
{
  /// Residues per key at 5 bits each in a 128-bit word.
  ///
  /// A uint64 holds only 12 residues (5*13 = 65 > 64) while the realised tag
  /// length reaches seed + 2*extension = 15 at the documented budget, so 64-bit
  /// packing would silently collide on exactly the long tags that carry the most
  /// information. Over-length tags are rejected, never truncated.
  constexpr int MAX_FILTER_LEN = 25;

  /// Sentinel for ambiguity codes; no tag can contain it.
  constexpr char AMBIG = '#';

  /// A 125-bit k-mer key as two 64-bit halves.
  ///
  /// Was __uint128_t, which is a GCC/Clang extension MSVC does not provide -- so
  /// the whole tool failed to compile there, for one type in one file. Two
  /// uint64_t are portable, the same 16 bytes, and the only operations needed
  /// are shift-in-from-the-right and equality.
  struct Kmer128
  {
    uint64_t hi = 0, lo = 0;

    /// Shift left 5 bits and insert a residue code in the low bits.
    void push5(uint64_t code) noexcept
    {
      hi = (hi << 5) | (lo >> 59);   // the 5 bits leaving lo enter hi
      lo = (lo << 5) | (code & 0x1full);
    }
    static Kmer128 invalid() noexcept { return {~0ull, ~0ull}; }

    bool operator==(const Kmer128& o) const noexcept { return hi == o.hi && lo == o.lo; }
    bool operator!=(const Kmer128& o) const noexcept { return !(*this == o); }
  };

  struct Kmer128Hash
  {
    size_t operator()(const Kmer128& v) const noexcept
    {
      auto mix = [](uint64_t x) {
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
      };
      return static_cast<size_t>(mix(v.lo) ^ (mix(v.hi) * 3));
    }
  };

  class FastaFilter
  {
  public:
    enum class Hit { None, Forward, Reverse };

    /// @param both_orientations also match the reversed tag.
    ///   Tags are stored N->C under a y-ion assumption, so b-derived tags (~57%
    ///   of output) read reversed relative to the sequence. Forward-only would
    ///   discard most of them. Costs a factor of two in chance rate -- 0.26
    ///   residues of effective length -- and buys 2.3x recall.
    explicit FastaFilter(bool both_orientations = true) : both_(both_orientations) {}

    /// Adopt sequences parsed by OpenMS.
    ///
    /// I is folded to L unconditionally and there is deliberately no flag: the
    /// tagger traverses on canonical L because I and L are isobaric, so an
    /// emitted L carries no information about which residue it was. A "literal"
    /// match would not be strict, it would silently reject every match landing on
    /// an isoleucine.
    bool load(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries, std::string* err = nullptr);

    /// Interchangeable residue/pair rules, derived from the residue masses at
    /// @p tol rather than hardcoded.
    ///
    /// A tag reading N where the sequence spells GG is a mass-correct reading of
    /// the true peptide that exact string matching rejects -- a recall loss that
    /// is sequence-correlated, not random. At 0.04 Da (20 ppm over two peaks at
    /// m/z 1000) the derived set is N=GG, Q=GA/AG, R=GV/VG, W=AD/GE/SV, K=GA.
    /// No residue equals a sum of three, so one level of collapse suffices.
    void deriveCollapses(double tol);
    size_t collapseRules() const { return rules_.size(); }

    /// Below this length the filter is close to the identity function, so it is
    /// derived from database size unless set explicitly.
    void setMinLen(int n) { min_len_ = n; min_len_auto_ = false; }
    int  minLen() const { return min_len_ > 0 ? min_len_ : autoMinLen(); }
    bool minLenAuto() const { return min_len_auto_; }

    /// Smallest length whose expected chance-match rate is below 5%:
    /// k = ceil(log_A(20*c*M)), with A the effective alphabet and c=2 for both
    /// orientations. Gives 4 for a 392-residue protein, 8 for the human proteome.
    int autoMinLen() const;

    /// Expected fraction of random length-k tags that match, so an observed count
    /// is never mistaken for signal.
    double chanceRate(int k) const;

    /// Build indices for every length that can occur. Call once, after load()
    /// and any deriveCollapses(); the filter is immutable and thread-safe after.
    ///
    /// A separate full index is materialised per length, so the footprint scales
    /// with (database residues x number of lengths) -- roughly 2.4 GB for a 3 M
    /// residue database at seed 6 with extension 6, and far more for a proteome.
    /// Above @p max_bytes this THROWS rather than returning a status: every
    /// caller would otherwise have to remember to check, and a filter that
    /// silently built nothing would reject 100% of tags while looking healthy.
    ///
    /// @throws OpenMS::Exception::InvalidValue if the projection exceeds max_bytes
    void build(int min_k, int max_k, size_t max_bytes = size_t(4) << 30);

    /// Projected index footprint in bytes for the given length range.
    size_t projectBytes(int min_k, int max_k) const;

    /// Keys actually indexed, valid after build().
    size_t indexedKeys() const;

    Hit match(const std::string& tag) const;

    size_t sequenceCount() const { return seqs_.size(); }
    size_t residueCount() const { return residues_; }

  private:
    /// Amino acids are far from uniform: with I/L folded, sum(f^2) = 0.068 over
    /// the human proteome, so the effective alphabet is 1/0.068 = 14.7, not 19.
    /// Using 19 understates the chance rate and picks a floor one residue short.
    static constexpr double EFF_ALPHABET = 14.7;

    struct Collapse { char a, b, one; };

    static Kmer128 encode(const char* s, int n);
    bool contains(const std::string& t) const;
    void emitReadings(const std::string& seq, size_t i, int k, std::string& cur,
                      std::unordered_set<Kmer128, Kmer128Hash>& out, size_t& budget) const;

    bool both_;
    int  min_len_ = 0;
    bool min_len_auto_ = true;
    size_t residues_ = 0;
    std::vector<std::string> seqs_;
    std::vector<Collapse> rules_;
    std::vector<std::unordered_set<Kmer128, Kmer128Hash>> sets_;  ///< indexed by length
    std::vector<char> built_;
  };
}
