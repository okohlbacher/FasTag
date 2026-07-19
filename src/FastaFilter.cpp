// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "FastaFilter.h"

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace OpenMS;

namespace FasTag
{
  namespace
  {
    /// Upper-case, fold I to L, mark ambiguity codes. Returns 0 to skip.
    ///
    /// Ambiguity codes must not become matchable residues: an indexed X would
    /// match a tag spelling X, and B/Z/J each stand for two residues.
    inline char norm(char c)
    {
      if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
      if (c < 'A' || c > 'Z') return 0;
      if (c == 'X' || c == 'B' || c == 'Z' || c == 'J' || c == 'U' || c == 'O') return AMBIG;
      return c == 'I' ? 'L' : c;
    }
  }

  __uint128_t FastaFilter::encode(const char* s, int n)
  {
    __uint128_t v = 0;
    for (int i = 0; i < n; ++i)
    {
      const char c = s[i];
      if (c < 'A' || c > 'Z') return ~static_cast<__uint128_t>(0);
      v = (v << 5) | static_cast<__uint128_t>(c - 'A');
    }
    return v;
  }

  bool FastaFilter::load(const std::vector<FASTAFile::FASTAEntry>& entries, std::string* err)
  {
    for (const auto& e : entries)
    {
      std::string cur;
      cur.reserve(e.sequence.size());
      for (char c : e.sequence) { const char n = norm(c); if (n) cur.push_back(n); }
      if (cur.empty()) continue;
      residues_ += cur.size();
      seqs_.push_back(std::move(cur));
    }
    if (seqs_.empty()) { if (err) *err = "no usable sequence"; return false; }
    return true;
  }

  void FastaFilter::deriveCollapses(double tol)
  {
    rules_.clear();
    std::vector<std::pair<char, double>> R;
    for (const Residue* r : ResidueDB::getInstance()->getResidues("Natural19WithoutI"))
      R.emplace_back(r->getOneLetterCode()[0], r->getMonoWeight(Residue::Internal));

    for (const auto& one : R)
      for (const auto& a : R)
        for (const auto& b : R)
          if (std::fabs(one.second - (a.second + b.second)) <= tol)
            rules_.push_back({a.first, b.first, one.first});
  }

  int FastaFilter::autoMinLen() const
  {
    const double c = both_ ? 2.0 : 1.0;
    if (residues_ == 0) return 1;
    const int k = static_cast<int>(
        std::ceil(std::log(20.0 * c * static_cast<double>(residues_)) / std::log(EFF_ALPHABET)));
    return std::max(1, std::min(k, MAX_FILTER_LEN));
  }

  double FastaFilter::chanceRate(int k) const
  {
    if (k < 1 || k > MAX_FILTER_LEN) return 0.0;
    const double c = both_ ? 2.0 : 1.0;
    double space = 1.0;
    for (int i = 0; i < k; ++i) space *= EFF_ALPHABET;
    return 1.0 - std::exp(-c * static_cast<double>(residues_) / space);
  }

  /// Emit every length-k reading from position i, applying collapses. Depth is k
  /// and branching at most 2 per step, but real sequences give ~1; capped anyway.
  void FastaFilter::emitReadings(const std::string& seq, size_t i, int k, std::string& cur,
                                 std::unordered_set<__uint128_t, Kmer128Hash>& out,
                                 size_t& budget) const
  {
    if (static_cast<int>(cur.size()) == k) { out.insert(encode(cur.data(), k)); return; }
    if (i >= seq.size() || budget == 0 || seq[i] == AMBIG) return;

    cur.push_back(seq[i]);
    emitReadings(seq, i + 1, k, cur, out, budget);
    cur.pop_back();

    if (i + 1 < seq.size() && seq[i + 1] != AMBIG)
      for (const auto& r : rules_)
        if (r.a == seq[i] && r.b == seq[i + 1])
        {
          if (budget == 0) break;
          --budget;
          cur.push_back(r.one);
          emitReadings(seq, i + 2, k, cur, out, budget);
          cur.pop_back();
        }
  }

  size_t FastaFilter::projectBytes(int min_k, int max_k) const
  {
    // Distinct k-mers are bounded by both the residue count and the alphabet
    // space; short k saturates the alphabet, long k saturates the sequence.
    // 48 bytes/key is measured steady-state for unordered_set<__uint128_t>.
    min_k = std::max(min_k, minLen());
    max_k = std::min(max_k, MAX_FILTER_LEN);
    double total = 0;
    for (int k = min_k; k <= max_k; ++k)
    {
      double space = 1.0;
      for (int i = 0; i < k && space < 4e18; ++i) space *= 19.0;
      const double n = std::min(static_cast<double>(residues_), space);
      total += n * (rules_.empty() ? 1.0 : 1.6);   // collapse readings inflate the set
    }
    return static_cast<size_t>(total * 48.0);
  }

  size_t FastaFilter::indexedKeys() const
  {
    size_t n = 0;
    for (const auto& s : sets_) n += s.size();
    return n;
  }

  void FastaFilter::build(int min_k, int max_k, size_t max_bytes)
  {
    const size_t projected = projectBytes(min_k, max_k);
    if (projected > max_bytes)
    {
      throw Exception::InvalidValue(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "FASTA index would need about " + String(projected / (1024 * 1024)) +
          " MB, above the " + String(max_bytes / (1024 * 1024)) + " MB limit. "
          "Reduce 'extension', raise 'min_filter_length', or use a smaller database.",
          String(projected));
    }
    const int floor_k = minLen();
    min_k = std::max(min_k, floor_k);
    max_k = std::min(max_k, MAX_FILTER_LEN);
    if (max_k < min_k)
    {
      // No tag this run can produce is long enough to be worth matching, so the
      // filter would reject every tag while reporting "0 hits" as if that were a
      // result. Refuse loudly instead: silence here is indistinguishable from a
      // genuine absence of evidence.
      throw Exception::InvalidValue(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "No reportable tag can reach the minimum filter length. This database of " +
          String(residues_) + " residues needs tags of at least " + String(floor_k) +
          " residues to beat chance, but the longest tag this run can produce is " +
          String(max_k) + ". Raise 'tag_length' or 'extension', or set "
          "'min_filter_length' explicitly to accept a higher chance-match rate.",
          String(floor_k));
    }
    sets_.resize(static_cast<size_t>(max_k) + 1);
    built_.resize(static_cast<size_t>(max_k) + 1, 0);

    for (int k = min_k; k <= max_k; ++k)
    {
      auto& s = sets_[static_cast<size_t>(k)];
      for (const auto& seq : seqs_)
      {
        if (static_cast<int>(seq.size()) < k) continue;
        for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size(); ++i)
        {
          if (rules_.empty())
          {
            if (std::memchr(seq.data() + i, AMBIG, static_cast<size_t>(k))) continue;
            s.insert(encode(seq.data() + i, k));
          }
          else
          {
            std::string cur;
            cur.reserve(static_cast<size_t>(k));
            size_t budget = 64;
            emitReadings(seq, i, k, cur, s, budget);
          }
        }
      }
      built_[static_cast<size_t>(k)] = 1;
    }
  }

  bool FastaFilter::contains(const std::string& t) const
  {
    const int k = static_cast<int>(t.size());
    if (k < 1 || k > MAX_FILTER_LEN) return false;
    if (static_cast<size_t>(k) >= sets_.size() || !built_[static_cast<size_t>(k)]) return false;
    const __uint128_t e = encode(t.data(), k);
    if (e == ~static_cast<__uint128_t>(0)) return false;
    return sets_[static_cast<size_t>(k)].count(e) > 0;
  }

  FastaFilter::Hit FastaFilter::match(const std::string& tag) const
  {
    if (static_cast<int>(tag.size()) < minLen()) return Hit::None;
    std::string t = tag;
    for (char& c : t) if (c == 'I') c = 'L';
    if (contains(t)) return Hit::Forward;
    if (both_)
    {
      std::reverse(t.begin(), t.end());
      if (contains(t)) return Hit::Reverse;
    }
    return Hit::None;
  }
}
