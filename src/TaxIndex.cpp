// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TaxIndex.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace FasTag
{
  const std::vector<uint32_t> TaxIndex::empty_;

  std::string TaxIndex::fold(const std::string& s)
  {
    std::string out = s;
    for (char& c : out) if (c == 'I') c = 'L';
    return out;
  }

  void TaxIndex::build(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries,
                       const std::vector<uint32_t>& taxids, int k)
  {
    k_ = k;
    index_.clear();
    postings_.clear();
    if (k <= 0) return;

    for (size_t e = 0; e < entries.size(); ++e)
    {
      const uint32_t tax = e < taxids.size() ? taxids[e] : 0;
      if (tax == 0) continue;
      const std::string seq = fold(entries[e].sequence);
      if (static_cast<int>(seq.size()) < k) continue;

      // Distinct k-mers of THIS protein, so one long protein cannot inflate a
      // taxon's presence for a k-mer it only carries once. Collected per protein
      // then merged into the taxon set below.
      for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size(); ++i)
      {
        // A k-mer containing a non-standard residue (X/B/Z/U/*) is ambiguous and
        // would match spuriously; skip it.
        bool ok = true;
        for (int j = 0; j < k; ++j)
        {
          const char c = seq[i + static_cast<size_t>(j)];
          if (c < 'A' || c > 'Y' || c == 'X' || c == 'B' || c == 'J' || c == 'O' ||
              c == 'U' || c == 'Z') { ok = false; break; }
        }
        if (!ok) continue;

        std::vector<uint32_t>& set = index_[seq.substr(i, static_cast<size_t>(k))];
        // Insert tax keeping the set sorted and unique; postings count only the
        // first time a taxon is added to a given k-mer.
        auto it = std::lower_bound(set.begin(), set.end(), tax);
        if (it == set.end() || *it != tax)
        {
          set.insert(it, tax);
          ++postings_[tax];
        }
      }
    }
  }

  const std::vector<uint32_t>& TaxIndex::lookup(const std::string& kmer) const
  {
    auto it = index_.find(kmer);
    return it == index_.end() ? empty_ : it->second;
  }

  uint64_t TaxIndex::postings(uint32_t taxid) const
  {
    auto it = postings_.find(taxid);
    return it == postings_.end() ? 0 : it->second;
  }

  std::vector<uint32_t> TaxIndex::taxa() const
  {
    std::vector<uint32_t> out;
    out.reserve(postings_.size());
    for (const auto& kv : postings_) out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
  }

  // --- Serialization -------------------------------------------------------
  //
  // Little-endian, fixed-width. Layout:
  //   magic "FTXI"  (4 bytes)
  //   version u32   (=1)
  //   k u32
  //   n_kmers u64
  //   repeated n_kmers times:
  //     klen u8 | k bytes of the k-mer | setlen u32 | setlen * u32 taxids
  // postings_ is recomputed on load rather than stored (cheap, avoids drift).
  namespace
  {
    template <class T> void wr(std::ostream& o, T v) { o.write(reinterpret_cast<const char*>(&v), sizeof v); }
    template <class T> bool rd(std::istream& i, T& v) { return static_cast<bool>(i.read(reinterpret_cast<char*>(&v), sizeof v)); }
  }

  bool TaxIndex::save(const std::string& path, std::string* err) const
  {
    std::ofstream o(path, std::ios::binary);
    if (!o) { if (err) *err = "cannot open " + path + " for writing"; return false; }
    o.write("FTXI", 4);
    wr<uint32_t>(o, 1);
    wr<uint32_t>(o, static_cast<uint32_t>(k_));
    wr<uint64_t>(o, static_cast<uint64_t>(index_.size()));
    for (const auto& kv : index_)
    {
      wr<uint8_t>(o, static_cast<uint8_t>(kv.first.size()));
      o.write(kv.first.data(), static_cast<std::streamsize>(kv.first.size()));
      wr<uint32_t>(o, static_cast<uint32_t>(kv.second.size()));
      for (uint32_t t : kv.second) wr<uint32_t>(o, t);
    }
    return static_cast<bool>(o);
  }

  bool TaxIndex::load(const std::string& path, std::string* err)
  {
    std::ifstream i(path, std::ios::binary);
    if (!i) { if (err) *err = "cannot open " + path; return false; }
    char magic[4];
    if (!i.read(magic, 4) || std::memcmp(magic, "FTXI", 4) != 0)
    { if (err) *err = "not a FasTag taxonomy index"; return false; }
    uint32_t ver = 0, k = 0; uint64_t n = 0;
    if (!rd(i, ver) || ver != 1) { if (err) *err = "unsupported index version"; return false; }
    if (!rd(i, k) || !rd(i, n)) { if (err) *err = "truncated index header"; return false; }
    k_ = static_cast<int>(k);
    index_.clear(); postings_.clear();
    index_.reserve(n);
    for (uint64_t e = 0; e < n; ++e)
    {
      uint8_t klen = 0;
      if (!rd(i, klen)) { if (err) *err = "truncated index"; return false; }
      std::string key(klen, '\0');
      if (!i.read(&key[0], klen)) { if (err) *err = "truncated index key"; return false; }
      uint32_t setlen = 0;
      if (!rd(i, setlen)) { if (err) *err = "truncated index set length"; return false; }
      std::vector<uint32_t> set(setlen);
      for (uint32_t s = 0; s < setlen; ++s)
        if (!rd(i, set[s])) { if (err) *err = "truncated index set"; return false; }
      for (uint32_t t : set) ++postings_[t];
      index_.emplace(std::move(key), std::move(set));
    }
    return true;
  }
}
