// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TaxIndex.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>

#ifdef _WIN32
#  include <windows.h>
#else
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <unistd.h>
#endif

namespace FasTag
{
  namespace
  {
    /// The folded alphabet, exactly 19 letters: the 20 standard residues with I
    /// folded onto L. Anything else (X/B/J/O/U/Z, digits, '*') is ambiguous and
    /// would match spuriously, so it has no code and kills the k-mer -- the same
    /// rejection v1's character loop performed, expressed once.
    constexpr char ALPHABET[] = "ACDEFGHKLMNPQRSTVWY";
    constexpr int RADIX = 19;

    const int8_t* codeTable()
    {
      // Function-local static with a lambda initializer: C++11 guarantees the
      // initialization is thread-safe and runs exactly once. The previous
      // `static bool init` guard was a data race -- a v2 load never touches this
      // table, so the FIRST lookups happen straight from the OpenMP loop, with
      // multiple threads writing `t` and `init` at once.
      static const std::array<int8_t, 256> table = [] {
        std::array<int8_t, 256> t;
        t.fill(-1);
        for (int i = 0; i < RADIX; ++i) t[static_cast<uint8_t>(ALPHABET[i])] = static_cast<int8_t>(i);
        t[static_cast<uint8_t>('I')] = t[static_cast<uint8_t>('L')];  // fold, defensively
        return t;
      }();
      return table.data();
    }

    inline uint64_t rdKey(const uint8_t* p, int width)
    {
      uint64_t v = 0;
      for (int i = width - 1; i >= 0; --i) v = (v << 8) | p[i];
      return v;
    }
    inline void wrKey(uint8_t* p, int width, uint64_t v)
    {
      for (int i = 0; i < width; ++i) { p[i] = static_cast<uint8_t>(v & 0xFF); v >>= 8; }
    }
    inline size_t align8(size_t n) { return (n + 7u) & ~static_cast<size_t>(7); }

    template <class T> void wr(std::ostream& o, T v)
    { o.write(reinterpret_cast<const char*>(&v), sizeof v); }
    template <class T> bool rd(std::istream& i, T& v)
    { return static_cast<bool>(i.read(reinterpret_cast<char*>(&v), sizeof v)); }
  }

  int TaxIndex::keyBytes(int k)
  {
    // k is capped at MAX_K, where 19^k still fits in uint64 with room to spare.
    // Beyond that the code overflows (silent collisions) and the width formula
    // truncates -- both corrupt the index. buildtaxdb enforces the same cap.
    if (k < 1) return 1;
    if (k > MAX_K) k = MAX_K;
    uint64_t max = 1;
    for (int i = 0; i < k; ++i) max *= RADIX;   // 19^MAX_K << 2^64
    int b = 1;
    while (b < 8 && ((max - 1) >> (8 * b)) != 0) ++b;
    return b;
  }

  uint64_t TaxIndex::encode(const char* s, int k)
  {
    const int8_t* t = codeTable();
    uint64_t v = 0;
    for (int i = 0; i < k; ++i)
    {
      const int8_t c = t[static_cast<uint8_t>(s[i])];
      if (c < 0) return UINT64_MAX;
      v = v * RADIX + static_cast<uint64_t>(c);
    }
    return v;
  }

  std::string TaxIndex::fold(const std::string& s)
  {
    std::string out = s;
    for (char& c : out) if (c == 'I') c = 'L';
    return out;
  }

  TaxIndex::~TaxIndex() { reset(); }

  void TaxIndex::reset()
  {
    if (map_base_ != nullptr)
    {
#ifdef _WIN32
      UnmapViewOfFile(map_base_);
      if (map_handle_ != nullptr) CloseHandle(static_cast<HANDLE>(map_handle_));
      if (file_handle_ != nullptr) CloseHandle(static_cast<HANDLE>(file_handle_));
      map_handle_ = nullptr; file_handle_ = nullptr;
#else
      munmap(map_base_, map_len_);
      if (map_fd_ >= 0) ::close(map_fd_);
      map_fd_ = -1;
#endif
      map_base_ = nullptr;
      map_len_ = 0;
    }
    keys_ = nullptr; post_ = nullptr; offsets_ = nullptr; taxa_ = nullptr; taxpost_ = nullptr;
    n_kmers_ = 0; n_taxa_ = 0; key_bytes_ = 0; tax_bytes_ = 0; k_ = 0;
    index_.clear(); postings_.clear(); legacy_ = false;
  }

  void TaxIndex::build(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries,
                       const std::vector<uint32_t>& taxids, int k)
  {
    reset();
    k_ = k;
    legacy_ = true;   // build fills the maps; save() packs them into the v2 image
    if (k < 1 || k > MAX_K) { k_ = 0; return; }

    for (size_t e = 0; e < entries.size(); ++e)
    {
      const uint32_t tax = e < taxids.size() ? taxids[e] : 0;
      if (tax == 0) continue;
      const std::string seq = fold(entries[e].sequence);
      if (static_cast<int>(seq.size()) < k) continue;

      for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size(); ++i)
      {
        if (encode(seq.data() + i, k) == UINT64_MAX) continue;

        std::vector<uint32_t>& set = index_[seq.substr(i, static_cast<size_t>(k))];
        // Keep the set sorted and unique; postings count the FIRST time a taxon
        // is added to a k-mer, so they measure breadth, not abundance.
        auto it = std::lower_bound(set.begin(), set.end(), tax);
        if (it == set.end() || *it != tax)
        {
          set.insert(it, tax);
          ++postings_[tax];
        }
      }
    }
  }

  uint64_t TaxIndex::nKmers() const
  {
    return legacy_ ? static_cast<uint64_t>(index_.size()) : n_kmers_;
  }

  std::vector<uint32_t> TaxIndex::taxa() const
  {
    if (!legacy_ && taxa_ != nullptr) return std::vector<uint32_t>(taxa_, taxa_ + n_taxa_);
    std::vector<uint32_t> out;
    out.reserve(postings_.size());
    for (const auto& kv : postings_) out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
  }

  uint64_t TaxIndex::postings(uint32_t taxid) const
  {
    if (!legacy_ && taxa_ != nullptr)
    {
      const uint32_t* e = taxa_ + n_taxa_;
      const uint32_t* it = std::lower_bound(taxa_, e, taxid);
      if (it == e || *it != taxid) return 0;
      return taxpost_[it - taxa_];
    }
    auto it = postings_.find(taxid);
    return it == postings_.end() ? 0 : it->second;
  }

  void TaxIndex::lookup(const std::string& kmer, std::vector<uint32_t>& out) const
  {
    out.clear();
    if (k_ <= 0 || static_cast<int>(kmer.size()) != k_) return;

    if (legacy_)
    {
      // fold() here too: the mapped path folds I->L via encode(), so an
      // I-containing query must behave the same before a save/load round-trip.
      auto it = index_.find(fold(kmer));
      if (it != index_.end()) out = it->second;
      return;
    }
    if (keys_ == nullptr || n_kmers_ == 0) return;

    const uint64_t want = encode(kmer.data(), k_);
    if (want == UINT64_MAX) return;

    // Binary search the packed, sorted key array in place -- no deserialization.
    uint64_t lo = 0, hi = n_kmers_;
    while (lo < hi)
    {
      const uint64_t mid = lo + (hi - lo) / 2;
      const uint64_t got = rdKey(keys_ + mid * static_cast<size_t>(key_bytes_), key_bytes_);
      if (got < want) { lo = mid + 1; continue; }
      if (got > want) { hi = mid; continue; }

      const uint32_t b = offsets_[mid], e = offsets_[mid + 1];
      out.reserve(e - b);
      for (uint32_t p = b; p < e; ++p)
      {
        uint32_t idx = 0;
        for (int j = tax_bytes_ - 1; j >= 0; --j)
          idx = (idx << 8)
              | post_[static_cast<size_t>(p) * static_cast<size_t>(tax_bytes_) + static_cast<size_t>(j)];
        if (idx < n_taxa_) out.push_back(taxa_[idx]);
      }
      return;
    }
  }

  // --- Serialization -------------------------------------------------------
  //
  // v2 ("FTX2"), little-endian, every section 8-byte aligned so the mapped
  // image is usable in place:
  //
  //   0  magic "FTX2" | 4 u32 version=2 | 8 u32 k | 12 u32 n_taxa
  //   16 u64 n_kmers  | 24 u64 n_postings | 32 u8 key_bytes | 33 u8 tax_bytes
  //   34 u8 reserved[6]
  //   40 u32 taxa[n_taxa]             ascending
  //      u64 taxon_postings[n_taxa]   STORED, not recomputed: scanning the whole
  //                                   postings array on every load is exactly the
  //                                   cost mapping the file exists to avoid
  //      keys[n_kmers * key_bytes]    base-19 codes, ascending
  //      u32 offsets[n_kmers + 1]     into postings
  //      postings[n_post * tax_bytes] INDICES into taxa[], ascending per k-mer
  bool TaxIndex::save(const std::string& path, std::string* err) const
  {
    if (k_ < 1 || k_ > MAX_K)
    { if (err) *err = "refusing to save an index with no valid k"; return false; }
    const std::vector<uint32_t> tx = taxa();
    if (tx.empty())
    { if (err) *err = "refusing to save an empty index"; return false; }
    const uint32_t n_taxa = static_cast<uint32_t>(tx.size());
    const int kb = keyBytes(k_);
    // Width follows the largest INDEX stored (n_taxa-1), not the count:
    // 256 taxa still index 0..255, one byte.
    const uint32_t max_idx = n_taxa == 0 ? 0 : n_taxa - 1;
    const int tb = max_idx <= 0xFFu ? 1 : (max_idx <= 0xFFFFu ? 2 : 4);

    std::vector<std::pair<uint64_t, std::vector<uint32_t>>> rows;
    if (legacy_)
    {
      rows.reserve(index_.size());
      for (const auto& kv : index_)
      {
        const uint64_t code = encode(kv.first.data(), k_);
        if (code == UINT64_MAX) continue;
        std::vector<uint32_t> idx;
        idx.reserve(kv.second.size());
        for (uint32_t t : kv.second)
        {
          auto it = std::lower_bound(tx.begin(), tx.end(), t);
          if (it != tx.end() && *it == t) idx.push_back(static_cast<uint32_t>(it - tx.begin()));
        }
        rows.emplace_back(code, std::move(idx));
      }
      std::sort(rows.begin(), rows.end(),
                [](const std::pair<uint64_t, std::vector<uint32_t>>& a,
                   const std::pair<uint64_t, std::vector<uint32_t>>& b) { return a.first < b.first; });
    }
    else
    {
      // Re-saving an already-mapped index: its keys are already sorted.
      rows.reserve(static_cast<size_t>(n_kmers_));
      for (uint64_t i = 0; i < n_kmers_; ++i)
      {
        std::vector<uint32_t> idx;
        for (uint32_t p = offsets_[i]; p < offsets_[i + 1]; ++p)
        {
          uint32_t v = 0;
          for (int j = tax_bytes_ - 1; j >= 0; --j)
            v = (v << 8)
              | post_[static_cast<size_t>(p) * static_cast<size_t>(tax_bytes_) + static_cast<size_t>(j)];
          idx.push_back(v);
        }
        rows.emplace_back(rdKey(keys_ + i * static_cast<size_t>(key_bytes_), key_bytes_), std::move(idx));
      }
    }

    uint64_t n_post = 0;
    for (const auto& r : rows) n_post += r.second.size();
    if (n_post > 0xFFFFFFFFull)
    { if (err) *err = "index has too many postings for the 32-bit offset format"; return false; }
    std::vector<uint64_t> tpost(n_taxa, 0);
    for (const auto& r : rows) for (uint32_t i : r.second) if (i < n_taxa) ++tpost[i];

    std::ofstream o(path, std::ios::binary);
    if (!o) { if (err) *err = "cannot open " + path + " for writing"; return false; }

    o.write("FTX2", 4);
    wr<uint32_t>(o, 2);
    wr<uint32_t>(o, static_cast<uint32_t>(k_));
    wr<uint32_t>(o, n_taxa);
    wr<uint64_t>(o, static_cast<uint64_t>(rows.size()));
    wr<uint64_t>(o, n_post);
    wr<uint8_t>(o, static_cast<uint8_t>(kb));
    wr<uint8_t>(o, static_cast<uint8_t>(tb));
    for (int i = 0; i < 6; ++i) wr<uint8_t>(o, 0);

    // Track the position ourselves: a failed stream returns tellp()==-1, which
    // as size_t is SIZE_MAX and would spin the pad loop forever.
    size_t pos = 40;
    auto put = [&](const void* d, size_t n) { o.write(static_cast<const char*>(d), static_cast<std::streamsize>(n)); pos += n; };
    auto pad = [&]() { while ((pos % 8u) != 0 && o) { char z = 0; o.write(&z, 1); ++pos; } };

    for (uint32_t t : tx) put(&t, 4);
    pad();
    for (uint64_t v : tpost) put(&v, 8);
    pad();

    std::vector<uint8_t> kbuf(static_cast<size_t>(kb));
    for (const auto& r : rows) { wrKey(kbuf.data(), kb, r.first); put(kbuf.data(), static_cast<size_t>(kb)); }
    pad();

    uint32_t run = 0;
    for (const auto& r : rows) { put(&run, 4); run += static_cast<uint32_t>(r.second.size()); }
    put(&run, 4);
    pad();

    std::vector<uint8_t> pbuf(static_cast<size_t>(tb));
    for (const auto& r : rows)
      for (uint32_t i : r.second)
      { uint32_t v = i; for (int j = 0; j < tb; ++j) { pbuf[j] = static_cast<uint8_t>(v & 0xFF); v >>= 8; } put(pbuf.data(), static_cast<size_t>(tb)); }

    o.flush();
    return static_cast<bool>(o);
  }

  bool TaxIndex::load(const std::string& path, std::string* err)
  {
    reset();

    char magic[4] = {0, 0, 0, 0};
    {
      std::ifstream probe(path, std::ios::binary);
      if (!probe) { if (err) *err = "cannot open " + path; return false; }
      probe.read(magic, 4);
    }

    // ---- v1: parse into the maps, so pre-v2 indexes still load --------------
    if (std::memcmp(magic, "FTXI", 4) == 0)
    {
      std::ifstream i(path, std::ios::binary);
      i.seekg(4);
      uint32_t ver = 0, k = 0; uint64_t n = 0;
      if (!rd(i, ver) || ver != 1) { if (err) *err = "unsupported index version"; return false; }
      if (!rd(i, k) || !rd(i, n)) { if (err) *err = "truncated index header"; return false; }
      k_ = static_cast<int>(k);
      legacy_ = true;
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

    if (std::memcmp(magic, "FTX2", 4) != 0)
    { if (err) *err = "not a FasTag taxonomy index"; return false; }

    // ---- v2: map it and point into the image -------------------------------
#ifdef _WIN32
    HANDLE fh = CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (fh == INVALID_HANDLE_VALUE) { if (err) *err = "cannot open " + path; return false; }
    LARGE_INTEGER sz;
    if (GetFileSizeEx(fh, &sz) == 0)
    { CloseHandle(fh); if (err) *err = "cannot size " + path; return false; }
    HANDLE mh = CreateFileMappingA(fh, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (mh == nullptr)
    { CloseHandle(fh); if (err) *err = "cannot map " + path; return false; }
    void* base = MapViewOfFile(mh, FILE_MAP_READ, 0, 0, 0);
    if (base == nullptr)
    { CloseHandle(mh); CloseHandle(fh); if (err) *err = "cannot map " + path; return false; }
    file_handle_ = fh;
    map_handle_ = mh;
    map_len_ = static_cast<size_t>(sz.QuadPart);
#else
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { if (err) *err = "cannot open " + path; return false; }
    struct stat st;
    if (fstat(fd, &st) != 0) { ::close(fd); if (err) *err = "cannot size " + path; return false; }
    void* base = mmap(nullptr, static_cast<size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) { ::close(fd); if (err) *err = "cannot map " + path; return false; }
    map_fd_ = fd;
    map_len_ = static_cast<size_t>(st.st_size);
#endif
    map_base_ = base;

    if (map_len_ < 40) { reset(); if (err) *err = "index is truncated"; return false; }
    const uint8_t* p = static_cast<const uint8_t*>(base);
    auto u32at = [p](size_t o) { uint32_t v; std::memcpy(&v, p + o, 4); return v; };
    auto u64at = [p](size_t o) { uint64_t v; std::memcpy(&v, p + o, 8); return v; };

    // Re-check the magic in the MAPPED bytes, not just the earlier probe read:
    // the file could have been replaced between the two opens.
    if (std::memcmp(p, "FTX2", 4) != 0) { reset(); if (err) *err = "not a FasTag taxonomy index"; return false; }
    if (u32at(4) != 2) { reset(); if (err) *err = "unsupported index version"; return false; }
    k_ = static_cast<int>(u32at(8));
    n_taxa_ = u32at(12);
    n_kmers_ = u64at(16);
    const uint64_t n_post = u64at(24);
    key_bytes_ = p[32];
    tax_bytes_ = p[33];

    if (key_bytes_ <= 0 || key_bytes_ > 8 || tax_bytes_ <= 0 || tax_bytes_ > 4
        || k_ <= 0 || k_ > MAX_K
        || key_bytes_ != keyBytes(k_)
        || n_kmers_ == 0 || n_taxa_ == 0)
    { reset(); if (err) *err = "index header is malformed"; return false; }

    // Every section is bounds-checked BEFORE its pointer is formed, with
    // subtraction rather than addition so attacker-controlled counts cannot wrap
    // a 64-bit offset past the end and back inside the mapping. This matters in
    // practice, not just against a hostile file: the ~1 GB index is a downloaded
    // asset, so a truncated or partial download is a realistic corruption.
    // `avail` is what remains after `off`; a section of `count` elements of
    // `size` bytes fits iff count <= avail/size.
    size_t off = 40;
    auto fits = [&](uint64_t count, size_t size) -> bool {
      if (off > map_len_) return false;
      const size_t avail = map_len_ - off;
      return size == 0 || count <= avail / size;
    };
    auto bump = [&](uint64_t count, size_t size) { off = align8(off + static_cast<size_t>(count) * size); };

    if (!fits(n_taxa_, 4)) { reset(); if (err) *err = "index is truncated (taxa)"; return false; }
    taxa_ = reinterpret_cast<const uint32_t*>(p + off);
    bump(n_taxa_, 4);
    if (!fits(n_taxa_, 8)) { reset(); if (err) *err = "index is truncated (taxpost)"; return false; }
    taxpost_ = reinterpret_cast<const uint64_t*>(p + off);
    bump(n_taxa_, 8);
    if (!fits(n_kmers_, static_cast<size_t>(key_bytes_))) { reset(); if (err) *err = "index is truncated (keys)"; return false; }
    keys_ = p + off;
    bump(n_kmers_, static_cast<size_t>(key_bytes_));
    if (!fits(n_kmers_ + 1, 4)) { reset(); if (err) *err = "index is truncated (offsets)"; return false; }
    offsets_ = reinterpret_cast<const uint32_t*>(p + off);
    bump(n_kmers_ + 1, 4);
    if (!fits(n_post, static_cast<size_t>(tax_bytes_))) { reset(); if (err) *err = "index is truncated (postings)"; return false; }
    post_ = p + off;

    // The offsets array must be a valid CSR structure, or lookup() would read
    // arbitrary bytes of the postings section (or far past it). Checked once here
    // so the hot path can trust it.
    if (offsets_[0] != 0 || offsets_[n_kmers_] != n_post)
    { reset(); if (err) *err = "index offsets are malformed"; return false; }
    for (uint64_t i = 1; i <= n_kmers_; ++i)
      if (offsets_[i] < offsets_[i - 1])
      { reset(); if (err) *err = "index offsets are not monotonic"; return false; }

    legacy_ = false;
    return true;
  }
}
