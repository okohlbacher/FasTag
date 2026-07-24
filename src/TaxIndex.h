// Taxonomic k-mer index for tag-driven species detection.
//
// Maps each amino-acid k-mer (I/L folded) occurring in a curated reference
// protein set to the SET of taxa whose proteins contain it. Taxon sets are
// kept, not collapsed to a per-k-mer LCA: collapsing at index time is
// irreversible and one misannotated protein would drag a k-mer to the root,
// where keeping the set lets outliers be down-weighted and the
// lowest-common-ancestor be computed at report time.
//
// ON-DISK FORMAT AND MEMORY
//
// v1 stored `unordered_map<string, vector<uint32_t>>` and rebuilt it on load.
// That cost ~90 bytes per k-mer in RAM against 16.9 on disk -- a map node, an
// SSO string, and a heap-allocated vector holding on average 1.2 taxids. A
// 50-taxon index measured 2.13 GB on disk and 10.2 GB resident, which put the
// feature out of reach on an ordinary laptop.
//
// v2 is a flat, packed, memory-MAPPED image:
//   * k-mers are packed base-19 (the folded alphabet has exactly 19 letters),
//     so k=7 needs 30 bits -> 4 bytes, against 7 bytes of characters.
//   * postings store an INDEX into the taxon table, not a taxid: 50 taxa need
//     one byte where a uint32 was used.
//   * nothing is parsed at load. The file is mapped and the arrays are used in
//     place, so resident memory is file-backed and EVICTABLE under pressure
//     rather than anonymous heap the OS cannot reclaim.
//
// v1 files still load (via the legacy in-memory path) so existing indexes keep
// working; save() always writes v2.
//
// ACCEPTED LIMITATIONS (deliberate, after an adversarial review):
//   * Little-endian, 64-bit hosts only. macOS/Windows/Linux on x86-64/arm64 are
//     all LE; a big-endian reader would see byte-swapped numeric fields, and a
//     32-bit host cannot map a >4 GB index. Both are non-targets.
//   * load() bounds-checks every section against the file length and validates
//     the offset array, so a truncated or partial download is rejected, not
//     crashed on. It does NOT fully validate SEMANTICS (keys sorted, taxon
//     indices in range) -- a corrupt-but-in-bounds file yields wrong results,
//     not a fault. The index is a trusted, checksummed release asset.
//   * Reading uint32/uint64 arrays straight from the mapping is technically
//     object-lifetime UB in C++17, but the sections are 8-byte aligned and the
//     map base is page-aligned, so it is well-defined on every target compiler.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#pragma once

#include <OpenMS/FORMAT/FASTAFile.h>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace FasTag
{
  /// k-mer -> taxon set index over a curated reference protein database.
  ///
  /// Thread-safe for concurrent const reads after build()/load().
  class TaxIndex
  {
  public:
    TaxIndex() = default;
    ~TaxIndex();
    TaxIndex(const TaxIndex&) = delete;             ///< owns a mapping
    TaxIndex& operator=(const TaxIndex&) = delete;

    /// Build from protein entries whose taxid is supplied in parallel
    /// (taxids[i] is the taxon of entries[i]; 0 skips the entry). Only k-mers of
    /// the single canonical length @p k are indexed, I is folded to L. A taxon is
    /// recorded once per distinct k-mer regardless of how many of its proteins or
    /// positions carry it, so posting counts reflect taxon breadth, not abundance.
    void build(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries,
               const std::vector<uint32_t>& taxids, int k);

    int k() const { return k_; }

    /// Taxa whose proteins contain @p kmer (I/L folded, length k), ascending.
    ///
    /// Fills @p out rather than returning a reference: the mapped image stores
    /// taxon indices, so there is no vector<uint32_t> in the file to point at.
    /// Callers reuse one buffer instead of allocating per lookup.
    void lookup(const std::string& kmer, std::vector<uint32_t>& out) const;

    /// Background for the null "taxon T absent": the number of indexed k-mers
    /// whose taxon set contains T (its posting count). p_T = postings(T)/nKmers().
    uint64_t postings(uint32_t taxid) const;
    uint64_t nKmers() const;

    /// Every taxon that appears anywhere in the index, ascending.
    std::vector<uint32_t> taxa() const;

    /// Compact binary (de)serialization. save() writes v2; load() accepts v1
    /// and v2 and returns false on any other magic or version.
    bool save(const std::string& path, std::string* err = nullptr) const;
    bool load(const std::string& path, std::string* err = nullptr);

    static std::string fold(const std::string& s);  ///< I -> L, uppercase passthrough

    /// Pack a folded k-mer into its base-19 code, or UINT64_MAX if it contains a
    /// residue outside the 19-letter alphabet (X/B/J/O/U/Z and anything else).
    static uint64_t encode(const char* s, int k);
    /// Bytes needed to store a base-19 code of length @p k.
    static int keyBytes(int k);

    /// k must be in [1, MAX_K]: 19^MAX_K fits in a uint64 code with margin.
    /// Higher k overflows the packing and silently corrupts the index.
    static constexpr int MAX_K = 12;

  private:
    void reset();

    int k_ = 0;

    // ---- v2: the mapped image (all pointers are into map_base_) -------------
    void* map_base_ = nullptr;
    size_t map_len_ = 0;
#ifdef _WIN32
    void* map_handle_ = nullptr;
    void* file_handle_ = nullptr;
#else
    int map_fd_ = -1;
#endif
    const uint8_t* keys_ = nullptr;       ///< n_kmers_ * key_bytes_, sorted
    const uint32_t* offsets_ = nullptr;   ///< n_kmers_ + 1
    const uint8_t* post_ = nullptr;       ///< n_postings * tax_bytes_
    const uint32_t* taxa_ = nullptr;      ///< n_taxa_, ascending
    const uint64_t* taxpost_ = nullptr;   ///< n_taxa_ posting counts
    uint64_t n_kmers_ = 0;
    uint32_t n_taxa_ = 0;
    int key_bytes_ = 0;
    int tax_bytes_ = 0;

    // ---- v1 legacy / build scratch -----------------------------------------
    std::unordered_map<std::string, std::vector<uint32_t>> index_;
    std::unordered_map<uint32_t, uint64_t> postings_;
    bool legacy_ = false;
  };
}
