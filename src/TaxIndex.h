// Taxonomic k-mer index for tag-driven species detection.
//
// Maps each amino-acid k-mer (I/L folded) occurring in a reduced, curated
// reference protein set to the SET of taxa whose proteins contain it. Taxon
// sets are kept, not collapsed to a per-k-mer LCA: collapsing at index time is
// irreversible and one misannotated protein would drag a k-mer to the root,
// where a reduced curated database still lets outliers be down-weighted and the
// lowest-common-ancestor be computed at report time (see the adversarial review
// in the design notes). At the reduced-taxa scale (hundreds to low thousands of
// genera) a plain sorted vector<uint32_t> per k-mer is compact enough that no
// roaring-bitmap dependency is warranted.
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
  /// k-mer -> taxon set index over a reduced reference protein database.
  ///
  /// Thread-safe for concurrent const reads after build()/load().
  class TaxIndex
  {
  public:
    /// Build from protein entries whose taxid is supplied in parallel
    /// (taxids[i] is the taxon of entries[i]; 0 skips the entry). Only k-mers of
    /// the single canonical length @p k are indexed, I is folded to L. A taxon is
    /// recorded once per distinct k-mer regardless of how many of its proteins or
    /// positions carry it, so posting counts reflect taxon breadth, not abundance.
    void build(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries,
               const std::vector<uint32_t>& taxids, int k);

    int k() const { return k_; }

    /// Taxa whose proteins contain @p kmer (I/L folded, length k). Empty if none.
    /// The returned vector is sorted ascending and de-duplicated.
    const std::vector<uint32_t>& lookup(const std::string& kmer) const;

    /// Background for the null "taxon T absent": the number of indexed k-mers
    /// whose taxon set contains T (its posting count). p_T = postings(T)/nKmers().
    uint64_t postings(uint32_t taxid) const;
    uint64_t nKmers() const { return static_cast<uint64_t>(index_.size()); }

    /// Every taxon that appears anywhere in the index, ascending.
    std::vector<uint32_t> taxa() const;

    /// Compact binary (de)serialization. Format is versioned; load() returns
    /// false on a version or magic mismatch.
    bool save(const std::string& path, std::string* err = nullptr) const;
    bool load(const std::string& path, std::string* err = nullptr);

    static std::string fold(const std::string& s);  ///< I -> L, uppercase passthrough

  private:
    int k_ = 0;
    std::unordered_map<std::string, std::vector<uint32_t>> index_;  ///< kmer -> sorted taxid set
    std::unordered_map<uint32_t, uint64_t> postings_;               ///< taxid -> posting count
    static const std::vector<uint32_t> empty_;
  };
}
