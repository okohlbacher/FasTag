// Deterministic random subsampling of spectra.
//
// For a taxon call (and quick previews) a run of hundreds of thousands of
// spectra does not need every spectrum: a uniformly random subset gives the
// same lowest-common-ancestor at a fraction of the cost, and repeated draws
// with different seeds give a bootstrap distribution over the call. Selection
// is seeded so a run is exactly reproducible, and the mask is indexed by the
// spectrum's global position so downstream output order is unaffected.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

namespace FasTag
{
  /// A 0/1 mask over [0, n_total): mask[i] != 0 iff spectrum i is selected.
  using SampleMask = std::vector<char>;

  /// Select @p n_want of @p n_total indices uniformly at random without
  /// replacement, deterministically from @p seed. Returns a mask of length
  /// n_total. n_want >= n_total selects everything; n_want == 0 selects nothing.
  ///
  /// Partial Fisher-Yates: shuffle only the first n_want positions of an index
  /// array, so it is O(n_total) time and memory and draws a uniform subset (each
  /// n_want-subset equally likely). mt19937 keeps the draw identical across
  /// platforms for a given seed.
  inline SampleMask sampleByCount(size_t n_total, size_t n_want, uint32_t seed)
  {
    SampleMask mask(n_total, 0);
    if (n_total == 0 || n_want == 0) return mask;
    if (n_want >= n_total) { std::fill(mask.begin(), mask.end(), char(1)); return mask; }

    std::vector<size_t> idx(n_total);
    std::iota(idx.begin(), idx.end(), size_t(0));
    std::mt19937 rng(seed);
    for (size_t i = 0; i < n_want; ++i)
    {
      std::uniform_int_distribution<size_t> pick(i, n_total - 1);
      std::swap(idx[i], idx[pick(rng)]);
      mask[idx[i]] = 1;
    }
    return mask;
  }

  /// Select a @p fraction (0..1] of @p n_total indices, rounding to nearest,
  /// at least one when fraction > 0 and n_total > 0.
  inline SampleMask sampleByFraction(size_t n_total, double fraction, uint32_t seed)
  {
    if (fraction <= 0.0 || n_total == 0) return SampleMask(n_total, 0);
    if (fraction >= 1.0) return SampleMask(n_total, char(1));
    size_t n_want = static_cast<size_t>(fraction * static_cast<double>(n_total) + 0.5);
    if (n_want == 0) n_want = 1;
    return sampleByCount(n_total, n_want, seed);
  }
}
