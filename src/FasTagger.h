// FasTag -- partial sequence tags from peptide MS/MS spectra.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
//
// Reimplements the DirecTag algorithm (Tabb et al., J Proteome Res 2008, 7:3838)
// from its publication. No code was copied from the reference implementation.
//
// Chemistry, spectrum handling, peak filtering and statistics come from OpenMS
// and Boost. The only mathematics kept in-house is the rank-sum null in
// ranksum.h, which neither provides.
#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <FasTag/ranksum.h>

#include <string>
#include <vector>

namespace FasTag
{
  /// Tagging parameters. Every field is set by the TOPP tool; there are no
  /// hidden knobs and no compatibility modes.
  struct Param
  {
    int    tag_length     = 3;      ///< seed length in residues
    int    max_extension  = 0;      ///< residues appended per terminus; 0 = off
    double frag_tol       = 20.0;   ///< fragment tolerance
    bool   tol_ppm        = true;
    double complement_tol = 0.02;   ///< always absolute; complements span the whole range
    double precursor_tol  = 1.5;
    size_t max_peak_count = 100;
    int    max_tag_count  = 50;     ///< 0 = unlimited
    double max_evalue     = 20.0;   ///< 0 = no cutoff
    unsigned seed         = 20080717u;  ///< the m/z-fidelity null is Monte Carlo
    int    mzfidelity_samples = 10000;
  };

  /// One tag. Flanking masses assume the low-m/z peak is a y ion, as DirecTag
  /// does; a tag read off the b series therefore reads reversed relative to the
  /// peptide and its flanks carry a one-water offset.
  struct Tag
  {
    std::string seq;            ///< N->C
    double nterm_mass = 0;      ///< residue mass N-terminal to the tag
    double cterm_mass = 0;      ///< residue mass C-terminal to the tag
    int    charge = 1;          ///< fragment charge the tag was read at
    double p_intensity = 1, p_mzfidelity = 1, p_complement = 1;
    double evalue = 1;          ///< combined p-value x tags enumerated
    bool   extended = false;
    double low_mz = 0;          ///< lowest peak; also the output tie-break
  };

  /// Precomputed nulls, built once and shared read-only.
  ///
  /// Eager construction is what lets several spectra be tagged concurrently: a
  /// lazily-filled cache would be written from every thread. Everything needed is
  /// known from Param, since the realised tag length spans
  /// [tag_length+1, tag_length+1+2*max_extension] peaks.
  class Tables
  {
  public:
    explicit Tables(const Param& p);

    /// P(rank sum <= observed) for a tag of @p tag_peaks peaks among @p npeaks.
    /// Conservative (1.0) outside the built range rather than out-of-bounds.
    double intensityP(int tag_peaks, int npeaks, int ranksum) const;

    /// P(SSE <= observed) for the m/z-fidelity statistic.
    double mzFidelityP(int tag_peaks, double sse) const;

  private:
    int k_min_, k_max_;
    size_t n_max_;
    std::vector<std::vector<std::vector<float>>> ranksum_;  ///< [k][npeaks] -> CDF
    std::vector<std::vector<double>> mz_null_;              ///< [k] -> sorted SSE samples
  };

  /// Infer sequence tags from one centroided MS/MS spectrum.
  ///
  /// @param spec       fragment spectrum; must be sorted by m/z
  /// @param precursor_mz  precursor m/z
  /// @param charge     precursor charge; <= 0 is treated as 2
  /// @return tags ordered by E-value, best first
  std::vector<Tag> tagSpectrum(const OpenMS::MSSpectrum& spec, double precursor_mz,
                               int charge, const Param& p, const Tables& tables);
}
