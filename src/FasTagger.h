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

#include <memory>
#include <string>
#include <vector>

namespace FasTag
{
  struct Alphabet;  ///< residue/pair tables built from Param.mods; defined in the .cpp

  /// One residue-specific modification, resolved to a mass by the TOPP tool
  /// (from OpenMS ModificationsDB) so this library needs no modification
  /// database of its own.
  ///
  /// Terminal modifications are NOT represented here: a fixed N- or C-terminal
  /// mod shifts the whole precursor mass and is therefore absorbed into the
  /// reported flanking masses automatically, with no effect on which internal
  /// edges the graph has. Only residue-specific mods change the alphabet.
  struct ModSpec
  {
    char        residue;   ///< the residue this applies to (e.g. 'C', 'S', 'K')
    double      delta;     ///< monoisotopic mass shift
    std::string name;      ///< display name, e.g. "Phospho", "Carbamidomethyl"
    bool        variable;  ///< true: adds a modified alternative to the residue;
                           ///< false: shifts the residue's mass unconditionally
  };

  /// Tagging parameters. Every field is set by the TOPP tool; there are no
  /// hidden knobs and no compatibility modes.
  struct Param
  {
    int    tag_length     = 3;      ///< seed length in residues
    int    max_extension  = 0;      ///< residues appended per terminus; 0 = off
    int    max_gaps       = 0;      ///< gap edges allowed per tag; 0 = off, 1 = max
    /// Collapse isotope clusters to their monoisotopic peak before selecting
    /// peaks, and move multiply-charged fragments onto the singly-charged scale.
    bool   deisotope      = false;
    int    deisotope_max_charge = 3;
    double frag_tol       = 20.0;   ///< fragment tolerance
    bool   tol_ppm        = true;
    double complement_tol = 0.02;   ///< always absolute; complements span the whole range
    double precursor_tol  = 1.5;
    size_t max_peak_count = 100;   ///< hard ceiling; also bounds the null tables
    /// Peaks kept per 100 Da window; 0 selects the strongest max_peak_count
    /// overall instead. Windowed selection makes the effective budget scale with
    /// the spectrum's own m/z range rather than being a constant.
    int    peaks_per_window = 0;
    int    max_tag_count  = 50;     ///< 0 = unlimited
    double max_evalue     = 20.0;   ///< 0 = no cutoff
    /// Multiplier applied to a gapped tag's E-value for RANKING. 1.0 disables it.
    ///
    /// Gapped tags are systematically over-scored against contiguous ones: a gap
    /// edge is chosen as the best fit from a 190-entry two-residue table where an
    /// ordinary edge tests 19 residues, and it spends one fewer peak to spell one
    /// more residue. Left uncorrected, gapped tags took 95% of rank-1 slots while
    /// being 3.6x less likely to be right.
    ///
    /// The multiplicity argument gives |pairs|/|residues| = 10 as a LOWER bound.
    /// Measurement says the real over-scoring is larger: rank-1 accuracy keeps
    /// improving to ~100 and is flat above it, consistently across all four
    /// benchmark profiles and tag lengths 3-5 (bench/benchmark.cpp). 100 is
    /// therefore calibrated, not derived -- and calibrated on synthetic spectra,
    /// which is the weakest part of the claim. It is a knob for that reason.
    ///
    /// Ranking only: the E-value cutoff is applied to the UNCORRECTED value, so
    /// this never removes a tag, only reorders. See tagSpectrum().
    double gap_penalty    = 100.0;
    unsigned seed         = 20080717u;  ///< the m/z-fidelity null is Monte Carlo
    int    mzfidelity_samples = 10000;
    /// Residue-specific modifications. Fixed mods shift a residue's mass;
    /// variable mods add a modified alternative that the graph may also use.
    /// Empty means the unmodified 19-residue alphabet. Resolved by the TOPP
    /// tool; the default it passes is a single fixed Carbamidomethyl (C).
    std::vector<ModSpec> mods;
  };

  /// One tag. Flanking masses assume the low-m/z peak is a y ion, as DirecTag
  /// does; a tag read off the b series therefore reads reversed relative to the
  /// peptide and its flanks carry a one-water offset.
  struct Tag
  {
    std::string seq;            ///< N->C. Modified residues are written inline as
                                ///< X[ModName] (e.g. "GS[Phospho]TK"), so seq is
                                ///< not one char per residue -- use n_res.
    size_t n_res = 0;           ///< residue count (bracket annotations excluded)
    double nterm_mass = 0;      ///< residue mass N-terminal to the tag
    double cterm_mass = 0;      ///< residue mass C-terminal to the tag
    int    charge = 1;          ///< fragment charge the tag was read at
    double p_intensity = 1, p_mzfidelity = 1, p_complement = 1;
    double evalue = 1;          ///< combined p-value x tags enumerated
    /// Per-edge confidence in [0,1], summarised. Each edge (one residue, or two
    /// for a gap) scores m/z fit x endpoint intensity; min_conf is the weakest
    /// edge (the residue most likely wrong), mean_conf the average. The E-value
    /// is the primary, validated tag confidence; these localise WHERE a tag is
    /// weak, which the single E-value cannot.
    float  min_conf = 1.0f;
    float  mean_conf = 1.0f;
    bool   extended = false;
    bool   gapped = false;      ///< crossed a missing peak; two residues are a
                                ///< mass-only inference, not two observed steps
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

    /// The residue/pair alphabet for this run, built once from Param.mods.
    const Alphabet& alphabet() const { return *alpha_; }

  private:
    int k_min_, k_max_;
    size_t n_max_;
    std::vector<std::vector<std::vector<float>>> ranksum_;  ///< [k][npeaks] -> CDF
    std::vector<std::vector<double>> mz_null_;              ///< [k] -> sorted SSE samples
    std::shared_ptr<const Alphabet> alpha_;                 ///< residue + pair tables
  };

  /// Infer sequence tags from one centroided MS/MS spectrum.
  ///
  /// @param spec       fragment spectrum; must be sorted by m/z
  /// @param precursor_mz  precursor m/z
  /// @param charge     precursor charge; <= 0 is treated as 2
  /// @return tags ordered by E-value, best first
  std::vector<Tag> tagSpectrum(const OpenMS::MSSpectrum& spec, double precursor_mz,
                               int charge, const Param& p, const Tables& tables);

  /// Strip inline modification annotations, leaving one base residue per
  /// position: "GS[Phospho]TK[TMT6plex]" -> "GSTK". Used for FASTA matching,
  /// which is against unmodified protein sequences.
  inline std::string baseSequence(const std::string& annotated)
  {
    std::string out;
    out.reserve(annotated.size());
    int depth = 0;
    for (char c : annotated)
    {
      if (c == '[') { ++depth; continue; }
      if (c == ']') { if (depth > 0) --depth; continue; }
      if (depth == 0) out.push_back(c);
    }
    return out;
  }
}
