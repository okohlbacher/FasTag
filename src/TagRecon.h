// Tag reconciliation: place a sequence tag onto database peptides by its
// flanking masses, and localize the single mass gap (modification or mutation)
// that the flanks imply. The database-search counterpart to tagging, after
// DirecTag/TagRecon (Dasari et al., J Proteome Res 2010, 9:1716).
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
  /// One reconciliation of a tag against one database peptide placement.
  struct Reconciliation
  {
    std::string protein;    ///< accession of the protein the peptide came from
    std::string peptide;    ///< the digested peptide (I/L as in the database)
    size_t      pos = 0;    ///< 0-based residue index where the tag starts in peptide
    bool        reversed = false;  ///< tag matched the peptide reading C->N (b-derived)
    /// Flank mass agreement. Both true = exact placement (no mass gap).
    bool        nterm_match = false;
    bool        cterm_match = false;
    double      delta_mass = 0.0;  ///< spectrum flank - database flank on the mismatched side; 0 if exact
    /// Inclusive residue range (in peptide coordinates) the mass gap localizes to
    /// -- the residues between the tag and the mismatched terminus. Empty (lo>hi)
    /// when exact.
    int         region_lo = 0, region_hi = -1;
    /// Stage B: what the mass gap most likely IS. Empty when exact. Otherwise one
    /// of "mod:Name@X" (a known modification on a residue X present in the
    /// region), "sub:X->Y" (a single residue X in the region replaced by Y), or
    /// "?" (delta explained by neither within tolerance). The delta on its own is
    /// a number; this is the interpretation a search would act on.
    std::string delta_interp;
  };

  /// A candidate variable modification for stage-B interpretation: its
  /// monoisotopic mass shift and the residues it can sit on ("" = any). Supplied
  /// by the caller (the TOPP tool resolves names via OpenMS ModificationsDB), so
  /// this library needs no modification database of its own.
  struct ModCandidate
  {
    std::string name;
    double      delta = 0.0;
    std::string residues;  ///< e.g. "STY" for Phospho; "" matches any residue
  };

  /// Reconciles tags against a tryptically digested protein database.
  ///
  /// Holds the digested peptides and, for every tag-length k-mer, the peptide
  /// placements it occurs at, so a tag can be located and its database flanking
  /// masses computed and compared to the spectrum's. I is folded to L, exactly
  /// as the tagger and FastaFilter do, since a spectrum cannot distinguish them.
  class TagReconciler
  {
  public:
    /// @param frag_tol,tol_ppm  flank-mass agreement tolerance (per flank).
    ///   ppm is evaluated at the flank mass itself.
    /// @param both_orientations also try the tag reversed (b-ion reading).
    TagReconciler(double frag_tol, bool tol_ppm, bool both_orientations = true)
      : frag_tol_(frag_tol), tol_ppm_(tol_ppm), both_(both_orientations) {}

    /// Digest the FASTA entries (trypsin, no cut before proline, up to
    /// @p missed_cleavages) and index every @p k-mer of the resulting peptides.
    /// Residue masses come from OpenMS with the given fixed modifications applied
    /// (variable mods are not part of the database; they are what reconciliation
    /// discovers as mass gaps). Call once; reconcile() is const/thread-safe after.
    void build(const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries, int k,
               int missed_cleavages,
               const std::vector<std::pair<char, double>>& fixed_mods);

    /// All placements of @p tag consistent with the flanking masses, at most one
    /// flank mismatch each. @p tag is base residues (bracket annotations already
    /// stripped), N->C, at the k used in build().
    std::vector<Reconciliation> reconcile(const std::string& tag, double nterm_mass,
                                          double cterm_mass) const;

    /// Candidate variable modifications for stage-B delta interpretation. Set
    /// before reconcile(); read-only after. Empty means only substitutions and
    /// "?" are ever reported.
    void setModCandidates(std::vector<ModCandidate> mods) { mods_ = std::move(mods); }

    size_t peptideCount() const { return peptides_.size(); }
    size_t indexedKmers() const { return index_.size(); }

  private:
    struct Peptide
    {
      std::string seq;              ///< I/L folded
      std::string protein;          ///< source accession
      std::vector<double> prefix;   ///< prefix[i] = residue mass of seq[0..i); size seq.len()+1
    };
    struct Placement { uint32_t pep; uint32_t pos; };

    double tolAt(double m) const { return tol_ppm_ ? m * frag_tol_ * 1e-6 : frag_tol_; }
    std::string fold(const std::string& s) const;
    void tryPlace(const Peptide& pep, size_t pos, const std::string& tag, bool reversed,
                  double nterm_mass, double cterm_mass, std::vector<Reconciliation>& out) const;
    /// Best explanation of @p delta over the residues in @p region: a candidate
    /// mod on a present residue, a single substitution of a present residue, or
    /// "?" if neither fits within tolerance. Prefers whichever fits with the
    /// smaller mass error; on a near-tie a modification wins (it is the more
    /// common explanation).
    std::string interpretDelta_(double delta, const std::string& region) const;

    double frag_tol_;
    bool   tol_ppm_;
    bool   both_;
    int    k_ = 0;
    double residue_masses_[128] = {0};  ///< by ASCII code, fixed mods folded in
    std::vector<ModCandidate> mods_;    ///< candidate variable mods for stage B
    std::vector<Peptide> peptides_;
    std::unordered_map<std::string, std::vector<Placement>> index_;  ///< k-mer -> placements
  };
}
