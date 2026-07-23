// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TagRecon.h"

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <algorithm>
#include <cmath>

using namespace OpenMS;

namespace FasTag
{
  std::string TagReconciler::fold(const std::string& s) const
  {
    std::string out = s;
    for (char& c : out) if (c == 'I') c = 'L';
    return out;
  }

  void TagReconciler::build(const std::vector<FASTAFile::FASTAEntry>& entries, int k,
                            int missed_cleavages,
                            const std::vector<std::pair<char, double>>& fixed_mods)
  {
    k_ = k;

    // Internal residue masses from OpenMS, with fixed mods folded into the
    // residue they modify -- the database peptides carry fixed mods, so their
    // flanking masses must too. Variable mods are deliberately NOT applied here:
    // they are exactly what reconciliation is meant to discover as a mass gap.
    for (const Residue* r : ResidueDB::getInstance()->getResidues("Natural19WithoutI"))
    {
      char c = r->getOneLetterCode()[0];
      residue_masses_[static_cast<unsigned char>(c)] = r->getMonoWeight(Residue::Internal);
    }
    residue_masses_[static_cast<unsigned char>('I')] = residue_masses_[static_cast<unsigned char>('L')];
    for (const auto& m : fixed_mods)
      if (m.first >= 0) residue_masses_[static_cast<unsigned char>(m.first)] += m.second;

    ProteaseDigestion dig;
    dig.setEnzyme("Trypsin");
    dig.setMissedCleavages(missed_cleavages);

    for (const auto& e : entries)
    {
      std::vector<AASequence> parts;
      dig.digest(AASequence::fromString(e.sequence), parts);
      for (const AASequence& part : parts)
      {
        const std::string seq = fold(part.toUnmodifiedString());
        if (static_cast<int>(seq.size()) < k_) continue;
        // A residue outside the canonical 20 (X/B/Z/U/*) makes the prefix masses
        // meaningless; skip the whole peptide rather than mis-mass it.
        bool ok = true;
        for (char c : seq) if (residue_masses_[static_cast<unsigned char>(c)] <= 0.0) { ok = false; break; }
        if (!ok) continue;

        const uint32_t pid = static_cast<uint32_t>(peptides_.size());
        Peptide pep;
        pep.seq = seq;
        pep.protein = e.identifier;
        pep.prefix.resize(seq.size() + 1);
        pep.prefix[0] = 0.0;
        for (size_t i = 0; i < seq.size(); ++i)
          pep.prefix[i + 1] = pep.prefix[i] + residue_masses_[static_cast<unsigned char>(seq[i])];
        peptides_.push_back(std::move(pep));

        const std::string& s = peptides_.back().seq;
        for (size_t i = 0; i + static_cast<size_t>(k_) <= s.size(); ++i)
          index_[s.substr(i, static_cast<size_t>(k_))].push_back({pid, static_cast<uint32_t>(i)});
      }
    }
  }

  void TagReconciler::tryPlace(const Peptide& pep, size_t pos, const std::string& tag,
                               bool reversed, double nterm_mass, double cterm_mass,
                               std::vector<Reconciliation>& out) const
  {
    const size_t k = tag.size();
    const size_t L = pep.seq.size();
    const double db_n = pep.prefix[pos];               // residues N-terminal to the tag
    const double db_c = pep.prefix[L] - pep.prefix[pos + k];  // residues C-terminal to the tag

    // Under the y-ion assumption the tag's nterm_mass is the N-side flank and
    // cterm_mass the C-side. A reversed (b-derived) placement reads the peptide
    // C->N, so the two spectrum flanks swap relative to the peptide.
    const double want_n = reversed ? cterm_mass : nterm_mass;
    const double want_c = reversed ? nterm_mass : cterm_mass;

    const bool n_ok = std::fabs(db_n - want_n) <= tolAt(std::max(db_n, want_n));
    const bool c_ok = std::fabs(db_c - want_c) <= tolAt(std::max(db_c, want_c));

    // Only ONE flank may carry a mass gap: a peptide with mismatches on both
    // sides of the tag is under-constrained and, in TagRecon's own words, hurts
    // both speed and accuracy. Reject those.
    if (!n_ok && !c_ok) return;

    Reconciliation r;
    r.protein = pep.protein;
    r.peptide = pep.seq;
    r.pos = pos;
    r.reversed = reversed;
    r.nterm_match = n_ok;
    r.cterm_match = c_ok;
    if (n_ok && c_ok) { r.delta_mass = 0.0; r.region_lo = 0; r.region_hi = -1; }
    else if (!n_ok)   // gap is on the N-side residues [0, pos)
    {
      r.delta_mass = want_n - db_n;
      r.region_lo = 0;
      r.region_hi = static_cast<int>(pos) - 1;
    }
    else              // gap is on the C-side residues [pos+k, L)
    {
      r.delta_mass = want_c - db_c;
      r.region_lo = static_cast<int>(pos + k);
      r.region_hi = static_cast<int>(L) - 1;
    }
    out.push_back(std::move(r));
  }

  std::vector<Reconciliation> TagReconciler::reconcile(const std::string& tag_in,
                                                       double nterm_mass, double cterm_mass) const
  {
    std::vector<Reconciliation> out;
    if (static_cast<int>(tag_in.size()) != k_) return out;  // index is per exact k

    const std::string tag = fold(tag_in);
    auto it = index_.find(tag);
    if (it != index_.end())
      for (const Placement& pl : it->second)
        tryPlace(peptides_[pl.pep], pl.pos, tag, /*reversed=*/false, nterm_mass, cterm_mass, out);

    if (both_)
    {
      std::string rev(tag.rbegin(), tag.rend());
      if (rev != tag)
      {
        auto rit = index_.find(rev);
        if (rit != index_.end())
          for (const Placement& pl : rit->second)
            tryPlace(peptides_[pl.pep], pl.pos, rev, /*reversed=*/true, nterm_mass, cterm_mass, out);
      }
    }
    return out;
  }
}
