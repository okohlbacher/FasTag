// Stage-A tag reconciliation: place a tag on a database peptide by flanking
// masses and localize a single mass gap. Dataset-independent -- the "spectrum"
// flank masses are computed from the peptide itself, so a correct placement
// must reconcile with delta 0, and an injected gap must localize.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TagRecon.h"

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace FasTag;
using namespace OpenMS;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    if (!ok) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
  }

  double resMass(char c)
  {
    return ResidueDB::getInstance()->getResidue(String(c))->getMonoWeight(Residue::Internal);
  }
  double sumMass(const std::string& s)
  {
    double m = 0; for (char c : s) m += resMass(c); return m;
  }

  TagReconciler makeRecon(const std::string& protein_seq)
  {
    std::vector<FASTAFile::FASTAEntry> entries;
    entries.emplace_back("PROT1", "test protein", protein_seq);
    TagReconciler r(0.02, /*ppm=*/false, /*both=*/true);
    r.build(entries, /*k=*/4, /*missed=*/0, /*fixed_mods=*/{});
    return r;
  }
}

int main()
{
  // A tryptic peptide (ends in K, no internal K/R): trypsin yields it whole.
  // Flanked by tryptic boundaries so the digest produces it cleanly.
  const std::string pep = "SAMPLEVGATSDGK";
  TagReconciler r = makeRecon("MSTVWYAAR" + pep + "GTANLEFDSK");
  check(r.peptideCount() >= 1, "digestion produced peptides");

  // --- Exact placement: an internal 4-mer with its true flanking masses ---
  {
    const size_t pos = 4;                // tag = pep[4..8)
    const std::string tag = pep.substr(pos, 4);      // "LEVG"
    const double nflank = sumMass(pep.substr(0, pos));           // S A M P
    const double cflank = sumMass(pep.substr(pos + 4));          // A T S D G K
    auto res = r.reconcile(tag, nflank, cflank);
    bool exact = false;
    for (const auto& x : res)
      if (x.peptide == pep && x.pos == pos && x.nterm_match && x.cterm_match &&
          std::fabs(x.delta_mass) < 1e-6) exact = true;
    check(exact, "exact placement of an internal tag reconciles with delta 0");
  }

  // --- Modification gap on the C-side: +79.96633 (phospho) ---
  {
    const size_t pos = 4;
    const std::string tag = pep.substr(pos, 4);
    const double nflank = sumMass(pep.substr(0, pos));
    const double cflank = sumMass(pep.substr(pos + 4)) + 79.96633;  // a phospho somewhere C-side
    auto res = r.reconcile(tag, nflank, cflank);
    bool localized = false;
    for (const auto& x : res)
      if (x.peptide == pep && x.pos == pos && x.nterm_match && !x.cterm_match &&
          std::fabs(x.delta_mass - 79.96633) < 1e-3 &&
          x.region_lo == static_cast<int>(pos + 4) && x.region_hi == static_cast<int>(pep.size()) - 1)
        localized = true;
    check(localized, "a C-side +79.966 gap localizes to the C-flank region with the right mass");
  }

  // --- Modification gap on the N-side ---
  {
    const size_t pos = 6;
    const std::string tag = pep.substr(pos, 4);
    const double nflank = sumMass(pep.substr(0, pos)) + 42.01057;   // acetyl somewhere N-side
    const double cflank = sumMass(pep.substr(pos + 4));
    auto res = r.reconcile(tag, nflank, cflank);
    bool localized = false;
    for (const auto& x : res)
      if (x.peptide == pep && x.pos == pos && !x.nterm_match && x.cterm_match &&
          std::fabs(x.delta_mass - 42.01057) < 1e-3 &&
          x.region_lo == 0 && x.region_hi == static_cast<int>(pos) - 1)
        localized = true;
    check(localized, "an N-side +42.011 gap localizes to the N-flank region");
  }

  // --- Both flanks wrong: rejected (only one gap allowed) ---
  {
    const size_t pos = 5;
    const std::string tag = pep.substr(pos, 4);
    const double nflank = sumMass(pep.substr(0, pos)) + 16.0;
    const double cflank = sumMass(pep.substr(pos + 4)) + 16.0;
    auto res = r.reconcile(tag, nflank, cflank);
    bool any_this_pep = false;
    for (const auto& x : res) if (x.peptide == pep && x.pos == pos) any_this_pep = true;
    check(!any_this_pep, "a tag with both flanks off is rejected");
  }

  // --- Reversed (b-derived) placement: reversed tag, swapped flanks ---
  {
    const size_t pos = 4;
    const std::string tag = pep.substr(pos, 4);
    const std::string rev(tag.rbegin(), tag.rend());
    // A b-derived reading stores the tag reversed; its nterm/cterm flanks are the
    // peptide's C/N flanks respectively.
    const double nflank_spec = sumMass(pep.substr(pos + 4));   // becomes the tag's "nterm"
    const double cflank_spec = sumMass(pep.substr(0, pos));    // becomes the tag's "cterm"
    auto res = r.reconcile(rev, nflank_spec, cflank_spec);
    bool ok = false;
    for (const auto& x : res)
      if (x.peptide == pep && x.pos == pos && x.reversed && x.nterm_match && x.cterm_match)
        ok = true;
    check(ok, "a reversed tag reconciles with swapped flanks");
  }

  if (failures == 0) std::cout << "tagrecon_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
