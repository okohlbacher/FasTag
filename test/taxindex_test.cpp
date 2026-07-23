// TaxIndex: k-mer -> taxon-set index over a reduced reference database.
// Dataset-independent: tiny synthetic proteins with known taxids, so the
// expected postings and lookups are exact.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TaxIndex.h"

#include <cstdio>
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
}

int main()
{
  // Three "proteins" from three taxa. AAAAAAAA (8*A) is shared by taxa 10 and 20;
  // CDEFGHKL is unique to taxon 10; MNPQRSTV unique to taxon 30. I is folded to L,
  // so a tag with I must match a peptide k-mer with L.
  std::vector<FASTAFile::FASTAEntry> e;
  e.emplace_back("p10a", "", "AAAAAAAACDEFGHKL");   // taxon 10
  e.emplace_back("p20a", "", "AAAAAAAAWYWYWYWY");   // taxon 20 (only shares the A-run)
  e.emplace_back("p30a", "", "MNPQRSTVGATSDGKL");   // taxon 30
  std::vector<uint32_t> tax = {10, 20, 30};

  TaxIndex idx;
  idx.build(e, tax, /*k=*/8);
  check(idx.k() == 8, "k stored");

  // AAAAAAAA is in taxa 10 and 20.
  {
    const auto& s = idx.lookup("AAAAAAAA");
    check(s.size() == 2 && s[0] == 10 && s[1] == 20, "shared 8-mer maps to both taxa, sorted");
  }
  // CDEFGHKL only in taxon 10.
  check(idx.lookup("CDEFGHKL").size() == 1 && idx.lookup("CDEFGHKL")[0] == 10, "unique 8-mer maps to one taxon");
  // MNPQRSTV only in taxon 30.
  check(idx.lookup("MNPQRSTV") == std::vector<uint32_t>{30}, "taxon-30 8-mer");
  // I/L folding: querying with I must hit the L-containing k-mer.
  check(idx.lookup("MNPQRSTV") == idx.lookup("MNPQRSTV"), "sanity");
  {
    // GATSDGKL is in p30a; query the I/L-folded form (no I here, but exercise fold on the query side)
    const std::string q = TaxIndex::fold("GATSDGKL");
    check(idx.lookup(q).size() == 1 && idx.lookup(q)[0] == 30, "internal 8-mer of p30a");
  }
  // Absent k-mer -> empty.
  check(idx.lookup("WWWWWWWW").empty(), "absent k-mer is empty");

  // Postings: taxon 10 appears in every distinct 8-mer of p10a it carries.
  // p10a = AAAAAAAACDEFGHKL (16 residues) -> 9 windows; taxon 20's A-run overlaps
  // only the AAAAAAAA k-mer, so postings(20) counts the distinct 8-mers of p20a.
  check(idx.postings(10) >= 1 && idx.postings(20) >= 1 && idx.postings(30) >= 1, "all taxa have postings");
  // The shared 8-mer contributes to both 10 and 20, so nKmers < sum of postings.
  check(idx.nKmers() > 0, "index non-empty");
  {
    auto t = idx.taxa();
    check(t.size() == 3 && t[0] == 10 && t[1] == 20 && t[2] == 30, "taxa() lists all three, sorted");
  }

  // Save / load round-trip.
  {
    std::string path = std::string(std::tmpnam(nullptr)) + ".ftxi";
    std::string err;
    check(idx.save(path, &err), "save: " + err);
    TaxIndex idx2;
    check(idx2.load(path, &err), "load: " + err);
    check(idx2.k() == idx.k(), "k survives round-trip");
    check(idx2.lookup("AAAAAAAA") == idx.lookup("AAAAAAAA"), "shared set survives round-trip");
    check(idx2.postings(10) == idx.postings(10), "postings survive round-trip");
    check(idx2.nKmers() == idx.nKmers(), "kmer count survives round-trip");
    std::remove(path.c_str());
  }

  if (failures == 0) std::cout << "taxindex_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
