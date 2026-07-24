// Build a FasTag taxonomic k-mer index from a directory of per-taxon FASTA files.
//
// Each input file is named <taxid>.fasta and holds one organism's proteins; the
// taxid is taken from the filename, so no header parsing is needed. All proteins
// are pooled and handed to TaxIndex::build, which records, for every I/L-folded
// k-mer, the SET of taxa whose proteins contain it (breadth, not abundance). The
// result is a compact index the classifier loads with -taxdb.
//
// Deliberately a reduced reference set (one representative proteome per
// genus/family), not complete NR: short peptide k-mers cannot resolve species,
// and a curated set keeps the LCA honest at genus/family level while staying
// small enough to index in bounded RAM.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include <OpenMS/FORMAT/FASTAFile.h>

#include "TaxIndex.h"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace OpenMS;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: buildtaxdb <proteomes_dir> <out.taxdb> [k=7]\n"
                 "  <proteomes_dir> holds files named <taxid>.fasta\n";
    return 2;
  }
  const std::string dir = argv[1];
  const std::string out = argv[2];
  const int k = argc > 3 ? std::atoi(argv[3]) : 7;
  if (k < 3 || k > 30) { std::cerr << "k out of range (3..30)\n"; return 2; }

  std::vector<FASTAFile::FASTAEntry> entries;
  std::vector<uint32_t> taxids;
  size_t files = 0;
  for (const auto& e : fs::directory_iterator(dir))
  {
    if (!e.is_regular_file()) continue;
    const std::string path = e.path().string();
    const std::string stem = e.path().stem().string();  // "<taxid>"
    if (e.path().extension() != ".fasta") continue;
    char* end = nullptr;
    const unsigned long tx = std::strtoul(stem.c_str(), &end, 10);
    if (end == stem.c_str() || tx == 0)
    {
      std::cerr << "  skip (filename is not <taxid>.fasta): " << path << "\n";
      continue;
    }
    std::vector<FASTAFile::FASTAEntry> part;
    try { FASTAFile().load(path, part); }
    catch (const std::exception& ex) { std::cerr << "  FASTA load failed: " << path << ": " << ex.what() << "\n"; continue; }
    for (auto& fe : part) { entries.push_back(std::move(fe)); taxids.push_back(static_cast<uint32_t>(tx)); }
    ++files;
    std::cerr << "  taxid " << tx << ": " << part.size() << " proteins (" << stem << ")\n";
  }
  if (entries.empty()) { std::cerr << "no proteins loaded from " << dir << "\n"; return 1; }

  std::cerr << "building k=" << k << " index over " << entries.size()
            << " proteins from " << files << " taxa...\n";
  FasTag::TaxIndex idx;
  idx.build(entries, taxids, k);
  std::cerr << "  " << idx.nKmers() << " distinct k-mers, " << idx.taxa().size() << " taxa\n";

  std::string err;
  if (!idx.save(out, &err)) { std::cerr << "save failed: " << err << "\n"; return 1; }
  std::cerr << "wrote " << out << "\n";
  return 0;
}
