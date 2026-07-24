// TaxIndex: k-mer -> taxon-set index over a reduced reference database.
// Dataset-independent: tiny synthetic proteins with known taxids, so the
// expected postings and lookups are exact.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT

#include "TaxIndex.h"

#include <cstdint>

#include <cstdio>
#include <fstream>
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

  // lookup() fills a caller buffer: the v2 image stores taxon INDICES, so there
  // is no vector in the file to hand back a reference to.
  std::vector<uint32_t> s;

  // AAAAAAAA is in taxa 10 and 20.
  idx.lookup("AAAAAAAA", s);
  check(s.size() == 2 && s[0] == 10 && s[1] == 20, "shared 8-mer maps to both taxa, sorted");
  // CDEFGHKL only in taxon 10.
  idx.lookup("CDEFGHKL", s);
  check(s.size() == 1 && s[0] == 10, "unique 8-mer maps to one taxon");
  // MNPQRSTV only in taxon 30.
  idx.lookup("MNPQRSTV", s);
  check(s == std::vector<uint32_t>{30}, "taxon-30 8-mer");
  {
    const std::string q = TaxIndex::fold("GATSDGKL");
    idx.lookup(q, s);
    check(s.size() == 1 && s[0] == 30, "internal 8-mer of p30a");
  }
  // Absent k-mer -> empty.
  idx.lookup("WWWWWWWW", s);
  check(s.empty(), "absent k-mer is empty");
  // A k-mer containing a residue outside the folded 19-letter alphabet has no
  // code at all and must not match anything.
  idx.lookup("AAAAAAAX", s);
  check(s.empty(), "k-mer with an ambiguous residue is empty");
  // Wrong length is rejected rather than silently truncated.
  idx.lookup("AAAA", s);
  check(s.empty(), "wrong-length query is empty");

  check(idx.postings(10) >= 1 && idx.postings(20) >= 1 && idx.postings(30) >= 1, "all taxa have postings");
  check(idx.nKmers() > 0, "index non-empty");
  {
    auto t = idx.taxa();
    check(t.size() == 3 && t[0] == 10 && t[1] == 20 && t[2] == 30, "taxa() lists all three, sorted");
  }

  // base-19 packing: distinct k-mers must get distinct codes, and anything with
  // an out-of-alphabet residue must be rejected.
  check(TaxIndex::encode("AAAAAAAA", 8) != TaxIndex::encode("AAAAAAAC", 8), "codes are distinct");
  check(TaxIndex::encode("AAAAAAAX", 8) == UINT64_MAX, "ambiguous residue has no code");
  check(TaxIndex::encode("AAAAAAAI", 8) == TaxIndex::encode("AAAAAAAL", 8), "I folds onto L in the code");
  check(TaxIndex::keyBytes(7) == 4 && TaxIndex::keyBytes(8) == 5, "key width follows 19^k");

  // Save / load round-trip -- this is the v2 packed+mapped path.
  {
    std::string path = std::string(std::tmpnam(nullptr)) + ".ftx2";
    std::string err;
    check(idx.save(path, &err), "save: " + err);
    TaxIndex idx2;
    check(idx2.load(path, &err), "load: " + err);
    check(idx2.k() == idx.k(), "k survives round-trip");

    std::vector<uint32_t> a, b;
    idx.lookup("AAAAAAAA", a); idx2.lookup("AAAAAAAA", b);
    check(a == b, "shared set survives round-trip");
    idx.lookup("MNPQRSTV", a); idx2.lookup("MNPQRSTV", b);
    check(a == b, "singleton set survives round-trip");
    idx2.lookup("WWWWWWWW", b);
    check(b.empty(), "absent k-mer still absent after round-trip");
    check(idx2.postings(10) == idx.postings(10), "postings survive round-trip");
    check(idx2.postings(20) == idx.postings(20), "postings(20) survives round-trip");
    check(idx2.nKmers() == idx.nKmers(), "kmer count survives round-trip");
    check(idx2.taxa() == idx.taxa(), "taxa survive round-trip");

    // Re-saving a MAPPED index must reproduce it exactly (the non-legacy save path).
    std::string path2 = std::string(std::tmpnam(nullptr)) + ".ftx2";
    check(idx2.save(path2, &err), "re-save mapped: " + err);
    TaxIndex idx3;
    check(idx3.load(path2, &err), "re-load: " + err);
    idx3.lookup("AAAAAAAA", b);
    check(a != b || true, "");
    idx.lookup("AAAAAAAA", a);
    check(a == b, "mapped re-save round-trips");
    check(idx3.nKmers() == idx.nKmers(), "mapped re-save keeps kmer count");
    std::remove(path.c_str());
    std::remove(path2.c_str());
  }

  // Corruption rejection: a v2 file must fail to LOAD, not crash, when
  // truncated or given a hostile header (the case an adversarial review found
  // could read past the mapping). Realistic because the ~1 GB index is a
  // downloaded asset, so a partial download is ordinary corruption.
  {
    std::string good = std::string(std::tmpnam(nullptr)) + ".ftx2";
    std::string err;
    check(idx.save(good, &err), "save for corruption test: " + err);

    // Read the good file's bytes.
    std::ifstream gi(good, std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(gi)), std::istreambuf_iterator<char>());
    gi.close();

    auto writeAndTryLoad = [&](const std::string& data, const char* what) {
      std::string p = std::string(std::tmpnam(nullptr)) + ".ftx2";
      { std::ofstream o(p, std::ios::binary); o.write(data.data(), static_cast<std::streamsize>(data.size())); }
      TaxIndex bad;
      std::string e;
      check(!bad.load(p, &e), std::string("rejects ") + what);
      std::vector<uint32_t> out;
      bad.lookup("AAAAAAAA", out);  // must be safe on a failed load, not a crash
      check(out.empty(), std::string("no lookup after rejected ") + what);
      std::remove(p.c_str());
    };

    if (bytes.size() > 100) writeAndTryLoad(bytes.substr(0, bytes.size() / 2), "a truncated file");
    writeAndTryLoad(bytes.substr(0, 40), "a header-only file");
    writeAndTryLoad("FTX2\x02\x00\x00\x00", "a 8-byte stub");
    writeAndTryLoad("not an index at all", "a non-index file");
    // Hostile header: huge n_kmers must be rejected by the bounds check, not wrap.
    std::string hostile = bytes.substr(0, 40);
    for (int i = 16; i < 24; ++i) hostile[i] = static_cast<char>(0xFF);  // n_kmers = UINT64_MAX
    writeAndTryLoad(hostile, "a hostile huge-count header");

    std::remove(good.c_str());
  }

  if (failures == 0) std::cout << "taxindex_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
