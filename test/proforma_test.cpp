// Checks the ProForma tag rendering (F13).
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#include "Proforma.h"

#include <iostream>
#include <string>

using FasTag::toProforma;

static int failures = 0;
static void check(bool ok, const std::string& what)
{
  if (!ok) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
}

int main()
{
  // Both flanks non-trivial: full [+n]-SEQ-[+c].
  check(toProforma("WPNVVNP", 86.0027, 156.0947) == "[+86.0027]-WPNVVNP-[+156.0947]",
        "both flanks");

  // A flank at the terminus (~0) is omitted rather than written as [+0.0000].
  check(toProforma("GPS", 0.0, 419.2481) == "GPS-[+419.2481]", "n-flank omitted");
  check(toProforma("GPS", 525.1847, 0.0) == "[+525.1847]-GPS", "c-flank omitted");
  check(toProforma("GPS", 0.0, 0.0) == "GPS", "both flanks omitted");

  // Inline residue mods pass through verbatim -- FasTag's X[Name] IS ProForma.
  check(toProforma("GS[Phospho]P", 100.0, 200.0) == "[+100.0000]-GS[Phospho]P-[+200.0000]",
        "inline mod preserved");

  // Sign: a (rounded) negative flank must read [-..], never [+-..].
  check(toProforma("AAA", -12.3456, 50.0) == "[-12.3456]-AAA-[+50.0000]", "negative flank sign");

  // The 0.5 Da threshold: just under is dropped, just over is kept.
  check(toProforma("AAA", 0.49, 50.0) == "AAA-[+50.0000]", "sub-threshold flank dropped");
  check(toProforma("AAA", 0.51, 50.0) == "[+0.5100]-AAA-[+50.0000]", "supra-threshold flank kept");

  if (failures == 0) std::cout << "proforma_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
