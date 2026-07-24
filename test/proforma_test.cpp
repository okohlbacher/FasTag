// Checks the ProForma tag rendering (F13).
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#include "Proforma.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using FasTag::toProforma;

static int failures = 0;
static void check(bool ok, const std::string& what)
{
  if (!ok) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
}

int main()
{
  // Both flanks; note L is rendered J (FasTag folds I->L, so L means "I or L").
  check(toProforma("WPNVVNP", 86.0027, 156.0947) == "[+86.0027]-WPNVVNP-[+156.0947]",
        "both flanks, no L");
  check(toProforma("GLPSL", 100.0, 200.0) == "[+100.0000]-GJPSJ-[+200.0000]", "L rendered as J");

  // A flank at the terminus (~0) is omitted, not written as [+0.0000].
  check(toProforma("GPS", 0.0, 419.2481) == "GPS-[+419.2481]", "n-flank omitted");
  check(toProforma("GPS", 525.1847, 0.0) == "[+525.1847]-GPS", "c-flank omitted");
  check(toProforma("GPS", 0.0, 0.0) == "GPS", "both flanks omitted");

  // Global fixed-mod prefix (the default Carbamidomethyl C case).
  check(toProforma("GCS", 100.0, 200.0, "<[Carbamidomethyl]@C>")
            == "<[Carbamidomethyl]@C>[+100.0000]-GCS-[+200.0000]",
        "global fixed-mod prefix");

  // Inline variable mods pass through; the L INSIDE a mod token is not folded,
  // only residue-letter L is.
  check(toProforma("GS[Phospho]P", 100.0, 200.0) == "[+100.0000]-GS[Phospho]P-[+200.0000]",
        "inline mod preserved");
  check(toProforma("L[Foo]", 0.0, 0.0) == "J[Foo]", "residue L folded, mod name untouched");

  // Sign and threshold.
  check(toProforma("AAA", -12.3456, 50.0) == "[-12.3456]-AAA-[+50.0000]", "negative flank sign");
  check(toProforma("AAA", 0.49, 50.0) == "AAA-[+50.0000]", "sub-threshold flank dropped");
  check(toProforma("AAA", 0.51, 50.0) == "[+0.5100]-AAA-[+50.0000]", "supra-threshold flank kept");

  // Non-finite / absurd masses are skipped, never emitted as [+inf] or garbage.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  check(toProforma("AAA", inf, 50.0) == "AAA-[+50.0000]", "infinite flank skipped");
  check(toProforma("AAA", nan, 50.0) == "AAA-[+50.0000]", "NaN flank skipped");
  check(toProforma("AAA", 1e13, 50.0) == "AAA-[+50.0000]", "absurd flank skipped");

  if (failures == 0) std::cout << "proforma_test: all checks passed\n";
  return failures == 0 ? 0 : 1;
}
