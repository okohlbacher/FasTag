// ProForma 2.0 rendering of a sequence tag, for interoperability (F13).
//
// A FasTag tag is an internal stretch of residues with a known mass to each
// peptide terminus but UNKNOWN flanking composition. ProForma has no first-class
// "mass gap of unknown sequence at a terminus", so the closest standard-conformant
// notation is a terminal mass tag: `[+<nterm>]-RESIDUES-[+<cterm>]`. Downstream
// tools that speak ProForma (e.g. for USI construction or proteoform display)
// can then consume a tag without a bespoke parser.
//
// The residue string is passed through verbatim: FasTag already writes variable
// modifications inline as `X[Name]`, which is exactly ProForma residue-mod
// syntax, so `S[Phospho]` needs no translation.
//
// This is deliberately a lossless-of-what-we-know rendering, not a claim that the
// tag is a fully characterised proteoform. The header documents the terminal-mass
// convention so a reader is not misled.
//
// Copyright (c) 2026 Oliver Kohlbacher and contributors
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstdio>
#include <string>

namespace FasTag
{
  /// ProForma string for a tag: `[+<nterm>]-<residues>-[+<cterm>]`.
  ///
  /// A terminal mass is emitted only when it is non-trivial (> 0.5 Da): a tag
  /// that sits at the very peptide terminus has essentially no flank there, and
  /// `[+0.0000]-` would be noise. Masses are signed (`[+..]` / `[-..]`) as
  /// ProForma requires; FasTag's flanks are positive in practice, but a rounded
  /// negative is still written correctly rather than as a malformed `[+-1]`.
  inline std::string toProforma(const std::string& residues, double nterm, double cterm)
  {
    auto massTag = [](double m) {
      char b[48];
      std::snprintf(b, sizeof b, "[%+.4f]", m);   // %+ gives the leading + or -
      return std::string(b);
    };
    std::string out;
    if (std::fabs(nterm) > 0.5) out += massTag(nterm) + "-";
    out += residues;
    if (std::fabs(cterm) > 0.5) out += "-" + massTag(cterm);
    return out;
  }
}
