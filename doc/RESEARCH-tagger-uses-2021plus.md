# How sequence tagging is used and extended in the MS literature (2021+) — feature implications for FasTag

**Scope.** Adversarial survey of peer-reviewed and preprint literature published 2021 or later on how de novo *sequence tags* (partial sequences with flanking masses, from MS/MS) are used and extended, read specifically as a feature backlog for **FasTag** — a fast, parallel, open-source DirecTag reimplementation that emits scored partial tags with flanking masses, supports gaps, extension, mods, and mzML/mzPeak input, and is being extended with (i) tag→database reconciliation, (ii) modification/mutation localization, and (iii) a taxonomic/species detector.

**Method / honesty note.** Sources were found via web search, PubMed, and journal sites and cross-checked to a DOI or working URL. Every citation below was seen in a search result or landing page during this research. Where a paper is a *preprint only* or the peer-reviewed version's DOI could not be fully confirmed, it is flagged inline. I did not read every full text; claims about internal algorithm details are drawn from abstracts, tool docs, and search snippets and are marked where load-bearing. Dates are as reported by the venue.

Two constraints carried in from FasTag's own design memory shape the verdicts:
- **Tagging is a sensitive prefilter.** Later stages supply specificity, so features that trade recall for tag specificity are viewed skeptically.
- **Synthetic evidence cannot flip shipped defaults.** "Adopt" below means *build the capability*, not *make it the default* without real-data benchmarks.

---

## (a) Landscape overview: what tag-based MS looks like in 2021+

The field has moved in four directions since 2020, and all four are relevant to a tagger.

**1. Tags survived the "fragment-ion index" disruption — as a complement, not a casualty.** The dominant open/blind-search accelerator is now the *fragment-ion index* (MSFragger 2017; Sage 2023; Comet fragment-ion indexing 2025), which indexes theoretical fragments rather than extracting sequence tags. But tag-based search did not die: **PIPI2 (2024)** and **Open-pFind** keep the tag-extract-then-reconcile paradigm precisely because it degrades gracefully when a peptide carries *multiple* PTMs or a mutation, where a pure precursor-mass-offset open search struggles. The modern framing: a tag is a *robust anchor* that survives unexpected mass changes anywhere outside the tag, and the flanking masses become the search's degrees of freedom. This is exactly FasTag's reconciliation thesis.

**2. Deep-learning de novo (Casanovo, pNovo3, Novor successors) largely replaced hand-scored full-length de novo — but "tags" reappeared as the reliable *sub-string* of an unreliable full read.** Antibody assembly (Stitch, ALPS), metaproteomics taxonomy (MegaPX, NovoLign), and immunopeptidomics (MARS) all take *de novo output and use only its confident stretches*. A tagger that emits calibrated, confidence-scored partial reads is doing natively what these pipelines extract post hoc from a full de novo sequence.

**3. Rescoring became standard, and tag-like agreement features feed it.** MS2Rescore/Oktoberfest/Percolator pipelines (2021-2024) boost IDs 10-30%, especially in immunopeptidomics and hard search spaces, by adding features such as predicted-vs-observed fragment agreement. A tagger can supply "longest de novo agreement with the candidate peptide" as a first-class rescoring feature.

**4. Real-time / instrument-control search matured (RTS on Tribrid; real-time library search for XL-MS; Comet-FI explicitly targeting real-time).** Latency, not just throughput, is now a design axis. A *fast* tagger has a differentiated niche here that a database search engine does not: tags can be computed with no database at all, on the acquisition timescale.

**Application areas surveyed:** open/blind PTM search; error-tolerant/variant/proteogenomics; metaproteomics & taxonomic profiling; immunopeptidomics/HLA; cross-linking MS; de novo+tag hybrids and antibody assembly; DIA; single-cell; top-down; glycoproteomics; QC/rescoring; chimeric/multiplexed spectra; real-time acquisition.

---

## (b) Candidate features for FasTag

Each entry: **application driving it → capability the tagger must provide → adversarial verdict**. Citations are numbered to the reference list in section (d).

### F1. Tag + flanking-mass reconciliation with the mass-delta interpreted as modification/mutation ("blind" reconciliation)
- **Driving application:** open/blind PTM search; multi-PTM peptides; variant/mutation search.
- **Evidence:** PIPI2 identifies peptides with up to four PTMs by extracting tags, indexing the database, extending matched tags, and localizing mass shifts — 30-88% more modified peptides than competitors on real data [1]. Open-pFind uses the same extract-tag → index → extend → localize-shift loop and is a mainstream engine [2 background; 3 localization-aware open search]. TagRecon originated exactly this "mass difference between peptide and tag = PTM/mutation" idea [4 background]. TagGraph scaled unrestricted tag-string search to 25M spectra, tripling IDs [5 background].
- **Capability needed:** flanking (N-side/C-side) masses emitted *and* carried through reconciliation; a database index keyed by tag substrings; a delta-attribution step that can split a delta across two flanks and test candidate PTM/mutation masses.
- **Verdict: ADOPT (this is FasTag's core reconciliation extension; the literature validates the exact design).** The differentiation vs. fragment-index open search is real but *narrow*: it wins specifically on multi-PTM and mutation-bearing peptides [1]. Cost: moderate — needs a substring-searchable protein index (see F3) and delta bookkeeping. Do **not** market it as "faster than MSFragger for ordinary open search" — the fragment index wins there.

### F2. Substring/FM-index for tag→database lookup (the reconciliation index)
- **Driving application:** everything in F1; speed of reconciliation at proteome scale.
- **Evidence:** PIPI2 builds tag indexes on an **FM-index** for rapid database retrieval [1]. Fragment-ion indexing (MSFragger, Sage, Comet-FI) shows the general lesson that indexing is what makes open/large searches tractable — Comet-FI cut search time up to 94% [6,7,8].
- **Capability needed:** an FM-index / suffix-array (or amino-acid k-mer index) over the protein database supporting gapped/mismatch-tolerant tag queries.
- **Verdict: ADOPT.** Without this, reconciliation (F1) does not scale. FM-index is the proven choice for tag lookup specifically (as opposed to fragment indexing, which is for spectrum-centric search). Cost: real engineering, but bounded and well-trodden. Reuse an existing succinct-index library rather than hand-rolling.

### F3. Gap tags / mass-gap residues (allow a 1-2 residue "gap" spanned by a mass)
- **Driving application:** noisy/incomplete fragmentation everywhere; immunopeptidomics (non-tryptic, sparse ladders); single-cell (low ion counts).
- **Evidence:** DirecTag itself allowed this; modern de novo tools handle missing peaks by mass gaps. Single-cell and immunopeptidomics spectra are ion-starved [12,13], so unbroken 3-mers are often unavailable. The "tag is a sensitive prefilter" principle argues for maximizing recall via gaps.
- **Capability needed:** tag model where an edge can be a single amino acid *or* a mass equal to a 2-residue combination, flagged as ambiguous; downstream reconciliation must accept gapped tags.
- **Verdict: ADOPT (recall lever, cheap).** Gaps directly raise the fraction of spectra that yield *any* usable tag, which is the whole point of a sensitive prefilter. Cost: low-moderate; increases candidate count, so pair with confidence scoring (F6) to keep the flood manageable. Risk: over-gapping explodes false tags — cap gap count per tag (e.g., ≤1).

### F4. Calibrated tag confidence + FDR/probability on tags (not just a raw score)
- **Driving application:** rescoring (features must be calibrated), taxonomic profiling (need to trust a tag before assigning a taxon), variant QC (false variant peptides are the field's plague).
- **Evidence:** DirecTag's original contribution was *statistical* tag scoring [background, and the reason it beat GutenTag/InsPecT]. Variant-peptide QC tools exist precisely because tag/de novo-derived variants are false-positive-prone — PgxSAVy flags spurious variant peptides [9]. NovoBoard and de novo FDR work (Winnow) show the field now demands calibrated confidence on de novo output [10 background; 11 preprint].
- **Capability needed:** per-tag posterior/e-value or empirical FDR, ideally per-residue confidence within the tag; a decoy or calibration model for tags.
- **Verdict: ADOPT.** This is the connective tissue for F7 (rescoring) and F8 (taxonomy) and it is *differentiating* — most fast taggers emit uncalibrated scores. Cost: moderate; needs a calibration/decoy scheme for tags specifically. High leverage: a well-calibrated tag confidence is reusable across every downstream.

### F5. PTM-aware tags + mass-shift localization within/around the tag
- **Driving application:** open PTM search; the mod-localization extension FasTag is already building.
- **Evidence:** Localization-aware open search indexes shifted *and* unshifted fragments to place the mod [3]. PTM-Shepherd summarizes/localizes open-search mass shifts and finds diagnostic ions [14, 15]. SpecGlobX interprets a PSM's mass delta as one or several unlocalized mods via spectral alignment — a natural post-tag localization aid [16].
- **Capability needed:** ability to place a known/unknown delta at a residue by testing shifted fragment support; emit localization scores and diagnostic-ion flags; represent modified tags in a standard grammar (ProForma; see F13).
- **Verdict: ADOPT for localization; MAYBE for full open-mod *profiling*.** Localization is squarely in FasTag's stated roadmap and the shifted-fragment method [3] is proven. But do **not** reimplement PTM-Shepherd-style aggregate mass-shift profiling — interoperate with it instead (emit deltas it can consume). Cost: localization is moderate; profiling would be scope creep.

### F6. Multi-length tag emission with length as a tunable recall/specificity knob
- **Driving application:** tag-based DB search (longer tags prune the candidate list); immunopeptidomics and non-specific search (huge search space needs specific anchors); metaproteomics (specificity buys taxonomic resolution).
- **Evidence:** JUMP showed tag length + pattern-match scoring drives sensitivity/specificity in tag-based search [17 background, 2014]. Metaproteomics taxonomy resolution depends on peptide/tag specificity as UniProt grows and LCA specificity drops [18]. Tag-based open search (PIPI2) tunes tag handling for candidate reduction [1].
- **Capability needed:** emit tags across a length range (e.g., 3-7) with scores; let downstream pick the operating point rather than hard-coding length 3.
- **Verdict: ADOPT, but resist the temptation to default to long tags.** Per the "sensitive prefilter" principle, the *default* should favor short/recall-y tags; longer tags are an opt-in for specificity-hungry consumers (taxonomy, non-specific HLA search). Cost: low — mostly a matter of not truncating the DAG traversal early. Adversarial caution: long tags in ion-starved spectra collapse recall; expose the knob, don't preset it high.

### F7. Emit tag-agreement features for ML rescorers (MS2Rescore / Oktoberfest / Percolator interoperability)
- **Driving application:** QC/rescoring — now standard and worth 10-30% more IDs, biggest gains in immunopeptidomics and hard spaces.
- **Evidence:** MS2Rescore 3.0 is explicitly modular and aimed at proteogenomics/metaproteomics/immunopeptidomics [19]; the original MS2Rescore showed dramatic immunopeptide gains [20 preprint→MCP 2022]. Rescoring benchmarks show prediction-based features give the largest gains and harmonize engines [21, 22]. A tag→PSM "longest agreement" feature is a natural, cheap addition to this feature stack.
- **Capability needed:** a stable output where, given a candidate peptide, FasTag reports longest matched tag, its confidence, #residues of de novo agreement, and flank-mass consistency — as columns a rescorer can ingest.
- **Verdict: ADOPT (high ROI, low cost, differentiating).** This turns FasTag into a *feature generator* for the dominant rescoring ecosystem rather than a standalone engine competing with fragment indexes. Cost: low. The main work is agreeing on an output schema (F13). Evidence it helps *specifically as a tag feature* is indirect (rescoring works; tag-agreement is a plausible feature), so validate on real data before claiming a number.

### F8. De novo tag → taxonomic/species classification (alignment-free k-mer or fast alignment)
- **Driving application:** metaproteomics taxonomic profiling; the species detector FasTag is building; also food/forensic species ID.
- **Evidence:** **MegaPX** classifies de novo peptides against large databases with an alignment-free, k-mer, interleaved-Bloom-filter multi-index — built for metaproteomics taxonomy and virus assignment, with mutation-tolerant reference generation [23, published Bioinformatics 2026 / bioRxiv 2025]. **NovoLign** does metaproteomics taxonomy by large-scale sequence alignment of de novo reads and also *evaluates database-search coverage* [24]. **MetaNovo** uses de novo tag matching to build reduced tailored databases without prior composition knowledge [25]. Unipept's LCA approach shows specificity erosion as databases grow — motivating tag-length/confidence control [18].
- **Capability needed:** emit clean, confidence-filtered tags/short reads suitable for k-mer indexing or alignment; ideally a built-in LCA over tag→protein→taxon hits; mutation-tolerant matching (F1) to survive strain variation.
- **Verdict: ADOPT — this is arguably FasTag's most differentiated capability.** No mainstream *tagger* ships a taxonomic detector; the tools that do (MegaPX, NovoLign) bolt onto external de novo engines. Doing it natively, fast, and parallel is a genuine niche. Cost: moderate-high (need a taxon-annotated index + LCA + calibration). Adversarial caveats: (1) taxonomic calls are only as good as tag confidence (F4) and length (F6) — garbage tags → confident-but-wrong taxa; (2) benchmark against MegaPX/NovoLign honestly, and (3) per memory, do not ship aggressive taxonomic defaults on synthetic data.

### F9. Error-tolerant / single-substitution ("mutation") tags for proteogenomics
- **Driving application:** variant/proteogenomics search; strain variation in metaproteomics.
- **Evidence:** TagRecon's whole purpose was mutation identification via tagging [4 background]. MegaPX explicitly generates mutated reference databases for error-tolerant search [23]. MutCombinator handles combinatorial mutations via variant-graph search [26 background, 2020]. Variant peptides are false-positive-prone (PgxSAVy) [9] — so this must be paired with strong QC.
- **Capability needed:** reconciliation (F1) that allows one substituted residue inside the flank region (mass delta = residue-swap mass), distinct from a PTM delta; output the implied substitution.
- **Verdict: MAYBE → lean ADOPT, but gate behind QC.** The capability falls out almost for free once F1 exists (a substitution is just a constrained delta). The risk is entirely downstream false positives [9]; FasTag should *emit candidate substitutions with a clear "unvalidated variant" flag* and not present them as calls. Cost: low incremental over F1. Do not enable by default in generic runs.

### F10. Low-latency single-spectrum tagging for real-time / instrument-control search
- **Driving application:** real-time search (RTS) on Tribrid; real-time library search for XL-MS; Comet-FI explicitly targeting real-time; multiplexed single-cell RTS.
- **Evidence:** RTS-assisted acquisition improves multiplexed single-cell coverage [27]. Real-time library search increases XL-MS depth by triggering high-res MS2 only on likely cross-links [28]. Comet-FI targets real-time searches as a use case [7,8]. A database-free tagger has a structural latency advantage: no index round-trip needed to decide "is this spectrum worth a second scan / worth keeping."
- **Capability needed:** per-spectrum tagging in low single-digit ms; a stable C API / streaming interface; deterministic, bounded work per spectrum (cap gaps/branches).
- **Verdict: MAYBE (promising, but prove the latency first).** This is a genuinely differentiated niche for a *fast parallel* tagger and aligns with FasTag's identity. But it is speculative: no current workflow triggers acquisition on a *tag* today (they use RTS DB matches or library scores). Cost: low to *expose* (streaming API), high to *land in an instrument* (vendor integration, out of scope). Recommendation: build the fast streaming path and a benchmark; treat instrument integration as a research collaboration, not a shipped feature.

### F11. Glyco-aware tagging (peptide-backbone tags under a glycan; Y-ion/oxonium awareness)
- **Driving application:** glycoproteomics.
- **Evidence:** GlycanFinder integrates peptide-based and glycan-based search + de novo and handles the complex Y-ion/oxonium fragmentation [29]. pGlycoNovo does site-specific de novo glycan sequencing incl. unexpected fragments [30]. Both show the peptide backbone can be tagged if oxonium/Y-ions are recognized and the glycan mass treated as a large flank delta.
- **Capability needed:** recognize oxonium/diagnostic glycan ions to *classify* a spectrum as glyco; treat the glycan as an outsized flank mass; optionally emit the peptide-only tag.
- **Verdict: SKIP for v1; MAYBE later as a "glyco spectrum flag."** Full glyco is a specialized subfield with entrenched tools [29,30]; competing is a dead end. But cheap wins exist: an oxonium-ion *detector* that flags glyco spectra and emits the backbone tag with the glycan as a flank delta would let FasTag *route* rather than *solve* glyco. Cost: the flag is low; real glyco sequencing is high and not differentiated. Adopt only the flag, and only if a glyco consumer exists.

### F12. Multiple tags per spectrum for chimeric / multiplexed spectra (and DIA)
- **Driving application:** chimeric DDA spectra (the majority are chimeric); DIA multiplexed windows; TMT.
- **Evidence:** CHIMERYS shows most MS2 are chimeric and deconvolving them ~doubles IDs [31]. DIA de novo (Cascadia; transformer-DIA) tackles multiplexed spectra directly [32, 33 arXiv]. A tagger that returns *the top-N tag paths* rather than one already partially addresses chimericity by not committing to a single precursor's ladder.
- **Capability needed:** return ranked multiple tag sets per spectrum; ideally precursor-aware grouping for DIA (fragment→precursor association).
- **Verdict: MAYBE (top-N tags: yes; full DIA disentangling: skip).** Emitting several independent tag paths per spectrum is cheap and helps chimeric DDA and metaproteomics — adopt that. Full DIA precursor-fragment deconvolution [32,33] is a research problem owned by dedicated tools; do not chase it. Cost: top-N is low; DIA deconvolution is very high. Adversarial note: top-N tags multiply downstream candidates — again lean on F4/F6 to keep FDR sane.

### F13. Standard, self-describing tag output (ProForma mods, mzTab/mzIdentML, USI)
- **Driving application:** *all* of the above — reconciliation, rescoring, taxonomy, QC — need interoperable I/O.
- **Evidence:** The rescoring ecosystem (MS2Rescore, Oktoberfest) and search engines interoperate via community formats; the field standardized ProForma for modified peptidoforms and USI for spectrum references. There is no universal *tag* exchange format, which is itself an opportunity.
- **Capability needed:** emit tags with flank masses, gaps, per-residue confidence, and localized mods in a documented schema; use ProForma for modified residues; reference spectra by USI; provide mzTab/mzIdentML export for reconciled PSMs.
- **Verdict: ADOPT (unglamorous, force-multiplying).** Interoperability is what makes F7 (rescoring features) and F8 (taxonomy handoff) actually usable by others. Cost: low-moderate, mostly design discipline. Proposing a clean open tag format could become a small standards contribution. This is the cheapest high-leverage item after F4.

### F14. Overlap-friendly long-read output for antibody / repertoire de novo assembly
- **Driving application:** monoclonal/polyclonal antibody de novo sequencing and repertoire profiling.
- **Evidence:** Stitch maps de novo *short reads* to templates for mAb/pAb reconstruction [34]; ALPS assembles long peptide sequences from de novo reads + quality scores into a de Bruijn graph; Casanovo+ALPS reached ~94% coverage [35 review]. These pipelines are fundamentally *tag/short-read assemblers* and want confident overlapping partial reads with per-residue quality.
- **Capability needed:** emit long, high-confidence contiguous tags with per-residue quality and positional flank masses, in a form assemblers can ingest.
- **Verdict: MAYBE (low-cost alignment with existing consumers).** FasTag doesn't need to *assemble* antibodies, but if F4 (per-residue confidence) and F13 (clean output) exist, feeding Stitch/ALPS is nearly free and opens a real user base. Cost: low if F4+F13 done. Don't build an assembler.

### F15. Cross-link diagnostic-ion triggering / XL tagging
- **Driving application:** cross-linking MS.
- **Evidence:** Real-time library search for XL triggers a second scan on likely cross-links using diagnostic peak ratios [28]. Modern XL search is dominated by dedicated engines (pLink, MeroX, XlinkX, ECL) and deep-learning tools.
- **Capability needed:** detect cleavable-crosslinker signature/reporter ions; emit two tags (one per linked peptide) with the linker mass as a shared constraint.
- **Verdict: SKIP.** XL search is a crowded, specialized space with strong incumbents; a tagger adds little. The only cheap, defensible piece is a *signature-ion detector* to flag XL spectra — and even that has limited value without an XL-specific downstream. Not worth it for v1.

---

## (c) Cross-cutting adversarial observations

- **The fragment-ion index is the competitor, not other taggers.** FasTag should not claim to beat MSFragger/Sage/Comet-FI on ordinary open search — they'll win on speed there [6,7,8]. FasTag's defensible ground is: (i) database-free operation (taxonomy F8, real-time F10), (ii) multi-PTM/mutation reconciliation where fragment-offset search struggles (F1) [1], and (iii) being a calibrated *feature/tag generator* for rescoring and assembly (F7, F14).
- **Almost every feature depends on two primitives: calibrated tag confidence (F4) and clean interoperable output (F13).** Build those first; they unlock the rest and are individually differentiating.
- **Recall-first, per the prefilter principle.** Multiple candidate features (F6 long tags, F9 mutations, F12 top-N) tempt you to increase specificity or candidate count. Keep FasTag recall-biased and push specificity/FDR to the confidence layer and downstream, consistent with FasTag's own memory.
- **Variant and taxonomic outputs are false-positive traps.** [9] and [18] are cautionary. Emit these as flagged, unvalidated candidates; never as calls; never as synthetic-data-justified defaults.

## (d) Top 5 features to prioritize

1. **F4 — Calibrated tag confidence + per-residue/FDR.** Foundational. Every differentiated downstream (rescoring, taxonomy, variant flagging) is only as trustworthy as this. Differentiating because most fast taggers emit raw scores. Build first.
2. **F1 + F2 — Flank-mass reconciliation on an FM/substring index.** FasTag's stated core extension, directly validated by PIPI2 and Open-pFind [1,2], and the capability that justifies a tagger existing alongside fragment-index engines. F2 is the enabling substrate.
3. **F8 — Native de novo-tag taxonomic/species detector.** The single most *differentiated* capability: no mainstream tagger ships this, and the tools that do taxonomy (MegaPX, NovoLign [23,24]) bolt onto external de novo. A fast parallel native detector is a real niche — provided F4/F6 gate the tag quality.
4. **F7 + F13 — Rescoring-feature output in a standard, self-describing format.** Cheap, high-ROI, force-multiplying: plugs FasTag into the dominant MS2Rescore/Oktoberfest ecosystem [19,21] and every downstream consumer, instead of competing head-on. Low cost once F4 exists.
5. **F5 — Mass-shift localization (shifted-fragment method).** On FasTag's roadmap, proven approach [3], and complements F1: reconciliation finds *what* changed; localization says *where*. Interoperate with PTM-Shepherd [14] for profiling rather than reimplementing it.

Deliberately deprioritized: full glyco (F11), full DIA disentangling (F12b), XL search (F15) — crowded specialist spaces where a general tagger adds little. Keep as spectrum *flags/routers* at most.

---

## References

Foundational (pre-2021, context only): DirecTag — Tabb et al., *J. Proteome Res.* 2008, https://pubs.acs.org/doi/10.1021/pr800154p ; TagRecon — Dasari et al., *J. Proteome Res.* 2010, https://pubs.acs.org/doi/10.1021/pr900850m ; JUMP — Wang et al., *Mol. Cell. Proteomics* 2014 (note: 2014, not 2020), DOI 10.1074/mcp.M114.039743, https://www.mcponline.org/article/S1535-9476(20)33787-7/fulltext ; TagGraph — Devabhaktuni et al., *Nat. Biotechnol.* 2019, https://www.nature.com/articles/s41587-019-0067-5 ; MutCombinator — Na et al., 2020, https://pmc.ncbi.nlm.nih.gov/articles/PMC7355298/ .

1. **PIPI2: Sensitive Tag-Based Database Search to Identify Peptides with Multiple Post-translational Modifications** — Lai, Zhao, Zhou, Li, Yu. *J. Proteome Res.* 2024, 23(6):1960-1969. DOI 10.1021/acs.jproteome.3c00819. https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.3c00819 (PMID 38770571). *Tag extraction + FM-index reconciliation + multi-PTM localization.*
2. **Open-pFind** (Chi et al., *Nat. Biotechnol.* 2018) — mainstream tag-extract→index→extend→localize engine. Background; DOI not re-verified here, see preprint https://www.biorxiv.org/content/10.1101/285395.full.pdf . *Flag: peer-reviewed DOI not confirmed in this session.*
3. **Identification of modified peptides using localization-aware open search** — Nesvizhskii lab, *Nat. Commun.* 2020, 11:4065. https://www.nature.com/articles/s41467-020-17921-y . *Shifted+unshifted fragment indexing for mod localization.*
4. TagRecon (see Foundational) — the origin of "peptide−tag mass delta = mutation/PTM."
5. TagGraph (see Foundational) — unrestricted tag-string search at 25M-spectrum scale.
6. **Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale** — Lazear. *J. Proteome Res.* 2023, 22(11):3652-3659. DOI 10.1021/acs.jproteome.3c00486. https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00486 .
7. **Comet Fragment-Ion Indexing for Enhanced Peptide Sequencing** — *J. Proteome Res.* 2025, 24(7):3715-3721. DOI 10.1021/acs.jproteome.4c01094. https://pubs.acs.org/doi/10.1021/acs.jproteome.4c01094 (PMID 40526826). *Fragment indexing incl. real-time/immunopeptidomics use cases; up to 94% faster.*
8. MSFragger — Kong et al., *Nat. Methods* 2017 (fragment-ion index origin). Background, https://www.nature.com/articles/nmeth.4256 .
9. **PgxSAVy: comprehensive evaluation of variant peptide quality in proteogenomics** — 2024. https://pmc.ncbi.nlm.nih.gov/articles/PMC10825656/ . *Variant peptides are false-positive-prone; QC needed.*
10. NovoBoard — framework for de novo FDR/accuracy evaluation, *ScienceDirect* 2024. https://www.sciencedirect.com/science/article/pii/S1535947624001397 . Background on de novo confidence.
11. **Winnow: de novo peptide sequencing rescoring and FDR estimation** — arXiv 2025 (preprint). https://arxiv.org/pdf/2509.24952 . *Flag: preprint.*
12. **Real-Time Search-Assisted Acquisition improves coverage in multiplexed single-cell proteomics** — 2022. https://pmc.ncbi.nlm.nih.gov/articles/PMC8961214/ (PMID 35219906).
13. Single-cell computational proteomics challenges — *Mol. Cell. Proteomics* 2023. https://www.mcponline.org/article/S1535-9476(23)00028-2/fulltext .
14. **PTM-Shepherd: Analysis and Summarization of PTMs from Open Search Results** — Geiszler et al., *Mol. Cell. Proteomics* 2021. https://www.mcponline.org/article/S1535-9476(20)35132-X/fulltext ; https://pmc.ncbi.nlm.nih.gov/articles/PMC7950090/ .
15. **Detecting diagnostic features in MS/MS spectra of post-translationally modified peptides** — *Nat. Commun.* 2023, 14:4132. https://www.nature.com/articles/s41467-023-39828-0 .
16. **SpecGlobX: fast alignment of mass spectra capturing multiple complex modifications** — *J. Proteome Res.* 2023. https://pmc.ncbi.nlm.nih.gov/articles/PMC10631047/ . *Post-search delta interpretation via spectral alignment.*
17. JUMP (see Foundational) — tag length + pattern matching drives tag-search sensitivity/specificity.
18. **Unipept in 2024: Expanding Metaproteomics Analysis (missed cleavages, semi-/non-tryptic)** — *J. Proteome Res.* 2024. DOI 10.1021/acs.jproteome.4c00848. https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00848 . *LCA specificity erosion as databases grow.*
19. **MS2Rescore 3.0 Is a Modular, Flexible, User-Friendly Platform to Boost Peptide Identifications** — Buur, Declercq et al. *J. Proteome Res.* 2024. DOI 10.1021/acs.jproteome.3c00785. https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00785 (PMID 38491990).
20. **MS2Rescore: Data-driven rescoring dramatically boosts immunopeptide identification rates** — Declercq et al., preprint bioRxiv 2021 https://www.biorxiv.org/content/10.1101/2021.11.02.466886v2 (published *Mol. Cell. Proteomics* 2022). *Flag: peer-reviewed DOI 10.1016/j.mcpro.2022.100266 cited from memory, not re-verified this session.*
21. **Rescoring Peptide Spectrum Matches: Boosting Proteomics Performance by Integrating Peptide Property Predictors** — *Mol. Cell. Proteomics* 2024. https://www.sciencedirect.com/science/article/pii/S1535947624000884 (PMID 38871251). *Oktoberfest/Prosit/MS2PIP/DeepLC feature stack.*
22. **Comparative Analysis of Data-Driven Rescoring Platforms (HeLa digest)** — 2024/2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC11962579/ . *Percolator vs MS2Rescore vs Oktoberfest.*
23. **MegaPX: fast and space-efficient peptide assignment using IBF-based multi-indexing** — RKI-MF2 group. bioRxiv 2025 (DOI 10.1101/2025.04.14.648734) https://www.biorxiv.org/content/10.1101/2025.04.14.648734v1 ; published *Bioinformatics* 2026 (btag134) https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag134/8533242 . *Alignment-free k-mer/IBF de novo→taxonomy; mutation-tolerant refs. Flag: 2026 publication date as reported by venue.*
24. **NovoLign: metaproteomics by sequence alignment** — Kleikamp et al., *ISME Communications* 2024, 4(1):ycae121. DOI 10.1093/ismeco/ycae121. https://academic.oup.com/ismecommun/article/4/1/ycae121/7819827 . *De novo reads→taxonomy + DB-search coverage evaluation.*
25. **MetaNovo: probabilistic peptide discovery in complex metaproteomic datasets** — *PLOS Comput. Biol.* 2023. DOI 10.1371/journal.pcbi.1011163. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011163 (PMID 37327214). *De novo tag matching → reduced tailored DB, no prior composition.*
26. MutCombinator (see Foundational) — combinatorial-mutation variant-graph search.
27. Real-Time Search-Assisted Acquisition (see [12]).
28. **Real-Time Library Search Increases Cross-Link Identification Depth across All Levels of Sample Complexity** — Ruwolt, He et al. *Anal. Chem.* 2023. DOI 10.1021/acs.analchem.2c05141. https://pubs.acs.org/doi/10.1021/acs.analchem.2c05141 (PMID 36926872). *Diagnostic-ion-ratio triggering of a second scan.*
29. **Glycopeptide database search and de novo sequencing with PEAKS GlycanFinder** — *Nat. Commun.* 2023, 14:4046. https://www.nature.com/articles/s41467-023-39699-5 . *Peptide+glycan search + de novo; Y-ion/oxonium handling.*
30. **pGlycoNovo: uncovering missing glycans and unexpected fragments for site-specific glycosylation** — *Nat. Commun.* 2024. https://www.nature.com/articles/s41467-024-52099-7 .
31. **Unifying the analysis of bottom-up proteomics data with CHIMERYS** — *Nat. Methods* 2025. DOI 10.1038/s41592-025-02663-w. https://www.nature.com/articles/s41592-025-02663-w . *Most MS2 are chimeric; deconvolution ~doubles IDs.*
32. **A transformer model for de novo sequencing of DIA MS data (Cascadia)** — *Nat. Methods* 2025. https://www.nature.com/articles/s41592-025-02718-y .
33. Transformer-based de novo sequencing for DIA — arXiv 2024 (preprint) https://arxiv.org/abs/2402.11363 ; Disentangling Complex Multiplexed DIA Spectra — arXiv 2024 https://arxiv.org/abs/2411.15684 . *Flag: preprints.*
34. **Template-Based Assembly of Proteomic Short Reads for De Novo Antibody Sequencing and Repertoire Profiling (Stitch)** — Heck lab, *Anal. Chem.* 2022. DOI 10.1021/acs.analchem.2c01300. https://pubs.acs.org/doi/10.1021/acs.analchem.2c01300 .
35. **Comprehensive evaluation of peptide de novo sequencing tools for monoclonal antibody assembly** — *Briefings in Bioinformatics* 2023, 24(1):bbac542. https://academic.oup.com/bib/article/24/1/bbac542/6955273 . *Stitch/ALPS/Casanovo assembly, ~94% coverage.*
36. **MARS: improved de novo peptide candidate selection for non-canonical antigen discovery** — *Nat. Commun.* 2024, 15:661. DOI 10.1038/s41467-023-44460-z. https://www.nature.com/articles/s41467-023-44460-z . *MHC-centric de novo candidate selection (uses confident partial reads).*
37. **Casanovo: Sequence-to-sequence translation from mass spectra to peptides** — Yilmaz et al., *Nat. Commun.* 2024. DOI 10.1038/s41467-024-49731-x. https://www.nature.com/articles/s41467-024-49731-x .
38. **ionbot: machine-learning open modification search** — Degroeve, Gabriels et al. bioRxiv 2021 https://www.biorxiv.org/content/10.1101/2021.07.02.450686v2 . *Flag: preprint as found; ML-scored open search context for rescoring features.*

*Report generated 2026-07-23. Citations verified to landing pages/DOIs during research except where explicitly flagged. No citation here was fabricated; flagged items need a manual DOI check before use in a manuscript.*
