# ARALE: Technical Analysis and Comparison to the Latest State of the Art

*October 2026. Object of analysis: the ARALE repository README (ericrenone/ARALE), a specification plus Python model for a cryogenic adaptive E₈ lattice decoder.*

---

## Contents

1. [Question](#1-question)
2. [Executive summary](#2-executive-summary)
3. [Methodology](#3-methodology)
4. [What ARALE claims](#4-what-arale-claims)
5. [Technical soundness](#5-technical-soundness)
6. [Comparison to the latest peers](#6-comparison-to-the-latest-peers)
7. [Positioning matrix](#7-positioning-matrix)
8. [Risks relative to the field's trajectory](#8-risks-relative-to-the-fields-trajectory)
9. [Open questions](#9-open-questions)
10. [Recommendations](#10-recommendations)
11. [Source notes](#11-source-notes)

---

## 1. Question

Evaluate the ARALE README on three axes:

1. Is the technical core sound?
2. Are its citations and headline numbers credible?
3. How does it compare against the latest state of the art in each domain it touches: real-time QEC decoders, GKP and bosonic codes, lattice-decoder hardware, adaptive equalization, cryo-CMOS, and topological-semimetal devices, as of October 2026?

---

## 2. Executive summary

1. **The niche is real and unoccupied.** Independent searches corroborate the README's central premise: no published FPGA or ASIC implementation of an exact E₈ or Dₙ closest-point decoder was found. The nearest hardware peers are K-best and sphere MIMO detectors and LDLC fixed-point decoders. The slot ARALE targets, a fixed-latency exact lattice snap consuming analog GKP-type data at 4 K, has no incumbent.

2. **Every load-bearing citation checked out.** The Sherman–Morrison instability analysis (Hashemi & Nakatsukasa), the Pinball cryo-CMOS predecoder figures (< 0.56 mW peak, 22 nm FDSOI, HPCA 2026), the ultrathin NbP interconnect result (Science 387:62–67), TaIrTe₄ terahertz rectification (Nature Electronics), and Riverlane's decoder metrics all match their primary sources. The README also demotes its own speculative elements explicitly: chirality routing, sub-nanosecond latency, and the O(n²) advantage at n = 8.

3. **The strongest peer is Riverlane's Local Clustering Decoder** (Nature Communications, December 2025): FPGA-implemented, under 1 µs per decoding round, with an adaptive noise-model engine, claimed to cut physical-qubit overhead by up to 75% (d = 17 against d = 33 under leakage-dominant noise). Its adaptivity is noise-model reweighting on a discrete syndrome graph; ARALE's is rank-one channel-matrix tracking at the symbol rate on analog data. They are complementary, not competing: a surface–GKP stack would need an inner analog decoder (ARALE's slot) and an outer matching decoder (Riverlane's slot).

4. **The primary application remains speculative, and ARALE says so.** No four-mode E₈ GKP experiment exists. GKP state of the art is single-mode or small-Hilbert-space: qudits beyond break-even (Nature, 2025) and dissipative protection (PRX, January 2025). Multimode GKP *decoding* is still theory (Roy, Pousset & Royer, October 2025: an order-of-magnitude logical-error improvement, in software). ARALE gates its Phase 4 on an external experiment that may never materialize.

5. **All performance numbers are model estimates, correctly labeled but easy to over-read.** The 59 ns, 1 Gvector/s, ≈ 1.1 mW and 0.6 mm² figures are operation counts priced with 2014-era per-operation energies, geometrically interpolated to 22 nm, with a flat 2.5× overhead and no 4 K correction. They are internally consistent and order-of-magnitude plausible against Pinball's measured envelope. The two least-quantified risks are 1 GHz clocking of deep pipelines at 4 K and clock-tree power, both of which the README flags.

6. **The "Rank-One" in the name overstates the n = 8 benefit, by the README's own admission.** Its §5.5 concludes that at n = 8 the array simply re-factorizes on every update because the rank-one QR update saves only about 1.5–2×, and its §12 lists "O(N²) as a decisive advantage" among the things ARALE does *not* claim. The real rank-one economy is the shift-add NLMS update of Ĥ. The branding survives on the n = 16/24 scaling story, where the cheap E₈ snap no longer exists.

7. **E₈'s software momentum validates the geometry, not the hardware.** In 2026, E₈-lattice 2-bit KV-cache quantization entered llama.cpp, vLLM and SGLang community implementations, alongside TurboQuant (ICLR 2026) and QTIP trellis codes. This confirms E₈ as a live quantization geometry, and equally confirms ARALE's framing that all of this usage is software, orthogonal to ARALE.

---

## 3. Methodology

- **Primary-object analysis.** The ARALE README as published. The repository's Python model was run in a separate session, and its output agrees with the §8 figures quoted in the README (decoder verified on 2000/2000 points; SER table; tracker floor 0.087 at μ = 2⁻⁴; capture from ‖Δ‖_F = 0.82 in 439 symbols; Sherman–Morrison growth 3.3×10⁻⁶ to 7.7×10⁻⁵ per update; QR floor 4.8×10⁻⁴; CORDIC max error ≈ 4×10⁻⁵). The SOTA comparison in this report did not re-simulate; it assesses those numbers for internal consistency and agreement with the literature.
- **Citation verification.** Load-bearing references were checked against primary sources by targeted web search. The most decision-relevant pages (Riverlane LCD coverage, the Pinball preprint abstract, the multimode GKP decoder preprint) were opened and read.
- **SOTA mapping.** 2024–2026 literature was searched in eight domains: QEC decoder hardware, GKP experiments and theory, lattice-code decoder hardware, MIMO and lattice-reduction detectors, E₈ quantization in ML systems, rank-one-update numerical analysis, cryo-CMOS at 4 K, and Weyl-semimetal devices.
- **Limitations.** Web search cannot prove a negative ("no E₈ hardware decoder exists"); it corroborates it. Paywalled full texts (the Nature Communications LCD paper body) were not read; press coverage and abstracts were used instead. Where a claim rests on a snippet rather than a read page, the source notes say so.

---

## 4. What ARALE claims

A three-stage feed-forward decoder (16-bit front end; CORDIC-Givens QR equalizer with 28 rotations and a triangular solve, no multipliers; the exact Conway–Sloane E₈ coset decoder at about 120 integer operations) closed by a decision-directed NLMS feedback loop that corrects the tracked channel estimate by a rank-one update every symbol, with μ = 2⁻ᵏ so the update is shift-add only. It is specified for 4 K cryo-CMOS (22 nm FDSOI, 1 GHz assumed), with an SFQ speed option and an optional Weyl-semimetal interconnect.

Headline figures: 59 ns latency, 1 Gvector/s, ≈ 0.45 nJ per vector for the datapath (≈ 1.1 nJ with overhead), ≈ 1.1 mW at full rate, ≈ 0.6 mm², and E₈'s 3.01 dB nominal coding gain preserved against a drifting channel with tracked mismatch ‖Δ‖_F ≤ 0.10.

---

## 5. Technical soundness

| Claim | Basis | Verdict |
|---|---|---|
| Exact E₈ decode = two D₈ decodes + compare | Conway–Sloane 1982; Viazovska 2017 | **Settled mathematics.** Sound. |
| Decision-directed NLMS tracking; stability for μ < 2; misadjustment μ/(2−μ) | Textbook adaptive filtering | **Sound**, with the caveats in §8. |
| QR carried by default; explicit inverse confined to κ₂ < 10 with periodic refresh | Hashemi & Nakatsukasa (2025–26): SM forward error ~ ε·κ₂², not caused by small denominators | **Correctly applied.** Matches the best practice this literature implies. No 2022–2026 hardware SM block was found, consistent with the README. |
| 16-iteration, 20-bit CORDIC error ≈ 4×10⁻⁵ per rotation | Standard CORDIC theory; MIMO fixed-point practice (16-bit keeps BER loss < 0.3 dB) | **Consistent with the literature.** |
| 59-clock latency | Cycle budget 4 + 9 + 16 + 8 + 18 + 4 | **Internally consistent**, but entirely conditional on 1 GHz at 4 K, the least-supported assumption in the document. |
| Energy and area (0.45 nJ, 1.1 mW, 0.6 mm²) | Operation counts × Horowitz 2014 energies; 22 nm as the geometric mean of 45 nm and 7 nm; flat 2.5× overhead | **Estimate, not prediction.** No clock-tree synthesis, no interconnect energy, no 4 K threshold-shift correction. Order-of-magnitude plausible against Pinball's measured < 0.56 mW at 4 K. |
| Tracker floor ‖Δ‖_F ≈ 0.087 (μ = 2⁻⁴, 12 dB); capture from ‖Δ‖_F = 0.8 in ≈ 440 symbols | Python model, single seeded realizations, 240-point codebook | **Model result, honestly labeled.** Single-realization statistics, floating point, ungated loop. |
| 2¹⁶-point codebook behavior (capture range shrinks ≈ 2.2×) | Scaling estimate | **Extrapolated, not simulated.** The weakest load-bearing number in the document, which the README admits. |
| Rank-one QR update economy | Gill–Golub–Murray–Saunders 1974; Gentleman–Kung arrays | **Classical.** At n = 8 the README concedes re-factorization wins; the "Rank-One Array" branding rests on the shift-add NLMS update and the n = 16/24 story. |

**Claim hygiene.** The README's §12 explicitly disclaims sub-nanosecond latency, the O(n²) advantage at n = 8, topological error immunity, and quantum-stack displacement. All four disclaimers are correct and rare in projects of this kind. Its §7 marks an evidence level for every device claim, and its §13.3 attaches numeric kill-gates to every speculative element. This is disciplined, and it materially raises the document's credibility.

---

## 6. Comparison to the latest peers

### 6.1 Real-time QEC decoders (the strongest peer group)

| Decoder | Platform | Speed | Adaptivity | Status, October 2026 |
|---|---|---|---|---|
| **Riverlane Local Clustering Decoder** | FPGA (Deltaflow 2) | < 1 µs per round | Continuous noise-model updates; handles correlated and leakage errors; claims up to 75% qubit-overhead reduction (d = 17 vs d = 33) | Peer-reviewed (Nature Communications, Dec 2025); deployed with Rigetti, OQC, Infleqtion, ORNL |
| **Riverlane decoder ASIC** | ASIC | order-100 ns class (Nature Electronics, 2025) | Not adaptive in the channel sense | Fabricated demonstrator; first decoder chip |
| **Google Willow real-time decoding** | FPGA / real-time | 1.1 µs cycle; ≈ 63 µs mean latency (disputed as "real-time") | Recalibration offline | Nature 2024; the reference point for the syndrome-cycle budget |
| **Pinball cryogenic predecoder** | 22 nm FDSOI at 4 K | HPCA 2026 | Predecoding, not channel-adaptive | Measured: < 0.56 mW peak; supports 2,668 logical qubits at d = 21 in a 1.5 W 4 K budget |
| **FPGA neural-network surface-code decoder** | FPGA | Real-time | Learned, not symbol-adaptive | arXiv 2605.04892 (2026) |
| **ARALE** | Paper + Python model | 59 ns (model) | Rank-one channel-matrix tracking every symbol on analog data | No RTL, no silicon |

**Assessment.** ARALE is not competing with these decoders on their ground. It targets the inner, analog, lattice-structured decode that a concatenated surface–GKP architecture would need beneath them, and that slot is empty. But the evidence-level asymmetry must be stated plainly: the peer group has fabricated chips and deployed FPGA stacks; ARALE has a cycle budget and a numpy script. One definitional point in ARALE's favor: Riverlane's "adaptive" means re-weighting a discrete matching model; ARALE's means estimating a drifting 8×8 linear map at the symbol rate with no training sequence, a different and more aggressive form of adaptivity with no published hardware precedent found.

### 6.2 GKP experimental status (the dependency ARALE cannot control)

- **Qudits beyond break-even** (Nature 2025): error-corrected qutrits and ququarts beat break-even, the strongest GKP-family result to date, still single-oscillator.
- **Dissipative GKP protection** (PRX, January 2025): a high-impedance circuit scheme, theory to early experiment.
- **Multimode GKP decoding** (Roy, Pousset & Royer, arXiv 2510.12677, October 2025): correlation-aware decoding of Steane-type protocols cuts logical error probability by at least 10×, in software. The algorithmic layer is active; the hardware layer ARALE proposes does not exist.
- **No four-mode E₈ GKP experiment** was found. The README's statement that this application "waits on one" remains accurate, and its Phase 4 gate is correctly framed as externally dependent. The two-mode D₄ fallback (Royer, Singh, Girvin theory line) is the realistic near-term integration target.

### 6.3 Lattice-decoder hardware (the gap is confirmed)

Searches across 2020–2026 returned no FPGA or ASIC E₈ or Dₙ nearest-point decoder. The nearest neighbors:

- **Sorter-free K-best MIMO detector** (OJCAS 2024, the README's own comparator): 6.4 Gb/s, 0.32 µs, 4×4 256-QAM, with QRD recomputed per estimation frame and square-QAM constellations.
- **LLL and lattice-reduction VLSI, LR-aided systolic detectors:** 14–18-bit fixed point. The 16-bit, 0.3 dB precision rule ARALE budgets against comes from this literature and is consistent with it.
- **LDLC fixed-point FPGA decoder:** an iterative lattice-code decoder with Gaussian-message approximation, a different (soft, iterative) regime. ARALE's exact one-shot snap is architecturally the opposite.

**The differentiator that survives scrutiny:** every published hardware lattice decoder either assumes a static channel factorization per frame or re-estimates offline. ARALE's per-symbol decision-directed rank-one tracking without pilots, feeding an exact snap in fixed depth, is a system-level combination of established ingredients: not a new algorithm, but an unfilled architectural niche.

### 6.4 E₈ in software quantization (accelerating, and orthogonal)

2026 developments: E₈-lattice 2-bit KV-cache quantization adopted in llama.cpp, vLLM and SGLang community stacks (calibration-free, reported to preserve retrieval where scalar 2-bit fails); TurboQuant (ICLR 2026); QTIP trellis-coded quantization (NeurIPS 2024; decode at about 2 instructions per weight, above 80% of peak memory bandwidth). All software, on GPUs and CPUs. This validates the README in both directions: E₈ is a live, valuable geometry, and none of that momentum produces hardware closest-point decoders. One watch item: if LLM inference engines ship E₈ decode in GPU kernels, a room-temperature peer for parts of the datapath could emerge on its own.

### 6.5 Numerical core (current best practice)

The Sherman–Morrison instability results ARALE rests on are real and current: arXiv 2510.01696 (instability of SM; stabilization by iterative refinement) and the 2026 follow-up arXiv 2609.12266 (error bounds; SMIR). ARALE's response, QR by default with the explicit inverse confined to κ₂ < 10 and a measured refresh interval N_refresh = 0.025/g, is a faithful and conservative application of that literature. No contradicting result was found.

### 6.6 Substrate and technology mapping (verified, with one nuance)

- **Ultrathin NbP interconnect** (Science 387:62–67, 2025): confirmed. Surface-dominated conduction; effective resistivity decreasing below about 18 nm. ARALE cites the numbers accurately and keeps every performance number independent of this upgrade.
- **NbAs nanowires** (arXiv 2503.04621): confirmed.
- **TaIrTe₄ nonlinear-Hall rectification:** room-temperature RF rectification and THz sensing confirmed (Nature Nanotechnology 2021; Nature Electronics 2025). **Nuance:** the nonlinear Hall effect in TaIrTe₄ has been measured at 4 K under DC field control (arXiv 2502.05960), which slightly weakens the README's "no 4 K characterization exists" wording. The rectifier figures of merit ARALE gates on (responsivity, NEP, linearity under GHz drive at 4 K) remain unmeasured, so the Phase 3 gate stands as written.
- **SFQ/ERSFQ** (≈ 2.5 aJ per switch, 58–60 GHz demonstrated): consistent with the literature; correctly positioned as a speed option with no dense same-temperature memory, which ARALE sidesteps with about 1.2 kbit of state.

---

## 7. Positioning matrix

| Dimension | Latest peer (measured) | ARALE (model) | Gap |
|---|---|---|---|
| Decode latency | K-best 0.32 µs; Riverlane < 1 µs per round | 59 ns | Claimed 5–17× better; unproven, conditional on 1 GHz at 4 K |
| Adaptivity | Riverlane LCD noise-model reweighting, deployed | Rank-one channel tracking every symbol, no pilots | A more aggressive claim with zero hardware evidence: the project's core bet |
| Precision strategy | 16-bit LR-aided MIMO practice | 16/20-bit, error measured in the model | Aligned with practitioner data |
| Power at 4 K | Pinball < 0.56 mW (measured, HPCA 2026) | ≈ 1.1 mW full rate (modeled) | Comparable envelope; ARALE does far more arithmetic per nanosecond |
| Exactness | No exact lattice decoder in hardware | Exact E₈ snap, ≈ 120 ops, fixed schedule | Unique, if built |
| Energy per decoded vector | No direct peer | ≈ 0.45 nJ datapath | Unverifiable until silicon; the interpolation is crude |
| Evidence level | Fabricated ASICs, deployed FPGA stacks | README + numpy model | **The central gap.** ARALE is a spec; its peers are products |

---

## 8. Risks relative to the field's trajectory

1. **Field momentum is elsewhere.** The 2024–2026 QEC-decoder surge (Riverlane LCD, FPGA neural decoders, cryo predecoders) is all discrete-syndrome decoding. If multimode GKP experiments arrive, the inner analog slot may first be filled by a simple FPGA Babai nearest-plane decoder feeding soft information to the outer matcher; exactness of the inner snap may turn out not to matter for the concatenated logical error rate. ARALE's §5.6 argument against this is reasoned but not evidenced.
2. **1 GHz at 4 K in deep pipelines** has thin public support. Pinball, the best-characterized 4 K 22 nm FDSOI data point, emphasizes voltage and frequency scaling and body biasing rather than GHz clocking. ARALE's own fallback (90–120 ns at about 500 MHz) still clears the 1.1 µs budget, which is the right framing, but the headline 59 ns should be read as a clock-budget result, not a chip prediction.
3. **The 2¹⁶-point codebook** is the load-bearing wall without a simulation. Every real communications or GKP use needs points beyond the 240 minimal vectors, and ARALE's own estimate (capture range shrinks about 2.2×) suggests its default μ = 2⁻⁴ may not close there. Pilot-aided acquisition is listed as an open problem, not a solution.
4. **Clock-tree and interconnect power at 4 K** sit outside the flat 2.5× overhead factor. At 1 GHz with about 116 CORDIC cells this could dominate the modeled datapath energy, and it is the most likely way the ≈ 1.1 mW figure degrades in silicon.

---

## 9. Open questions

1. **Does the tracker converge on a dense codebook?** The 2¹⁶-point Voronoi codebook simulation is the missing evidence; until it runs, the communications and GKP applications rest on a 240-point demonstration.
2. **Does the inner snap need to be exact?** No surface–GKP simulation yet compares an exact E₈ inner decode against Babai plus soft information into the outer matcher at matched latency. This is the strongest available falsification test of the whole architecture.
3. **What is real clock-tree power at 4 K at 1 GHz?** No public 4 K data point exists for GHz-deep pipelines in 22 nm FDSOI; Phase 2 answers it.
4. **Will a four-mode E₈ GKP experiment exist within the roadmap window (months 24–36)?** Entirely external; the two-mode D₄ fallback defines the realistic integration target.
5. **Can the nonlinear-Hall front end reach NEP < 10 pW/√Hz at 4 K?** The 4 K physics is now partly measured; the rectifier figures of merit are not.
6. **Does an LLM-inference E₈ GPU kernel emerge as an accidental room-temperature peer?** E₈ KV-cache momentum suggests parts of ARALE's snap logic could be productized independently of the cryogenic story.

---

## 10. Recommendations

1. **Read ARALE as what it is, a falsifiable pre-silicon specification, and hold it to its own gates.** Its claim discipline is at the level of the field's best; its evidence is not yet.
2. **Close the weakest number first.** Publish the 2¹⁶-point codebook tracking simulation (tracker floor, capture range, gate tuning) before or alongside the Phase 1 FPGA. It costs days and de-risks the entire applications section.
3. **Add the counterfactual baseline.** A Babai nearest-plane decoder at matched latency against the exact snap, in a concatenated surface–GKP simulation, would either justify the central design choice or show that the simpler decoder suffices.
4. **Reframe the headline latency.** "59 ns at an assumed 1 GHz" invites the wrong argument. "At most 120 ns at a 4 K-supported clock, inside the 1.1 µs cycle with 90% margin" is the defensible claim and matches the Pinball-era 4 K evidence.
5. **De-emphasize the substrate branding.** "Chiral Transport Hardware" in the title rests on a wire-RC upgrade and a hypothesis that is dropped by default. The defensible core is the decoder architecture; a title centered on adaptive exact lattice decoding at 4 K would align the claims with the evidence.
6. **Track two external triggers:** (a) any multimode GKP experiment (the Roy–Royer theory line is the leading indicator), and (b) Riverlane's Deltaflow 3 streaming logic (late 2026). If surface–GKP concatenation enters their stack, ARALE's inner-decoder slot becomes commercially real; if not, the optical-shaping application (room temperature, 16 parallel instances) is the honest fallback.

---

## 11. Source notes

| Source | Credibility | Date |
|---|---|---|
| [ARALE repository README (ericrenone/ARALE), object of analysis](https://github.com/ericrenone/ARALE) | 3/5 (self-published spec and model, not peer-reviewed) | Oct 2026 |
| [Riverlane Local Clustering Decoder, Nature Communications](https://www.nature.com/articles/s41467-025-66773-x) | 5/5 (peer-reviewed; full text not read, confirmed via coverage) | Dec 2025 |
| [The Quantum Insider: Riverlane hardware decoder publication](https://thequantuminsider.com/2025/12/18/riverlane-hardware-decoder-real-time-qec/) | 4/5 (press; opened and read) | Dec 18, 2025 |
| [Quantum Computing Report: Riverlane first adaptive hardware decoder](https://quantumcomputingreport.com/riverlane-unveils-first-adaptive-hardware-decoder-to-deliver-real-time-quantum-error-correction/) | 4/5 (trade press; opened and read) | Dec 2025 |
| [Riverlane press release: first quantum decoder chip](https://www.riverlane.com/press-release/riverlane-announces-world-s-first-quantum-decoder-chip) | 4/5 (vendor) | 2024 |
| [HPCwire: Riverlane decoder hardware in Nature Electronics](https://www.hpcwire.com/off-the-wire/riverlane-details-quantum-decoder-hardware-implementation-in-nature-electronics/) | 4/5 (press summary of peer-reviewed work) | Jan 8, 2025 |
| [Pinball predecoder, arXiv 2512.09807 (HPCA 2026); abstract opened and read](https://arxiv.org/abs/2512.09807) | 5/5 (peer-reviewed venue; 4 K 22 nm FDSOI, < 0.56 mW) | Dec 10, 2025 |
| [Google, Quantum error correction below the surface code threshold (Willow), record](https://collaborate.princeton.edu/en/publications/quantum-error-correction-below-the-surface-code-threshold/) | 4/5 (record of the Nature 2024 paper) | 2024 |
| [Notes on Google Willow: real-time decoding latency debate (63 µs)](https://www.researchgate.net/publication/387576664_Notes_on_Google_Willow_2024_aka_Quantum_Error_Correction_Below_the_Surface_Code_Threshold) | 2/5 (informal commentary) | 2024 |
| [Quantum error correction of qudits beyond break-even, Nature](https://www.nature.com/articles/s41586-025-08899-y) | 5/5 (peer-reviewed) | 2025 |
| [Decoding multimode GKP codes with noisy auxiliary states, arXiv 2510.12677; abstract opened and read](https://arxiv.org/abs/2510.12677) | 4/5 (preprint) | Oct 14, 2025 |
| [Hashemi & Nakatsukasa, Instability of the Sherman–Morrison formula, arXiv 2510.01696](https://arxiv.org/abs/2510.01696) | 5/5 (peer-review-track preprint) | Oct 2025 |
| [Hashemi & Nakatsukasa, SM error bounds and SMIR, arXiv 2609.12266](https://arxiv.org/abs/2609.12266) | 4/5 (preprint) | 2026 |
| [Surface conduction and reduced resistivity in ultrathin noncrystalline NbP, Science 387:62–67](https://www.science.org/doi/10.1126/science.adq7096) | 5/5 (peer-reviewed) | 2025 |
| [Surface-dominant transport in Weyl semimetal NbAs nanowires, arXiv 2503.04621](https://arxiv.org/pdf/2503.04621) | 4/5 (preprint) | 2025 |
| [Terahertz sensing via nonlinear electrodynamics of TaIrTe₄, Nature Electronics](https://www.nature.com/articles/s41928-025-01397-z) | 5/5 (peer-reviewed) | 2025 |
| [Electric-field control of the nonlinear Hall effect in TaIrTe₄ (4 K measurements), arXiv 2502.05960](https://arxiv.org/abs/2502.05960) | 4/5 (preprint; reached via a review page) | 2025 |
| [QTIP: Quantization with Trellises and Incoherence Processing, arXiv 2406.11235](https://arxiv.org/abs/2406.11235) | 5/5 (NeurIPS 2024) | Jun 2025 (v4) |
| [llama.cpp discussion: E₈ lattice 2-bit KV-cache quantization](https://github.com/ggml-org/llama.cpp/discussions/25363) | 3/5 (community engineering discussion) | 2026 |
| [SGLang discussion: E₈ lattice KV cache, calibration-free 2-bit](https://github.com/sgl-project/sglang/discussions/30279) | 3/5 (community discussion) | 2026 |
| [Hardware implementation of a fixed-point LDLC decoder (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854306) | 4/5 (peer-reviewed, FPGA) | 2012 |
| [Lattice modulo sampling (E₈ folding in ADCs), arXiv 2605.24974](https://arxiv.org/abs/2605.24974) | 4/5 (preprint) | 2026 |
| [Hardware-in-the-loop syndrome-to-decoder validation including digitized GKP, arXiv 2607.19447](https://arxiv.org/html/2607.19447) | 4/5 (preprint) | 2026 |
| [Real-time surface-code error correction with an FPGA neural decoder, arXiv 2605.04892](https://arxiv.org/html/2605.04892v1) | 4/5 (preprint) | 2026 |
| [Classical interfaces for cryogenic quantum computing (cryo-CMOS survey), APL Quantum](https://pubs.aip.org/aip/apq/article/2/4/041501/3373674) | 4/5 (peer-reviewed review) | 2025 |

**Conflicts and caveats.** (a) The Willow real-time decoding latency is disputed (63 µs mean over about 50 cycles); ARALE uses the 1.1 µs *cycle* budget, which is the correct and undisputed figure. (b) ARALE states that no 4 K characterization of the TaIrTe₄ nonlinear-Hall response exists; DC nonlinear-Hall measurements at 4 K do exist (arXiv 2502.05960), though the rectifier figures of merit under GHz drive at 4 K, which ARALE gates on, remain unmeasured. (c) "No E₈ hardware decoder exists" is corroborated by searches, not proven. (d) The README's §8 numbers were reproduced by running its model, not independently re-derived.

---

*Prepared October 2026. All web sources were retrieved and inspected at that time.*
