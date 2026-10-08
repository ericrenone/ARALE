# ARALE

**Adaptive Rank-One Array Lattice Engine for E₈ Decoding on Chiral Transport Hardware**


---

<img width="380" height="526" alt="image" src="https://github.com/user-attachments/assets/5918b7da-4cc2-4638-b1a1-9bb21aedf5fc" />




---

ARALE is a cryogenic closest-point decoder for the E₈ lattice. It recovers an E₈ point from an 8-dimensional observation that has passed through an unknown, slowly drifting linear channel, and it keeps its channel estimate current at the symbol rate without a training sequence. It is specified for a 4 K cryo-CMOS implementation (with an SFQ speed option and an optional topological-semimetal interconnect upgrade) and is accompanied by a self-contained Python model, `arale_model.py`, that produces every simulated number in this README.

This file is complete in itself. It states the problem, the mathematics, the architecture, the numerical results with the exact commands that reproduce them, the physical technology mapping with its evidence levels, the target applications, the limits of the claims, a validation roadmap, and the references.

---

## Contents

1. [Summary and headline numbers](#1-summary-and-headline-numbers)
2. [Problem statement](#2-problem-statement)
3. [Background facts the design rests on](#3-background-facts-the-design-rests-on)
4. [System overview](#4-system-overview)
5. [Mathematical core](#5-mathematical-core)
6. [Hardware architecture](#6-hardware-architecture)
7. [Physical substrate and technology mapping](#7-physical-substrate-and-technology-mapping)
8. [Numerical results (Python model)](#8-numerical-results-python-model)
9. [Performance model: latency, throughput, energy, area](#9-performance-model-latency-throughput-energy-area)
10. [Comparison with alternative approaches](#10-comparison-with-alternative-approaches)
11. [Target applications](#11-target-applications)
12. [What ARALE does not claim](#12-what-arale-does-not-claim)
13. [Open problems, risks, and validation roadmap](#13-open-problems-risks-and-validation-roadmap)
14. [Repository and reproduction](#14-repository-and-reproduction)
15. [References](#15-references)

---

## 1. Summary and headline numbers

ARALE joins three components that are each individually established:

1. **Decision-directed rank-one channel tracking.** The estimate Ĥ of the 8×8 channel matrix is corrected after every decoded vector by a normalized-LMS update Ĥ ← Ĥ + μ·e·x̂ᵀ/‖x̂‖², a rank-one correction. With μ = 2⁻ᵏ the update contains only shifts and adds.
2. **A multiplier-free triangular systolic array of CORDIC rotators.** The equalizer is carried in orthogonal (QR) form and applied as 28 planar Givens rotations followed by a CORDIC triangular solve. No multiplier appears in the data path.
3. **The exact Conway–Sloane coset decoder.** The equalized vector is snapped to the nearest E₈ point by two D₈ decodes and one distance comparison, about 120 integer operations with a fixed schedule.

Four design decisions follow from the literature and are kept throughout:

- **A fixed E₈ lattice never needs a matrix inverse.** The coset decoder is already exact and linear-time. Inversion earns its place only when the channel skews the lattice.
- **The equalizer is carried in QR form by default.** Sherman–Morrison updates of an explicit inverse have forward error that scales as ε·κ₂² and are reserved for well-conditioned channels (κ₂ < 10) with periodic re-factorization.
- **Quantum phase estimation and HHL are not comparators.** They estimate an eigenphase or return a quantum state; neither returns an argmin lattice vector. The relevant comparators are software coset decoders, FPGA sphere/K-best detectors, and real-time quantum-error-correction decoders.
- **Weyl-semimetal "topological protection" protects the existence of band nodes, not a transport coefficient.** No claim of immunity to scattering or to arithmetic error is made. Weyl materials are used only for what they demonstrably do.

### Headline results

All performance figures are model results from operation counts and published per-operation energies, not silicon measurements. Error-rate, tracking, and numerical results come from the Python model in this repository.

| Quantity | Value | Where it comes from |
|---|---|---|
| Pipeline latency, cryo-CMOS at 1 GHz | **59 clocks = 59 ns** (4 + 9 + 16 + 8 + 18 + 4) | §6.5, `budget` |
| Throughput, full-rate configuration | 1 vector per clock = **1 Gvector/s** (8 Gsample/s) | §9.2 |
| Throughput, reduced configuration | 1 vector per 3 clocks = 0.33 Gvector/s | §9.2 |
| Datapath energy per decoded vector (22 nm, interpolated) | ≈ **0.45 nJ**; ≈ 1.1 nJ with clock/register overhead | §9.3 |
| Power at full rate / at QEC duty (100 Mvector/s) | ≈ 1.1 mW / ≈ 0.11 mW | §9.3 |
| Area (22 nm estimate) | full-rate ≈ 0.6 mm² (≈ 1.8 MGE); reduced ≈ 0.3 mm² (≈ 0.85 MGE) | §9.4 |
| E₈ nominal coding gain over Z⁸ | 3.01 dB, preserved by the tracked decoder | §8.2 |
| Mismatch tolerance of the snap (SER, 100,000 trials per cell) | ≤ 0.05 dB loss at ‖Δ‖_F = 0.10; ≈ 0.12–0.16 dB at 0.20 | §8.2 |
| Tracker floor, μ = 2⁻⁴, SNR 12 dB, 65-symbol loop delay | ‖Δ‖_F = 0.087 | §8.3 |
| Recovery from a 10% single-line gain step, μ = 2⁻⁴ | 1/e in 105 symbols, 0 decision errors | §8.3 |
| Capture range | re-acquires from 50% full-rank error in 382 symbols, from 80% in 439, with ≤ 0.5% SER | §8.4 |
| CORDIC error (16 iterations, 20-bit) | max 4.22×10⁻⁵, RMS 1.44×10⁻⁵ per rotation | §8.6 |
| QR-mode representation error, 10⁴ updates, κ₂ = 10 | worst ‖R⁻¹QᵀH − I‖_F = 4.82×10⁻⁴ (no growth) | §8.5 |
| Sherman–Morrison growth per update, 16-bit fraction | 3.3×10⁻⁶ (κ₂ = 2), 1.5×10⁻⁵ (κ₂ = 10), 7.7×10⁻⁵ (κ₂ = 30) | §8.5 |
| State | Ĥ 1,152 b; R and Qᵀy 648 b; y delay line ≈ 9 kb; pipeline ≈ 2 kb | §6.6 |
| Primary application | multi-mode GKP closest-point decoding inside a 1.1 µs syndrome cycle | §11.1 |

The one-sentence defensible claim: *a 4 K decoder that holds E₈'s 3 dB coding gain on a drifting channel, re-estimating the channel every symbol, at a rate three to four orders of magnitude faster than software re-calibration, in a power envelope about a hundred times below the 4 K cooling budget.*

---

## 2. Problem statement

The problem is **adaptive closest-point decoding**: given

```
y = H x + w ,      y ∈ R⁸,   x ∈ α·E₈,   w ~ N(0, σ² I₈),
```

where H is an unknown, invertible, slowly drifting 8×8 matrix and α is the constellation scale, return the lattice point x̂ that minimizes the decoding error rate, within a latency budget below one microsecond, with the estimate of H kept current without stalling the data path.

E₈ is the right lattice to target because its properties are settled. It is the densest sphere packing in eight dimensions (Viazovska 2017, density π⁴/384 ≈ 0.2537), it has the maximal kissing number 240, it is even and unimodular, and its normalized second moment G(E₈) ≈ 0.0717 lies 0.88 dB from the sphere bound, against 1.53 dB for scalar quantization. These facts give E₈ a nominal coding gain of 3.01 dB over Z⁸ and a granular gain of about 0.65 dB, which is why it recurs in 2024–2026 work on LLM weight quantization (QuIP#'s E8P codebook), coherent-optical shaping (Chalmers), ML-KEM message encoding, and bosonic quantum codes.

E₈ is unusual among dense lattices in that its nearest-point problem is cheap. Conway and Sloane showed that E₈ = D₈ ∪ (D₈ + ½·**1**), so decoding reduces to two D₈ decodes (round every coordinate; fix parity by re-rounding the worst coordinate the other way) and one distance comparison, roughly 15 operations per coordinate and about 120 per vector. A fixed E₈ lattice therefore needs no matrix inverse, no sphere search, and no systolic array.

The matrix machinery becomes necessary exactly when H ≠ I. A physical channel (eight parallel analog lines, eight quadrature components of four optical modes, eight homodyne outcomes from four bosonic modes) introduces gain imbalance, crosstalk, phase rotation, and offset, and these drift with temperature, bias, and aging. The received constellation is then the skewed lattice H·E₈ whose Voronoi cells are no longer those of E₈. Two decoders are available: equalize with Ĥ⁻¹ and run the coset decoder on Ĥ⁻¹y (zero-forcing), or run Babai's nearest-plane algorithm directly on the basis ĤM. Both need Ĥ, and Ĥ must be tracked.

Three existing approaches each fail one requirement:

- **Software re-calibration** (CPU, periodic least-squares re-estimation of H, then an O(n³) inverse) has millisecond control-loop latency, three orders of magnitude slower than the 1.1 µs syndrome cycle of a superconducting error-correction experiment and far too slow for a 10¹⁰-symbol/s optical link.
- **FPGA MIMO detectors** (K-best, sphere) reach 0.3 µs latency and several Gb/s at 4×4 with 256-QAM, but assume a QR factorization recomputed per channel-estimation frame (47–220 cycles) and target square QAM constellations rather than lattice cosets.
- **Quantum linear-algebra subroutines** (QPE, HHL) produce a quantum state rather than an argmin vector, carry κ-to-κ² complexity dependence, and have no demonstrated use for closest-vector decoding.

ARALE fills the gap: a decoder whose channel estimate is updated by rank-one corrections at the symbol rate, whose equalizer is applied by multiplier-free rotations in a fixed-depth pipeline, and whose final snap is the exact coset decoder. The design questions are which representation of Ĥ to carry (explicit inverse, R factor, or R⁻¹), what word length suffices, what the loop's stability margin is, and which cryogenic technology can host it within the 4 K thermal budget of a dilution refrigerator (roughly 0.7–1.5 W residual at the 4 K stage of a current large system).

---

## 3. Background facts the design rests on

| Area | Result | Numbers | Source |
|---|---|---|---|
| E₈ structure | densest 8-D packing, proven | density π⁴/384 ≈ 0.2537; kissing number 240 | Viazovska 2017 |
| E₈ decoding | coset decoder via two D₈ decodes | ≈ 15 ops/coordinate, ≈ 120 ops/vector | NestQuant 2025; Conway–Sloane |
| E₈ quantizer | normalized second moment | G = 0.0717, 0.88 dB from sphere bound | HyperQuant 2026 |
| E₈ in LLMs | E8P 2-bit codebook, Llama-2-70B | PPL 4.16 (E8P) vs 4.45 (8-D k-means); FP16 3.12 | QuIP# 2024 |
| E₈ in optics | 8-D Voronoi constellations, coherent fiber | up to 1.84 dB SNR gain over Gray-QAM BICM | Chalmers 2025 |
| E₈ in PQC | 2E₈ ∩ Z₄⁸ message code for Kyber-1024 | decryption-failure rate 2⁻¹⁷⁴ → 2⁻²⁸⁶ | arXiv 2601.08452 |
| GKP lattices | multimode GKP decoding is a closest-point problem | D₄ distance √(2π) ≈ 2.51; exact CVP exponential in modes | Lin, Chamberland, Noh 2023 |
| QEC decoder budget | real-time surface-code decoding | 1.1 µs cycle; 63 µs mean latency; d = 17 in < 1 µs/round on FPGA | Google 2024; Riverlane 2024 |
| Cryo predecoder | 22 nm FDSOI at 4 K | < 0.56 mW peak; fits 1.5 W 4 K budget | Pinball 2025 |
| MIMO detector HW | sorter-free K-best, 4×4 256-QAM, 90 nm | 6.4 Gb/s, 0.32 µs latency | IEEE OJCAS 2024 |
| CORDIC QRD | 4×4 systolic CORDIC, 15–16-bit | 568 MHz, max error 6.9×10⁻⁴ at 16 bits | Muñoz & Hormigo |
| CORDIC precision | guard bits for n iterations | ≈ log₂ n extra fraction bits | MathWorks guide |
| Sherman–Morrison | instability analysis | forward error ≈ ε·κ₂(A)²; one refinement step restores stability | Hashemi & Nakatsukasa 2025 |
| Energy per op | 45 nm CMOS | INT32 add 0.1 pJ, INT32 mult 3.1 pJ | Horowitz 2014 |
| SFQ logic | RSFQ/ERSFQ switching energy | ≈ 2.5 aJ per junction switch; 58–60 GHz demonstrated | Holmes et al. |
| Weyl v_F | TaAs quantum oscillations | 3.0–7.6×10⁵ m/s; mobility 1.2–1.9×10⁴ cm²/Vs | arXiv 1603.08846 |
| Topological interconnect | NbP noncrystalline films, 400 °C sputter | ≈ 34 µΩ·cm at 1.5 nm vs ≈ 100 for metals | Science 2024 |
| Topological interconnect | NbAs nanowires, 40 nm | 9.7 µΩ·cm, 3–4× below bulk | Science 2026 |
| Chiral-anomaly device | PtSe₂ field-effect chirality device | ON/OFF ≈ 10³ at 6–9 T, cryogenic | arXiv 2103.00279 |
| Fluxonium | two-qubit gate and readout | CNOT 99.94% in 60 ns; readout 97.8% in 260 ns | arXiv 2407.15783, 2501.16691 |
| Fridge budget | 4 K residual cooling | 0.7–1.5 W (XLD class); 12–14 µW at 20 mK | arXiv 2608.00990; Bluefors |

### 3.1 E₈ decoding and lattice hardware

The E₈ decoder is a solved problem in software and an unbuilt one in hardware. No FPGA or ASIC implementation of the Conway–Sloane E₈ decoder, of Babai nearest-plane as a stand-alone block, or of a Leech decoder with reported numbers surfaced in 2020–2026 searches. The nearest hardware is lattice-reduction-aided MIMO detection (a 4×4 LLL core in 130 nm at 333 MHz, 14 cycles per matrix) and systolic lattice-reduction detectors with word lengths of 14–18 bits. Fixed-point studies of lattice-reduction-aided soft MIMO show that a 16-bit total word length keeps the BER loss under 0.3 dB and that 14 bits costs about 1.3 dB. The lesson for ARALE: the snap stage is a known quantity and the precision budget should be 16–18 bits, not 32.

Application pull comes from three directions. LLM quantization adopted E₈ (E8P, NestQuant, HyperQuant) and then moved to trellis-coded quantization at effective dimension 256 (QTIP), because codebook cost caps lattice vector quantization near eight dimensions; that is a software, not a real-time, use. Coherent optics uses E₈ and Λ₁₆/Λ₂₄ as shaping lattices around a Zⁿ coding lattice, with experimental 1.7 dB OSNR gains. And the one place E₈ enters standardized post-quantum cryptography is as a message encoder in Kyber-style schemes, lowering decryption-failure rate; the underlying ML-KEM arithmetic is NTT-based polynomial multiplication, and profiling on OpenTitan puts 60–77% of cycles in hashing. A closest-point engine does not accelerate ML-KEM, ML-DSA, or HQC. Falcon's fast Fourier sampling is the only standardized algorithm with a nearest-plane relative, and it must be constant-time and exactly Gaussian, which excludes heuristic adaptive hardware.

### 3.2 Bosonic codes and the decoder latency budget

GKP codes are symplectic lattices in 2n-dimensional phase space, and decoding Gaussian shift noise is the closest-point problem in the symplectic dual lattice (Conrad, Eisert, Arzani 2022; Lin, Chamberland, Noh 2023). Two-mode D₄ codes beat square codes by a factor 2¹ᐟ⁴ in distance; four-mode E₈ codes are the natural next step, and a 2025 construction from cryptographic lattices confirms that E₈ and the Leech lattice give the strongest low-dimensional GKP distances. Experiments have reached beyond-break-even single-mode GKP (gain 2.27, Yale 2023) and GKP qudits (gain 1.87, 2024); no four-mode E₈ GKP experiment has been reported.

The decoder budget is set by the surface-code experiments: 1.1 µs per cycle, with FPGA decoders under 1 µs per round at d = 17 and a 12 nm ASIC at 240 ns and 8 mW (Riverlane). A cryo-CMOS predecoder at 4 K drawing under 0.56 mW has been characterized. An E₈ GKP inner decoder would have to run inside that budget and that power envelope.

### 3.3 Cryogenic logic

Cryo-CMOS at 4 K is mature enough to carry control logic: IBM's 14 nm transmon controller draws 23 mW per qubit, and 2025 surveys place typical designs at 1–10 mW per qubit. RSFQ/ERSFQ logic switches at about 2.5 aJ per junction and runs at 58–60 GHz in demonstrated processors, with AQFP at about 15 aJ per operation; the cryocooler penalty is 350–3500 W at room temperature per watt at 4.2 K. Fluxonium is a qubit, not a memory. Classical cryogenic memory remains an open problem (nanowire nTron arrays reach 2.6 Mb/cm² at 1.3 K with BER 10⁻⁵). ARALE's state, an 8×8 matrix at 18 bits, about 1.2 kbit, is small enough to live in registers.

### 3.4 CORDIC, systolic QRD, and rank-one updates

CORDIC's scale factor K = 1.6468 is compensated once at the end; precision is about one bit per iteration plus log₂ n guard bits. Recent FPGA work gives 6-cycle radix-4 latency at 500 MHz and scaling-free six-stage pipelines at 216 MHz with RMS error 6.5×10⁻⁵. The Gentleman–Kung triangular array computes R in A = QR with boundary cells in vectoring mode and internal cells in rotation mode; McWhirter's extension extracts the least-squares residual without back-substitution, and inverse-QR-RLS propagates R⁻¹ so weights emerge directly (a 2017 FPGA design reaches 52-cycle latency and 3.92 MS/s). The RLS covariance recursion is a Sherman–Morrison update, and its known failure modes (loss of symmetry, covariance wind-up without persistent excitation) are documented.

Most important for this design: Hashemi and Nakatsukasa (2025, 2026) show that Sherman–Morrison's instability comes from ill-conditioning of A, not from a small denominator; that forward error scales as ε·κ₂(A)²; and that one step of iterative refinement or a bordered-system reformulation restores backward stability. No 2022–2026 hardware paper implements Sherman–Morrison as a named block. In-memory "one-step" matrix inverters (memristor crossbars) reach under 20 ns settling but at about 6–10-bit precision, below the 16 bits a lattice snap needs.

---

## 4. System overview

ARALE is a three-stage feed-forward decoder closed by one feedback loop. The forward path equalizes, snaps, and emits; the feedback path uses the emitted decision to correct the equalizer by a rank-one update. Nothing in the forward path waits for the feedback path, so decoding latency is fixed and the channel estimate lags by a bounded number of symbols.

```
                     ┌──────────────────────── feedback (rank-one update) ───────────────────────┐
                     │                                                                           │
                     │   e = y(t−D) − Ĥ x̂         Ĥ ← Ĥ + μ e x̂ᵀ / ‖x̂‖²   (gated)             │
                     │                                                                           │
   y (8 lines)       ▼                                                                           │
  ─────────► [A] Front end ──► [B] Equalizer ──────────► z ──────► [C] E₈ snap ──────► x̂ ──────┴──► out
             sample, align     QR array + solve           ≈ x      two D₈ decodes                  (point + two
             16-bit            z = R⁻¹ Qᵀ y                         + compare                      distances)
                                    ▲
                                    │  R (and Qᵀ via rotation history) from Ĥ = QR
                                    └─────────── Ĥ (64 × 18 b) ◄──────────────────────────────────
```

The tracker closes the loop with the decoder's own output, so the channel estimate is corrected at the symbol rate without a training sequence.

**Stage A, front end.** Eight analog lines are sampled to 16-bit fixed point. In a quantum-error-correction setting these are the real and imaginary quadratures of four homodyne measurements; in an optical link they are four complex symbols after carrier recovery; in a memory controller they are eight cell read-outs. The front end's only job is to deliver y with a known, bounded delay. It is the natural home for the one Weyl-semimetal element established at device level: a zero-bias nonlinear-Hall rectifier (TaIrTe₄ class) can serve as a passive envelope detector for the input lines without a bias current.

**Stage B, equalizer.** ARALE carries the channel estimate as a QR factorization Ĥ = QR rather than an explicit inverse. The equalized vector is z = R⁻¹Qᵀy. Qᵀ is applied as 28 planar Givens rotations (n(n−1)/2 for n = 8) in a triangular systolic array of CORDIC rotators, and the triangular solve with R is done by linear-mode CORDIC operations along the diagonal. No multiplier appears in the data path. The explicit-inverse alternative (carry Ĥ⁻¹, update by Sherman–Morrison) is retained as a lower-latency mode for channels with condition number below about 10.

**Stage C, lattice snap.** The equalized z is quantized to the nearest point of E₈ by the Conway–Sloane coset decoder: round z and z − ½·**1** coordinate-wise, correct the parity of each by re-rounding the coordinate with the largest rounding error, and keep the candidate closer to z. This is exact and costs about 120 integer operations with no data-dependent branching beyond two compares, so it fits in a fixed-depth block. The output x̂ is an E₈ point in canonical coordinates; information bits are read off by the coset encoder in reverse.

**Feedback, channel tracker.** The residual e = y − Ĥx̂ is zero when both the estimate and the decision are right. When it is not, the normalized least-mean-squares correction Ĥ ← Ĥ + μ·e·x̂ᵀ/‖x̂‖² is a rank-one update, and the corresponding update to a QR factorization costs O(n²) via two sweeps of Givens rotations, the same CORDIC primitive the forward path uses. This is the central economy of the framework: the equalizer's representation (rotations) and its update mechanism (rotations) are the same hardware.

### 4.1 What each stage guarantees

| Stage | Guarantee | Condition |
|---|---|---|
| Front end | y delivered with fixed delay D₀, 16-bit resolution | ADC ENOB ≥ 12 bits; input within ± full scale |
| Equalizer | z = R⁻¹Qᵀy with relative error below 10⁻³ | κ₂(Ĥ) ≤ 100; 16 CORDIC iterations + 4 guard bits |
| Lattice snap | x̂ = argmin over E₈ of ‖z − x‖, exactly | z finite; no condition on noise |
| Tracker | ‖Ĥ − H‖_F converges to a floor set by μ and the noise | decision-error rate below ~10%; persistent excitation of x |

The guarantees compose: when the tracker's floor keeps ‖Ĥ⁻¹H − I‖_F below about 0.10, the skew of the equalized constellation costs essentially nothing against the exact decoder (§8.2), and the decision-error rate stays low enough for the tracker's own convergence condition. Above that floor the loop can lose lock; §8.4 and §8.7 give capture range and the re-acquisition procedure.

### 4.2 Technology mapping at a glance

ARALE's logic is technology-agnostic. Three mappings are specified in order of readiness:

1. **Cryo-CMOS at 4 K (baseline).** 22 nm FDSOI characterized at 4 K, as in the Pinball predecoder, at 1 GHz. Buildable today; sets the numbers in §9.
2. **SFQ logic at 4 K (speed option).** RSFQ bit-serial arithmetic at 20–60 GHz with about 2.5 aJ per junction switch. CORDIC's shift-add structure maps well to bit-serial SFQ; the cost is roughly 10× larger area per bit of state and no dense same-temperature memory, which ARALE avoids by keeping only ~1.2 kbit of state.
3. **Topological-semimetal interconnect (materials upgrade).** The array's dense local wiring is replaced with sputtered NbP (400 °C, BEOL-compatible) or CoSi lines below 5 nm, where resistivity falls with thinning instead of rising. This changes wire RC, not logic, and is the only use of Weyl materials the 2024–2026 device literature supports.

The chirality-routed deserializer is a hypothesis with a defined experiment (§7.4). It is not required for ARALE and no performance claim depends on it.

---

## 5. Mathematical core

The decoder is zero-forcing equalization followed by exact E₈ quantization, with the equalizer's factorization maintained by rank-one, decision-directed updates.

### 5.1 The E₈ lattice and its coset decoder

In the even coordinate system,

```
E₈ = { x ∈ Z⁸ ∪ (Z+½)⁸ : Σᵢ xᵢ ≡ 0 (mod 2) } = D₈ ∪ (D₈ + ½·1),
```

where D₈ is the set of integer vectors with even coordinate sum. E₈ is unimodular (det M = 1 for any generator M). Its 240 minimal vectors have squared norm 2: 112 of the form (±1, ±1, 0⁶) and 128 of the form (±½)⁸ with an even number of minus signs. Its packing radius is ρ = 1/√2, and its nominal coding gain over Z⁸ is γ = d²_min / (det M)^{2/8} = 2, i.e. 3.01 dB.

A convenient lower-triangular generator (rows are basis vectors):

```
M = [  2    0    0    0    0    0    0    0  ]
    [ −1    1    0    0    0    0    0    0  ]
    [  0   −1    1    0    0    0    0    0  ]
    [  0    0   −1    1    0    0    0    0  ]
    [  0    0    0   −1    1    0    0    0  ]
    [  0    0    0    0   −1    1    0    0  ]
    [  0    0    0    0    0   −1    1    0  ]
    [ ½    ½    ½    ½    ½    ½    ½    ½  ]
```

with lattice points x = Mᵀb for b ∈ Z⁸. The generator is needed only for encoding (bits to lattice point) and for reading bits back from a decoded point; the decoder never inverts it. The model verifies `det M = 1.000000` and that the Gram matrix has an even diagonal, `[4, 2, 2, 2, 2, 2, 2, 2]` (all even), as an even lattice requires.

The closest point Q_E₈(z) for z ∈ R⁸ is found as follows. Define Q_D₈(z): round each coordinate to the nearest integer; if the sum is odd, find the coordinate whose rounding error is largest in magnitude and round it the other way. Then

```
Q_E₈(z) = argmin { ‖z − c‖ : c ∈ { Q_D₈(z),  Q_D₈(z − ½·1) + ½·1 } }.
```

```python
def dec_D8(z):
    r = rint(z)
    if int(r.sum()) % 2 != 0:                 # odd parity: re-round the worst coordinate
        err = z - r
        k = argmax(abs(err))
        r[k] += 1.0 if err[k] > 0 else -1.0
    return r

def dec_E8(z):
    a = dec_D8(z)                             # integer coset
    b = dec_D8(z - 0.5) + 0.5                 # half-integer coset
    return a if norm(z - a) <= norm(z - b) else b
```

This is exact, not approximate: the two cosets partition E₈, and Q_D₈ is exact for D₈. The model confirms it against brute-force search over a neighbourhood of lattice points: **agreement on 2000 of 2000 test points** near the origin.

**Operation count per vector.** 16 roundings, 16 subtractions for rounding errors, two 8-way argmax searches, two parity checks, two conditional corrections, two squared-distance accumulations (16 squarings of values below 1, or 16 table lookups at 8-bit precision), and one compare: about 120 operations, all integer or short fixed-point, on a fixed data-independent schedule. In hardware this is a 4–5 stage combinational block.

### 5.2 Channel model and decoding rule

The observation is y = Hx + w with x ∈ αE₈, w ~ N(0, σ²I₈), H invertible and slowly varying. Maximum-likelihood decoding is the closest-vector problem in the skewed lattice H·αE₈, which has no fast decoder in general. ARALE computes

```
x̂ = α · Q_E₈( (1/α) · Ĥ⁻¹ y ),
```

with Ĥ the tracked estimate of H. Write Ĥ⁻¹H = I + Δ for the residual mismatch; the equalized vector is z = x + Δx + Ĥ⁻¹w. Two losses relative to ML arise.

*Noise enhancement.* The noise covariance after equalization is σ²Ĥ⁻¹Ĥ⁻ᵀ, whose largest eigenvalue is σ²/σ²_min(Ĥ). For a well-conditioned channel (κ₂(H) ≤ 3) this costs at most 20·log₁₀κ₂ ≈ 9.5 dB in the worst direction but under 1 dB on average when H is a small perturbation of a scaled orthogonal matrix, which is the physical case (gain imbalance and small crosstalk).

*Mismatch.* The term Δx is a signal-dependent offset of norm at most ‖Δ‖₂‖x‖. For the 240 minimal vectors ‖x‖ = α√2 and the packing radius is α/√2, so a mismatch ‖Δ‖₂ = 0.05 consumes 10% of the decision margin in the worst direction, a loss of 20·log₁₀(1/0.9) ≈ 0.9 dB, and about 0.3 dB averaged over the constellation. That is the worst-case direction. Averaged over isotropic mismatch and the constellation, the Monte Carlo of §8.2 finds a loss of at most about 0.05 dB at ‖Δ‖_F = 0.10 and about 0.12–0.16 dB at 0.20. The design budget is therefore ‖Δ‖_F ≤ 0.10 for isotropic mismatch; the 0.05 directional bound is retained for the structured case of a single-line gain error, which is the physical drift mode.

When H is a scalar multiple of an orthogonal matrix, zero-forcing is exactly ML, since Ĥ⁻¹ then preserves Euclidean distances up to scale. The decoder is near-ML precisely in the regime it is designed for, and its degradation is graceful and computable as κ₂(H) grows.

### 5.3 Rank-one drift model and decision-directed tracking

Physical drift is low-rank per event. A gain change on line i changes row i of H (rank one); a crosstalk change between lines i and j changes two entries (rank at most two); a common phase rotation of one complex pair is a 2×2 block (rank two). ARALE therefore models H(t) as H₀ plus a sequence of rank-one increments and tracks it by normalized least-mean-squares driven by the decoder's own decisions:

```
e_t = y_t − Ĥ_t x̂_t ,        Ĥ_{t+1} = Ĥ_t + μ · e_t x̂_tᵀ / ‖x̂_t‖² .
```

This is a rank-one update with u = μe_t/‖x̂_t‖², v = x̂_t. Standard NLMS results apply: the update is stable for 0 < μ < 2, the misadjustment is about μ/(2 − μ), and the tracking time constant is about ‖x̂‖²/(μ·λ_avg) symbols, where λ_avg is the mean eigenvalue of the input correlation. With μ = 2⁻⁴ the misadjustment is 3.2%. Because x̂ ∈ αE₈ has coordinates in α·{0, ±½, ±1, ±3/2, ±2} for any codebook of practical size, the outer product e·x̂ᵀ is computed by shifts and adds alone, and μ = 2⁻ᵏ is a shift: the tracker contains no multiplier.

Decision-directed tracking is self-consistent only while decision errors are rare. A wrong x̂ injects an update in a wrong direction of size μ‖e‖/‖x̂‖. For symbol error rates p ≲ 10⁻² the noise share of the update energy is below the NLMS misadjustment and the loop converges; above p ≈ 10⁻¹ it does not, and the tracker must fall back to a pilot-aided or decision-gated mode (§8.7).

### 5.4 Explicit-inverse representation: Sherman–Morrison

If the decoder carries Ĥ⁻¹ directly, the rank-one update of Ĥ induces

```
Ĥ⁻¹_{t+1} = Ĥ⁻¹_t − (Ĥ⁻¹_t u vᵀ Ĥ⁻¹_t) / (1 + vᵀ Ĥ⁻¹_t u),     u = μ e_t/‖x̂_t‖²,   v = x̂_t.
```

Cost per update at n = 8: two matrix-vector products (128 multiply-adds), one inner product (8), one division, one outer product (64), and 64 subtractions; about 270 operations, of which 192 are general multiplies that cannot be reduced to shifts because Ĥ⁻¹ is dense. Applying the equalizer is 64 more multiply-adds per symbol. The denominator 1 + μ·x̂ᵀĤ⁻¹e/‖x̂‖² is close to 1 whenever μ is small and Ĥ is well-conditioned, so the division is a one-step Newton correction from 1.

The stability result that governs this mode: for a solve with a rank-one-updated matrix computed by Sherman–Morrison in precision ε, the forward error is of order ε·κ₂(Ĥ)² and the residual of order ε·κ₂(Ĥ), regardless of how well-conditioned Ĥ + uvᵀ is, and this is *not* caused by a small denominator. The model measures the actual growth in §8.5; it is linear in the number of updates and, over the tested range, sub-quadratic in κ₂. ARALE schedules a full re-derivation of Ĥ⁻¹ from the QR path every N_refresh updates, with N_refresh set from the measured κ₂.

### 5.5 Orthogonal representation: QR and its rank-one update

The stable alternative is to carry Ĥ = QR and never form an inverse. Equalization is z = R⁻¹(Qᵀy): apply Qᵀ as a product of planar rotations, then back-substitute. The Givens rotation in the (i, j) plane by angle θ is

```
G(i,j,θ) = I + (cos θ − 1)(eᵢeᵢᵀ + eⱼeⱼᵀ) − sin θ (eᵢeⱼᵀ − eⱼeᵢᵀ).
```

A QR factorization of an 8×8 matrix by Givens rotations uses 28 of them, one per sub-diagonal entry, in an order the triangular systolic array realizes with 8 boundary cells and 28 internal cells. Each rotation is one CORDIC operation: in vectoring mode the boundary cell drives a pair (a, b) to (√(a²+b²), 0) and records the micro-rotation direction bits; in rotation mode the internal cells replay those bits on their own pairs. The angle is never computed, stored, or multiplied; what passes between cells is a 16-bit direction vector.

The rank-one update of a QR factorization (Gill, Golub, Murray, Saunders 1974) proceeds in two sweeps. With Ĥ + uvᵀ = Q(R + ũvᵀ), ũ = Qᵀu: a first sweep of n − 1 rotations from the bottom reduces ũ to ‖ũ‖e₁ and turns R into an upper Hessenberg matrix; adding ‖ũ‖e₁vᵀ keeps it Hessenberg; a second sweep of n − 1 rotations restores upper-triangular form. Total: 2(n − 1) = 14 vectoring operations, each applied across one row pair of R (8 entries) and one column pair of Q (8 entries), about 14 × 16 = 224 rotation-mode operations. Every step is orthogonal, so the update is backward stable and error after k updates grows only linearly in kε, not in κ².

**What the rank-one structure buys at n = 8.** A full Givens QR of Ĥ from scratch costs 28 vectoring and Σₖ₌₁⁷ k² = 140 rotation operations. The rank-one update costs 14 vectoring and about 224 rotation operations if Q is maintained explicitly, or 14 vectoring and 112 if only R is maintained and Qᵀy is formed by replaying the rotation history. The asymptotic advantage O(n²) versus O(n³) is a factor of about 1.5–2 at this dimension, not an order of magnitude. The design decision follows: at n = 8 the array simply re-factorizes Ĥ on every update, which is O(n) latency on n²/2 cells, and the rank-one structure is exploited where it is unambiguous: in the shift-add formation of the update to Ĥ itself and in the explicit-inverse mode. The framework scales to n = 16 and n = 24, where the rank-one QR update becomes decisively cheaper than re-factorization, with unchanged array geometry (§9.6).

### 5.6 Why not Babai nearest-plane on the skewed basis

An alternative is to run Babai's nearest-plane algorithm directly on B = ĤM with its QR factorization, avoiding the equalization step. For a generic lattice this is the method of choice. For E₈ it is worse than equalize-then-snap: nearest-plane on a non-orthogonal basis only approximates the closest point, with an approximation factor that grows with the basis's orthogonality defect, whereas the coset decoder is exact. Since Ĥ is close to a scaled orthogonal matrix in the operating regime, equalization costs little and preserves exactness of the snap. Nearest-plane is retained as the **fallback decoder during initial acquisition**, when Ĥ is poor and no decisions can be trusted.

### 5.7 Precision

CORDIC with N iterations resolves N bits of angle, loses about log₂N bits to accumulated rounding, and scales the result by K = ∏ₖ√(1 + 2⁻²ᵏ) → 1.646760, compensated once by a shift-add multiplication by 1/K = 0.60725. For 16-bit data and N = 16 iterations, a 20-bit internal datapath (16 + 4 guard bits) keeps the rotation error near 4×10⁻⁵ per operation (measured in §8.6). Fixed-point studies of QR-based lattice-reduction MIMO detectors confirm that a 16-bit total word length holds the BER loss below 0.3 dB and that 14 bits costs about 1.3 dB. ARALE adopts **16-bit input, 20-bit internal, 18-bit stored Ĥ, and a 16-iteration CORDIC.** A 32-bit datapath would double area and latency for a precision no stage of the decoder can use.

---

## 6. Hardware architecture

The engine is three tiers in one 4 K package: a **substrate tier** that moves the eight input lines and the array's local wiring, an **adaptive tier** that holds Ĥ and forms its rank-one update, and a **logic tier** that factorizes, equalizes, and snaps. The logic tier is the bulk of the area and all of the latency.

The array is the Gentleman–Kung triangular QRD array with y appended as a ninth column, so one pass through the array yields both R and Qᵀy. The solve row and the snap block complete the decode without any multiplier in the data path.

```
          column →   1      2      3     ...     8      y
  row 1             [B]───►[I]───►[I]───►...───►[I]───►[Y]
                     │      │      │             │      │
  row 2                    [B]───►[I]───►...───►[I]───►[Y]
                            │      │             │      │
  row 3                           [B]───►...───►[I]───►[Y]
   ...                                   ...
  row 8                                         [B]───►[Y]
                                                  │      │
                                                  ▼      ▼
                       solve row (linear-mode CORDIC, right → left)  ──►  z ──► snap ──► x̂

  B = boundary cell (vectoring), I = internal cell (rotation), Y = y-column cell (rotation)
  8 B + 28 I + 8 Y = 44 array cells.  Direction bits flow rightward; partial results flow downward.
```

### 6.1 Cell microarchitecture

**Boundary cell (8 instances).** A 16-iteration circular-mode CORDIC in vectoring mode with a 20-bit datapath. Input: the pair (r_ii, a), where r_ii is the stored diagonal and a is the incoming element. Output: the new diagonal √(r_ii² + a²)/K (K-compensation folded into a final shift-add stage) and a 16-bit direction vector d ∈ {−1, +1}¹⁶ broadcast rightward. The cell is a 16-stage pipeline: a new pair enters every clock, and the direction bits emerge 16 clocks later. Area per cell is three 20-bit adders per stage × 16 stages plus pipeline registers, comparable to the 790-LUT reference pipeline reported for a Q1.15 16-stage CORDIC on a Zynq-7000.

**Internal cell (28 instances) and y-column cell (8 instances).** The same 16-stage shift-add pipeline in rotation mode, consuming the direction vector from the left and applying the identical micro-rotation sequence to its own pair (r_ij, b). Because the direction bits are known one stage ahead, the internal cell has no decision logic at all; its critical path is one 20-bit addition. The y-column cells are internal cells whose stored value is the running Qᵀy component.

**Solve row (8 cells per row).** Back-substitution z_i = (q_i − Σ_{j>i} r_ij z_j)/r_ii, executed right to left by linear-mode CORDIC cells: the multiply-accumulate is a linear-mode rotation and the division is a linear-mode vectoring, both shift-add. Latency is 8 cells × (16 + 2) clocks when unpipelined. ARALE pipelines two vectors through the row by interleaving even and odd symbols, giving one vector per 9 clocks per row, and three rows in parallel sustain one vector per 3 clocks (reduced configuration). For the full-rate one-vector-per-clock target the row is replicated nine times (72 cells), the single largest area item after the array.

**Snap block.** A four-stage combinational pipeline: (1) round z and z − ½·**1**; (2) compute rounding errors and parities; (3) find the max-error coordinate in each candidate and apply the parity correction; (4) compute the two squared distances (16 squarers on 8-bit residuals, implemented as lookup) and select. No feedback, no multiplier, one vector per clock. The two squared distances are exported alongside x̂ because their difference is the likelihood information a concatenated outer decoder needs.

### 6.2 Adaptive tier

The adaptive tier holds Ĥ as 64 words of 18 bits (1,152 bits) in flip-flops, and a delay line of D₀ + D₁ ≈ 70 samples of y (8 × 16 × 70 ≈ 9 kbit) so that the residual e_t = y_t − Ĥx̂_t can be formed when x̂_t emerges. Forming Ĥx̂_t for an E₈ point whose coordinates lie in α·{0, ±½, ±1, ±3/2, ±2} is a sum of shifted columns of Ĥ: 64 shift-add operations. The outer-product update μ·e·x̂ᵀ/‖x̂‖² is likewise 64 shift-adds plus one normalization by ‖x̂‖², which takes only four distinct values for the minimal-vector codebook and is a table lookup. An update gate (§8.7) decides whether the update is applied or discarded.

In explicit-inverse mode the same tier holds Ĥ⁻¹ and the Sherman–Morrison update needs 192 general multiplies, implemented as linear-mode CORDIC multiplies on a time-shared 16-cell bank taking 12 clocks per update. The equalizer in this mode is a dense 8×8 matrix-vector product, which the same bank executes at 4 clocks per vector. This mode is for channels with κ₂(Ĥ) < 10 and update rates below one per 20 symbols; outside that envelope the QR mode is used.

### 6.3 Substrate tier

The substrate tier is the package, the input lines, and the on-chip interconnect. Baseline: 22 nm FDSOI cryo-CMOS, where 4 K characterization exists and threshold shifts are modeled. Two upgrades are specified:

- **Topological-semimetal local interconnect.** Sputtered noncrystalline NbP at 400 °C (BEOL-compatible) or CMOS-compatible CoSi, used for the array's dense local mesh where wire pitch is below 20 nm. The measured resistivity of NbP films falls to about 34 µΩ·cm at 1.5 nm, against roughly 100 µΩ·cm for Cu at that thickness, so mesh RC delay is reduced by about 3× in the thinnest wires.
- **Zero-bias front-end rectification** by a nonlinear-Hall element (TaIrTe₄ class, demonstrated at 0.1–10 THz and 0.3 A/W at room temperature) for envelope detection of the eight input lines without bias current, which removes eight DC bias loads from the 4 K budget.

### 6.4 Explicit-inverse mode versus QR mode

| Property | QR mode (default) | Explicit-inverse mode |
|---|---|---|
| Stored representation | R and Qᵀy via rotations | Ĥ⁻¹ (64 × 18 b) |
| Condition-number envelope | κ₂ ≤ 100 | κ₂ < 10 |
| Error growth per update | none (rounding floor only) | linear in updates; ≈ 3×10⁻⁶ to 8×10⁻⁵ per update (κ₂ = 2 to 30) |
| Latency | 59 clocks | 4 + 16 + 4 = 24 clocks |
| Throughput | 1 vector/clock (full rate) | 0.25 vector/clock (16-cell bank) |
| Refresh | not needed | re-derive from QR every N_refresh updates |
| Hardware | CORDIC array | time-shared CORDIC multiply bank |

### 6.5 Dataflow timing

At a 1 GHz clock the forward path is:

| Stage | Clocks | Note |
|---|---|---|
| Front-end sample and align | 4 | ADC pipeline, skew registers |
| Array fill (9 columns skewed) | 9 | one column per clock |
| Boundary CORDIC depth | 16 | 16 iterations |
| Propagation across 8 rows | 8 | diagonal wavefront |
| Solve row | 18 | pipelined, interleaved |
| Snap block | 4 | combinational pipeline |
| **Total latency** | **59** | **≈ 59 ns at 1 GHz** |

Throughput is one 8-vector per clock in QR mode once the pipeline is full, provided the solve row is replicated; with a single solve row throughput drops to one vector per 3 clocks. The feedback loop adds D₁ = 59 clocks of decision delay plus 6 clocks to form e_t and the update, so the channel estimate used for symbol t was formed from decisions up to symbol t − 65. At μ = 2⁻⁴ the tracker's time constant is 16 symbols, and the normalized loop gain μ·D = 65/16 × 2⁻⁴ ≈ 0.25 stays well below 1; the model (§8.3) confirms the loop is stable at every tested μ with this delay.

### 6.6 Precision and state budget

| Signal | Width | Reason |
|---|---|---|
| Input y | 16 bits | ADC ENOB 12–14 plus headroom |
| CORDIC datapath | 20 bits | 16 + log₂(16) guard bits |
| Stored R, Ĥ | 18 bits | 2 bits over input to absorb update rounding |
| Direction vector | 16 bits | one bit per iteration |
| Solve row output z | 16 bits | matched to snap input |
| Snap internals | 8-bit residuals | distances need only the fractional part |

Total state: 1,152 bits of Ĥ, 648 bits of R and Qᵀy in the array, ≈ 9 kbit of y delay line, and ≈ 2 kbit of pipeline registers. This is well under the capacity of flip-flops in a cryo-CMOS block, so **no cryogenic memory technology is required**.

---

## 7. Physical substrate and technology mapping

Weyl semimetals enter ARALE in three roles with three different evidence levels: as a low-resistivity nanoscale interconnect (**established** in 2024–2026 device experiments), as a zero-bias rectifying front end (**established at room temperature, untested at 4 K**), and as a chirality-routed deserializer (**a hypothesis with no device demonstration**). This section states the physics each role relies on and marks the boundary of what is known.

### 7.1 Transport physics that is established

A Weyl semimetal has pairs of band-touching points of opposite chirality in momentum space. Each node is a monopole of Berry curvature and can be removed only by annihilation with its partner, which is the precise sense in which the nodes are topologically protected. Linearly dispersing bulk states near the nodes have Fermi velocities of 3–8×10⁵ m/s in the TaAs family (quantum-oscillation measurements on TaAs give 2.96×10⁵ and 7.55×10⁵ m/s for the two pockets), two to four times below graphene's 1.1–2.5×10⁶ m/s. The surface carries open Fermi arcs connecting the projections of opposite-chirality nodes. Under parallel electric and magnetic fields, charge is pumped between nodes at a rate proportional to **E·B** and relaxed by internode scattering. The resulting negative longitudinal magnetoresistance is the transport signature of the chiral anomaly, but in high-mobility crystals it is routinely mimicked by current jetting, so only measurements that vary contact geometry and exclude the artifact are reliable; a 2025 TaRhTe₄ study that did so found a negative longitudinal magnetoresistance of only about 0.2% at 10 T.

Three quantitative facts constrain any device claim:

1. **Band velocity is not drift velocity.** At the ~10³ cm²/Vs mobilities of epitaxial TaAs films and the roughly 10 cm²/Vs surface mobilities measured in sputtered NbP films, drift velocity in practical fields is far below v_F, and signal propagation across a 1 µm device is set by RC, not by v_F.
2. **Backscattering suppression is narrow.** It applies to intranode, chirality-conserving processes on a clean surface. It does not suppress internode scattering, phonon scattering, or scattering from grain boundaries and interfaces, which dominate in thin films.
3. **Node protection places no bound on resistivity, noise, or arithmetic error.** ARALE therefore makes no claim of immunity to scattering or to computational error from the substrate; its error budget is the fixed-point analysis of §8, unchanged by the material.

### 7.2 Role 1, interconnect: established

The strongest device result in the 2024–2026 literature is that the resistivity of nanoscale topological-semimetal conductors *decreases* as they shrink, the opposite of copper. Sputtered noncrystalline NbP films grown at 400 °C on SiO₂/Si with a 1.4–4 nm Nb seed show about 34 µΩ·cm at 1.5 nm, below 45 µΩ·cm under 3 nm, against 60–70 µΩ·cm for bulk NbP and roughly 100 µΩ·cm for conventional metals at the same thickness; from 20 nm to 5 nm the sheet resistance of Cu rises 10–100× while NbP's rises less than 2× (Science 2024). Single-crystal NbAs nanowires at 40 nm show 9.7 ± 1.6 µΩ·cm, three to four times below bulk, beating Co and Ru at that size but not yet 10 nm Cu, with the crossover projected near 12 nm (Science 2026). CoSi is the CMOS-compatible member of the family, with films to 5 nm and wafer-scale nanowires reported. Surface-state conduction is the proposed mechanism, supported by first-principles transport showing that surface states carry about 76% of the ballistic conductance in 2 nm NbAs slabs.

ARALE uses this where it helps: the local mesh of the 44-cell array and solve rows, whose wire pitch is set by cell density and falls below 20 nm in a 22 nm node. The gain is a wire-RC reduction of order 3× in the thinnest wires, which lowers the array's wire-dominated critical path and gives the 1 GHz target margin (or a path to 1.5–2 GHz). The limits are also clear: NbP films are noncrystalline with low surface mobility, so the benefit is purely resistive; the sputtering route is BEOL-compatible but nanomolded NbAs is not, and arsenic handling is a fabrication cost. **No part of ARALE's performance model depends on this upgrade;** every number in §9 is stated without it.

### 7.3 Role 2, zero-bias rectifying front end: established at 300 K, open at 4 K

The Berry-curvature dipole of a low-symmetry Weyl semimetal produces a second-order (nonlinear) Hall current under an AC drive with no bias and no magnetic field. TaIrTe₄ demonstrates room-temperature RF rectification and a terahertz detector with responsivity 0.3 A/W (18 A/W with gate-tuned correlations), noise-equivalent power near 1 pW/√Hz, and picosecond intrinsic response across 0.1–10 THz; NbIrTe₄ rectifies 20–820 GHz without junction, bias, or field. ARALE's front end can use such an element as a passive envelope detector per input line, eliminating eight DC bias currents from the 4 K stage, where each milliwatt is a material fraction of the 0.7–1.5 W budget.

What is open: the nonlinear Hall response in TaIrTe₄ changes sign near 175 K, and no 4 K characterization of responsivity, noise, or linearity for this use exists. A conventional cryo-CMOS front end is the fallback and costs nothing in the model.

### 7.4 Role 3, chirality-routed deserializer: hypothesis

A tempting extension is to route time-multiplexed input pulses into eight parallel channels by chiral-anomaly pumping, replacing a shift register. The elements that would make this work exist individually: nonlocal valley-polarization transport over 2–7 µm has been measured in gated PtSe₂ at cryogenic temperature with ON/OFF ratios near 10³, but only at 6–9 T; room-temperature nonlocal chiral charge pumping has been reported in Cd₃As₂; and a field-effect chirality device has performed AND/OR logic with gate voltage and magnetic field as inputs. No experiment has routed a signal by temporal sequence into spatially distinct channels, no demonstration exists below 1 T, and a multi-tesla field at the 4 K stage is incompatible with superconducting qubits in the same cryostat.

ARALE therefore does not depend on this mechanism. It is recorded as a hypothesis with a falsifiable test: a four-terminal PtSe₂ or TaRhTe₄ device with two chirality-selective contacts, driven by a 100 MHz pulse train under **E**∥**B** at 1 T and 4 K, should show channel-selective nonlocal voltage with a **contrast ratio above 10 and a settling time below 10 ns**. Failure of either number rules it out as a deserializer at ARALE's symbol rate. Until that experiment succeeds, deserialization is done by a cryo-CMOS shift register, which costs 8 × 16 flip-flops and no performance.

### 7.5 Fabrication stack and evidence levels

| Layer | Material and process | Temperature | Status |
|---|---|---|---|
| Logic and adaptive tiers | 22 nm FDSOI cryo-CMOS | standard BEOL | characterized at 4 K |
| Speed option | Nb/AlOₓ/Nb RSFQ, ERSFQ bias | 4.2 K operation | demonstrated at 58–60 GHz |
| Local interconnect upgrade | sputtered NbP on Nb seed, SiNₓ cap | 400 °C | BEOL-compatible, demonstrated |
| Local interconnect alternative | CoSi films and nanowires | Si-compatible | demonstrated, numbers not public |
| Front-end rectifier | exfoliated or CVD TaIrTe₄ flakes | room-temperature transfer | demonstrated at 300 K only |
| Deserializer hypothesis | gated PtSe₂ or TaRhTe₄ Hall bar | cryogenic, 1 T | no device |

Epitaxial TaAs (MBE on GaAs at 640–650 °C, 16% lattice mismatch, 10²⁰ cm⁻³ unintentional doping) and nanomolded NbAs are excluded from the stack: neither is BEOL-compatible and neither offers a transport advantage over sputtered NbP for a resistive interconnect.

---

## 8. Numerical results (Python model)

All results in this section are produced by `arale_model.py` (numpy only). Each subsection names the command that reproduces it. Four error sources act on the decoder: CORDIC rounding in the array, growth of representation error in the tracked channel estimate, the tracker's own noise floor and lag, and decision errors fed back into the tracker. Each is bounded and measured below.

**What the model is, and is not.** The decoder, the SER experiments, the CORDIC test, and the Sherman–Morrison and Givens-QR experiments model the stated fixed-point word lengths explicitly (rounding to 2⁻¹⁶ or 2⁻²⁰ after each operation). The closed-loop tracker simulation (§8.3–8.4) uses floating-point arithmetic for H, Ĥ and the NLMS update, with the 65-symbol decision delay and the 240-point minimal-vector codebook, so it isolates the loop dynamics from word-length effects, which §8.5–8.6 quantify separately. The update gate and re-acquisition fallback of §8.7 are design policy; the tracker simulation runs the plain decision-directed loop without them, so its capture-range results are conservative with respect to the gated hardware. All tracking runs use a fixed seed and are single realizations; the SER experiments report binomial standard errors.

**SNR convention.** SNR is per dimension with signal energy Es = 2 (the squared norm of a minimal vector): σ² = Es / (8·10^(SNR/10)), i.e. SNR/dim = (Es/8)/σ² with Es/8 = 0.25.

### 8.1 Decoder verification (`python3 -I arale_model.py decoder`)

```
[decoder] agrees with brute force on 2000/2000 points near the origin
[decoder] generator determinant = 1.000000 (unimodular expected: 1)
[decoder] Gram diagonal = [4. 2. 2. 2. 2. 2. 2. 2.]  (all even: True)
```

Test points are drawn as z ~ N(0, 0.25²·I₈); for these the closest lattice point is the origin or one of the 240 minimal vectors, so brute force over those 241 candidates is exhaustive. The coset decoder agrees on every point. The generator has determinant 1 (unimodular) and an even Gram diagonal (even lattice).

### 8.2 Symbol error rate and mismatch tolerance (`ser`, `ser-precise`)

The experiment draws x uniformly from the 240 minimal vectors, forms z = x + Δx + w with Δ an isotropic Gaussian matrix scaled to a prescribed Frobenius norm ‖Δ‖_F and w white Gaussian, snaps with the E₈ decoder, and counts errors. It isolates the mismatch term of §5.2; noise enhancement from a non-unitary H is exercised in the tracking experiments, where the channel has κ₂ = 2.

**20,000 trials per cell** (`python3 -I arale_model.py ser`):

| SNR/dim (dB) | exact | ‖Δ‖_F = 0.05 | ‖Δ‖_F = 0.10 | ‖Δ‖_F = 0.20 |
|---|---|---|---|---|
| 6 | 0.2003 | 0.2001 | 0.2043 | 0.2150 |
| 7 | 0.0897 | 0.0898 | 0.0917 | 0.0963 |
| 8 | 0.0266 | 0.0296 | 0.0316 | 0.0328 |
| 9 | 0.0062 | 0.0070 | 0.0057 | 0.0084 |
| 10 | 0.0009 | 0.0010 | 0.0008 | 0.0015 |

At 20,000 trials the binomial standard error is 0.0028 at 6 dB, 0.0020 at 7 dB, 0.0011 at 8 dB, 0.0005 at 9 dB, and 0.0002 at 10 dB, so differences between neighbouring columns are mostly sampling noise. The 8 and 9 dB rows are therefore repeated at five times the sample size.

**100,000 trials per cell, ± one standard error** (`python3 -I arale_model.py ser-precise`):

| SNR/dim (dB) | exact | ‖Δ‖_F = 0.05 | ‖Δ‖_F = 0.10 | ‖Δ‖_F = 0.20 |
|---|---|---|---|---|
| 8 | 0.0277 ± 0.0005 | 0.0284 ± 0.0005 | 0.0299 ± 0.0005 | 0.0334 ± 0.0006 |
| 9 | 0.0058 ± 0.0002 | 0.0064 ± 0.0003 | 0.0063 ± 0.0002 | 0.0079 ± 0.0003 |

Converting SER ratios to equivalent SNR loss uses the local slope of the exact curve (a factor of about 4.7 per dB near 8 dB and about 7 per dB near 9 dB). The results:

- **‖Δ‖_F = 0.05:** no loss distinguishable from sampling noise at either SNR.
- **‖Δ‖_F = 0.10:** SER ratio 1.08 at 8 dB and 1.09 at 9 dB, equivalent to about **0.05 dB or less**. The 8 dB difference is roughly three standard errors; the 9 dB difference is under two. The loss is real but small.
- **‖Δ‖_F = 0.20:** SER ratio 1.21 at 8 dB and 1.36 at 9 dB, equivalent to **about 0.12–0.16 dB**.

Across independent runs the single-cell values at a given nominal setting differ by one to two standard errors, which is why the 100,000-trial rows are the reference. Reading the curve against the E₈ geometry: at 10 dB the exact error rate is 8×10⁻⁴ and no errors occur in 20,000 trials at 12 dB in the tracking runs; the 3.01 dB nominal coding gain of E₈ over Z⁸ is preserved by a decoder whose tracked mismatch stays under 0.10.

The design budget of ‖Δ‖_F ≤ 0.10 for isotropic mismatch follows. The 0.05 directional bound of §5.2 is retained for the structured single-line gain error, which projects onto the decision margin more efficiently than a random Δ.

### 8.3 Tracker noise floor, lag, and stability with delay (`track`, `track-snr`)

The NLMS tracker has a noise floor set by μ and the noise power and a lag set by 1/μ. The channel has κ₂ = 2; the loop includes a 65-symbol decision delay (the forward latency of §6.5 plus update formation); SNR is 12 dB; at symbol 1000 a 10% gain step is applied to one line of H (all columns of a single row).

`python3 -I arale_model.py track` (first block):

| μ | Steady ‖Δ‖_F floor | Peak after the 10% step | 1/e recovery (symbols) | Decision errors during recovery |
|---|---|---|---|---|
| 2⁻³ | 0.126 | 0.181 | 87 | 0 |
| 2⁻⁴ | 0.087 | 0.138 | 105 | 0 |
| 2⁻⁵ | 0.060 | 0.119 | 264 | 0 |
| 2⁻⁶ | 0.040 | 0.110 | 385 | 0 |

The floor is the median of ‖Ĥ⁻¹H − I‖_F over symbols 500–1000, the peak is the maximum after the step, and the recovery time is the first symbol at which the mismatch falls below base + (peak − base)/e. The floor scales as √μ as NLMS theory predicts: the ratios between successive rows are 1.45, 1.45 and 1.50 against √2 = 1.41.

**Floor versus SNR** at μ = 2⁻⁴ (`python3 -I arale_model.py track-snr`, mean of ‖Δ‖_F over symbols 1500–3000):

| SNR/dim (dB) | Floor ‖Δ‖_F | Decision errors in 3,000 symbols |
|---|---|---|
| 8 | 0.147 | 5 |
| 10 | 0.116 | 0 |
| 12 | 0.092 | 0 |
| 14 | 0.073 | 0 |

The floor falls as σ, as expected, and the loop stays stable at every SNR tested. (The 12 dB entry, 0.092, is a mean over a later window than the 0.087 median in the table above; the two statistics agree to within 6%.) At 8 dB the floor of 0.147 exceeds the 0.10 budget and the few decision errors visible there are the first sign of the regime in which the update gate and pilot fallback matter; §13.2 lists operating above 9 dB as the residual mitigation.

**Stability with delay.** The loop is stable with the 65-symbol delay at every μ tested because the normalized loop gain stays below 1: with μ = 2⁻⁴ the open-loop time constant is 16 symbols, and D = 65 gives μ·D/τ ≈ 0.25.

**Default operating point.** ARALE's default is μ = 2⁻⁴ at the symbol rate, which holds the floor near 0.09 at 12 dB and recovers from a 10% step in about 100 symbols (≈ 100 ns at 1 Gsymbol/s) with no decision errors. A slower μ = 2⁻⁶ is selected when the drift rate is known to be low, trading a ≈ 4× longer recovery for a floor of 0.04. Because the step size is a shift, switching between them is free. For the 2¹⁶-point Voronoi codebook, whose vectors reach squared norm 10, the mismatch term scales with ‖x‖ (√10 against √2), so the same absolute floor costs more: roughly 0.3 dB at μ = 2⁻⁴ and under 0.1 dB at μ = 2⁻⁶ (a scaling estimate, not a simulation), and the slower step is the default for that codebook.

### 8.4 Capture range (`track`, second block)

Starting from an erroneous Ĥ, with the decoder running and the channel static, μ = 2⁻⁴, SNR 12 dB, 240-point codebook:

| Initial error in Ĥ | Initial ‖Δ‖_F | Symbols to ‖Δ‖_F < 0.10 | SER in first 200 symbols |
|---|---|---|---|
| Single-line gain error of 30% | 0.248 | 268 | 0.000 |
| Single-line gain error of 120% | 0.586 | 676 | 0.000 |
| Rotation error of 20° in one coordinate pair | 0.491 | 431 | 0.000 |
| Rotation error of 45° in one coordinate pair | 1.082 | 780 | 0.075 |
| Random full-rank error, ‖·‖_F = 0.5 | 0.503 | 382 | 0.000 |
| Random full-rank error, ‖·‖_F = 0.8 | 0.822 | 439 | 0.005 |

The decision-directed loop's capture range is wide, so the first response to lost lock is to do nothing different. Single-line errors are the easiest case because E₈'s parity structure corrects any distortion confined to one coordinate: a D₈ decode that rounds a shrunken coordinate to zero fails the parity check and re-rounds exactly that coordinate. **Widening μ does not speed re-acquisition**; it raises the floor as √μ, and at μ = 2⁻³ the floor (0.126) already exceeds the 0.10 target. Even a 45° rotation error in one coordinate pair (initial ‖Δ‖_F = 1.08) is recovered, at the price of a 7.5% symbol error rate during the first 200 symbols.

Two caveats bound these numbers. They are for the 240-point codebook at 12 dB; for a larger codebook such as a 2¹⁶-point Voronoi codebook with vectors of squared norm up to 10, the offset Δx scales with ‖x‖ and the capture range is estimated to shrink by about 2.2×. And the channel is static during re-acquisition.

### 8.5 Representation error: Sherman–Morrison versus QR (`sm`)

**Explicit inverse.** The model stores Ĥ⁻¹ with 16 fraction bits and applies 300 rank-one updates (u with i.i.d. N(0, 0.02²) entries, norm about 0.06; v a random minimal vector), rounding every intermediate product and the final result to 2⁻¹⁶, and measures ‖Ĥ⁻¹H − I‖_F against the exactly updated H:

| κ₂(Ĥ) | After 30 updates | After 300 updates | Growth per update | Updates to reach 0.05 (extrapolated) | Refresh interval N_refresh = 0.025/g |
|---|---|---|---|---|---|
| 2 | 0.0005 | 0.0014 | 3.3×10⁻⁶ | ≈ 15,000 | ≈ 7,600 |
| 10 | 0.0014 | 0.0053 | 1.5×10⁻⁵ | ≈ 3,300 | ≈ 1,700 |
| 30 | 0.0033 | 0.0240 | 7.7×10⁻⁵ | ≈ 650 | ≈ 320 |

The growth rate g is (e₃₀₀ − e₃₀)/270 and the extrapolations assume that linear growth continues. Over this range growth is linear in the update count and between linear and roughly κ₂^1.5 in the condition number, well below the κ₂² worst case of the forward-error bound, because the updates are small and well scaled. The bound is a worst case; the measurement is the typical case.

The refresh policy is set from this table: the explicit inverse is re-derived from a fresh QR factorization every N_refresh = 0.025/g updates, where g is interpolated from the measured κ₂, which keeps the representation error at half the 0.05 budget. For κ₂ ≤ 10 that is a refresh at least every ≈ 1,700 updates, each costing one array pass (59 clocks), an overhead of about 3.5%. For κ₂ > 30 explicit-inverse mode is disabled.

**Orthogonal representation.** Re-factorizing by Givens QR with every intermediate rounded to 20 fraction bits, over 10⁴ rank-one updates at κ₂ = 10 (evaluated every 500 updates and at the end), the model gives

```
worst |R⁻¹ Qᵀ H − I|_F = 4.82e-04
```

which is the CORDIC rounding floor, not a growth: every update is a product of orthogonal rotations in 20-bit arithmetic, and the stored R inherits only the per-operation rounding.

### 8.6 CORDIC rounding (`cordic`)

```
[cordic] 16 iterations, 20-bit: max err=4.22e-05 rms=1.44e-05 K=1.646760
```

Over 5,000 random angles in (−π/2, π/2) and random unit-box input pairs, a 16-iteration circular CORDIC with a 20-bit datapath and a single final scaling by 1/K, K = 1.646760, has a maximum error of 4.22×10⁻⁵ and an RMS error of 1.44×10⁻⁵. The deepest path through the array is 8 boundary cells and 8 internal cells in series along the wavefront, so the worst-case accumulated rotation error is 16 × 4.22×10⁻⁵ ≈ 6.8×10⁻⁴ and the RMS accumulation is √16 × 1.44×10⁻⁵ ≈ 5.8×10⁻⁵. The solve row adds 8 linear-mode operations with the same per-operation bound. The equalized vector z therefore carries a numerical error below 10⁻³ relative (worst case), about 1% of the mismatch budget and two orders of magnitude below the noise at any SNR where the decoder is useful. The whole-pipeline figure from §8.5 (4.82×10⁻⁴ over 10⁴ updates) is consistent with this bound. Sixteen-bit word length matches the practitioner data for QR-based MIMO detectors (0.3 dB loss at 16 bits, 1.3 dB at 14).

### 8.7 Update gate and re-synchronization policy (design policy)

The following rules are specified for the RTL; they are not exercised by the Python tracker, which runs ungated.

**Denominator check.** The Sherman–Morrison denominator 1 + vᵀĤ⁻¹u equals 1 + μ·x̂ᵀĤ⁻¹e/‖x̂‖². Since ‖Ĥ⁻¹e‖ is bounded by condition-scaled noise and μ ≤ 2⁻³, it lies in (0.9, 1.1) under operating conditions; a value outside (0.5, 1.5) indicates a gross decision error or a channel discontinuity and triggers the gate.

**Update gate.** Three tests are applied before an update is committed: (1) the residual norm ‖e‖ is below 3σ_e, with σ_e from a running average; (2) the decoded point is inside the active codebook; (3) the two coset candidates in the snap differed in distance by more than a confidence margin, which excludes decisions near a Voronoi facet. An update failing any test is discarded, which costs the tracker one symbol of lag and nothing else. The design expectation is that at symbol error rates below 10⁻² a small fraction (a few percent) of updates is gated.

**Lock-loss detection.** Lock is declared lost when the running mean of ‖e‖ exceeds 2σ_e for 64 consecutive symbols or when more than 25% of updates are gated over a 256-symbol window. The first response is to continue unchanged, because the capture range in §8.4 is wide.

**Fallback to acquisition.** Only if the gated fraction stays above 25% for 1,024 symbols (expected to require an initial error beyond ‖Δ‖_F ≈ 1 or SNR below 8 dB) does the engine fall back. With pilots available it estimates H by least squares from 64 pilot symbols (one array pass per column of the normal equations, about 600 clocks). Without pilots it replaces the snap with Babai nearest-plane on the current ĤM basis, which tolerates larger mismatch at the cost of exactness, until the gated fraction falls below 10%, then returns to the coset decoder.

---

## 9. Performance model: latency, throughput, energy, area

Every figure below follows from a stated assumption; none is measured on silicon. The model counts operations from the architecture of §6, prices them with published per-operation energies, and takes error rates from §8. The counts are generated by `python3 -I arale_model.py budget`.

### 9.1 Assumptions

| Quantity | Value | Basis |
|---|---|---|
| Clock, cryo-CMOS baseline | 1 GHz | 22 nm FDSOI at 4 K; array critical path is one 20-bit add plus wire |
| CORDIC operation | 16 stages × 3 adds, 20-bit | §6.1 |
| CORDIC operations per vector, QR mode | 268 | 28 vectoring + 140 rotation + 64 y-column + 36 solve-row |
| Snap and tracker operations per vector | 256 | 120 snap + 136 tracker shift-adds |
| Energy per 20-bit add, 45 nm | 0.063 pJ | 0.1 pJ for a 32-bit add, scaled by width |
| Energy per 20-bit add, 7 nm | 0.019 pJ | 0.03 pJ for a 32-bit add, scaled by width |
| Overhead factor (clock tree, registers, control) | 2.5× | typical for deeply pipelined datapaths |
| Cells, full-rate configuration | 44 array + 72 solve = 116 | solve row replicated 9× for 1 vector/clock |
| Cells, reduced configuration | 44 array + 8 solve = 52 | 1 vector per 3 clocks |
| Gate equivalents per CORDIC cell | ≈ 15 kGE | 48 × 20-bit adders (≈ 9.6 kGE) + 960 pipeline flip-flops (≈ 5.8 kGE) |
| Codebook | 240 minimal vectors (7.9 bits) or 2¹⁶-point Voronoi codebook (16 bits) | §8 uses the former |

### 9.2 Latency and throughput

| Configuration | Latency | Throughput | Equivalent |
|---|---|---|---|
| Cryo-CMOS, full-rate | 59 clocks = 59 ns | 1 vector/clock = 1 Gvector/s | 8 Gsample/s; 7.9 Gb/s (240-point) or 16 Gb/s (2¹⁶-point) |
| Cryo-CMOS, reduced | 59 ns | 0.33 Gvector/s | 2.7 Gsample/s |
| Explicit-inverse mode, 16-cell CORDIC bank | 4 + 16 + 4 = 24 clocks = 24 ns | 0.25 Gvector/s | for κ₂ < 10 only |
| RSFQ bit-serial, 30 GHz (design estimate) | ≈ 150 ns | 1 vector per 22 clocks ≈ 1.4 Gvector/s | higher throughput, longer latency than cryo-CMOS at n = 8 |

The 59 ns figure is the number to compare with the 1.1 µs syndrome cycle of a superconducting surface-code experiment and with the 0.32 µs latency of the 90 nm K-best MIMO detector. **Latency and throughput interval are distinct:** no pipelined fixed-point decoder with a 16-iteration CORDIC and an 8-deep array has sub-nanosecond latency, but a throughput interval of 1 ns per vector is achievable.

The feedback loop adds no latency to the data path. Its lag of 65 symbols sets the drift bandwidth the engine can follow: with μ = 2⁻⁴ the tracker recovers from a 10% step in about 100 symbols (§8.3). At 1 Gsymbol/s that corresponds to tracking fractional channel changes of order 10⁵ per second, far above any thermal or aging mechanism; the loop is limited by decision quality, not by speed.

### 9.3 Energy and power

| Item | 45 nm energies | 7 nm energies | 22 nm interpolated |
|---|---|---|---|
| Datapath energy per vector, QR mode | 820 pJ | 246 pJ | ≈ 450 pJ |
| With 2.5× overhead | 2.05 nJ | 0.61 nJ | ≈ 1.1 nJ |
| Power at 1 Gvector/s | 2.05 mW | 0.61 mW | ≈ 1.1 mW |
| Power at 100 Mvector/s (QEC duty) | 0.21 mW | 0.06 mW | ≈ 0.11 mW |
| Share of a 0.7 W 4 K budget at full rate | 0.29% | 0.09% | 0.16% |

The 22 nm column is the geometric mean of the 45 nm and 7 nm values (√(820 × 246) ≈ 450 pJ). These use room-temperature CMOS energies; 4 K operation lowers leakage and wire resistance but raises threshold voltage, and **no 4 K correction is applied.** For scale, the cryo-CMOS predecoder characterized at 4 K draws under 0.56 mW and a 12 nm ASIC surface-code decoder draws 8 mW; ARALE at QEC duty sits roughly an order of magnitude below both.

In RSFQ the switching energy is about 19.3 pJ per vector (268 CORDIC operations × 48 adds × 20 bits × 30 junctions × 2.5 aJ). Static bias current, which ERSFQ removes, would otherwise dominate.

The topological-interconnect upgrade of §7.2 changes none of these energies. Its effect is on the array's wire delay, which at a 1 GHz target is a timing margin, not a power item; it is a path to a 1.5–2 GHz clock in the same node.

### 9.4 Area

| Configuration | Cells | CORDIC-cell logic | Total logic (with snap and tracker) | Estimated area, 22 nm |
|---|---|---|---|---|
| Full-rate | 116 CORDIC cells | 1.74 MGE | ≈ 1.8 MGE | ≈ 0.6 mm² |
| Reduced | 52 CORDIC cells | 0.78 MGE | ≈ 0.85 MGE | ≈ 0.3 mm² |

The CORDIC-cell figures are 116 × 15 kGE and 52 × 15 kGE. The totals add a few percent for the snap block and tracker (an estimate, not a synthesis result), and the area conversion assumes about 3 MGE/mm² at 22 nm. For comparison, a 4×4 CORDIC-Givens MIMO detector in 90 nm used 266 kGE, and FalconSign in 28 nm used 0.71 mm². The reduced configuration is the recommended first silicon.

### 9.5 Decoding error rate

From §8.2 (240-point codebook): SER ≈ 2.8×10⁻² at 8 dB per-dimension SNR, 5.8×10⁻³ at 9 dB, 8×10⁻⁴ at 10 dB, and no errors in 3,000 tracked symbols at 12 dB. Mismatch up to ‖Δ‖_F = 0.10 costs about 0.05 dB or less; 0.20 costs about 0.12–0.16 dB. The tracker holds the mismatch at 0.087 (μ = 2⁻⁴) or 0.040 (μ = 2⁻⁶) at 12 dB, so the tracked decoder performs as the exact decoder to within about 0.05 dB. The E₈ nominal coding gain of 3.0 dB over Z⁸ is preserved.

### 9.6 Scaling to n = 16 and n = 24

The array geometry is unchanged. The cell count grows as n(n+1)/2 + n (the triangle plus the y column): 152 cells at n = 16 and 324 at n = 24. Latency grows by about 2n clocks, to roughly 75 ns and 91 ns. The snap stage changes: the Barnes–Wall lattice Λ₁₆ and the Leech lattice Λ₂₄ have no two-coset decoder, and the best Leech decoders cost 2,700–56,000 operations, so the snap becomes the critical path and needs its own pipelined block. At these dimensions the rank-one QR update becomes decisively cheaper than re-factorization (2(n−1) rotations against n(n−1)/2), which is where the adaptive tier's rank-one structure pays in full.

---

## 10. Comparison with alternative approaches

The right comparators are decoders that take a real vector and return a lattice point under a latency budget. Quantum phase estimation is not one of them: it estimates the eigenphase of a unitary from an eigenstate, its inverse QFT costs Θ((log N)²) gates on log N qubits, and its output is a quantum state that must be measured, collapsing it to one basis value. No construction maps closest-vector decoding onto an eigenphase. Quantum linear solvers (HHL and successors) return a state |x⟩ from which only observables, not the argmin vector, can be read, at cost that scales as κ to κ² in the condition number. Quantum sieving for lattice problems at cryptographic dimensions is estimated to need 10¹³ physical qubits for no net speedup over a single classical core. The quantum row appears in the table to record this, not because those methods compete.

| Approach | What it decodes | Latency | Throughput | Adaptivity | Precision | Evidence |
|---|---|---|---|---|---|---|
| Conway–Sloane in software (CPU) | fixed E₈, exact | ≈ 50–100 ns per vector (120 ops, one core) | ≈ 10–20 Mvector/s per core | re-estimate H and re-invert in O(n³) at ms cadence | float | NestQuant op count |
| FPGA K-best / sphere MIMO detector | QAM symbols on 4×4 channel, near-ML | 0.32 µs | 6.4 Gb/s | QRD per channel-estimate frame, 47–220 cycles | 14–18 bit | OJCAS 2024 |
| Lattice-reduction-aided systolic detector | QAM on 4×4 channel | ≈ 80 cycles per channel update | 249 MHz | LLL per frame | 18/13 fixed | Wang, Biglieri, Yao |
| In-memory one-step inverter (memristor crossbar) | linear solve | < 20 ns settle | continuous | re-program conductances | 6–10 bit | JETCAS 2026 |
| Real-time surface-code decoder (FPGA) | syndrome graph, MWPM / clustering | < 1 µs per round | 1 round per 1.1 µs | n/a | integer | Riverlane |
| Quantum phase estimation / HHL | eigenphase of a unitary; a state \|x⟩ | not a decoder | n/a | n/a | n/a | n/a |
| **ARALE** | E₈ on a drifting 8×8 channel, exact snap after ZF | **59 ns** (model) | **1 Gvector/s** (model) | rank-one update every symbol, ≈ 100-symbol recovery | 16/20 bit | §8–9, Monte Carlo |

Three distinctions matter.

- **Against software,** ARALE's gain is not the snap, which is cheap in software too, but the co-location of the tracker with the equalizer in a fixed-depth pipeline: software re-estimates H on a control-loop cadence of milliseconds, ARALE on a symbol cadence of nanoseconds.
- **Against FPGA MIMO detectors,** ARALE decodes a lattice coset rather than a QAM constellation, and it never recomputes a factorization as a separate step; the detectors' QRD preprocessing is ARALE's steady state.
- **Against in-memory solvers,** ARALE's 16-bit fixed point is the precision a lattice snap needs. Analog inverters at 6–10 bits would place the equalized vector within 2–3% of a decision boundary on a 2¹⁶-point codebook, which the 0.20-mismatch results of §8.2 show is not free.

---

## 11. Target applications

The applications are ordered by how directly the literature supports them. The first has a documented latency budget and a documented lattice structure; the last is out of scope.

### 11.1 Multi-mode GKP decoding in superconducting quantum error correction (primary)

A Gottesman–Kitaev–Preskill code on n bosonic modes is a symplectic lattice in R²ⁿ, and correcting Gaussian shift noise is a closest-point problem in its dual. Two-mode D₄ codes are established theory; four-mode codes on E₈ are the natural next construction, and 2025 work confirms E₈ and the Leech lattice give the strongest low-dimensional GKP distances. The syndrome is an 8-vector of homodyne outcomes. The "channel" is the measurement chain (amplifier gain imbalance, quadrature crosstalk, local-oscillator phase drift), and it drifts on the timescale of minutes to hours. The decoder must return the shift estimate within the 1.1 µs cycle of a surface-code experiment and, for a concatenated surface–GKP code, hand soft information to the outer matching decoder.

ARALE's 59 ns latency leaves about 95% of the cycle for the outer decoder; its tracker absorbs measurement-chain drift without pausing the experiment for recalibration; and its ≈ 0.1 mW at QEC duty fits beside a cryo-CMOS predecoder in the 4 K budget. The output interface is the two squared distances from the snap stage, which provide the likelihood ratio the outer decoder needs. **No four-mode E₈ GKP experiment exists yet; this application waits on one.**

### 11.2 Coherent optical shaping with 8-dimensional Voronoi constellations

Chalmers' 2025 work shows 8-D Voronoi constellations with an E₈ shaping lattice gain up to 1.84 dB over Gray-QAM BICM and 0.99 dB over QAM MLCM, with 1.7 dB experimental OSNR gain over 80 km, targeting 800 Gb/s to 1.25 Tb/s. The receiver's channel is four complex symbols after carrier recovery; residual IQ imbalance, polarization rotation, and skew form an 8×8 real matrix that drifts with temperature and fiber state. This is a room-temperature design (the lattice logic is technology-agnostic; only the substrate tier is cryogenic), at 1 Gvector/s per instance and 16 instances in parallel for a 100 Gbaud link. The coset decoder replaces the 144 squared-distance terms per 2-D symbol that the Chalmers MLCM decoder computes, and the tracker replaces the periodic pilot-based MIMO equalizer update.

### 11.3 Lattice message encoding in ML-KEM-style KEMs

The one place E₈ enters standardized post-quantum cryptography is as a message encoder. Replacing coordinate-wise rounding in Kyber-1024 decryption with a 2E₈ ∩ Z₄⁸ block code lowers the decryption-failure rate from 2⁻¹⁷⁴ to 2⁻²⁸⁶ and permits tighter ciphertext compression; Barnes–Wall-16 and Leech-24 encoders give similar gains with 20–24% ciphertext reduction. Here the "channel" is the LWE error distribution, which is static and known, so ARALE's tracker is idle and only its snap stage is used, at 59 ns against the 11–20 µs of a full ML-KEM operation in 5 nm. The contribution is small (the decode is under 1% of ML-KEM's cycles), but the constant-time property of the coset decoder is a security asset. ARALE does not accelerate the NTT, hashing, or sampling that dominate ML-KEM, ML-DSA, and Falcon, and **the adaptive tier must be disabled in any cryptographic use,** since a decoder whose behavior depends on past ciphertexts is a side channel.

### 11.4 Multi-level cryogenic memory read-out (design study)

E₈ has been proposed for multi-level flash (8 cells per lattice point, 1.6–1.8 dB gain over BCH with Gray-PAM at WER 10⁻⁶), and the same construction applies to any multi-level cell array whose read-out drifts. Superconducting nanowire memories at 1.3 K currently report BER 10⁻⁵ in 4×4 arrays with multi-flux-quantum storage. A controller that reads eight cells as one E₈ point and tracks per-cell gain and offset drift with ARALE's rank-one loop would convert the lattice's 3 dB into margin against retention loss. No cryogenic multi-level memory dense enough to need this exists today.

### 11.5 Out of scope: deep-space optical links

Deep-space communication is sometimes proposed as a target for multidimensional lattice decoders. The deployed system, NASA's Deep Space Optical Communications demonstration on Psyche, is a photon-starved link using 16- and 32-ary pulse-position modulation with serially concatenated PPM codes decoded on an FPGA and received by a superconducting nanowire detector array; it reached 267 Mb/s at 55 million km and 6–8 Mb/s at 386 million km. That is a Poisson timing channel, not a coherent Euclidean one, and its losses are scintillation fades, cloud, and pointing, which are erasures that only aperture, station arraying, interleaving and forward error correction address. A multidimensional lattice decoder has no role in it, and "atmospheric distortion" on such a link is not a linear transform an equalizer can invert. **The application is out of scope.**

---

## 12. What ARALE does not claim

Four claims sometimes made for engines of this kind are not made here, because the model does not support them:

- **Sub-nanosecond latency.** The figure is 59 ns, with a 1 ns throughput interval.
- **O(N²) complexity as a decisive advantage over O(N³).** At n = 8 the factor is about 1.5–2 (§5.5).
- **Error immunity from a topological substrate.** The substrate provides interconnect resistivity and nothing else (§7.1).
- **A replacement for quantum coherence stacks.** No quantum stack performs this task (§10).

Further limits:

- **Generality.** The exact two-coset snap is specific to E₈ and the D_n family. Λ₁₆ and Λ₂₄ need their own snap block (2,700–56,000 operations for Leech).
- **Measured silicon.** Nothing here is measured on silicon. Energies are published room-temperature CMOS figures with no 4 K correction; latency is a cycle count at an assumed 1 GHz clock.
- **Weyl hardware.** No performance number depends on any Weyl-semimetal element.
- **Cryptography.** The tracker must be disabled in any cryptographic use.

What the model does support is narrower and defensible: a 4 K decoder that holds a 3 dB E₈ coding gain on a drifting channel at a symbol rate three to four orders of magnitude faster than software re-calibration, in a power envelope a hundred times below the 4 K budget.

---

## 13. Open problems, risks, and validation roadmap

The framework's logic is specified to the level at which an FPGA prototype can be written; its physical claims are specified to the level at which each can be falsified by one experiment.

### 13.1 Open problems (in the order they would stop the project)

1. **Decision-directed tracking on dense codebooks.** All loop results in §8 use the 240-point minimal-vector codebook. For a 2¹⁶-point Voronoi codebook the offset Δx grows with ‖x‖ and the capture range is estimated to shrink by about 2.2×; the model must be rerun with that codebook and the gate thresholds re-tuned. Pilot-aided acquisition may be needed at the start of every session.
2. **Soft output for concatenated codes.** The snap stage produces two squared distances, giving a hard decision and one likelihood ratio. A surface–GKP outer decoder benefits from the full set of nearby-point distances; extending the snap to the k nearest E₈ points (k = 2 to 8) without a search is open, though the coset structure suggests enumerating the 16 coordinate re-roundings.
3. **Noise correlation.** The model assumes white noise after the front end. Correlated noise (a colored measurement chain) makes zero-forcing suboptimal and would call for a whitening stage before the array, a second triangular factorization of the noise covariance, at the cost of doubling the array.
4. **4 K behavior of 22 nm FDSOI at 1 GHz.** Threshold shifts, kink effects, and interconnect behavior at 4 K are characterized for control logic at hundreds of MHz; 1 GHz deep pipelines at 4 K have fewer public data points. The Phase 2 chip answers this.
5. **Thermal and electrical coupling to qubits.** 0.1–1 mW at the 4 K stage is small, but the engine's switching noise couples to the readout chain that produces y. Isolation and filtering between the engine and the homodyne lines are a packaging problem not addressed here.
6. **Nonlinear Hall rectification at 4 K.** TaIrTe₄'s Berry-curvature dipole changes sign near 175 K; its responsivity, linearity and noise at 4 K for envelope detection are unmeasured. Phase 3 measures them, and a conventional cryo-CMOS front end is the fallback.
7. **Chirality routing.** The deserializer hypothesis requires a field-free or low-field chirality-selective channel with nanosecond settling, which no experiment has shown. It is in Phase 3 with an explicit drop criterion.

### 13.2 Risks

| Risk | Effect if realized | Mitigation | Residual |
|---|---|---|---|
| Decision errors destabilize the tracker at low SNR | loss of lock below 8 dB | gate and pilot fallback (§8.7) | operate above 9 dB |
| CORDIC rounding exceeds budget in the solve row | mismatch floor rises | 2 more guard bits, 22-bit datapath | 10% area |
| 1 GHz not reached at 4 K in 22 nm | latency 90–120 ns | still 10× inside the 1.1 µs budget | none for QEC use |
| NbP interconnect not integrable in the chosen fab | lose the 1.5–2 GHz path | baseline Cu wiring at 1 GHz | none for the model |
| No four-mode E₈ GKP experiment materializes | primary application absent | optical shaping (§11.2) is room temperature and does not need it | the cryogenic case is then unmotivated |
| Chirality routing fails its test | none | cryo-CMOS shift register | none |

### 13.3 Validation roadmap (36 months)

Phases 1 and 2 are sequential and carry the engine. Phase 3 runs beside them and can fail without stopping the project. Phase 4 depends on an external experiment and is the first point at which the cryogenic design is tested against its primary application.

**Phase 1, FPGA prototype (months 0–9).** Implement the reduced configuration (44 array cells, one solve row, snap, tracker) in 16-bit fixed point on a mid-range FPGA at 200–300 MHz. Drive it with a software channel emulator that applies rank-one drift and additive noise. *Gate:* measured latency of 59 clocks; SER curves within Monte Carlo resolution of §8.2; tracker floor and recovery times within 20% of the §8.3 table. *Deliverable:* RTL, testbench, and the measured curves.

**Phase 2, cryo-CMOS test chip (months 6–24).** Port the RTL to 22 nm FDSOI, reduced configuration, 0.3 mm² target. Characterize at 300 K and at 4 K in a dilution refrigerator's 4 K stage. *Gate:* 1 GHz operation at 4 K, power below 1 mW at full rate, and lock held against a live drifting analog front end for 10⁴ seconds without a pilot. *Deliverable:* the chip, its 4 K characterization, and the measured power.

**Phase 3, materials (months 0–18, in parallel).** Three experiments with numeric gates.

- (a) NbP sputtered test structures at the array's wire pitch, measuring RC delay at 4 K against Cu controls. *Gate:* RC reduction of at least 2× below 5 nm.
- (b) TaIrTe₄ nonlinear-Hall rectifier at 4 K under a 1 GHz drive. *Gate:* responsivity above 0.1 A/W and NEP below 10 pW/√Hz at 4 K.
- (c) The chirality-routing test of §7.4: four-terminal PtSe₂ or TaRhTe₄ device, 100 MHz pulse train, 1 T, 4 K. *Gate:* channel contrast above 10 and settling under 10 ns.

A failed gate removes that element from the stack and nothing else changes.

**Phase 4, GKP integration (months 24–36).** Attach the Phase 2 chip to a four-mode bosonic experiment implementing an E₈ GKP code, if one exists by then, or to a two-mode D₄ experiment with the array configured for n = 4. *Gate:* shift estimates delivered to the outer decoder inside the 1.1 µs cycle, with a measured logical-error improvement over the experiment's existing software decoder. *Deliverable:* the integrated system and its logical error rate.

---

## 14. Repository and reproduction

```
README.md          this file
arale_model.py     the numerical model (numpy only)
```

`arale_model.py` contains the E₈ coset decoder (`dec_D8`, `dec_E8`), the 240 minimal vectors (`min_vectors`) and generator (`E8_GEN`), the SER experiments (`ser`, `run_ser`, `run_ser_precise`), the channel generator and decision-directed NLMS tracker with decision delay (`make_channel`, `track`, `run_track`, `run_track_snr`), the Sherman–Morrison and Givens-QR fixed-point experiments (`run_sm`), the CORDIC model (`cordic_rotate`, `run_cordic`), and the operation, latency, energy and area budget (`run_budget`).

```
python3 arale_model.py all            # every experiment (several minutes)
python3 arale_model.py decoder        # brute-force check of the E8 decoder          (§8.1)
python3 arale_model.py ser            # SER vs SNR, with and without mismatch, 20k trials/cell   (§8.2)
python3 arale_model.py ser-precise    # the 8 and 9 dB rows at 100k trials with standard errors  (§8.2)
python3 arale_model.py track          # tracker floor, step recovery, capture range  (§8.3, §8.4)
python3 arale_model.py track-snr      # tracker floor versus SNR at mu = 1/16        (§8.3)
python3 arale_model.py sm             # Sherman-Morrison vs QR error growth          (§8.5)
python3 arale_model.py cordic         # CORDIC rotation error                        (§8.6)
python3 arale_model.py budget         # operation, latency, energy, area counts      (§9)
```

Requirements: Python 3 and numpy. All random streams are seeded (`default_rng(1)` globally, and fixed seeds inside the tracker, Sherman–Morrison and QR experiments), so each command is reproducible when run on its own. Because the SER experiments share one global generator, running `all` consumes the stream in a different order from running a step alone and gives different but statistically equivalent single-cell SER values; the per-cell standard errors in §8.2 bound this variation.

Expected `budget` output:

```
[budget] CORDIC ops/vector = 268; other ops/vector = 256
[budget] 45 nm: datapath 820 pJ/vector; x2.5 overhead 2.05 nJ; 2.05 mW at 1 Gvector/s
[budget] 7 nm: datapath 246 pJ/vector; x2.5 overhead 0.61 nJ; 0.61 mW at 1 Gvector/s
[budget] 22 nm interpolated: ~450 pJ datapath, ~1.1 nJ with overhead
[budget] RSFQ switching only: 19.3 pJ/vector
[budget] latency: 4 + 9 + 16 + 8 + 18 + 4 = 59 clocks
[budget] full-rate: 116 CORDIC cells x 15 kGE = 1.74 MGE
[budget] reduced: 52 CORDIC cells x 15 kGE = 0.78 MGE
```

---

## 15. References

Pages opened for this work, grouped by topic. Numbers in the text are taken from these pages; where a figure was derived (a latency from a cycle count, an energy from TOPS/W) the text says so.

**E₈ lattice and decoding**
- Viazovska, M., *The sphere packing problem in dimension 8*, Annals of Mathematics 185 (2017); arXiv 1603.04246.
- Conway, J. H. and Sloane, N. J. A., *Fast quantizing and decoding algorithms for lattice quantizers and codes*, IEEE Trans. Inf. Theory 28 (1982); procedure as described in US patent 4,507,648; bibliography at neilsloane.com.
- E₈ structure, cosets, generator, minimal vectors: arXiv 1009.5764; Cioffi, Appendix B; E8 lattice, Wikipedia.
- NestQuant (operation count, normalized second moment): arXiv 2502.09720.
- QuIP# E8P codebook: arXiv 2402.04396; QTIP: arXiv 2406.11235; HyperQuant: arXiv 2606.23406.
- Modulo-ADC with E₈ folding, Voronoi-relevant vectors: arXiv 2605.24974.
- Babai nearest plane, lecture notes: Regev, NYU.
- E₆*/E₇* closest-point decoders: arXiv 2607.10885.

**Coded modulation and optics**
- Chalmers, Voronoi constellations with E₈/Λ₁₆/Λ₂₄ shaping, IEEE Trans. Commun. 73 (2025); ECOC 2021 experiment.

**Lattice-based PQC**
- E₈ message code for Kyber decryption failure: arXiv 2601.08452; BW₁₆/Leech encoders for Kyber: arXiv 2308.13981; Kyber with lattice quantizer: arXiv 2401.15534; 2-D codes for Kyber: ePrint 2024/1243.
- ML-KEM/ML-DSA profiling on OpenTitan: ePrint 2024/1192; Adams Bridge 5 nm accelerator: ePrint 2026/256; FalconSign: TCHES 2025; Falcon specification: falcon-sign.info; NIST PQC status: Moody 2025.
- Quantum algorithms for lattices: Regev factoring arXiv 2308.06572; Chen's LWE claim and bug note ePrint 2024/555; quantum sieving resource estimate arXiv 2410.13759.

**Quantum error correction and bosonic codes**
- GKP codes, lattice perspective: Conrad, Eisert, Arzani, Quantum 2022; arXiv 2109.14645.
- Closest-point decoding of multimode GKP: Lin, Chamberland, Noh 2023; D₄ grid codes: Royer, Singh, Girvin; GKP from SIS lattices: arXiv 2509.10183.
- GKP experiments: Sivak et al. 2023; GKP qudits 2024; ETH trapped-ion GKP; Error Correction Zoo.
- Real-time decoders: Google, arXiv 2408.13687; Riverlane Collision Clustering and Local Clustering; cryo-CMOS predecoder Pinball, arXiv 2512.09807.
- Quantum phase estimation and linear solvers: QPE and QFT (Wikipedia); QLS survey, arXiv 2411.02522.

**Fluxonium and cryogenic memory**
- arXiv 2407.15783 (CNOT 99.94%); arXiv 2501.16691 (tantalum fluxonium readout); arXiv 2507.14436; arXiv 2411.13437.
- Nanowire memory: arXiv 2503.22897; Buzzi et al., MIT; cold DRAM argument: Semiconductor Engineering.

**Cryogenic electronics**
- Dilution-refrigerator budgets: arXiv 2608.00990; Bluefors specifications and KIDE.
- Cryo-CMOS: IBM 14 nm controller; survey, arXiv 2511.13965.
- SFQ: Holmes et al., arXiv 1602.03546; Yoshikawa, EUCAS 2025.

**Weyl semimetals**
- Armitage, Mele, Vishwanath, Rev. Mod. Phys. 2018; TaAs discovery: Phys. Rev. X 5, 031013.
- Fermi velocities and mobilities: TaAs arXiv 1603.08846; NbP arXiv 1502.04361; graphene arXiv 1208.0567.
- Chiral anomaly and current jetting: dos Reis et al.; TaRhTe₄ arXiv 2502.18937; nonlocal valley transport proposal, Parameswaran et al.; disorder and localization arXiv 2402.14063.
- Interconnects: NbP films, Science 2024, arXiv 2409.17337; NbAs nanowires, Science 2026, arXiv 2503.04621; NbAs slab transport arXiv 2211.10426; CoSi (IBM).
- Chirality devices: PtSe₂ arXiv 2103.00279; Cd₃As₂ chiral pumping (ECNU); NbP/NbN transistor, Nat. Commun. 2024.
- Nonlinear Hall: TaIrTe₄ rectification, Nature Nanotechnology 2021; THz detector, Nat. Electron. 2025; NbIrTe₄ rectenna (SITP) 2026.
- Films: TaAs MBE (NREL), arXiv 2303.05469; Co₃Sn₂S₂ sputtering arXiv 2106.01843; Mn₃Sn arXiv 2207.12885.

**CORDIC, systolic arrays, rank-one updates**
- Meher et al., *50 years of CORDIC*; guard bits (MathWorks); FPGA error analysis arXiv 2308.01025; scaling-free CORDIC, TCAD 2026; radix-4 CORDIC, JSPS 2023.
- Gentleman & Kung 1981; McWhirter extension as described in Hsieh, NASA CR-191396; IQRD-RLS FPGA, Radioengineering 2017; CORDIC systolic QRD fixed point, Muñoz & Hormigo; CORDIC Givens MIMO detector, IEEE Access 2022; MVDR systolic QRD arXiv 2609.03137.
- Gill, Golub, Murray, Saunders, *Methods for modifying matrix factorizations*, Math. Comp. 1974.
- Sherman–Morrison stability: Hashemi & Nakatsukasa 2025 and 2026 follow-up (arXiv 2609.12266); Woodbury identity; RLS covariance update (SIFt-RLS).
- In-memory inversion: memristor MIMO precoding, JETCAS 2026; memristive linear algebra; in-memory eigenvector solver.
- MIMO detectors and fixed point: sorter-free K-best, OJCAS 2024; systolic lattice-reduction-aided detection; LLL VLSI (ETH); fixed-point lattice-reduction-aided soft MIMO (Rahman & Choi); LDLC decoder.
- Energy per operation: Horowitz 2014 as tabulated in SFU lecture and PokeBNN; IBM PCM chip arXiv 2212.02872.

**Deep-space optical communication**
- DSOC results: IEEE Photonics; The Register; JPL arraying poster; SNSPD receiver arXiv 2409.02356; SCPPM modem, NASA NTRS.

**Numerical model**
- The Monte Carlo and operation-count model used in §8 and §9 is the self-contained Python script `arale_model.py` (E₈ coset decoder verified against brute force; decision-directed NLMS tracker with decision delay; fixed-point Sherman–Morrison and Givens QR; CORDIC rotation error; energy, latency and area counts).
