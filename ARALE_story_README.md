# ARALE

### A decoder that keeps listening while the world drifts

---

## The shape that solved itself

There is a shape in eight dimensions called E₈. You cannot picture it, and that's fine. Mathematicians proved in 2016 that it packs spheres more tightly than any other arrangement in eight dimensions. Each sphere touches 240 neighbors. Nothing else in that space comes close.

Here is the odd part. You'd expect a shape this perfect to be a nightmare to work with. Finding the nearest point of a dense lattice is usually a brutally hard search problem. But E₈ is built from two simpler lattices stacked together, so the nearest-point question takes about 120 simple integer operations. No search is needed. A small piece of hardware can do it in a few nanoseconds, and it is always exactly right.

So the shape isn't the hard part. The hard part is that the world never hands you the shape cleanly.

## The problem nobody sees until it breaks

Imagine eight wires carrying a signal. Each wire amplifies a little differently. They leak into each other. A temperature change shifts them. A week later, it's different again.

By the time the signal reaches you, the perfect shape has been stretched and tilted. The clean answer is hiding inside a distorted picture, and the distortion keeps changing.

The usual fix is slow: stop, measure the distortion, do some heavy math, then resume. That cycle takes milliseconds. In a quantum computer, the clock ticks every 1.1 microseconds. Milliseconds are an eternity. The tool arrives long after the moment it was needed.

ARALE is built on a different idea. **Don't stop to recalibrate. Learn the distortion while you decode, one symbol at a time.**

## The trick: use your own answers as the teacher

Here's the loop.

1. A signal arrives.
2. ARALE undoes its best guess of the distortion.
3. It snaps the result to the nearest E₈ point. That's its answer.
4. It asks: *if the distortion really were what I think, would this answer have produced the signal I just saw?*
5. The gap between "what I expected" and "what I saw" is a small error. ARALE nudges its estimate to shrink it.

The nudge is tiny and simple, a single correction of the kind mathematicians call rank-one. It uses only shifts and additions. There are no multipliers anywhere in the data path.

This sounds circular. It's learning from its own answers. It works because the answers are almost always right, and a right answer is a free training example. The system doesn't need a teacher while it's mostly correct.

## What happens when it's wrong

We tested this the hard way: by breaking it on purpose in simulation.

- **A 10% jump in one line's gain.** ARALE recovered within about 100 symbols, with zero decision errors along the way.
- **A line off by 120%.** It recovered in 676 symbols with no errors in the first 200.
- **A scrambled estimate** (a random error of 0.5, then 0.8, on the whole 8×8 matrix). It found its way back in 382 and 439 symbols. At 0.8 it made roughly 1 mistake in 200 along the way.
- **A 45° twist in one pair of coordinates.** This one hurt. It recovered in 780 symbols, but it got about 7.5% of its first 200 answers wrong.

There's a reason single-line errors are easy. E₈ has a built-in parity check. When one coordinate is off, the check fails and points at exactly the guilty coordinate. The shape helps repair itself.

One surprise: turning up the learning rate does **not** make recovery faster. It just makes the system jumpier. The "floor", the small amount of error that remains when everything is steady, rises with the square root of the learning rate. In our tests the floor went from 0.040 to 0.126 as we sped up the learning. Calm and steady wins.

## How accurate is "good enough"?

E₈ gives a theoretical advantage of 3 dB over plain grid-based coding, which is the same as doubling the effective signal strength. The question is how much of that survives an imperfect channel estimate.

We ran 100,000 trials per setting. With a residual error of 0.10 in the estimate, the penalty was about 0.05 dB or less. At 0.20 it was about 0.12 to 0.16 dB. The tracker settles at around 0.09, so the 3 dB advantage is essentially intact.

## Why the engine is built from rotations

There are two ways to keep track of the channel. One stores its inverse and patches it as things change. That works, but the errors compound. We measured it: for a nicely behaved channel the error grows slowly, but for a poorly conditioned one it grew about 23 times faster per update (3.3×10⁻⁶ versus 7.7×10⁻⁵). You have to rebuild the inverse regularly or it drifts.

The other way stores the channel as a stack of **rotations**. Rotations are special: they never stretch anything, so they never amplify mistakes. After 10,000 updates, the error was still 0.00048, the same as the floor set by the arithmetic itself. It did not grow.

So ARALE uses rotations. The hardware is a triangle of 44 small units, each doing a rotation with shifts and adds only (a method called CORDIC). The same units that apply the correction also update it. One kind of part does two jobs.

## The numbers

| | |
|---|---|
| **Latency** | 59 clocks, or **59 ns** at 1 GHz |
| **Throughput** | one 8-number vector per clock (full build) |
| **Energy per vector** | about 0.45 nJ for the arithmetic, about 1.1 nJ all-in |
| **Power** | about 1.1 mW flat out, about 0.1 mW at quantum-computer duty |
| **Size** | about 0.3 mm² (reduced build), 0.6 mm² (full) |
| **Memory** | about 12 thousand bits. No special memory technology needed |

All of these are model estimates, not silicon measurements. Energy uses room-temperature figures with no correction for 4 K.

To put the speed in context: 59 nanoseconds leaves roughly 95% of a 1.1 µs quantum error-correction cycle for everything else. And a dilution refrigerator can spare around 0.7 to 1.5 watts at its 4 K stage. ARALE would use roughly 0.15% of that.

## Where this could matter

**Quantum error correction.** Some of the most promising ways to protect quantum information store it across several oscillator modes. Correcting errors in those codes is mathematically a nearest-lattice-point problem. A four-mode version would use E₈. The measurement chain drifts, and the answer is needed in under a microsecond. This is the primary target. No four-mode experiment exists yet, so the application is waiting on one.

**Fiber-optic links.** Shaping signals with E₈ has shown gains up to 1.84 dB in lab work. A receiver that tracks drift continuously could replace periodic recalibration. This version runs at room temperature.

**Cryptography.** E₈ can be used to encode messages in lattice-based encryption, cutting the chance of a decryption failure dramatically. The gain is real but small, and the adaptive part must be switched off there, because a decoder that changes with past inputs can leak information.

**Cryogenic memory.** Reading eight storage cells as a single E₈ point could turn the 3 dB advantage into reliability. No memory dense enough to need it exists today. That's a design study.

**Not a fit: deep-space laser links.** Those use pulse timing and photon counting. The losses are fades and dropouts, which a lattice decoder can't fix.

## The part we don't oversell

A good idea is made stronger by saying clearly what it isn't.

- **It isn't sub-nanosecond.** A 16-step rotation pipeline cannot be faster than the clocks it takes. The latency is 59 ns. What is one nanosecond is the *interval* between results.
- **It isn't a big algorithmic leap at this size.** The mathematical advantage of rank-one updates over full recalculation is about 1.5 to 2 times at eight dimensions. It becomes decisive at 16 or 24.
- **Exotic materials don't make it error-proof.** There is a family of materials, Weyl semimetals, whose electronic properties are protected by topology. What's protected is the existence of certain electronic states, not the accuracy of a calculation. We use them in one place that's proven: ultra-thin wires whose resistance falls as they shrink. That could shorten wire delay by about 3×. Nothing in the performance numbers depends on it.
- **One idea is pure speculation.** Routing signals by electron "handedness" instead of a shift register has been seen only with multi-tesla magnets, which would destroy a quantum computer sitting beside it. It's listed with a test and a way to drop it. If the test fails, nothing else changes.
- **It's not a quantum algorithm.** Quantum phase estimation finds eigenphases, and quantum linear solvers return a quantum state. Neither returns the nearest lattice point.
- **It's only for E₈ and its relatives.** Larger lattices like the 24-dimensional Leech lattice need a far heavier final step, thousands of operations rather than 120.

## How we'd find out if it works

1. **Months 0–9: FPGA prototype.** Build the reduced design, feed it a simulated drifting channel. Pass if it hits 59 clocks and reproduces the error curves.
2. **Months 6–24: cryo-CMOS test chip.** 22 nm, 0.3 mm². Pass if it runs at 1 GHz at 4 K, under 1 mW, and holds lock against a live drifting front end for 10,000 seconds without help.
3. **Months 0–18, in parallel: materials.** Three measurements, each with a number to hit or be dropped.
4. **Months 24–36: connect to a real bosonic experiment**, if one exists by then, or a two-mode version if not.

The biggest unanswered questions: does tracking still work with the much denser 65,536-point codebooks (the capture range is estimated to shrink about 2.2×)? How does a 1 GHz pipeline behave at 4 K? And can the engine's switching noise be kept out of the delicate quantum readout lines?

## Try it yourself

The whole model is one Python file and needs only numpy.

```
python3 arale_model.py all          # everything (a few minutes)
python3 arale_model.py decoder      # check the E₈ decoder against brute force
python3 arale_model.py ser          # error rates vs. noise and mismatch
python3 arale_model.py ser-precise  # 100,000-trial version with error bars
python3 arale_model.py track        # tracking, recovery, capture range
python3 arale_model.py track-snr    # tracking floor vs. noise
python3 arale_model.py sm           # inverse vs. rotation error growth
python3 arale_model.py cordic       # rotation arithmetic accuracy
python3 arale_model.py budget       # latency, energy, area
```

The tracking simulations use floating point and no update gating, so they show the behavior of the basic loop. Word-length effects are measured separately by the `sm` and `cordic` steps.

## The big picture

Most engineering treats a changing world as something to pause for. ARALE treats it as something to follow.

The pieces are old: a famous lattice from 1982, a rotation method from 1959, a learning rule from the 1960s. What's new is putting them in one pipeline small enough to sit in the coldest part of a quantum computer, quick enough to keep pace with it, and honest about where each claim stands.

Whether the first four-mode E₈ quantum code arrives in two years or ten, a decoder like this will be wanted when it does.
