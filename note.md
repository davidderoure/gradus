# A Two-Stage Perceptual Account of Interval Consonance

**David De Roure**  
*Working note, June 2026*

---

## Abstract

Euler's *Gradus Suavitatis* (1739) assigns a dissonance value to a musical interval p/q by the formula G(p/q) = 1 + Ω\*(p) + Ω\*(q), where Ω\*(n) = Σ eᵢ(pᵢ − 1) sums the weighted prime exponents of n. We show that the simpler asymmetric formula f(p/q) = p + Ω\*(q) achieves a higher rank correlation with human consonance ratings (ρ = 0.989 vs 0.979), and that this formula admits a natural two-stage perceptual interpretation grounded in harmonic series theory. We also show that Gradus is exactly equivalent to a weighted harmonic coincidence count with weights equal to Ω\*(n) itself, and that the simpler metric max(p,q) — equivalent to Galileo's pulse-coincidence model — achieves the same correlation as f(p/q).

---

## 1. Background

Euler's Gradus Suavitatis, catalogued in OEIS A275314, assigns to each musical interval a positive integer measuring its arithmetic dissonance. For a ratio p/q in lowest terms, the characteristic number is n = p·q, and:

> G(p/q) = 1 + Ω\*(p·q)

where Ω\* is the completely additive function Ω\*(n) = Σ eᵢ(pᵢ − 1). Since Ω\* is completely additive, G(p/q) = 1 + Ω\*(p) + Ω\*(q), treating numerator and denominator symmetrically. The weights (p − 1) for prime p encode Euler's aesthetic judgment that large primes represent greater "difficulty."

Euler's formula correctly ranks the 13 standard intervals of Western music against human consonance ratings (Krumhansl 1990) with Spearman ρ = 0.979, and is exactly tied only on the Major third / Major sixth pair and on a three-way tie at Gradus 8.

---

## 2. Gradus as Harmonic Coincidence

For two notes at frequency ratio p:q (coprime), exact harmonic coincidences occur at harmonics p, 2p, 3p, … of the lower note. For a fixed count M of coincidences, the weighted sum with weights w(n) = Ω\*(n) yields:

> Σₖ₌₁ᴹ Ω\*(k·p) = Ω\*(p)·M + C(M)

where C(M) = Σₖ₌₁ᴹ Ω\*(k) is a constant independent of p. This follows from the complete additivity of Ω\*: Ω\*(k·p) = Ω\*(k) + Ω\*(p). The symmetric score over both notes gives:

> Score = M · (G(p/q) − 1) + 2·C(M)

**Gradus is therefore exactly a weighted harmonic coincidence count**, with harmonics weighted by their own prime complexity. This confirms the physical intuition behind Euler's formula but reveals its self-referential structure: recovering Gradus from coincidences requires weights that are themselves Gradus-like.

Unweighted coincidence counting (flat weights) achieves ρ = 0.791 across all coprime ratios with p,q ≤ 16. The weight function log(n) achieves ρ = 0.795 without any prime structure, establishing a ceiling for physics-based models that do not invoke prime factorisation.

---

## 3. Galileo's Model and max(p,q)

Galileo (1638) described consonance in terms of pulse coincidences: two strings at ratio p:q produce coinciding pulses every p vibrations of the lower string. The "sweetness" of the interval grows with the frequency of coincidences, i.e. inversely with p. The metric max(p,q) = p is therefore:

1. **Galileo's pulse model** — vibrations between coincidences
2. **Harmonic series position** — the upper note is partial p of the implied fundamental (Rameau 1722)
3. **Farey index** — the smallest n such that q/p appears in the Farey sequence Fₙ (Cauchy 1816)
4. **First coincidence harmonic** — the position in the lower note's series where the two notes first meet

These four descriptions, discovered independently across four centuries, are mathematically identical. max(p,q) achieves ρ = 0.989 against human ratings, marginally exceeding Gradus, but cannot distinguish intervals with the same numerator (e.g. 5/4 and 5/3 both have max = 5).

---

## 4. An Improved Formula

Euler's formula is symmetric: Ω\*(p) + Ω\*(q). Human ratings, however, distinguish numerator and denominator — the bass note and the upper note play different perceptual roles. We find empirically that the asymmetric formula:

> **f(p/q) = p + Ω\*(q)**

achieves ρ = 0.989, matching max(p,q) and exceeding Gradus. The formula differs from Gradus only in replacing Ω\*(p) with p for the numerator. Optimising prime weights freely with six parameters (separate weights for primes 2, 3, 5 in numerator and denominator) achieves perfect rank correlation ρ = 1.000, confirming that the structure required to perfectly predict the human rankings is present in an asymmetric prime-weighted formula.

The formula f(p/q) = p + Ω\*(q) can be derived from first principles as the sum of two quantities:

| Component | Formula | Meaning |
|-----------|---------|---------|
| Upper note reach | p | Partial number of upper note in implied harmonic series |
| Bass context depth | Ω\*(q) | Prime complexity of lower note's position in series |

---

## 5. A Two-Stage Perceptual Model

The formula motivates a two-stage model of how the auditory system evaluates an interval:

**Stage 1 — Bass establishes harmonic context (cost Ω\*(q)).**  
The lower note is the q-th partial of an implied fundamental F. The brain extrapolates downward to F by "dividing out" the prime factors of q. Each factor of 2 (an octave) costs 1 unit; each factor of 3 costs 2; each factor of 5 costs 4. This is exactly Ω\*(q). Intervals whose bass note is a power of 2 (e.g. the octave, q=1; the major second bass, q=8=2³) establish the fundamental cheaply via octave equivalence. Intervals with an odd prime in q (the fourth, q=3; the minor third bass, q=5) require the brain to extrapolate through a genuinely new prime, incurring higher cost.

**Stage 2 — Upper note reaches into the series (cost p).**  
Once the fundamental is established, the upper note must be recognised as its p-th partial. Higher partials are quieter (amplitude ∝ 1/p in typical harmonic sounds), more narrowly tuned, more sensitive to mistuning, and have lower pitch salience (Terhardt 1979). The cost of confirming the upper note as a member of the harmonic series is therefore p itself — the direct measure of its height in the series.

**Total perceptual effort: f(p/q) = p + Ω\*(q).**

This account connects:
- **Galileo** (1638): pulse coincidences as physical mechanism
- **Rameau** (1722): the *corps sonore* (harmonic series) as the source of musical meaning
- **Euler** (1739): prime arithmetic as the language of harmonic complexity
- **Terhardt** (1979): virtual pitch and harmonic template matching as auditory mechanisms
- **Huron** (2001): statistical learning of harmonic patterns as the basis of consonance judgments

---

## 6. The Remaining Failure

No simple formula distinguishes the Major third (5/4, human rank 5) from the Major sixth (5/3, human rank 7). Both have p = 5 (upper note is P5) and Ω\*(q) = 2 (denominator is either 4 = 2² or 3, both with Ω\* = 2). The distinction requires knowing that q = 4 = 2² (a pure octave stack, "free" under octave equivalence) is simpler than q = 3 (the first genuinely new prime), yet Ω\*(4) = Ω\*(3) = 2 so Euler's weighting cannot see this. Resolving this tie requires either a richer formula or an appeal to learned cultural context (thirds occur far more frequently than sixths in tonal music).

---

## 7. Summary

| Formula | ρ vs human | Notes |
|---------|-----------|-------|
| Euler's Gradus: 1 + Ω\*(p) + Ω\*(q) | 0.979 | Symmetric; two ties unresolved |
| Galileo / max(p,q) | 0.989 | Simplest; same tie |
| **f(p/q) = p + Ω\*(q)** | **0.989** | Asymmetric; two-stage interpretation |
| Optimised asymmetric weights | 1.000 | Six parameters; no closed form |
| Roughness (Plomp-Levelt) | 0.709–0.863 | Register-dependent; fails on Tritone |

The formula f(p/q) = p + Ω\*(q) is marginally simpler than Gradus, slightly more accurate, and directly interpretable as a two-stage process of harmonic extrapolation. It suggests that Euler's prime-weighting is correct for the *denominator* but that the *numerator* is better described by the raw partial number — consistent with the auditory system's known sensitivity to partial height rather than prime structure in the upper voice.

---

*Scripts and data: `/Users/davidderoure/gradus/`*
