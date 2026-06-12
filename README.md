# Gradus

**A two-stage perceptual account of interval consonance**

*David De Roure, June 2026*

---

## Overview

Euler's *Gradus Suavitatis* (1739) measures the dissonance of a musical interval p/q by the formula:

> G(p/q) = 1 + Ω\*(p) + Ω\*(q)

where Ω\*(n) = Σ e_i(p_i − 1) is a weighted sum over the prime factorisation of n. This function is catalogued as [OEIS A275314](https://oeis.org/A275314).

This project investigates whether a simpler asymmetric formula:

> **f(p/q) = p + Ω\*(q)**

gives a better account of human consonance perception — and shows that it does, while admitting a natural two-stage perceptual interpretation.

---

## Key findings

- **Gradus is exactly a weighted harmonic coincidence count.** Using weights w(n) = Ω\*(n), the weighted sum of coincidences between two notes at ratio p:q equals M·(G(p/q) − 1) + constant, for any fixed count M. This follows from the complete additivity of Ω\*.

- **max(p,q) achieves Spearman ρ = 0.989** against human consonance ratings (Krumhansl 1990), beating Euler's Gradus (ρ = 0.979), and has four equivalent musical interpretations: Galileo's pulse-coincidence model (1638), harmonic series position (Rameau 1722), Farey index (Cauchy 1816), and first coincidence harmonic.

- **f(p/q) = p + Ω\*(q) also achieves ρ = 0.989**, matching max(p,q) while distinguishing more intervals. It differs from Gradus only in replacing Ω\*(p) with the raw partial number p.

- **Perfect ρ = 1.000** is achievable with an asymmetric six-parameter prime-weighted formula, confirming that the structure needed to perfectly predict human rankings is present in an asymmetric prime-weighted approach.

- **One stubborn tie remains:** the major third (5/4) and major sixth (5/3) both give f = 7, because Ω\*(4) = Ω\*(3) = 2. Resolving this requires either a richer formula or cultural context.

---

## Two-stage perceptual model

The formula f(p/q) = p + Ω\*(q) decomposes as two perceptual costs:

| Stage | Cost | Interpretation |
|-------|------|----------------|
| 1 — Bass establishes context | Ω\*(q) | Lower note is the q-th partial of an implied fundamental; prime complexity of extrapolating down to it |
| 2 — Upper note reaches into series | p | Upper note must be recognised as the p-th partial; higher partials are quieter and harder to confirm |

This connects Galileo (1638), Rameau (1722), Euler (1739), Terhardt (1979), and Huron (2001) within a single arithmetic framework.

---

## Files

| File | Description |
|------|-------------|
| `note.md` / `note.pdf` | Working note with full derivations and results |
| `harmonic_overlap.py` | Initial exploration of harmonic coincidence vs Gradus |
| `gradus_analysis.py` | Spearman correlations; Gradus ties and failures |
| `weighted_coincidence.py` / `weighted_coincidence2.py` | Proof that Gradus = weighted coincidence count |
| `perceptual_model.py` | Plomp–Levelt roughness model comparison |
| `gradus_failures.py` | Where Gradus fails; max(p,q) as tiebreaker |
| `maxpq_interpretation.py` | Four interpretations of max(p,q) |
| `formula_search.py` | Systematic search for improved formula; optimised weights |
| `table.py` / `table_plot.py` | Comparison table: Gradus, f(p/q), human ratings |
| `perceptual_account.py` / `perceptual_stem.py` | Two-stage model diagrams |
| `oeis_draft.md` | OEIS submission draft for triangle T(n,k) = n + Ω\*(k) |
| `oeis_generate.py` | Generate and verify triangle terms and b-file |
| `oeis_superparticular.py` | OEIS submission draft for superparticular subsequence a(n) = n + Ω\*(n−1) |
| `bfile_draft.txt` | OEIS b-file: triangle T(n,k), 774 terms |
| `bfile_superparticular.txt` | OEIS b-file: superparticular sequence, 300 terms |

---

## OEIS

Two entries are being prepared for the [OEIS](https://oeis.org):

- **Main triangle** T(n,k) = n + Ω\*(k), n ≥ 1, 1 ≤ k ≤ n, gcd(n,k) = 1  
  First 40 terms: `1, 2, 3, 4, 4, 6, 5, 6, 7, 7, 6, 10, 7, 8, 9, 9, 11, 10, 8, 10, 12, 14, 9, 10, 11, 13, 15, 12, 10, 12, 16, 14, 11, 12, 13, 13, 15, 14, 17, 14`

- **Superparticular subsequence** a(n) = n + Ω\*(n−1), n ≥ 2  
  First 40 terms: `2, 4, 6, 7, 10, 10, 14, 12, 14, 16, 22, 17, 26, 22, 22, 21, 34, 24, 38, 27, 30, 34, 46, 30, 34, 40, 34, 37, 58, 38`

The author is the creator of [A275314](https://oeis.org/A275314) (Euler's Gradus Suavitatis).

---

## References

- L. Euler, *Tentamen novae theoriae musicae*, St. Petersburg, 1739.
- G. Galilei, *Discorsi e Dimostrazioni Matematiche*, 1638.
- J.-P. Rameau, *Traité de l'Harmonie*, Paris, 1722.
- C. L. Krumhansl, *Cognitive Foundations of Musical Pitch*, Oxford University Press, 1990.
- E. Terhardt, Calculating virtual pitch, *Hearing Research* 1 (1979), 155–182.
- R. Plomp & W. J. M. Levelt, Tonal consonance and critical bandwidth, *JASA* 38 (1965), 548–560.
