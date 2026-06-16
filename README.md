# Gradus

**An Asymmetric Formula for Interval Consonance and its Relation to Harmonic Coincidence**

*David De Roure, June 2026*

Full derivations and results: [arXiv:2606.16412](https://arxiv.org/abs/2606.16412)

---

## Overview

Euler's *Gradus Suavitatis* (1739) measures the dissonance of a musical interval p/q by the formula:

> G(p/q) = 1 + Ω\*(p) + Ω\*(q)

where Ω\*(n) = Σ eᵢ(pᵢ − 1) is a weighted sum over the prime factorisation of n. This function is catalogued as [OEIS A275314](https://oeis.org/A275314).

This project investigates whether a simpler asymmetric formula:

> **f(p/q) = p + Ω\*(q)**

gives a better account of human consonance perception — and shows that it does, while admitting a natural two-stage perceptual interpretation.

---

## Key findings

- **Gradus is equivalent, under a discrete harmonic model, to a weighted harmonic coincidence count.** Using weights w(n) = Ω\*(n) and harmonics indexed by positive integers up to a fixed truncation level M, the weighted sum of coincidences between two notes at ratio p:q equals M·(G(p/q) − 1) + constant. This follows from the complete additivity of Ω\* and connects Gradus to Galileo's pulse-coincidence model (1638).

- **max(p,q) achieves Spearman ρ = 0.989** against human consonance ratings (Krumhansl 1990), beating Euler's Gradus (ρ = 0.979), and has four equivalent musical interpretations: Galileo's pulse-coincidence model (1638), harmonic series position (Rameau 1722), Farey index (Cauchy 1816), and first coincidence harmonic.

- **f(p/q) = p + Ω\*(q) also achieves ρ = 0.989**, matching max(p,q) while distinguishing more intervals. It differs from Gradus only in replacing Ω\*(p) with the raw partial number p.

- **One stubborn tie remains:** the major third (5/4) and major sixth (5/3) both give f = 7, because Ω\*(4) = Ω\*(3) = 2. Resolving this requires either a richer formula or cultural context.

- **The formula defines a coprime integer triangle** T(n,k) = n + Ω\*(k) whose rightmost diagonal gives the two-stage dissonance of the superparticular intervals n/(n-1). Both are submitted to the OEIS as A397104 and A397106.

---

## Two-stage perceptual hypothesis

The formula f(p/q) = p + Ω\*(q) decomposes as two perceptual costs, offered as a speculative hypothesis:

| Stage | Cost | Interpretation |
|-------|------|----------------|
| 1 — Bass establishes context | Ω\*(q) | Lower note is the q-th partial of an implied fundamental; prime complexity of extrapolating down to it (analytical proxy for subharmonic template matching) |
| 2 — Upper note reaches into series | p | Upper note must be recognised as the p-th partial; higher partials are quieter and harder to confirm |

This connects Galileo (1638), Rameau (1722), Euler (1739), Terhardt (1979), Parncutt (1988), and Huron (2001) within a single arithmetic framework.

---

## Files

| File | Description |
|------|-------------|
| `note2.tex` / `note2.pdf` | Working note with full derivations and results (current version) |
| `note.tex` / `note.pdf` | Earlier version of the working note |
| `harmonic_overlap.py` | Initial exploration of harmonic coincidence vs Gradus |
| `gradus_analysis.py` | Spearman correlations; Gradus ties and failures |
| `weighted_coincidence.py` / `weighted_coincidence2.py` | Proof that Gradus = weighted coincidence count |
| `perceptual_model.py` | Plomp–Levelt roughness model comparison |
| `gradus_failures.py` | Where Gradus fails; max(p,q) as tiebreaker |
| `maxpq_interpretation.py` | Four interpretations of max(p,q) |
| `formula_search.py` | Systematic search for improved formula |
| `table.py` / `table_plot.py` | Comparison table: Gradus, f(p/q), human ratings |
| `perceptual_account.py` / `perceptual_stem.py` | Two-stage model diagrams |
| `tolerance_model.py` | Partial-beating tolerance model (§9 of note) |
| `well_plot.py` | Consonance wells near just ratios (§5 of note) |
| `oeis_draft.md` | OEIS submission draft for triangle T(n,k) = n + Ω\*(k) |
| `oeis_generate.py` | Generate and verify triangle terms and b-file |
| `oeis_superparticular.py` | OEIS submission draft for superparticular subsequence a(n) = n + Ω\*(n−1) |
| `b397104.txt` | OEIS b-file: triangle T(n,k), 774 terms |
| `b397105.txt` | OEIS b-file: superparticular sequence, 300 terms |

---

## OEIS

Two entries have been submitted to the [OEIS](https://oeis.org) (under editorial review):

- **[A397104](https://oeis.org/A397104)** — Main triangle T(n,k) = n + Ω\*(k), for each row n, k ranges over integers with 1 ≤ k ≤ n and gcd(n,k) = 1; row length is φ(n) for n ≥ 2.  
  First 40 terms: `1, 2, 3, 4, 4, 6, 5, 6, 7, 7, 6, 10, 7, 8, 9, 9, 11, 10, 8, 10, 12, 14, 9, 10, 11, 13, 15, 12, 10, 12, 16, 14, 11, 12, 13, 13, 15, 14, 17, 14`

- **[A397106](https://oeis.org/A397106)** — Superparticular subsequence a(n) = n + Ω\*(n−1), n ≥ 2  
  First 40 terms: `2, 4, 6, 7, 10, 10, 14, 12, 14, 16, 22, 17, 26, 22, 22, 21, 34, 24, 38, 27, 30, 34, 46, 30, 34, 40, 34, 37, 58, 38`

The author is the creator of [A275314](https://oeis.org/A275314) (Euler's Gradus Suavitatis).

---

## References

- L. Euler, *Tentamen novae theoriae musicae*, St. Petersburg, 1739.
- G. Galilei, *Discorsi e Dimostrazioni Matematiche*, 1638.
- P. Barbieri, "Galileo's" coincidence theory of consonances, from Nicomachus to Sauveur, *Recercare* 13 (2001), 201–232.
- J.-P. Rameau, *Traité de l'Harmonie*, Paris, 1722.
- H. von Helmholtz, *On the Sensations of Tone*, 1863.
- R. Plomp & W. J. M. Levelt, Tonal consonance and critical bandwidth, *JASA* 38 (1965), 548–560.
- E. Terhardt, Calculating virtual pitch, *Hearing Research* 1 (1979), 155–182.
- R. Parncutt, Revision of Terhardt's psychoacoustical model of the root(s) of a musical chord, *Music Perception* 6 (1988), 65–93.
- C. L. Krumhansl, *Cognitive Foundations of Musical Pitch*, Oxford University Press, 1990.
- J. Tenney, *A History of 'Consonance' and 'Dissonance'*, 1988.
- D. Huron, Tone and voice, *Music Perception* 19 (2001), 1–64.
- I. Lahdelma & T. Eerola, Single chords convey distinct emotional qualities, *Psychology of Music* 44 (2016), 37–54.
- J. H. McDermott et al., Indifference to dissonance in native Amazonian, *Nature* 535 (2016), 547–550.
- D. L. Bowling, D. Purves & K. Z. Gill, Vocal similarity predicts chord attraction, *PNAS* 115 (2018), 216–221.
- P. Harrison & M. Pearce, Simultaneous consonance in music perception and composition, *Psychological Review* 127 (2020), 216–244.
- I. Lahdelma & T. Eerola, The anatomy of consonance/dissonance, *Music & Science* 4 (2021), 1–19.
