# Reconstruction record (October 2026)

The historical undergraduate source was unavailable in the supplied local
materials. All Python in this repository was newly written from the published
operators; no claim about historical code authorship or recovery is made.

1. The source is the 13-page published PRA article, DOI
   `10.1103/PhysRevA.102.022218`, pp. 1-13. Equations (1)-(2) and their S/R/G
   definitions were checked visually on p. 2. The defect multiplies the input
   ket site. Read products right to left.
2. Appendix A, p. 6, supplies d0(k) in (A2)/(A8). Only these analytic formulas
   and the explicitly stated rotations/gain matrices are reconstructed. The
   transfer-recursion solution (B1)-(B33) is not implemented or independently
   proved here. Numerical eigenstates replace it in the example plots.
3. Appendix D, pp. 10-12, Tables I-VI supplies numeric probability anchors.
   Five even positions and four columns per table were transcribed and
   visually checked. Odd entries are numerical roundoff and tested as zero.
   There is no optimization against these numbers; the angles and gamma
   are fixed to the article values, and the matching step only permutes modes.
4. Main Figs. 2-4 on pp. 3-4 establish the peak-shift check and gamma choices.
   Fig. 5 on p. 5 uses U1 gamma=.5 for its completely broken phase scan/time
   panel, while Fig. 4 uses .3. These distinct settings are kept separate.
   Figs. 6-9 and Appendix C on pp. 8-10 were also inspected; only the six
   phase-pi peak-shift cases are numerically checked here, not every scan.
5. The full periodic lattice contains 2N+1 sites. Although two-shift evolution
   preserves parity on the infinite line, odd periodic rings eventually mix
   it at the boundary. Evolution sizes are therefore chosen larger than the
   exact light cone; eigenstate tests use size and edge diagnostics.

The first exploratory localization filter (radius 8) omitted the broader U1
gamma=.3 side modes and included odd-sector modes. The final filter states its
radius 20 and even-sector criterion explicitly. The full matrix spectrum is
always available through `eigensystem`; this filter is a diagnostic, not a
general infinite-volume spectral classification.

All generated probabilities and figures are new numerical results. Agreement
with printed numerical anchors is a bounded reproduction claim for those
anchors; qualitative peak/delocalization checks do not establish that every
figure, every transient, or every analytic derivation has been reproduced.
