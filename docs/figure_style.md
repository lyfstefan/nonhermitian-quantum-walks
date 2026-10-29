# Figure conventions

Eigenstate probabilities are shown as compact three-dimensional layered bars.
Each layer contains one right eigenstate, the horizontal coordinate is the
physical site $n$, and the bar height is $P_n$. The depth coordinate is
categorical: the four layers have equal display spacing. The color legend
identifies the four eigenstates by $|\lambda_i|$, with $i=1,2,3,4$ ordered
from front to back, and gives their moduli. The depth axis has no tick labels;
its spacing is not a numerical scale for $|\lambda|$.

All cases use the same orthographic projection, elevation $22^\circ$,
azimuth $-55^\circ$, a coordinate-box aspect ratio of $1.4:1.2:1.8$,
and probability-axis range $0\leq P\leq1.08$.
Bars are centred on the original lattice sites; their width is 0.85 site
units and their thickness is 0.45 of the layer spacing. Values and
normalization are taken directly from the saved numerical arrays.

Each profile remains paired with the corresponding evolution heatmap.
The evolution arrays, site/time extents and color normalization $0\leq P\leq1$
are unchanged. The 3D viewport is 3.6 cm square; its projected coordinate box has a
comparable visual size to the 2.5 cm square heatmap frame. Each heatmap has its
own thin colorbar, 60% of the frame height and 4.5% of the frame width, bottom-aligned
with the frame and labeled above.

The overview has three rows and two case pairs per row. Each detailed figure
has two case pairs, stacked vertically. LaTeX and Computer Modern are used
throughout: axis labels 10 pt, ordinary ticks and model headings 9 pt, and
eigenvalue legends 7 pt. Axis lines are 1 pt. The
two nearby small moduli
in the $U_2$, $\gamma=0.1$ panel retain four decimal places so that they remain
distinguishable; the other displayed values use two.

PDF, SVG and 300 dpi PNG exports have identical geometry. The accompanying
typography JSON separates projected three-dimensional axes from square
heatmap frames and records the displayed profile values and evolution-data
checksums. [Figure validation](figure_validation.json) records the export
checks and unchanged scientific-input checksums.

## Dependencies

Use Matplotlib 3.10 or later and a LaTeX installation with `latex`, `dvipng`,
`type1cm`, `type1ec`, `amsmath` and `textcomp`. On Ubuntu:

```sh
sudo apt-get install texlive-latex-extra cm-super dvipng ghostscript
```

The preamble is `\usepackage{amsmath}\usepackage{textcomp}`. The legacy
`--font-path` argument remains accepted; TeX supplies the Computer Modern
fonts.
