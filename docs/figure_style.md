# Figure style v0

Figures use the user-supplied `publication_style.py` template without changes
to its style functions. Axis labels are 10 pt; ticks, legends and annotations
are 9 pt; frames are 1 pt. PNG exports are 300 dpi. Standalone PDF pages are
17.8 cm wide, with no bounding-box crop. SVG versions remain available.

The local release explicitly pins `lmroman10-regular.otf` (LM Roman 10 Regular);
PDF font resources were checked for `LMRoman10-Regular`. Math uses Matplotlib's
Computer Modern (`mathtext.fontset = cm`). Numerical arrays, fitted parameters,
probability normalization and heat-map color limits are unchanged. Captions
are kept in the project documentation; long panels have been rearranged.

For portable regeneration, use `--font-path /path/to/lmroman10-regular.otf`.
Without a supplied font file, the template uses its installed-font fallback
list. Some LM OpenType files report the family name `Latin Modern Roman`, so
the literal `LM Roman 10` name alone does not always find them. Pinning the
file avoids this mismatch. The font file itself is not bundled. Each figure's
`*_typography.json` records the exact file used and printed page width.

```sh
python -m examples.paper_cases --reuse-data --data-dir examples/output \
  --output docs --font-path /path/to/lmroman10-regular.otf
```

First generate the cached NPZ arrays with the normal example command if they
are absent. The overview contains six case rows; the three additional figures
separate unbroken, partially broken and completely broken bulk regimes.
