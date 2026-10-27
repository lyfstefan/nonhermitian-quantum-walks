"""Publication-template figures from calculated finite-ring quantum walks."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qwalk import Walk
from qwalk.publication_style import (
    apply_style,
    width,
    save,
    C,
    DASH,
    panel_colorbar,
    inset_legend,
    legend_value,
    INSET_LEGEND_FONT,
)
from validation.run import match
from validation.paper_anchors import TABLES


def draw_cases(cases, destination, font_path=None):
    """Render saved arrays without changing probability normalization or limits."""
    apply_style(font_path=font_path)
    overview = len(cases) > 2
    shape = (3, 4) if overview else (len(cases), 2)
    fig, axes = plt.subplots(
        *shape,
        figsize=(width(2), width(2) * (0.87 if overview else 1.06)),
        constrained_layout=True,
        squeeze=False,
    )
    eigen_axes = []
    evolution_axes = []
    for case_index, (model, gamma, data) in enumerate(cases):
        row = case_index // 2 if overview else case_index
        col = (case_index % 2) * 2 if overview else 0
        left, right = axes[row, col], axes[row, col + 1]
        eigen_axes.append(left)
        evolution_axes.append(right)
        sites, profiles = data["eigen_sites"], data["eigen_profiles"]
        use = abs(sites) <= 12
        labels = [legend_value(value) for value in data["moduli"]]
        if model == "U2" and np.isclose(gamma, 0.1):
            # Only this panel needs to distinguish its nearby small moduli.
            labels[2:] = [f"{value:.4f}" for value in data["moduli"][2:]]
        for index, (value, profile) in enumerate(zip(data["moduli"], profiles.T)):
            left.plot(
                sites[use],
                profile[use],
                marker=".",
                markersize=2.5,
                color=[C["blue"], C["red"], C["green"], C["purple"]][index],
                linestyle=DASH[index],
                label=f"Mode {index+1}" if overview else rf"${labels[index]}$",
            )
        left.set(
            title=rf"$U_{{{model[-1]}}}$, $\gamma={gamma}$",
            xlabel=r"$n$",
            ylabel=r"$P$",
        )
        left.set_xticks([-10, 0, 10])
        if not overview:
            inset_legend(
                left, title=r"$|\lambda|$", data_top=float(np.max(profiles[use]))
            )
        use = abs(data["sites"]) <= 30
        im = right.imshow(
            data["probability"][:, use],
            origin="lower",
            aspect="auto",
            extent=[-30.5, 30.5, -0.5, 300.5],
            cmap="magma",
            vmin=0,
            vmax=1,
        )
        right.set(title="", xlabel=r"$n$", ylabel=r"$T$")
        right.set_xticks([-20, 0, 20])
        panel_colorbar(fig, right, im, r"$P$")
    if overview:
        handles, labels = eigen_axes[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="outside upper center",
            ncol=4,
            frameon=False,
            handlelength=1.4,
            columnspacing=1.0,
        )
    save(fig, destination, legend_fontsize=8 if overview else INSET_LEGEND_FONT)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("examples/output"))
    parser.add_argument(
        "--reuse-data",
        action="store_true",
        help="render existing NPZ arrays without recalculation",
    )
    parser.add_argument(
        "--font-path", type=Path, help="legacy option; LaTeX selects Computer Modern"
    )
    parser.add_argument(
        "--data-dir", type=Path, help="optional directory of cached NPZ arrays"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data_dir = args.data_dir or args.output
    data_dir.mkdir(parents=True, exist_ok=True)
    cases = []
    for model, gamma, anchors in TABLES.values():
        path = data_dir / f"{model}_gamma_{gamma}.npz"
        if not args.reuse_data:
            w = Walk(model, gamma, 80)
            values, profiles, _ = match(w, np.asarray(anchors))
            large = Walk(model, gamma, 650)
            np.savez_compressed(
                path,
                sites=large.sites,
                probability=large.evolve(),
                moduli=abs(values),
                eigen_sites=w.sites,
                eigen_profiles=profiles,
            )
        with np.load(path) as stored:
            data = {name: stored[name].copy() for name in stored.files}
        cases.append((model, gamma, data))
    draw_cases(cases, args.output / "localization", font_path=args.font_path)
    for index, regime in enumerate(
        ["unbroken", "partially_broken", "completely_broken"]
    ):
        draw_cases(
            cases[2 * index : 2 * index + 2],
            args.output / f"localization_{regime}",
            font_path=args.font_path,
        )
    print(args.output / "localization.png")


if __name__ == "__main__":
    main()
