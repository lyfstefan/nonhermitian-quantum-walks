"""Publication-template figures from calculated finite-ring quantum walks."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
from qwalk import Walk
from qwalk.figures import draw_cases
from validation.run import match
from validation.paper_anchors import TABLES


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
