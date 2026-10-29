"""Compact layered eigenstate bars and unchanged normalized evolution maps."""

from io import BytesIO
from pathlib import Path
import hashlib
import json
import os
import tempfile

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
import numpy as np

from .publication_style import (
    C,
    apply_style,
    all_axes,
    panel_colorbar,
    pin_roman_face,
    unify_typography,
)

PROFILE_COLORS = [C[key] for key in ("blue", "red", "green", "purple")]
SITE_BAR_WIDTH = 0.85
LAYER_BAR_WIDTH = 0.45
CAMERA = (22, -55)
BOX_ASPECT = (1.4, 1.2, 1.8)


def profile_labels(model, gamma, moduli):
    labels = [f"{value:.2f}" for value in moduli]
    if model == "U2" and np.isclose(gamma, 0.1):
        labels[2:] = [f"{value:.4f}" for value in moduli[2:]]
    return labels


def draw_profile(ax, model, gamma, data):
    """Each layer is one right eigenstate; bar heights are the supplied P_n."""
    sites = np.asarray(data["eigen_sites"])
    profiles = np.asarray(data["eigen_profiles"])
    moduli = np.asarray(data["moduli"])
    if profiles.shape != (len(sites), 4) or moduli.shape != (4,):
        raise ValueError("Expected four eigenstate profiles and four moduli")
    if not np.isfinite(profiles).all() or np.any(profiles < 0):
        raise ValueError("Profile probabilities must be finite and nonnegative")
    use = np.abs(sites) <= 12
    n, values = sites[use], profiles[use]
    labels = profile_labels(model, gamma, moduli)
    ax.set_proj_type("ortho")
    for state, color in enumerate(PROFILE_COLORS):
        ax.bar3d(
            n - SITE_BAR_WIDTH / 2,
            np.full(len(n), state - LAYER_BAR_WIDTH / 2),
            np.zeros(len(n)),
            np.full(len(n), SITE_BAR_WIDTH),
            np.full(len(n), LAYER_BAR_WIDTH),
            values[:, state],
            color=color,
            edgecolor=color,
            linewidth=0.25,
            shade=True,
            zsort="average",
        )
        ax.plot(
            n, np.full(len(n), state), np.zeros(len(n)), color=color, lw=0.45, alpha=0.4
        )
    ax.set(xlim=(-12.8, 12.8), ylim=(-0.45, 3.5), zlim=(0, 1.08))
    ax.set_xticks([-10, 0, 10])
    value_labels = [
        rf"$|\lambda_{{{index+1}}}|={value}$" for index, value in enumerate(labels)
    ]
    ax.set_yticks([])
    ax.set_zticks([0, 0.5, 1])
    ax.set_xlabel(r"$n$", fontsize=10, labelpad=1)
    ax.set_ylabel("")
    ax.set_zlabel(r"$P$", fontsize=10, labelpad=1)
    ax.set_title(rf"$U_{{{model[-1]}}}$, $\gamma={gamma}$", fontsize=9, y=1.02, pad=0)
    ax.view_init(elev=CAMERA[0], azim=CAMERA[1])
    ax.zaxis.set_ticks_position("lower")
    ax.zaxis.set_label_position("lower")
    ax.zaxis.set_major_formatter(FuncFormatter(lambda value, pos: rf"${value:g}$"))
    ax.set_box_aspect(BOX_ASPECT)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("none")
        axis.line.set_color(C["black"])
        axis.line.set_linewidth(1.0)
        axis.set_tick_params(labelsize=9, pad=1)
    ax.yaxis.set_tick_params(labelsize=7)
    payload = {
        "kind": "layered_eigenstate_bars",
        "model": model,
        "gamma": float(gamma),
        "projection": "orthographic",
        "view_elevation_azimuth_degrees": list(CAMERA),
        "box_aspect": list(BOX_ASPECT),
        "depth_ticks_visible": False,
        "site_bar_width": SITE_BAR_WIDTH,
        "layer_bar_width": LAYER_BAR_WIDTH,
        "probability_limits": [0, 1.08],
        "depth_coordinate": "categorical state rows, equally spaced for display",
        "eigenvalue_moduli": moduli.tolist(),
        "displayed_modulus_labels": labels,
        "displayed_layer_labels": [],
        "displayed_value_labels": value_labels,
        "layer_indices_front_to_back": [1, 2, 3, 4],
        "sites": n.tolist(),
        "probabilities": values.tolist(),
    }
    ax._publication_profile = payload
    return payload


def build_case_figure(cases, font_path=None):
    """Pair each profile with its original evolution map in a fixed-size cell."""
    if not cases:
        raise ValueError("At least one quantum-walk case is required")
    apply_style(font_path=font_path)
    columns = 2 if len(cases) > 2 else 1
    rows = (len(cases) + columns - 1) // columns
    case_width, case_height = 10.4, 4.8
    page_width, page_height = columns * case_width, rows * case_height + 0.80
    fig = plt.figure(figsize=(page_width / 2.54, page_height / 2.54))
    evolution = []
    for index, (model, gamma, data) in enumerate(cases):
        row, column = divmod(index, columns)
        x, y = column * case_width, 0.65 + (rows - row - 1) * case_height
        profile = fig.add_axes(
            [
                (x + 1.60) / page_width,
                (y + 0.40) / page_height,
                3.6 / page_width,
                3.6 / page_height,
            ],
            projection="3d",
        )
        payload = draw_profile(profile, model, gamma, data)
        fig.legend(
            [Patch(facecolor=color) for color in PROFILE_COLORS],
            payload["displayed_value_labels"],
            loc="upper left",
            bbox_to_anchor=((x + 4.75) / page_width, (y + 4.15) / page_height),
            bbox_transform=fig.transFigure,
            frameon=False,
            fontsize=7,
            ncol=1,
            handlelength=0.6,
            handletextpad=0.3,
            borderaxespad=0,
            borderpad=0.1,
            labelspacing=0.15,
        )
        ax = fig.add_axes(
            [
                (x + 7.0) / page_width,
                (y + 1.05) / page_height,
                2.5 / page_width,
                2.5 / page_height,
            ],
        )
        use = np.abs(data["sites"]) <= 30
        probability = np.asarray(data["probability"])[:, use]
        im = ax.imshow(
            probability,
            origin="lower",
            aspect="auto",
            extent=[-30.5, 30.5, -0.5, 300.5],
            cmap="magma",
            vmin=0,
            vmax=1,
        )
        ax.set(xlabel=r"$n$", ylabel=r"$T$", xticks=[-20, 0, 20])
        ax.set_box_aspect(1)
        cb = panel_colorbar(fig, ax, im, r"$P$")
        cb.ax.set_title(r"$P$", loc="left", pad=6, fontsize=9)
        evolution.append(
            {
                "model": model,
                "gamma": float(gamma),
                "displayed_shape": list(probability.shape),
                "displayed_probability_sha256": hashlib.sha256(
                    np.ascontiguousarray(probability).tobytes()
                ).hexdigest(),
                "probability_color_limits": [0, 1],
                "site_limits": [-30.5, 30.5],
                "time_limits": [-0.5, 300.5],
            }
        )
    fig._publication_evolution = evolution
    return fig


def _atomic_write(path, content):
    with tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=".figure-", delete=False
    ) as stream:
        temporary = Path(stream.name)
        stream.write(content)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def save_case_figure(fig, destination):
    """Export mixed 3D/2D axes without applying a 2D aspect ratio to 3D bars."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    unify_typography(fig)
    for ax in all_axes(fig):
        if hasattr(ax, "_publication_profile"):
            ax.yaxis.label.set_fontsize(9)
            ax.yaxis.set_tick_params(labelsize=7)
            ax.zaxis.label.set_fontsize(10)
            ax.zaxis.set_tick_params(labelsize=9)
        else:
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)
            ax.tick_params(which="both", top=False, right=False, width=0.9)
            ax.tick_params(which="major", length=3.0)
        if hasattr(ax, "_colorbar"):
            ax._colorbar.outline.set_linewidth(0.9)
    for legend in fig.legends:
        for text in [*legend.get_texts(), legend.get_title()]:
            text.set_fontsize(7)
    pin_roman_face(fig)
    fig.canvas.draw()
    frames, projected, bars = [], [], []
    for extension, kwargs in [("pdf", {}), ("png", {"dpi": 300}), ("svg", {})]:
        buffer = BytesIO()
        fig.savefig(buffer, format=extension, **kwargs)
        if extension == "pdf":
            for ax in fig.axes:
                box = ax.get_position()
                geometry = {
                    "bounds_figure_fraction": list(box.extents),
                    "width_cm": float(box.width * fig.get_figwidth() * 2.54),
                    "height_cm": float(box.height * fig.get_figheight() * 2.54),
                }
                if hasattr(ax, "_publication_profile"):
                    projected.append({**ax._publication_profile, **geometry})
                else:
                    frames.append(geometry)
            for ax in all_axes(fig):
                if not hasattr(ax, "_publication_parent"):
                    continue
                parent = ax._publication_parent
                box, frame = ax.get_position(), parent.get_position()
                title_box = fig.transFigure.inverted().transform_bbox(
                    ax._left_title.get_window_extent(fig.canvas.get_renderer())
                )
                bars.append(
                    {
                        "bounds_figure_fraction": list(box.extents),
                        "parent_bounds_figure_fraction": list(frame.extents),
                        "height_fraction_of_frame": float(box.height / frame.height),
                        "width_fraction_of_frame": float(box.width / frame.width),
                        "bottom_offset_cm": float(
                            (box.y0 - frame.y0) * fig.get_figheight() * 2.54
                        ),
                        "label": ax._publication_label,
                        "label_position": "above",
                        "shared": False,
                        "label_bottom_clearance_cm": float(
                            (title_box.y0 - box.y1) * fig.get_figheight() * 2.54
                        ),
                        "scientific_exponent": ax._publication_exponent,
                    }
                )
        _atomic_write(
            destination.parent / (destination.name + "." + extension), buffer.getvalue()
        )
    metadata = {
        "template": "v1_compact_layered_bars_and_square_evolution",
        "text_engine": "latex",
        "text_usetex": True,
        "frame_pt": 1.0,
        "axis_labels_pt": 10,
        "ticks_pt": 9,
        "colorbar_ticks_pt": 9,
        "annotations_pt": 9,
        "depth_labels_pt": 7,
        "legends_pt": 7,
        "page_width_cm": float(fig.get_figwidth() * 2.54),
        "page_height_cm": float(fig.get_figheight() * 2.54),
        "font_family": "Computer Modern Roman",
        "math_font_family": "Computer Modern Math Italic",
        "target_frame_side_cm": 2.5,
        "data_frames": frames,
        "three_dimensional_axes": projected,
        "colorbars": bars,
        "evolution_maps": fig._publication_evolution,
        "latex_preamble": r"\usepackage{amsmath}\usepackage{textcomp}",
        "pinned_font_file": None,
    }
    (destination.parent / (destination.name + "_typography.json")).write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


def draw_cases(cases, destination, font_path=None):
    fig = build_case_figure(cases, font_path=font_path)
    try:
        save_case_figure(fig, destination)
    finally:
        plt.close(fig)
