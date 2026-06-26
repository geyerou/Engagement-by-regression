import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle
from nilearn import plotting

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 0.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def _centroids(img):
    data = np.asarray(img.get_fdata())
    labs = np.array([v for v in np.unique(data) if v > 0], dtype=int)
    rows = []
    for lab in labs:
        ijk = np.argwhere(data == lab)
        if ijk.size == 0:
            continue
        xyz = nib.affines.apply_affine(img.affine, ijk.mean(axis=0))
        rows.append((lab, ijk.shape[0], xyz))
    return rows


def _choose_parcels(img):
    rows = _centroids(img)
    specs = [
        ("posterior occipital", lambda xyz: xyz[1] < -55 and xyz[2] > -5),
        ("left temporal", lambda xyz: xyz[0] < -35 and -45 < xyz[1] < 15 and xyz[2] > -15),
        ("right frontal", lambda xyz: xyz[0] > 25 and xyz[1] > 0 and xyz[2] > 5),
    ]
    chosen = []
    used = set()
    for _, pred in specs:
        candidates = [r for r in rows if pred(r[2]) and r[0] not in used]
        if not candidates:
            candidates = [r for r in rows if r[0] not in used]
        lab, _, xyz = max(candidates, key=lambda r: r[1])
        chosen.append((lab, xyz))
        used.add(lab)
    return chosen


def _label_img(source_img, label):
    data = np.asarray(source_img.get_fdata())
    roi = (data == label).astype(np.int16)
    return nib.Nifti1Image(roi, source_img.affine, source_img.header)


def _wm_target(wm_img):
    data = np.asarray(wm_img.get_fdata()) > 0
    ijk_all = np.argwhere(data)
    xyz_all = nib.affines.apply_affine(wm_img.affine, ijk_all)
    desired = np.array([18.0, -22.0, 24.0])
    ix = int(np.argmin(np.sum((xyz_all - desired) ** 2, axis=1)))
    ijk = tuple(int(v) for v in ijk_all[ix])
    xyz = xyz_all[ix]
    single = np.zeros(data.shape, dtype=np.int16)
    single[ijk] = 1
    return ijk, xyz, nib.Nifti1Image(single, wm_img.affine, wm_img.header)


def _add_roi_display(fig, roi_img, bg_img, xyz, color, box, title):
    ax = fig.add_axes(box)
    display = plotting.plot_roi(
        roi_img,
        bg_img=bg_img,
        display_mode="ortho",
        cut_coords=tuple(float(v) for v in xyz),
        annotate=False,
        draw_cross=False,
        black_bg=False,
        axes=ax,
        cmap=ListedColormap([color]),
        alpha=0.90,
        dim=-0.20,
        colorbar=False,
        title="",
    )
    fig.text(box[0] - 0.025, box[1] + box[3] + 0.014, title,
             ha="left", va="bottom", fontsize=9.5, fontweight="bold", color="#202020")
    return display


def _draw_grid_inset(fig, box):
    ax = fig.add_axes(box)
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")
    for x in range(7):
        for y in range(7):
            face = "#f7e4e3" if abs(x - 3) + abs(y - 3) < 4 else "#f7f7f7"
            edge = "#d0d0d0"
            ax.add_patch(Rectangle((x, y), 1, 1, facecolor=face, edgecolor=edge, linewidth=0.45))
    ax.add_patch(Rectangle((3, 3), 1, 1, facecolor="#b2182b", edgecolor="#111111", linewidth=0.9))
    ax.text(3.5, -0.55, "single WM voxel", ha="center", va="top", fontsize=9, color="#222222")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--atlas", required=True)
    p.add_argument("--wm", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bg", default=None)
    p.add_argument("--dpi", type=int, default=600)
    args = p.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    atlas_img = nib.load(args.atlas)
    wm_img = nib.load(args.wm)
    bg_img = nib.load(args.bg) if args.bg else None

    parcels = _choose_parcels(atlas_img)
    colors = ["#2b83ba", "#3f9b59", "#6a51a3"]

    fig = plt.figure(figsize=(9.6, 3.45), dpi=args.dpi, facecolor="white")

    fig.text(0.030, 0.935, "Cortical GM parcels", fontsize=12, fontweight="bold",
             ha="left", va="top", color="#222222")
    y_boxes = [0.640, 0.405, 0.170]
    starts = []
    for i, ((lab, xyz), color, y0) in enumerate(zip(parcels, colors, y_boxes), start=1):
        _add_roi_display(
            fig,
            _label_img(atlas_img, lab),
            bg_img,
            xyz,
            color,
            [0.060, y0, 0.285, 0.178],
            f"GM{i}: Schaefer parcel {lab}",
        )
        starts.append((0.355, y0 + 0.088))

    fig.text(0.485, 0.935, "Nested ridge encoding", fontsize=12, fontweight="bold",
             ha="center", va="top", color="#222222")
    fig.text(0.190, 0.070, "...", fontsize=18, fontweight="bold",
             ha="center", va="center", color="#444444")
    fig.text(0.190, 0.034, "Schaefer-400 parcels", fontsize=9.5,
             ha="center", va="center", color="#333333")

    wm_ijk, wm_xyz, voxel_img = _wm_target(wm_img)
    ax = fig.add_axes([0.660, 0.315, 0.305, 0.470])
    display = plotting.plot_roi(
        wm_img,
        bg_img=bg_img,
        display_mode="ortho",
        cut_coords=tuple(float(v) for v in wm_xyz),
        annotate=False,
        draw_cross=False,
        black_bg=False,
        axes=ax,
        cmap=ListedColormap(["#f4a3a3"]),
        alpha=0.42,
        dim=-0.20,
        colorbar=False,
        title="",
    )
    display.add_overlay(voxel_img, threshold=0.5, cmap=ListedColormap(["#b2182b"]), alpha=1.0)
    fig.text(0.662, 0.805, f"WM95 target voxel vᵢ", fontsize=12, fontweight="bold",
             ha="left", va="top", color="#222222")
    fig.text(0.662, 0.276, f"MNI [{wm_xyz[0]:.0f}, {wm_xyz[1]:.0f}, {wm_xyz[2]:.0f}]",
             fontsize=8.5, ha="left", va="top", color="#444444")
    _draw_grid_inset(fig, [0.840, 0.115, 0.105, 0.210])

    for start, color in zip(starts, colors):
        arrow = FancyArrowPatch(
            start,
            (0.650, 0.540),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=1.15,
            color="#333333",
            alpha=0.78,
            connectionstyle="arc3,rad=0.12",
            zorder=1,
        )
        fig.patches.append(arrow)
    arrow = FancyArrowPatch(
        (0.560, 0.500),
        (0.650, 0.540),
        transform=fig.transFigure,
        arrowstyle="-|>",
        mutation_scale=11,
        linewidth=1.0,
        color="#333333",
        alpha=0.78,
        connectionstyle="arc3,rad=0.00",
    )
    fig.patches.append(arrow)

    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.025, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
