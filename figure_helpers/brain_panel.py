import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting, image
from matplotlib import text as mtext
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--title", default="")
    p.add_argument("--cmap", default="viridis")
    p.add_argument("--vmin", type=float, default=None)
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument("--threshold", default=None)
    p.add_argument("--kind", choices=["stat", "roi"], default="stat")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--symmetric", action="store_true")
    p.add_argument("--display", default="tiled")
    p.add_argument("--bg", default=None)
    p.add_argument("--cut", default="-8,-20,18")
    args = p.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    cmap_alias = {
        "coolwarm": "RdBu_r",
        "Reds": "Reds",
        "Purples": "Purples",
        "viridis": "viridis",
        "magma": "magma",
    }
    cmap = cmap_alias.get(args.cmap, args.cmap)

    img = image.load_img(args.input)
    bg = image.load_img(args.bg) if args.bg else None
    data = np.asarray(img.get_fdata(), dtype=float)
    vals = data[np.isfinite(data) & (data != 0)]
    vmin = args.vmin
    vmax = args.vmax
    if vals.size:
        if args.symmetric:
            lim = np.nanpercentile(np.abs(vals), 98.5) if vmax is None else abs(vmax)
            lim = float(lim) if np.isfinite(lim) and lim > 0 else 1.0
            vmin, vmax = -lim, lim
        elif args.kind == "stat" and vmax is None:
            lo = np.nanpercentile(vals, 1.0)
            hi = np.nanpercentile(vals, 99.0)
            vmin = float(lo) if vmin is None and lo < 0 else vmin
            vmax = float(hi) if np.isfinite(hi) and hi > 0 else vmax

    fig = plt.figure(figsize=(3.18, 2.05), dpi=args.dpi, facecolor="white")
    threshold = 1e-9 if args.threshold in (None, "", "none", "None") else float(args.threshold)
    common = dict(
        display_mode=args.display,
        cut_coords=tuple(float(x) for x in args.cut.split(",")),
        annotate=False,
        draw_cross=False,
        black_bg=False,
        figure=fig,
        colorbar=True,
        title="",
    )
    try:
        if args.kind == "roi":
            plotting.plot_stat_map(
                img,
                bg_img=bg,
                cmap=cmap,
                colorbar=False,
                threshold=0.5,
                vmin=0,
                vmax=1,
                dim=-0.35,
                **{k: v for k, v in common.items() if k != "colorbar"},
            )
        else:
            display = plotting.plot_stat_map(
                img,
                bg_img=bg,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                threshold=threshold,
                symmetric_cbar=args.symmetric,
                dim=-0.35,
                **common,
            )
    except Exception:
        plt.close(fig)
        fig = plt.figure(figsize=(3.18, 1.95), dpi=args.dpi, facecolor="white")
        common["figure"] = fig
        common["display_mode"] = "z"
        common["cut_coords"] = 7
        if args.kind == "roi":
            plotting.plot_stat_map(
                img,
                bg_img=bg,
                cmap=cmap,
                colorbar=False,
                threshold=0.5,
                vmin=0,
                vmax=1,
                dim=-0.35,
                **{k: v for k, v in common.items() if k != "colorbar"},
            )
        else:
            display = plotting.plot_stat_map(
                img,
                bg_img=bg,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                threshold=threshold,
                symmetric_cbar=args.symmetric,
                dim=-0.35,
                **common,
            )

    fig.patch.set_facecolor("white")
    # Panel titles are added as editable PowerPoint text boxes.
    # Nilearn's default colorbar is taller than the displayed brain panels.
    # Shorten narrow colorbar axes so the panel wastes less vertical space and
    # tick labels are not clipped during tight export.
    for ax in fig.axes:
        pos = ax.get_position()
        if pos.x0 > 0.80 and pos.width < 0.12 and pos.height > 0.35:
            new_h = min(0.40, pos.height)
            center = pos.y0 + pos.height / 2
            ax.set_position([min(pos.x0 + 0.025, 0.935), center - new_h / 2, pos.width * 0.72, new_h])
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("right")
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(labelright=True, labelleft=False, pad=2)
    for ax in fig.axes:
        ax.set_facecolor("white")
        ax.tick_params(labelsize=9, length=2, width=0.45, pad=2)
    for obj in fig.findobj(mtext.Text):
        obj.set_fontsize(10)
        obj.set_fontfamily("Arial")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.025, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
