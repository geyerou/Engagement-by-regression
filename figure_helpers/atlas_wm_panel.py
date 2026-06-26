import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from nilearn import image, plotting

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 0.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--atlas", required=True)
    p.add_argument("--wm", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--title", default="")
    p.add_argument("--bg", default=None)
    p.add_argument("--cut", default="-8,-20,18")
    p.add_argument("--dpi", type=int, default=600)
    args = p.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    atlas = image.load_img(args.atlas)
    wm = image.load_img(args.wm)
    bg = image.load_img(args.bg) if args.bg else None
    cut = tuple(float(x) for x in args.cut.split(","))

    fig = plt.figure(figsize=(3.35, 2.05), dpi=args.dpi, facecolor="white")
    display = plotting.plot_roi(
        atlas,
        bg_img=bg,
        display_mode="ortho",
        cut_coords=cut,
        annotate=False,
        draw_cross=False,
        black_bg=False,
        figure=fig,
        cmap="gist_ncar",
        alpha=0.78,
        dim=-0.25,
        colorbar=False,
        title="",
    )
    wm_cmap = ListedColormap(["#b2182b"])
    display.add_overlay(wm, threshold=0.5, cmap=wm_cmap, alpha=0.78)
    # Panel titles are added as editable PowerPoint text boxes.
    for ax in fig.axes:
        ax.set_facecolor("white")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.025, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
