import sys

from pathlib import Path
from io import BytesIO
from typing import Mapping

from PIL import Image
from matplotlib.axes import Axes

import pandas as pd
import numpy as np


class Tee():
    def __init__(self, streams):
        self.streams = streams
        self.orig_stdout = sys.stdout
        self.orig_stderr = sys.stderr

    def write(self, out):
        for stream in self.streams:
            stream.write(out)
    
    def __enter__(self):
        sys.stdout = self
        sys.stderr = self

    def __exit__(self, exc_type, exc_value, exc_tb):
        for stream in self.streams:
            stream.close()
        sys.stdout = self.orig_stdout
        sys.stderr = self.orig_stderr


def sort_adata_by_obs_names(adata):
    return adata[np.argsort(adata.obs_names), :]


def modify_index(df: pd.DataFrame, map_func):
    new_index = pd.Series(df.index.map(map_func), index=df.index)
    assert new_index.is_unique, "new index is not unique"
    df.index = new_index


def get_subdirs(p: Path):
    return [d for d in p.iterdir() if d.is_dir()]


def convert_rgba_to_rgb(img: Image):
    if img.mode == "RGBA":
        rgb = Image.new("RGB", img.size, (255,255,255)) 
        rgb.paste(img, mask=img.split()[3])
        return rgb
    else:
        return img


def ax2image(ax: Axes):
    buf = BytesIO()
    ax.get_figure().savefig(buf, transparent=False, dpi=600)
    buf.seek(0)
    img = Image.open(buf)
    img = convert_rgba_to_rgb(img)

    return img


def save_axes(axes: Mapping[str, Axes], out_fname: Path):
    figs = []
    for label, ax in axes.items():
        ax.set_title(label)
        figs.append(ax2image(ax))

    figs[0].save(
        out_fname,
        "PDF",
        dpi=(600, 600),
        save_all=True,
        append_images=figs[1:]
    )