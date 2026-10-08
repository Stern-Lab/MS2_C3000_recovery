"""
Whole-cell MS2 smFISH quantification: produces data/smFISH_<experiment>_per_cell.csv
(experiment 1: Figures 2D, S10; experiments 1 and 2: Table S3).

For each field in data/smFISH_<experiment>_image_manifest.csv with run_image = True:
  1. Cellpose-SAM segmentation (Cellpose 4.1.1) of the selected phase-contrast plane (segmentation_z_index).
  2. Removal of cells touching the image border; remaining labels renumbered 1..n.
  3. Mean projection of the RFP (TRITC) Z-stack; per-cell integrated signal, background-corrected with the
     median of a 0.20-0.50 um ring around each cell that excludes all segmented cells.
  4. whole_cell_signal_density = background-corrected integrated signal / cell area (pixels).

Fields excluded after segmentation QC are listed in the manifest with run_image = False and a reason in notes.

Environment: smFISH_env (envs/smFISH_env.yml)

Usage (from the repository root; image folders are set per experiment under smfish_nd2_dir in config.yaml):
    python preprocessing/smFISH_image_processing.py --experiment exp1
    python preprocessing/smFISH_image_processing.py --experiment exp2
Add --raw-mask-dir /path/to/masks to reuse existing Cellpose masks (<dir>/<image_name>/raw_mask.npy).
Masks saved by a previous run of this script in <out_dir>/smFISH/<experiment>/masks are reused automatically.
"""
import argparse
from pathlib import Path

import nd2
import numpy as np
import pandas as pd
import yaml
from scipy.ndimage import binary_dilation

CHANNELS = {"rfp": 0, "phase": 1, "dapi": 2}
FALLBACK_PIXEL_SIZE_UM = 0.064923
CELLPOSE_EVAL = {"diameter": 22, "flow_threshold": 0.4, "cellprob_threshold": 0.0}
BACKGROUND_RING_INNER_UM = 0.20
BACKGROUND_RING_OUTER_UM = 0.50
MIN_BACKGROUND_PIXELS = 50


def load_stack(path):
    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        if set(sizes) != {"Z", "C", "Y", "X"}:
            raise ValueError(f"{path.name}: expected Z/C/Y/X axes, found {sizes}")
        axes = list(sizes)
        stack = np.transpose(f.asarray(), [axes.index(k) for k in ("Z", "C", "Y", "X")])
        try:
            pixel_size = float(f.voxel_size().x)
            if not (np.isfinite(pixel_size) and pixel_size > 0):
                raise ValueError
        except (AttributeError, TypeError, ValueError):
            pixel_size = FALLBACK_PIXEL_SIZE_UM
    return stack, pixel_size


def segment(phase, model):
    return np.asarray(model.eval(phase.astype(np.float32, copy=False), **CELLPOSE_EVAL)[0], dtype=np.int32)


def remove_border_cells(raw):
    raw = raw.astype(np.int32, copy=False)
    border = set(map(int, np.unique(np.concatenate([raw[0], raw[-1], raw[:, 0], raw[:, -1]])))) - {0}
    clean = np.zeros_like(raw)
    new = 1
    for old in sorted(set(map(int, np.unique(raw))) - {0}):
        if old in border:
            continue
        clean[raw == old] = new
        new += 1
    return clean


def measure(mask, rfp_stack, pixel_size):
    projection = rfp_stack.astype(np.float64).mean(axis=0)
    all_cells = mask > 0
    inner = max(1, int(np.ceil(BACKGROUND_RING_INNER_UM / pixel_size)))
    outer = max(inner + 1, int(np.ceil(BACKGROUND_RING_OUTER_UM / pixel_size)))
    inner_struct = np.ones((2 * inner + 1, 2 * inner + 1), dtype=bool)
    outer_struct = np.ones((2 * outer + 1, 2 * outer + 1), dtype=bool)

    rows = []
    for label in np.unique(mask):
        if label == 0:
            continue
        cell = mask == label
        area = int(cell.sum())
        ring = binary_dilation(cell, structure=outer_struct) & ~binary_dilation(cell, structure=inner_struct) & ~all_cells
        n_bg = int(ring.sum())
        raw_integrated = float(projection[cell].sum(dtype=np.float64))
        background = float(np.median(projection[ring])) if n_bg else np.nan
        corrected = raw_integrated - background * area
        rows.append({
            "cell_id": int(label),
            "area_pixels": area,
            "area_um2": area * pixel_size ** 2,
            "rfp_mean_projection_raw_integrated": raw_integrated,
            "background_median": background,
            "background_ring_valid_pixels": n_bg,
            "background_ring_sufficient": n_bg >= MIN_BACKGROUND_PIXELS,
            "rfp_bg_corrected_integrated": corrected,
            "whole_cell_signal_density": corrected / area,
        })
    return pd.DataFrame(rows)


def main():
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load(open(root / "config.yaml"))

    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True, choices=sorted(config["smfish_nd2_dir"]))
    parser.add_argument("--raw-mask-dir", type=Path, default=None)
    args = parser.parse_args()
    args.nd2_dir = Path(config["smfish_nd2_dir"][args.experiment])
    args.manifest = root / "data" / f"smFISH_{args.experiment}_image_manifest.csv"
    args.out_dir = Path(config["out_dir"]) / "smFISH" / args.experiment

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest.run_image.astype(str).str.lower().isin(["true", "1"])]
    if manifest.image_name.duplicated().any():
        raise ValueError("duplicate image_name in manifest")

    mask_dir = args.out_dir / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    model = None
    tables = []

    for row in manifest.itertuples(index=False):
        stack, pixel_size = load_stack(args.nd2_dir / row.nd2_file)
        z = int(row.segmentation_z_index)

        saved = mask_dir / f"{row.image_name}_raw_mask.npy"
        if args.raw_mask_dir is not None:
            raw = np.load(args.raw_mask_dir / row.image_name / "raw_mask.npy", allow_pickle=False)
        elif saved.is_file():
            raw = np.load(saved, allow_pickle=False)
        else:
            if model is None:
                from cellpose import models
                model = models.CellposeModel(gpu=False)
            raw = segment(stack[z, CHANNELS["phase"]], model)

        if raw.shape != stack.shape[2:]:
            raise ValueError(f"{row.image_name}: mask {raw.shape} does not match image {stack.shape[2:]}")
        clean = remove_border_cells(raw)
        np.save(mask_dir / f"{row.image_name}_raw_mask.npy", raw)
        np.save(mask_dir / f"{row.image_name}_clean_mask.npy", clean)

        cells = measure(clean, stack[:, CHANNELS["rfp"]], pixel_size)
        cells.insert(0, "image_name", row.image_name)
        cells.insert(1, "sample_type", row.sample_type)
        cells.insert(2, "timepoint", row.timepoint)
        tables.append(cells)
        print(f"{row.image_name}: {len(np.unique(raw[raw > 0]))} segmented, {len(cells)} kept")

    per_cell = pd.concat(tables, ignore_index=True)
    out = args.out_dir / f"smFISH_{args.experiment}_per_cell.csv"
    per_cell.to_csv(out, index=False, float_format="%.8f")
    print(f"{len(per_cell)} cells from {per_cell.image_name.nunique()} fields written to {out}")


if __name__ == "__main__":
    main()