import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def safe_pearsonr(values_a, values_b):
    values_a = np.asarray(values_a, dtype=float)
    values_b = np.asarray(values_b, dtype=float)
    valid = ~(np.isnan(values_a) | np.isnan(values_b))
    values_a = values_a[valid]
    values_b = values_b[valid]

    if values_a.size < 2 or values_b.size < 2:
        return np.nan
    if np.nanstd(values_a) == 0 or np.nanstd(values_b) == 0:
        return np.nan

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return pearsonr(values_a, values_b)[0]
    except Exception:
        return np.nan


def load_csv(path, columns):
    df = pd.read_csv(path, header=None)
    if df.shape[1] != len(columns):
        raise ValueError(f"{path} has {df.shape[1]} columns, expected {len(columns)}")
    df.columns = columns
    return df


def save_corr_table(df, outpath, feature_cols=None, target_cols=("Et", "Class1", "Class2")):
    if feature_cols is None:
        feature_cols = [col for col in df.columns if col not in target_cols]

    corr_table = pd.DataFrame(index=feature_cols)
    for target in target_cols:
        corr_table[target] = [
            safe_pearsonr(df[feature_col], df[target])
            for feature_col in feature_cols
        ]

    height = max(4, 0.32 * len(feature_cols) + 1.0)
    fig, ax = plt.subplots(figsize=(7, height))
    ax.axis("off")
    table = ax.table(
        cellText=np.round(corr_table.values, 3),
        rowLabels=corr_table.index,
        colLabels=corr_table.columns,
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.25)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("k")
        cell.set_linewidth(0.5)
        if row == -1 or col == -1:
            cell.set_facecolor("#f0f0f0")
            cell.set_text_props(weight="bold")

    ax.set_facecolor("white")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def _edges_and_centers(values, max_bins):
    values = np.asarray(values, dtype=float)
    unique_values = np.unique(np.sort(values[~np.isnan(values)]))
    if unique_values.size == 0:
        raise ValueError("Cannot make bins from an empty column")

    if unique_values.size <= max_bins:
        if unique_values.size == 1:
            pad = 1e-6 if unique_values[0] == 0 else abs(unique_values[0]) * 1e-6
            edges = np.array([unique_values[0] - pad, unique_values[0] + pad])
            centers = unique_values
        else:
            diffs = np.diff(unique_values)
            left = unique_values[0] - diffs[0] / 2
            right = unique_values[-1] + diffs[-1] / 2
            mids = (unique_values[:-1] + unique_values[1:]) / 2
            edges = np.concatenate([[left], mids, [right]])
            centers = unique_values
    else:
        edges = np.quantile(values[~np.isnan(values)], np.linspace(0, 1, max_bins + 1))
        edges = np.unique(edges)
        if edges.size < 2:
            edges = np.array([np.nanmin(values), np.nanmax(values) + 1e-6])
        edges[-1] += 1e-12
        centers = 0.5 * (edges[:-1] + edges[1:])

    return edges, centers


def _binned_mean_grid(xvals, yvals, zvals, x_edges, y_edges):
    xvals = np.asarray(xvals, dtype=float)
    yvals = np.asarray(yvals, dtype=float)
    zvals = np.asarray(zvals, dtype=float)
    valid = ~(np.isnan(xvals) | np.isnan(yvals) | np.isnan(zvals))

    xvals = xvals[valid]
    yvals = yvals[valid]
    zvals = zvals[valid]

    xbin_count = len(x_edges) - 1
    ybin_count = len(y_edges) - 1
    sums = np.zeros((ybin_count, xbin_count), dtype=float)
    counts = np.zeros((ybin_count, xbin_count), dtype=int)

    xbins = np.searchsorted(x_edges, xvals, side="right") - 1
    ybins = np.searchsorted(y_edges, yvals, side="right") - 1

    in_range = (
        (xvals >= x_edges[0])
        & (xvals <= x_edges[-1])
        & (yvals >= y_edges[0])
        & (yvals <= y_edges[-1])
    )
    xbins = np.clip(xbins[in_range], 0, xbin_count - 1)
    ybins = np.clip(ybins[in_range], 0, ybin_count - 1)
    zvals = zvals[in_range]

    np.add.at(sums, (ybins, xbins), zvals)
    np.add.at(counts, (ybins, xbins), 1)

    heatmap = np.full((ybin_count, xbin_count), np.nan)
    np.divide(sums, counts, out=heatmap, where=counts > 0)
    return heatmap


def make_heatmap_mean(df, xcol, ycol, zcol, outpath, max_bins=40):
    xvals = df[xcol].to_numpy(dtype=float)
    yvals = df[ycol].to_numpy(dtype=float)
    zvals = df[zcol].to_numpy(dtype=float)

    x_edges, x_centers = _edges_and_centers(xvals, max_bins)
    y_edges, y_centers = _edges_and_centers(yvals, max_bins)
    heatmap = _binned_mean_grid(xvals, yvals, zvals, x_edges, y_edges)

    fig, ax = plt.subplots(figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, heatmap, shading="auto", cmap="viridis")
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(f"Mean {zcol} vs {xcol} & {ycol}")
    fig.colorbar(mesh, ax=ax, label=zcol)

    for edge in x_edges:
        ax.axvline(edge, color="k", linewidth=0.25, alpha=0.6)
    for edge in y_edges:
        ax.axhline(edge, color="k", linewidth=0.25, alpha=0.6)

    if len(x_centers) <= 20:
        ax.set_xticks(x_centers)
        ax.set_xticklabels([f"{value:.3g}" for value in x_centers], rotation=45, ha="right")
    else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    if len(y_centers) <= 20:
        ax.set_yticks(y_centers)
        ax.set_yticklabels([f"{value:.3g}" for value in y_centers])
    else:
        ax.yaxis.set_major_locator(plt.MaxNLocator(6))

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def make_strat_corr(df, xcol, ycol, var_corr, target, outpath, xbins=20, ybins=20, min_count=2):
    xvals = df[xcol].to_numpy(dtype=float)
    yvals = df[ycol].to_numpy(dtype=float)
    x_edges, x_centers = _edges_and_centers(xvals, xbins)
    y_edges, y_centers = _edges_and_centers(yvals, ybins)
    heatmap = np.full((len(y_centers), len(x_centers)), np.nan)

    for i in range(len(x_centers)):
        xmask = (xvals >= x_edges[i]) & (xvals < x_edges[i + 1])
        if i == len(x_centers) - 1:
            xmask = (xvals >= x_edges[i]) & (xvals <= x_edges[i + 1])
        for j in range(len(y_centers)):
            ymask = (yvals >= y_edges[j]) & (yvals < y_edges[j + 1])
            if j == len(y_centers) - 1:
                ymask = (yvals >= y_edges[j]) & (yvals <= y_edges[j + 1])
            mask = xmask & ymask
            if mask.sum() >= min_count:
                heatmap[j, i] = safe_pearsonr(df.loc[mask, var_corr], df.loc[mask, target])

    fig, ax = plt.subplots(figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, heatmap, shading="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(f"corr({var_corr},{target}) across {xcol} & {ycol}")
    fig.colorbar(mesh, ax=ax, label=f"corr({var_corr},{target})")

    for edge in x_edges:
        ax.axvline(edge, color="k", linewidth=0.25, alpha=0.5)
    for edge in y_edges:
        ax.axhline(edge, color="k", linewidth=0.25, alpha=0.5)

    ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def make_mean_scatter_matrix(df, knobs, zcol, outpath, max_bins=40):
    knob_count = len(knobs)
    if knob_count < 2:
        raise ValueError("Mean scatter matrix needs at least two knobs")

    fig_size = max(7, 1.45 * knob_count)
    fig, axes = plt.subplots(knob_count, knob_count, figsize=(fig_size, fig_size), squeeze=False)
    mesh_for_colorbar = None
    panel_grids = {}
    panel_values = []

    ranges = {}
    for knob in knobs:
        values = df[knob].to_numpy(dtype=float)
        valid_values = values[~np.isnan(values)]
        lower = valid_values.min()
        upper = valid_values.max()
        if lower == upper:
            pad = 1e-6 if lower == 0 else abs(lower) * 1e-6
            lower -= pad
            upper += pad
        ranges[knob] = (lower, upper)

    for ycol in knobs:
        for xcol in knobs:
            if xcol == ycol:
                continue
            xvals = df[xcol].to_numpy(dtype=float)
            yvals = df[ycol].to_numpy(dtype=float)
            zvals = df[zcol].to_numpy(dtype=float)
            x_edges, _ = _edges_and_centers(xvals, max_bins)
            y_edges, _ = _edges_and_centers(yvals, max_bins)
            heatmap = _binned_mean_grid(xvals, yvals, zvals, x_edges, y_edges)
            panel_grids[(xcol, ycol)] = (x_edges, y_edges, heatmap)
            finite_values = heatmap[np.isfinite(heatmap)]
            if finite_values.size:
                panel_values.append(finite_values)

    if not panel_values:
        raise ValueError(f"No finite binned {zcol} values found for {outpath}")

    finite_panel_values = np.concatenate(panel_values)
    zmin = np.nanmin(finite_panel_values)
    zmax = np.nanmax(finite_panel_values)
    if zmin == zmax:
        pad = 1e-6 if zmin == 0 else abs(zmin) * 1e-6
        zmin -= pad
        zmax += pad

    for row, ycol in enumerate(knobs):
        for col, xcol in enumerate(knobs):
            ax = axes[row, col]
            if row == col:
                ax.text(0.5, 0.5, xcol, ha="center", va="center", fontsize=9, transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                x_edges, y_edges, heatmap = panel_grids[(xcol, ycol)]
                mesh_for_colorbar = ax.pcolormesh(
                    x_edges,
                    y_edges,
                    heatmap,
                    shading="auto",
                    cmap="viridis",
                    vmin=zmin,
                    vmax=zmax,
                )
                ax.set_xlim(*ranges[xcol])
                ax.set_ylim(*ranges[ycol])
                for edge in x_edges:
                    ax.axvline(edge, color="k", linewidth=0.12, alpha=0.35)
                for edge in y_edges:
                    ax.axhline(edge, color="k", linewidth=0.12, alpha=0.35)

            if row == knob_count - 1:
                ax.set_xlabel(xcol, fontsize=8)
                ax.tick_params(axis="x", labelbottom=True, bottom=True)
            elif row == 0:
                ax.xaxis.set_ticks_position("top")
                ax.xaxis.set_label_position("top")
                ax.tick_params(axis="x", labeltop=True, top=True, labelbottom=False, bottom=False)
            else:
                ax.set_xticklabels([])
                ax.tick_params(axis="x", labelbottom=False, bottom=False)

            if col == 0:
                ax.set_ylabel(ycol, fontsize=8)
                ax.tick_params(axis="y", labelleft=True, left=True)
            elif col == knob_count - 1:
                ax.yaxis.set_ticks_position("right")
                ax.yaxis.set_label_position("right")
                ax.tick_params(axis="y", labelright=True, right=True, labelleft=False, left=False)
            else:
                ax.set_yticklabels([])
                ax.tick_params(axis="y", labelleft=False, left=False)
            ax.tick_params(labelsize=6, length=2)

    fig.suptitle(f"Mean {zcol} stratified by knob pairs", fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.84, 0.96))
    if mesh_for_colorbar is not None:
        cbar_ax = fig.add_axes([0.90, 0.12, 0.018, 0.76])
        cbar = fig.colorbar(mesh_for_colorbar, cax=cbar_ax)
        cbar.set_label(f"Mean {zcol}")

    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def generate_tables(config):
    ensure_dir(config["out_root"])
    written = []
    corr_table_func = config.get("save_corr_table", save_corr_table)

    for key, csv_path in config["datasets"].items():
        df = load_csv(csv_path, config["columns"])
        out_dir = os.path.join(config["out_root"], key)
        ensure_dir(out_dir)

        table_path = os.path.join(out_dir, f"{key}_correlation_table.png")
        written.append(corr_table_func(df, table_path, config.get("feature_cols")))

        for filename, xcol, ycol, zcol in config.get("heatmaps", []):
            written.append(make_heatmap_mean(df, xcol, ycol, zcol, os.path.join(out_dir, filename.format(key=key))))

        for filename, xcol, ycol, var_corr, target in config.get("stratcorrs", []):
            written.append(make_strat_corr(df, xcol, ycol, var_corr, target, os.path.join(out_dir, filename.format(key=key))))

        matrix_configs = config.get("mean_matrices")
        if matrix_configs is None and config.get("mean_matrix"):
            matrix_configs = [config["mean_matrix"]]

        for matrix_config in matrix_configs or []:
            written.append(
                make_mean_scatter_matrix(
                    df,
                    matrix_config["knobs"],
                    matrix_config.get("zcol", "Et"),
                    os.path.join(out_dir, matrix_config.get("filename", "{key}_MeanEt_knob_matrix.png").format(key=key)),
                    max_bins=matrix_config.get("max_bins", 40),
                )
            )

    print(f"Wrote {len(written)} files:")
    for path in written:
        print(path)
    return written
