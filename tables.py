# Generating requested PNGs for three datasets: P (nonnegative deltaK2), N (nonpositive deltaK2), U (union)
# Files: correlation table + 9 plots per dataset -> total 30 images (3 tables + 27 plots) + 3 stratified corr plots = 30? 
# User asked total of 30: for each dataset: 1 table + 9 plots = 10 -> 3*10 = 30. We'll produce exactly that.
import pandas as pd, numpy as np, matplotlib.pyplot as plt, os
from scipy.stats import pearsonr

# Helper functions
def ensure_dir(d): 
    os.makedirs(d, exist_ok=True)

def save_corr_table(df, outpath):
    corr_table = df.iloc[:,0:11].corrwith(df["Et"]).to_frame("Et")
    corr_table["Class1"] = df.iloc[:,0:11].corrwith(df["Class1"])
    corr_table["Class2"] = df.iloc[:,0:11].corrwith(df["Class2"])
    fig, ax = plt.subplots(figsize=(6,4))
    ax.axis('off')
    tbl = ax.table(cellText=np.round(corr_table.values,3),
                   rowLabels=corr_table.index, colLabels=corr_table.columns, loc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1,1.4)

    # add grid lines: set edgecolor and linewidth for all cells
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor('k')
        cell.set_linewidth(0.5)
        # emphasize header cells (column headers at row index -1, row labels at col index -1)
        if r == -1 or c == -1:
            cell.set_facecolor('#f0f0f0')
            cell.set_text_props(weight='bold')

    # ensure grid is visible on white backgrounds
    ax.set_facecolor('white')
    plt.savefig(outpath, dpi=300, bbox_inches='tight'); plt.close(fig)
    return outpath

def make_heatmap_mean(df, xcol, ycol, zcol, outpath, max_bins=40):
    xvals = df[xcol].values
    yvals = df[ycol].values
    zvals = df[zcol].values
    # choose bin edges: if unique counts small, use unique values; else quantile bins
    ux = np.unique(np.sort(xvals))
    uy = np.unique(np.sort(yvals))
    if len(ux) <= max_bins:
        x_edges = np.concatenate([ux, [ux[-1]+1e-6]])
        x_centers = ux
    else:
        x_edges = np.quantile(xvals, np.linspace(0,1,max_bins+1))
        x_centers = 0.5*(x_edges[:-1]+x_edges[1:])
    if len(uy) <= max_bins:
        y_edges = np.concatenate([uy, [uy[-1]+1e-6]])
        y_centers = uy
    else:
        y_edges = np.quantile(yvals, np.linspace(0,1,max_bins+1))
        y_centers = 0.5*(y_edges[:-1]+y_edges[1:])
    H = np.full((len(y_centers), len(x_centers)), np.nan)
    for i in range(len(x_centers)):
        for j in range(len(y_centers)):
            mask = (xvals >= x_edges[i]) & (xvals < x_edges[i+1]) & (yvals >= y_edges[j]) & (yvals < y_edges[j+1])
            if mask.sum() > 0:
                H[j,i] = np.nanmean(zvals[mask])
    # plot using pcolormesh then overlay grid lines
    fig, ax = plt.subplots(figsize=(6,5))
    pcm = ax.pcolormesh(x_edges, y_edges, np.nan_to_num(H, nan=np.nan), shading='auto', cmap='viridis')
    ax.set_xlabel(xcol); ax.set_ylabel(ycol); ax.set_title(f"Mean {zcol} vs {xcol} & {ycol}")
    plt.colorbar(pcm, ax=ax, label=zcol)

    # draw thin grid lines at bin edges for readability
    for xe in x_edges:
        ax.axvline(x=xe, color='k', linewidth=0.25, alpha=0.6)
    for ye in y_edges:
        ax.axhline(y=ye, color='k', linewidth=0.25, alpha=0.6)

    # set ticks to match actual data increments (use centers when many unique values)
    if len(x_centers) <= 20:
        ax.set_xticks(x_centers)
        ax.set_xticklabels([f"{v:.3g}" for v in x_centers], rotation=45, ha='right')
    else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    if len(y_centers) <= 20:
        ax.set_yticks(y_centers)
        ax.set_yticklabels([f"{v:.3g}" for v in y_centers])
    else:
        ax.yaxis.set_major_locator(plt.MaxNLocator(6))

    # tighten axis to match edges
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])

    plt.savefig(outpath, dpi=300, bbox_inches='tight'); plt.close(fig)
    return outpath

def make_strat_corr(df, xcol, ycol, var_corr, target, outpath=None, xbins=20, ybins=20, min_count=2):
    """Compute stratified Pearson correlation grid and optionally save a heatmap to outpath.
    Robust to constant inputs, small bins, and includes last-edge values."""
    xvals = df[xcol].values
    yvals = df[ycol].values

    # bin edges (quantile preferred, fallback to linspace)
    try:
        x_edges = np.quantile(xvals, np.linspace(0,1,xbins+1))
        y_edges = np.quantile(yvals, np.linspace(0,1,ybins+1))
    except Exception:
        x_edges = np.linspace(np.nanmin(xvals), np.nanmax(xvals), xbins+1)
        y_edges = np.linspace(np.nanmin(yvals), np.nanmax(yvals), ybins+1)

    # ensure monotonic unique edges and guarantee at least one bin
    x_edges = np.unique(x_edges)
    y_edges = np.unique(y_edges)
    if x_edges.size < 2:
        x_edges = np.array([np.nanmin(xvals), np.nanmax(xvals) + 1e-6])
    if y_edges.size < 2:
        y_edges = np.array([np.nanmin(yvals), np.nanmax(yvals) + 1e-6])

    # tiny pad to include max value in last bin
    x_edges[-1] = x_edges[-1] + 1e-12
    y_edges[-1] = y_edges[-1] + 1e-12

    x_centers = 0.5*(x_edges[:-1] + x_edges[1:])
    y_centers = 0.5*(y_edges[:-1] + y_edges[1:])
    H = np.full((len(y_centers), len(x_centers)), np.nan)

    for i in range(len(x_centers)):
        for j in range(len(y_centers)):
            # inclusive on last edge to avoid dropping max values
            if i == len(x_centers)-1:
                xmask = (xvals >= x_edges[i]) & (xvals <= x_edges[i+1])
            else:
                xmask = (xvals >= x_edges[i]) & (xvals < x_edges[i+1])
            if j == len(y_centers)-1:
                ymask = (yvals >= y_edges[j]) & (yvals <= y_edges[j+1])
            else:
                ymask = (yvals >= y_edges[j]) & (yvals < y_edges[j+1])

            mask = xmask & ymask
            if mask.sum() >= min_count:
                vals1 = df.loc[mask, var_corr].values
                vals2 = df.loc[mask, target].values
                # require at least two non-NaN and non-constant inputs
                if np.count_nonzero(~np.isnan(vals1)) < 2 or np.count_nonzero(~np.isnan(vals2)) < 2:
                    r = np.nan
                elif np.nanstd(vals1) == 0 or np.nanstd(vals2) == 0:
                    r = np.nan
                else:
                    try:
                        r = pearsonr(vals1, vals2)[0]
                    except Exception:
                        r = np.nan
                H[j,i] = r

    # Plot heatmap if outpath provided
    if outpath is not None:
        fig, ax = plt.subplots(figsize=(6,5))
        pcm = ax.pcolormesh(x_edges, y_edges, H, shading='auto', cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_xlabel(xcol); ax.set_ylabel(ycol)
        ax.set_title(f"corr({var_corr},{target}) across {xcol} & {ycol}")
        plt.colorbar(pcm, ax=ax, label=f"corr({var_corr},{target})")

        # overlay thin grid lines
        for xe in x_edges:
            ax.axvline(x=xe, color='k', linewidth=0.25, alpha=0.5)
        for ye in y_edges:
            ax.axhline(y=ye, color='k', linewidth=0.25, alpha=0.5)

        # reasonable tick selection
        if x_centers.size <= 20:
            ax.set_xticks(x_centers); ax.set_xticklabels([f"{v:.3g}" for v in x_centers], rotation=45, ha='right')
        else:
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
        if y_centers.size <= 20:
            ax.set_yticks(y_centers); ax.set_yticklabels([f"{v:.3g}" for v in y_centers])
        else:
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

        ax.set_xlim(x_edges[0], x_edges[-1])
        ax.set_ylim(y_edges[0], y_edges[-1])
        plt.savefig(outpath, dpi=300, bbox_inches='tight'); plt.close(fig)

    return H

def strat_corr_to_file(df, xcol, ycol, var_corr, target, outpath, xbins=20, ybins=20, min_count=2):
    """Non-interactive wrapper: compute stratified correlation and save heatmap to outpath."""
    try:
        H = make_strat_corr(df, xcol, ycol, var_corr, target, outpath=outpath, xbins=xbins, ybins=ybins, min_count=min_count)
        return outpath
    except Exception as e:
        print("Failed to generate stratified correlation:", e)
        return None

# Load datasets
df_U = pd.read_csv("Ms_res9U_U.csv", header=None)
df_N = pd.read_csv("Ms_res9U_N.csv", header=None)
df_P = pd.read_csv("Ms_res9U_P.csv", header=None)
# Assign column names per confirmation
cols = ["K1","K2","K3", "K4", "K5", "HD", "deltaD","deltaK","deltaC", "deltaG",
        "Et","Class1","Class2"]
df_U.columns = cols
df_N.columns = cols
df_P.columns = cols

datasets = {"U":df_U, "N":df_N, "P":df_P}

out_root = "Ms_outputsU"
ensure_dir(out_root)

results = {}
for key, df in datasets.items():
    out_dir = os.path.join(out_root, key)
    ensure_dir(out_dir)
    # Correlation table
    table_path = os.path.join(out_dir, f"{key}_correlation_table.png")
    results[key] = table_path

    # 9 plots as specified
    p1 = os.path.join(out_dir, f"{key}_Et_K1_deltaD.png"); make_heatmap_mean(df, "K1","deltaD","Et", p1)
    p2 = os.path.join(out_dir, f"{key}_Et_K1_deltaK.png"); make_heatmap_mean(df, "K1","deltaK","Et", p2)
    p3 = os.path.join(out_dir, f"{key}_Et_K1_deltaC.png"); make_heatmap_mean(df, "K1","deltaC","Et", p3)
    p4 = os.path.join(out_dir, f"{key}_Et_K1_deltaG.png"); make_heatmap_mean(df, "K1","deltaG","Et", p4)
    p5 = os.path.join(out_dir, f"{key}_Et_deltaD_deltaK.png"); make_heatmap_mean(df, "deltaD","deltaK","Et", p5)
    p6 = os.path.join(out_dir, f"{key}_Et_deltaC_deltaG.png"); make_heatmap_mean(df, "deltaC","deltaG","Et", p6)
    p7 = os.path.join(out_dir, f"{key}_Et_deltaD_deltaC.png"); make_heatmap_mean(df, "deltaD","deltaC","Et", p7)
    p8 = os.path.join(out_dir, f"{key}_Et_deltaK_deltaG.png"); make_heatmap_mean(df, "deltaK","deltaG","Et", p8)
    results[key] = [table_path, p1, p2, p3, p4, p5, p6, p7, p8]

    # Non-interactive stratified correlation: deltaD & deltaK, corr(Et, deltaD)
    strat_corr_to_file(df, "deltaD", "deltaK", "Et", "deltaG", os.path.join(out_dir, f"stratcorr_deltaD_deltaK_Et_vs_deltaG.png"))
    strat_corr_to_file(df, "deltaC", "deltaG", "Et", "deltaK", os.path.join(out_dir, f"stratcorr_deltaC_deltaG_Et_vs_deltaK.png"))
    strat_corr_to_file(df, "deltaK", "deltaG", "Et", "K1", os.path.join(out_dir, f"stratcorr_deltaK_deltaG_Et_vs_K1.png"))

# Flatten results and print paths
all_files = []
for k,v in results.items():
    if isinstance(v, list):
        all_files.extend(v)
    else:
        all_files.append(v)
all_files = sorted(all_files)
len(all_files), all_files[:5]
