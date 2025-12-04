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

def make_strat_corr(df, xcol, ycol, var_corr, target, outpath, xbins=20, ybins=20):
    xvals = df[xcol].values
    yvals = df[ycol].values
    # quantile bins to ensure coverage
    try:
        x_edges = np.unique(np.quantile(xvals, np.linspace(0,1,xbins+1)))
        y_edges = np.unique(np.quantile(yvals, np.linspace(0,1,ybins+1)))
    except Exception:
        x_edges = np.linspace(np.min(xvals), np.max(xvals), xbins+1)
        y_edges = np.linspace(np.min(yvals), np.max(yvals), ybins+1)
    # also track actual unique data increments
    ux = np.unique(np.sort(xvals))
    uy = np.unique(np.sort(yvals))
    x_centers = 0.5*(x_edges[:-1]+x_edges[1:])
    y_centers = 0.5*(y_edges[:-1]+y_edges[1:])
    H = np.full((len(y_centers), len(x_centers)), np.nan)
    for i in range(len(x_centers)):
        for j in range(len(y_centers)):
            mask = (xvals >= x_edges[i]) & (xvals < x_edges[i+1]) & (yvals >= y_edges[j]) & (yvals < y_edges[j+1])
            if mask.sum() >= 3:
                vals1 = df.loc[mask, var_corr].values
                vals2 = df.loc[mask, target].values
                if np.isnan(vals1).all() or np.isnan(vals2).all():
                    r = np.nan
                else:
                    r = pearsonr(vals1, vals2)[0]
                H[j,i] = r
    fig, ax = plt.subplots(figsize=(6,5))
    pcm = ax.pcolormesh(x_edges, y_edges, np.nan_to_num(H, nan=np.nan), shading='auto', cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xlabel(xcol); ax.set_ylabel(ycol); ax.set_title(f"corr({var_corr},{target}) across {xcol} & {ycol}")
    plt.colorbar(pcm, ax=ax, label=f"corr({var_corr},{target})")

    # overlay grid lines for cell boundaries
    for xe in x_edges:
        ax.axvline(x=xe, color='k', linewidth=0.25, alpha=0.6)
    for ye in y_edges:
        ax.axhline(y=ye, color='k', linewidth=0.25, alpha=0.6)

    # set ticks to match actual increments where reasonable
    if len(ux) <= max(20, len(x_edges)):
        # prefer showing the actual unique data increments if not too many
        # align ticks to data values (centers may equal ux when unique-based bins were used)
        tick_pos_x = ux if len(ux) <= 30 else x_centers
        ax.set_xticks(tick_pos_x)
        ax.set_xticklabels([f"{v:.3g}" for v in tick_pos_x], rotation=45, ha='right')
    else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(6))

    if len(uy) <= max(20, len(y_edges)):
        tick_pos_y = uy if len(uy) <= 30 else y_centers
        ax.set_yticks(tick_pos_y)
        ax.set_yticklabels([f"{v:.3g}" for v in tick_pos_y])
    else:
        ax.yaxis.set_major_locator(plt.MaxNLocator(6))

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])

    plt.savefig(outpath, dpi=300, bbox_inches='tight'); plt.close(fig)
    return outpath

# Load datasets
df_P = pd.read_csv("cons9U_sortedP.csv", header=None)
df_N = pd.read_csv("cons9U_sortedN.csv", header=None)
df_U = pd.read_csv("cons9U_sortedU.csv", header=None)
# Assign column names per confirmation
cols = ["K1","K2","K3","K4","K5","K6", "meanHD",
        "deltaD1","deltaK1","deltaHD","deltaD2","deltaK2",
        "Et","Class1","Class2"]
df_U.columns = cols
df_P.columns = cols 
df_N.columns = cols

datasets = {"P":df_P, "N":df_N, "U":df_U}

out_root = "outputs_30plotsU"
ensure_dir(out_root)

results = {}
for key, df in datasets.items():
    out_dir = os.path.join(out_root, key)
    ensure_dir(out_dir)
    # correlation table
    tbl_path = os.path.join(out_dir, f"{key}_correlation_table.png")
    save_corr_table(df, tbl_path)
    results[f"{key}_table"] = tbl_path
    # 9 plots as specified
    p1 = os.path.join(out_dir, f"{key}_Et_MeanHD_deltaHD.png"); make_heatmap_mean(df, "meanHD","deltaHD","Et", p1)
    p2 = os.path.join(out_dir, f"{key}_Et_deltaD1_deltaD2.png"); make_heatmap_mean(df, "deltaD1","deltaD2","Et", p2)
    p3 = os.path.join(out_dir, f"{key}_Et_deltaK1_deltaK2.png"); make_heatmap_mean(df, "deltaK1","deltaK2","Et", p3)
    p5a = os.path.join(out_dir, f"{key}_Et_deltaD1_K1.png"); make_heatmap_mean(df, "deltaD1","K1","Et", p5a)
    p5b = os.path.join(out_dir, f"{key}_Et_deltaD2_K4.png"); make_heatmap_mean(df, "deltaD2","K4","Et", p5b)
    p6 = os.path.join(out_dir, f"{key}_corr_deltaK1_vs_Et_K1_deltaD1.png"); make_strat_corr(df, "K1","deltaD1","deltaK1","Et", p6)
    p7 = os.path.join(out_dir, f"{key}_corr_deltaK2_vs_Et_K4_deltaD2.png"); make_strat_corr(df, "K4","deltaD2","deltaK2","Et", p7)
    p9a = os.path.join(out_dir, f"{key}_Et_deltaK1_deltaD1.png"); make_heatmap_mean(df, "deltaK1","deltaD1","Et", p9a)
    p9b = os.path.join(out_dir, f"{key}_Et_deltaK2_deltaD2.png"); make_heatmap_mean(df, "deltaK2","deltaD2","Et", p9b)
    results[key] = [tbl_path, p1,p2,p3,p5a,p5b,p6,p7,p9a,p9b]

# Flatten results and print paths
all_files = []
for k,v in results.items():
    if isinstance(v, list):
        all_files.extend(v)
    else:
        all_files.append(v)
all_files = sorted(all_files)
len(all_files), all_files[:5]
