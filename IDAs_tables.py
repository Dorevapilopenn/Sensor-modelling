from tables_common import generate_tables, save_corr_table


COLUMNS = ["K1", "K2", "K3", "HD", "deltaD", "deltaK", "Et", "Class1", "Class2"]

CONFIG = {
    "columns": COLUMNS,
    "feature_cols": [col for col in COLUMNS if col not in ("Et", "Class1", "Class2")],
    "datasets": {
        "U": "IDAs_res9U.csv",
    },
    "out_root": "IDAs_outputsU",
    "save_corr_table": save_corr_table,
    "heatmaps": [
        ("{key}_Et_K1_deltaD.png", "K1", "deltaD", "Et"),
        ("{key}_Et_K1_deltaK.png", "K1", "deltaK", "Et"),
        ("{key}_Et_deltad_deltak.png", "deltaK", "deltaD", "Et"),
    ],
    "stratcorrs": [
        ("stratcorr_deltaD_deltaK_Et_vs_HD.png", "deltaD", "deltaK", "Et", "HD"),
    ],
    "mean_matrix": {
        "knobs": ["HD", "deltaD", "deltaK"],
        "zcol": "Et",
        "filename": "{key}_MeanEt_knob_matrix.png",
    },
}


if __name__ == "__main__":
    generate_tables(CONFIG)
