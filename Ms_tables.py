from tables_common import generate_tables, save_corr_table


COLUMNS = [
    "K1",
    "K2",
    "K3",
    "K4",
    "K5",
    "HD",
    "deltaD",
    "deltaK",
    "deltaC",
    "deltaG",
    "Et",
    "Class1",
    "Class2",
]

CONFIG = {
    "columns": COLUMNS,
    "feature_cols": [col for col in COLUMNS if col not in ("Et", "Class1", "Class2")],
    "datasets": {
        "U": "Ms_res9U.csv",
    },
    "out_root": "Ms_outputsU",
    "save_corr_table": save_corr_table,
    "heatmaps": [
        ("{key}_Et_K1_deltaD.png", "K1", "deltaD", "Et"),
        ("{key}_Et_K1_deltaK.png", "K1", "deltaK", "Et"),
        ("{key}_Et_K1_deltaC.png", "K1", "deltaC", "Et"),
        ("{key}_Et_K1_deltaG.png", "K1", "deltaG", "Et"),
        ("{key}_Et_deltaD_deltaK.png", "deltaD", "deltaK", "Et"),
        ("{key}_Et_deltaC_deltaG.png", "deltaC", "deltaG", "Et"),
        ("{key}_Et_deltaD_deltaC.png", "deltaD", "deltaC", "Et"),
        ("{key}_Et_deltaK_deltaG.png", "deltaK", "deltaG", "Et"),
        ("{key}_Et_deltaD_deltaG.png", "deltaD", "deltaG", "Et"),
        ("{key}_Et_deltaK_deltaC.png", "deltaK", "deltaC", "Et"),
    ],
    "stratcorrs": [
        ("stratcorr_deltaD_deltaK_Et_vs_deltaG.png", "deltaD", "deltaK", "Et", "deltaG"),
        ("stratcorr_deltaC_deltaG_Et_vs_deltaK.png", "deltaC", "deltaG", "Et", "deltaK"),
        ("stratcorr_deltaC_deltaG_Et_vs_deltaD.png", "deltaC", "deltaG", "Et", "deltaD"),
        ("stratcorr_deltaD_deltaC_Et_vs_deltaG.png", "deltaD", "deltaC", "Et", "deltaG"),
        ("stratcorr_deltaG_deltaD_Et_vs_deltaC.png", "deltaG", "deltaD", "Et", "deltaC"),
        ("stratcorr_deltaC_deltaG_Et_vs_K1.png", "deltaC", "deltaG", "Et", "K1"),
        ("stratcorr_deltaK_deltaG_Et_vs_K1.png", "deltaK", "deltaG", "Et", "K1"),
        ("stratcorr_deltaD_deltaC_Et_vs_deltaK.png", "deltaD", "deltaC", "Et", "deltaK"),
        ("stratcorr_deltaD_deltaC_Et_vs_K1.png", "deltaD", "deltaC", "Et", "K1"),
    ],
    "mean_matrix": {
        "knobs": ["HD", "deltaD", "deltaK", "deltaC", "deltaG"],
        "zcol": "Et",
        "filename": "{key}_MeanEt_knob_matrix.png",
    },
}


if __name__ == "__main__":
    generate_tables(CONFIG)
