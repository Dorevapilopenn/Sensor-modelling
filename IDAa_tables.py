from tables_common import generate_tables, save_corr_table


COLUMNS = [
    "K1",
    "K2",
    "K3",
    "K4",
    "K5",
    "K6",
    "MeanHD",
    "deltaD1",
    "deltaK1",
    "deltaHD",
    "deltaD2",
    "deltaK2",
    "Et",
    "Class1",
    "Class2",
]

CONFIG = {
    "columns": COLUMNS,
    "feature_cols": [col for col in COLUMNS if col not in ("Et", "Class1", "Class2")],
    "datasets": {
        "P": "IDAa_res9U_P.csv",
        "N": "IDAa_res9U_N.csv",
        "U": "IDAa_res9U_sorted.csv",
    },
    "out_root": "IDAa_outputs",
    "save_corr_table": save_corr_table,
    "heatmaps": [
        ("{key}_Et_MeanHD_deltaHD.png", "MeanHD", "deltaHD", "Et"),
        ("{key}_Et_deltaD1_K1.png", "deltaD1", "K1", "Et"),
        ("{key}_Et_deltaD2_K4.png", "deltaD2", "K4", "Et"),
        ("{key}_Et_deltaD1_deltaD2.png", "deltaD1", "deltaD2", "Et"),
        ("{key}_Et_deltaK1_deltaD1.png", "deltaK1", "deltaD1", "Et"),
        ("{key}_Et_deltaK2_deltaD2.png", "deltaK2", "deltaD2", "Et"),
        ("{key}_Et_deltaK1_deltaK2.png", "deltaK1", "deltaK2", "Et"),
    ],
    "stratcorrs": [
        ("{key}_corr_deltaK1_vs_Et_K1_deltaD1.png", "K1", "deltaD1", "deltaK1", "Et"),
        ("{key}_corr_deltaK2_vs_Et_K4_deltaD2.png", "K4", "deltaD2", "deltaK2", "Et"),
    ],
    "mean_matrices": [
        {
            "knobs": ["K1", "deltaD1", "deltaK1"],
            "zcol": "Et",
            "filename": "{key}_MeanEt_component1_matrix.png",
        },
        {
            "knobs": ["K4", "deltaD2", "deltaK2"],
            "zcol": "Et",
            "filename": "{key}_MeanEt_component2_matrix.png",
        },
        {
            "knobs": ["MeanHD", "deltaHD", "deltaD1", "deltaD2", "deltaK1", "deltaK2"],
            "zcol": "Et",
            "filename": "{key}_MeanEt_corresponding_knobs_matrix.png",
        },
    ],
}


if __name__ == "__main__":
    generate_tables(CONFIG)
