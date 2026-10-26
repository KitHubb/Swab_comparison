from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
SOURCE_XLSX = (
    PROJECT_DIR.parent
    / "Paper"
    / "HN00182797_Original_Sample_QC_and_NGS_Library_QC.xlsx"
)
OUTPUT_CSV = SCRIPT_DIR / "input" / "concentration_metadata.csv"

REQUIRED_COLUMNS = [
    "Sample_Name",
    "Sample_Category",
    "Original_QC_Conc_ng_uL",
    "Library_QC_Conc_ng_uL",
    "Original_QC_Result",
    "Library_QC_Result",
    "Library_QC_Comment",
]


def main() -> None:
    data = pd.read_excel(SOURCE_XLSX, sheet_name="Sequencing_Metadata")
    missing = sorted(set(REQUIRED_COLUMNS) - set(data.columns))
    if missing:
        raise ValueError(f"Missing required Excel columns: {missing}")

    output = data.loc[:, REQUIRED_COLUMNS].copy()
    if output["Sample_Name"].isna().any():
        raise ValueError("Sample_Name contains missing values.")
    if output["Sample_Name"].duplicated().any():
        duplicates = output.loc[
            output["Sample_Name"].duplicated(keep=False), "Sample_Name"
        ].tolist()
        raise ValueError(f"Duplicate Sample_Name values: {duplicates}")

    concentration = pd.to_numeric(
        output["Original_QC_Conc_ng_uL"], errors="coerce"
    )
    if concentration.isna().any() or (concentration <= 0).any():
        invalid = output.loc[
            concentration.isna() | (concentration <= 0),
            ["Sample_Name", "Original_QC_Conc_ng_uL"],
        ]
        raise ValueError(
            "Original_QC_Conc_ng_uL must be positive for every sample:\n"
            + invalid.to_string(index=False)
        )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
    print(f"Wrote {len(output)} rows to {OUTPUT_CSV}")
    print(
        "Original_QC_Conc_ng_uL range: "
        f"{concentration.min():.3f} to {concentration.max():.3f} ng/uL"
    )


if __name__ == "__main__":
    main()
