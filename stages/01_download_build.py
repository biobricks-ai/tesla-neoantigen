#!/usr/bin/env python3
"""
Download and build TESLA (Tumor neoantigen SeLection Alliance) neoantigen
immunogenicity benchmark dataset.

Sources:
  - Table S4 (mmc4.xlsx): 608 neopeptides from Phase 1 (8 patients) with
    pMHC-multimer-validated immunogenicity outcome.
  - Table S7 (mmc7.xlsx): 310 neopeptides from Phase 2 (3 melanoma patients).

Both from:
  Wells et al., Cell 2020, "Key Parameters of Tumor Epitope Immunogenicity
  Revealed Through a Consortium Approach Improve Neoantigen Prediction"
  DOI: 10.1016/j.cell.2020.09.015
  License: CC-BY-NC-ND-4.0
"""

from pathlib import Path
import pandas as pd
import requests

PHASE1_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420311569-mmc4.xlsx"
PHASE2_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420311569-mmc7.xlsx"

DOWNLOAD_DIR = Path("download")
BRICK_DIR = Path("brick")


def download_file(url: str, dest: Path) -> Path:
    """Download a file from URL to dest path."""
    print(f"  Downloading {url} ...")
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    dest.write_bytes(response.content)
    print(f"  Saved {len(response.content):,} bytes to {dest}")
    return dest


def load_phase1(path: Path) -> pd.DataFrame:
    """
    Load TESLA Phase 1 peptide validation table (Table S4, mmc4.xlsx).
    Sheet: master-bindings-selected, 608 rows x 21 cols.
    Key columns: PMHC, PATIENT_ID, TISSUE_TYPE, MHC (HLA allele), ALT_EPI_SEQ
    (peptide), VALIDATED (immunogenicity outcome), plus prediction/quality features.
    """
    df = pd.read_excel(path, sheet_name="master-bindings-selected", engine="openpyxl")
    # Rename to canonical names
    df = df.rename(columns={
        "ALT_EPI_SEQ": "peptide",
        "MHC": "hla_allele",
        "VALIDATED": "immunogenic",
        "PATIENT_ID": "patient_id",
        "TISSUE_TYPE": "tissue_type",
        "PEP_LEN": "peptide_length",
        "MEASURED_BINDING_AFFINITY": "measured_binding_affinity_nM",
        "NETMHC_PAN_BINDING_AFFINITY": "netmhcpan_binding_affinity_nM",
        "TUMOR_ABUNDANCE": "tumor_abundance_tpm",
        "BINDING_STABILITY": "binding_stability_half_life_h",
        "FRAC_HYDROPHOBIC": "frac_hydrophobic",
        "AGRETOPICITY": "agretopicity",
        "FOREIGNNESS": "foreignness",
        "MUTATION_POSITION": "mutation_position",
        "NUMBER_PREDICTING": "num_teams_predicting",
        "TCR_FLOW_I": "tcr_flow_i",
        "TCR_FLOW_I_QUANT": "tcr_flow_i_quant",
        "TCR_NANOPARTICLE": "tcr_nanoparticle",
        "TCR_FLOW_II": "tcr_flow_ii",
        "TCR_FLOW_II_QUANT": "tcr_flow_ii_quant",
    })
    # Convert immunogenic to bool
    df["immunogenic"] = df["immunogenic"].astype(bool)
    df["phase"] = 1
    return df


def load_phase2(path: Path) -> pd.DataFrame:
    """
    Load TESLA Phase 2 peptide validation table (Table S7, mmc7.xlsx).
    Sheet: SM_Table_S7_PEPTIDE_VALIDATION_, 310 rows x 15 cols.
    Key columns: PMHC, PATIENT_ID, TISSUE_TYPE, ALT_EPI_SEQ (peptide),
    VALIDATED (immunogenicity), plus prediction/quality features.
    Note: Phase 2 uses a different pipeline (pMHCflurry-predicted affinities)
    and only TIL tissue from 3 additional melanoma patients.
    """
    df = pd.read_excel(path, sheet_name="SM_Table_S7_PEPTIDE_VALIDATION_", engine="openpyxl")
    df = df.rename(columns={
        "ALT_EPI_SEQ": "peptide",
        "VALIDATED": "immunogenic",
        "PATIENT_ID": "patient_id",
        "TISSUE_TYPE": "tissue_type",
        "PEP_LEN": "peptide_length",
        "PREDICTED_BINDING_AFFINITY": "pMHCflurry_binding_affinity_nM",
        "NETMHC_BINDING_AFFINITY": "netmhc_binding_affinity_nM",
        "TUMOR_ABUNDANCE": "tumor_abundance_tpm",
        "BINDING_STABILITY": "binding_stability_half_life_h",
        "AGRETOPICITY": "agretopicity",
        "FOREIGNNESS": "foreignness",
        "MUTATION_POSITION": "mutation_position",
        "TCR_FLOW_II": "tcr_flow_ii",
        "TCR_FLOW_II_QUANT": "tcr_flow_ii_quant",
    })
    # Extract HLA allele from PMHC field (format: "A*01:01_PEPTIDE")
    df["hla_allele"] = df["PMHC"].str.extract(r'^([^_]+)')
    df["immunogenic"] = df["immunogenic"].astype(bool)
    df["phase"] = 2
    return df


def main():
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    BRICK_DIR.mkdir(exist_ok=True)

    print("=" * 60)
    print("TESLA Neoantigen Immunogenicity Benchmark")
    print("Wells et al., Cell 2020, doi:10.1016/j.cell.2020.09.015")
    print("=" * 60)

    # Download Phase 1 (Table S4)
    phase1_path = DOWNLOAD_DIR / "mmc4_phase1_peptides.xlsx"
    download_file(PHASE1_URL, phase1_path)

    # Download Phase 2 (Table S7)
    phase2_path = DOWNLOAD_DIR / "mmc7_phase2_peptides.xlsx"
    download_file(PHASE2_URL, phase2_path)

    # Load and process Phase 1
    print("\nProcessing Phase 1 (Table S4) ...")
    df1 = load_phase1(phase1_path)
    print(f"  Phase 1: {len(df1)} peptides, {df1['immunogenic'].sum()} immunogenic")

    # Load and process Phase 2
    print("\nProcessing Phase 2 (Table S7) ...")
    df2 = load_phase2(phase2_path)
    print(f"  Phase 2: {len(df2)} peptides, {df2['immunogenic'].sum()} immunogenic")

    # Save Phase 1 parquet
    out1 = BRICK_DIR / "tesla_phase1_peptides.parquet"
    df1.to_parquet(out1, index=False)
    print(f"\nSaved Phase 1: {out1} ({len(df1)} rows, {len(df1.columns)} cols)")

    # Save Phase 2 parquet
    out2 = BRICK_DIR / "tesla_phase2_peptides.parquet"
    df2.to_parquet(out2, index=False)
    print(f"Saved Phase 2: {out2} ({len(df2)} rows, {len(df2.columns)} cols)")

    # Save combined parquet with common columns
    common_cols = ["PMHC", "phase", "patient_id", "tissue_type", "hla_allele",
                   "peptide", "peptide_length", "tumor_abundance_tpm",
                   "binding_stability_half_life_h", "agretopicity", "foreignness",
                   "mutation_position", "immunogenic"]

    df_all = pd.concat([
        df1[[c for c in common_cols if c in df1.columns]],
        df2[[c for c in common_cols if c in df2.columns]],
    ], ignore_index=True)

    out_all = BRICK_DIR / "tesla_neoantigen.parquet"
    df_all.to_parquet(out_all, index=False)
    total_immunogenic = df_all["immunogenic"].sum()
    print(f"Saved combined: {out_all} ({len(df_all)} rows, {len(df_all.columns)} cols)")
    print(f"  Immunogenic: {total_immunogenic} / {len(df_all)} "
          f"({100 * total_immunogenic / len(df_all):.1f}%)")

    # Summary
    print("\n" + "=" * 60)
    print("Output files:")
    print("=" * 60)
    for f in sorted(BRICK_DIR.glob("*.parquet")):
        df = pd.read_parquet(f)
        print(f"  {f.name}: {len(df):,} rows x {len(df.columns)} cols")
        print(f"    columns: {list(df.columns)}")
        if "immunogenic" in df.columns:
            print(f"    immunogenic: {df['immunogenic'].sum()} / {len(df)}")
        print(f"    sample:\n{df.head(2).to_string()}")
        print()


if __name__ == "__main__":
    main()
