# RNA-Sequence-Analysis

Transcriptomic and structural investigation of Parkinson's disease: identifying differentially expressed genes between control and disease samples, characterizing the associated protein classes, and following up on a candidate protein with structural analysis.

## Overview

Parkinson's disease (PD) is a neurodegenerative disorder marked by progressive loss of dopaminergic neurons in the substantia nigra. This project investigates the disorder from a transcriptomics angle: starting from public RNA-seq count data, it identifies genes that are differentially expressed between control and PD samples, characterizes what kind of proteins those genes encode, and uses that to motivate a structural follow-up on one candidate protein.

The full pipeline:
1. **Data acquisition** — control and PD sample count data from a public GEO dataset.
2. **Normalization** — per-sample z-score normalization of raw counts.
3. **Differential expression testing** — a per-gene t-test (control vs. disease) to call statistically significant differentially expressed (DE) genes.
4. **Visualization** — a volcano plot (fold change vs. significance) and a hierarchically clustered heatmap of the DE genes across samples.
5. **Functional characterization** — the most significantly changed genes are mapped to protein families/classes (via DAVID and PANTHER) to look for patterns in what kind of proteins are affected.
6. **Candidate hypothesis** — the protein class with the most disproportionate fold change is flagged as a candidate of interest in PD pathogenesis.
7. **Structural analysis** — a representative protein from that class is selected for structural characterization: fold, active site, and ligand-binding interactions, visualized from a solved crystal structure.
8. **Interaction network analysis** — the candidate protein's interaction network is examined to identify other proteins that may be functionally linked to it.

**Dataset:** [GEO accession GSE206308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206308) — "Genes critical for development and differentiation of dopaminergic neurons are downregulated in Parkinson's disease." The dataset profiles substantia nigra pars compacta transcriptomes for both a mouse MPTP-treatment model and human PD patients; this project analyzes only the human arm (3 control samples, 4 PD patient samples).

## Problem it addresses

The genes and pathways that drive dopaminergic neuron loss in PD are still not fully understood. Differential expression analysis of disease-vs-control transcriptomic data is a standard first step toward generating candidate genes and proteins that may be mechanistically involved — this project works through that process end to end, from raw counts to a structural hypothesis about one candidate protein.

## Repository structure

```
RNA-Sequence-Analysis/
├── README.md
├── requirements.txt          # Python dependencies
├── src/
│   └── rna_seq_ana.py        # Analysis script
├── data/
│   ├── README.md   # Expected input file layout 
│   ├── control/    # Control sample files (not included)
│   └── disease/    # Disease sample files (not included)
├── results/
│   └── genes_diff.csv        # Names of the 714 differentially expressed genes identified
└── docs/
    └── report.pdf            # Full write-up with figures, protein classification, and structural analysis
```

## Setup / installation

1. Clone the repository and create a Python environment (Python 3.9+ recommended).
2. Install dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Download the GSE206308 sample files from GEO and place them under `data/control/` and `data/disease/` as described in [`data/README.md`](data/README.md).

## How to run

From the project root:

```
python src/rna_seq_ana.py
```

The script reads the sample files from `data/control/` and `data/disease/`, and writes `genes_diff.csv`, `pvals_diff.csv`, and `log2fc_diff.csv` to the current working directory. The volcano plot and clustered heatmap are displayed as they are generated.

## Example usage

The core differential-expression logic, condensed:

```python
# Per-sample z-score normalization
data_norm = data_cd.iloc[:, 1:].apply(lambda x: (x - np.mean(x)) / np.std(x), axis=0)

# Per-gene t-test, control (first 3 columns) vs. disease (remaining columns)
diff_expr_genes = []
for i in data_norm.index:
    t, pval = stats.ttest_ind(data_norm.iloc[i, :3], data_norm.iloc[i, 3:])
    if pval < 0.05 and abs(t) > 1:
        diff_expr_genes.append(i)
```

## Results

- **Total genes analyzed:** 63,677
- **Differentially expressed genes (p < 0.05, |t| > 1):** 714 — 239 upregulated, 475 downregulated
- **Top upregulated gene:** RAD52 (RAD52 homolog, DNA repair protein), log2FC ≈ 36.65, p ≈ 0.00099
- **Top downregulated gene:** SLC16A4-AS1, log2FC ≈ -0.976, p ≈ 0.00463
- **Full list of DE gene names:** [`results/genes_diff.csv`](results/genes_diff.csv)

### Protein classification of the top differentially expressed genes

The 8 genes with the largest |log2FC| were mapped to PANTHER protein classes:

| Protein class | Genes |
|---|---|
| DNA metabolism protein | RAD52 |
| Transporter | SLC7A2, ABCC8 |
| Metabolite interconversion enzyme | HCCS, PNPLA4 |
| Defense/immunity protein | CEACAM21 |
| Unclassified | NDUFAF7, ST7 |

Although the transporter and metabolite-interconversion-enzyme classes had more genes represented, RAD52 — the sole DNA metabolism protein in the list — showed a disproportionately large fold change, motivating the hypothesis that DNA metabolism/repair proteins may play an underappreciated role in PD pathogenesis.

### Structural follow-up: RECQL4

RECQL4 (ATP-dependent DNA helicase Q4, RecQ subfamily, 649 aa) was selected as a related DNA-metabolism protein for structural characterization, using the solved crystal structure [PDB 5LST](https://www.rcsb.org/structure/5LST):

- Overall fold and domain architecture were visualized from the crystal structure.
- The ligand-binding pocket and interacting residues were examined to characterize the active site.
- An interaction-network analysis was performed to identify other proteins potentially linked to RECQL4's function.

The complete write-up — including the volcano plot, clustered heatmap, protein classification figures, structural visualizations, ligand-interaction detail, and interaction network — is in [`docs/report.pdf`](docs/report.pdf).

## How to reproduce

1. Download the GSE206308 sample files from GEO and place them under `data/control/` and `data/disease/` as described above.
2. Run `src/rna_seq_ana.py` as described in [How to run](#how-to-run).
3. Compare the resulting `genes_diff.csv` against [`results/genes_diff.csv`](results/genes_diff.csv) — with the same input files and thresholds, the DE gene list should match.
4. For the downstream protein classification, structural analysis, and network analysis steps, see the external tools referenced in `docs/report.pdf` (DAVID / Gene Accession Conversion Tool, PANTHER, RCSB PDB / molecular viewer for structure `5LST`).

## Assumptions and limitations

- The GSE206308 input files are not included in this repository and must be downloaded separately (see [Setup / installation](#setup--installation)).
- Only `genes_diff.csv` (the DE gene names) is committed to `results/`; `pvals_diff.csv` and `log2fc_diff.csv` are regenerated when the script is run but are not included in the repo.
- Statistical testing uses an unpaired t-test per gene without multiple-testing correction (e.g. FDR/Benjamini-Hochberg) — reasonable for exploratory candidate generation, but a limitation for rigorous DE calling at this scale (63k+ genes tested).
- Sample size is small (3 control, 4 disease replicates), which limits statistical power and generalizability of the findings.
- The protein classification, structural analysis, ligand-binding characterization, and interaction network analysis were performed using external web tools (DAVID, PANTHER, RCSB PDB, and a network-analysis web server) rather than automated in code — see `docs/report.pdf` for the full methodology and figures for those steps.

