# Temporal dynamics of RNA phage infection reveal bacterial recovery despite ongoing viral replication

Code to reproduce the figures and tables of:

> Harel N\*, Resisi Y\*, Meir M\*, Buchshtab N, Ghosh D, Golding I, Gophna U, Rosenthal A, Stern A. Temporal dynamics of RNA phage infection reveal bacterial recovery despite ongoing viral replication. 
>
> \* These authors contributed equally.

## Repository structure

```
.
├── LICENSE
├── config.example.yaml   # template for config.yaml (paths to data and outputs)
├── envs/                 # conda environments
├── data/                 # small input files not deposited elsewhere
├── preprocessing/        # scripts that produce some of the files in data/
└── notebooks/            # one notebook per dataset or experiment
```

## Data

Sequencing data are available from GEO SuperSeries [GSE349161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349161) and SRA BioProject [PRJNA1536727](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1536727). The *E. coli* C-3000 (ATCC 15597) and MS2 reference genomes are available from the [ATCC Genome Portal](https://genomes.atcc.org).

| Accession | Dataset | Files used by the notebooks |
|---|---|---|
| GSE348942 | Bulk RNA-seq time course | `bulk_RNAseq_htseq_raw_counts.csv` |
| GSE349085 | Isolate bulk RNA-seq | `isolates_RNAseq_htseq_raw_counts.csv` |
| GSE349086 | Probe-based scRNA-seq, 180 min | `INF180min_raw_feature_bc_matrix.h5` |
| GSE349090 | Probe-based scRNA-seq, 15, 30 and 45 min | `CTRL15min_raw_feature_bc_matrix.h5`, `INF15min_…`, `CTRL30min_…`, `INF30min_…`, `CTRL45min_…`, `INF45min_…` |
| PRJNA1536727 | Isolate whole-genome sequencing (Oxford Nanopore) | raw reads; processed to `data/WGS_plasmid_chromosome_coverage.csv` |

Download the GEO supplementary files, decompress them (`gunzip`), and place them in one folder per accession:

```
<data_dir>/
├── GSE348942/bulk_RNAseq_htseq_raw_counts.csv
├── GSE349085/isolates_RNAseq_htseq_raw_counts.csv
├── GSE349086/INF180min_raw_feature_bc_matrix.h5
└── GSE349090/CTRL15min_raw_feature_bc_matrix.h5, ...
```

GEO adds the sample accession to the names of sample-level files (e.g. `GSM…_CTRL15min_raw_feature_bc_matrix.h5`); the notebooks find these files with or without the prefix.

All other inputs (plate-reader, qPCR, plaque assay and smFISH measurements, sample metadata, gene annotation, and processed coverage tables) are in `data/`. Raw smFISH images are available from the authors upon request.

## Setup

### 1. Paths

Copy the template and set the paths:

```bash
cp config.example.yaml config.yaml
```

| Key | Meaning |
|---|---|
| `data_dir` | folder containing the GEO subfolders described above |
| `out_dir` | folder where figures, tables and intermediate results are written |
| `wgs_bam_dir` | only needed for `preprocessing/WGS_coverage.py` |
| `smfish_nd2_dir` | image folders of the two smFISH experiments (`exp1`, `exp2`); only needed for `preprocessing/smFISH_image_processing.py` |

### 2. Environments

Three conda environments are used by the notebooks; a fourth (`smFISH_env`) is only needed for the smFISH image-processing script. Create them and register the notebook environments as Jupyter kernels:

```bash
conda env create -f envs/python_env.yml
conda env create -f envs/bulk_R_env.yml
conda env create -f envs/single_cell_R_env.yml
conda env create -f envs/smFISH_env.yml

conda activate python_env
python -m ipykernel install --user --name python_env --display-name "Python (python_env)"

conda activate bulk_R_env
R -e 'IRkernel::installspec(name = "bulk_r_env", displayname = "R (bulk_R_env)")'

conda activate single_cell_R_env
R -e 'IRkernel::installspec(name = "single_cell_r_env", displayname = "R (single_cell_R_env)")'
```

Each notebook states its environment in its first cell and opens with the matching kernel. Installing [presto](https://github.com/immunogenomics/presto) in `single_cell_R_env` is optional; it speeds up cluster-marker computation.

## Running the notebooks

Run the notebooks from the `notebooks/` folder. They are independent, except that `bulk_RNAseq_1_DESeq2` must be run before `bulk_RNAseq_2_TableS1` and `bulk_RNAseq_3_supplementary`.

Each notebook writes its figures (PNG and PDF) to `out_dir/figures/`, tables to `out_dir/tables/`, and the data underlying each plotted panel to `out_dir/source_data/`.

| Notebook | Environment | Figures and tables |
|---|---|---|
| `CFU_PFU` | python_env | Fig 1A |
| `qPCR_MS2_strands` | python_env | Fig S1 |
| `bulk_RNAseq_1_DESeq2` | bulk_R_env | Fig 1B, 1D, 3A, 4A, S11; DESeq2 results |
| `bulk_RNAseq_2_TableS1` | python_env | Table S1 |
| `bulk_RNAseq_3_supplementary` | python_env | Fig S2, S3, S4, S14; gene lists for Tables S2, S4 |
| `scRNAseq_ProBacseq` | single_cell_R_env | Fig 2A, 2B, 2C, S5, S6, S7, S8, S9A, S9B |
| `smFISH` | python_env | Fig 2D, S10 (experiment 1); Table S3 (experiments 1 and 2) |
| `knockouts_growth_endpoint_titer` | python_env | Fig 3B, 3C |
| `WGS_plasmid_coverage` | python_env | Fig 4B |
| `isolates_RNAseq` | bulk_R_env | Fig 4C, 4D, S12A, S12B, S13 |

### Figure and table index

| Item | Notebook |
|---|---|
| Fig 1A | `CFU_PFU` |
| Fig 1B, 1D | `bulk_RNAseq_1_DESeq2` |
| Fig 2A–2C | `scRNAseq_ProBacseq` |
| Fig 2D | `smFISH` |
| Fig 3A | `bulk_RNAseq_1_DESeq2` |
| Fig 3B, 3C | `knockouts_growth_endpoint_titer` |
| Fig 4A | `bulk_RNAseq_1_DESeq2` |
| Fig 4B | `WGS_plasmid_coverage` |
| Fig 4C, 4D | `isolates_RNAseq` |
| Fig S1 | `qPCR_MS2_strands` |
| Fig S2, S3, S4 | `bulk_RNAseq_3_supplementary` |
| Fig S5–S9 | `scRNAseq_ProBacseq` |
| Fig S10 | `smFISH` |
| Fig S11 | `bulk_RNAseq_1_DESeq2` |
| Fig S12, S13 | `isolates_RNAseq` |
| Fig S14 | `bulk_RNAseq_3_supplementary` |
| Table S1 | `bulk_RNAseq_2_TableS1` |
| Tables S2, S4 | gene lists from `bulk_RNAseq_3_supplementary`; functional enrichment run on the STRING-DB website (see Methods) |
| Table S3 | `smFISH` |
| Table S5 | `preprocessing/WGS_lofreq.sh` and `preprocessing/WGS_entropy_filter.py` (see below) |
| Tables S6–S8 | not generated by code; S7 and S8 list the probe sets deposited in GSE349090 and GSE349086 |

Fig 1C and Fig 5 are schematics, and Fig 2E and 2F are representative microscopy images; they are not generated by code.

## Preprocessing

These scripts produce inputs in `data/` from alignment files that are not part of the repository. Their outputs are already included, so they only need to be run to reproduce those files.

| Script | Produces | Input |
|---|---|---|
| `bulk_MS2_coverage.py` | `data/bulk_MS2_strand_depth.csv.gz`, `data/bulk_mapped_reads.csv` (Fig S2) | strand-specific MS2 depth files and BAM files of the bulk RNA-seq samples |
| `WGS_coverage.py` | `data/WGS_plasmid_chromosome_coverage.csv` (Fig 4B) | minimap2 alignments of the PRJNA1536727 reads to the *E. coli* ATCC 15597 reference |
| `WGS_lofreq.sh`, `WGS_entropy_filter.py` | variant calls (Table S5) | the same alignments |
| `smFISH_image_processing.py` (environment `smFISH_env`) | `data/smFISH_exp1_per_cell.csv`, `data/smFISH_exp2_per_cell.csv` (Fig 2D, S10; Table S3) | ND2 image stacks of the two experiments, listed in `data/smFISH_exp1_image_manifest.csv` and `data/smFISH_exp2_image_manifest.csv`, which also record the focus plane used for segmentation and any excluded fields |

Read alignment and the other upstream processing steps are described in the Methods.

## License

This repository is released under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/); see [LICENSE](LICENSE). It may be used and adapted for non-commercial purposes with appropriate credit.

## Contact

Adi Stern (sternadi@tauex.tau.ac.il)
