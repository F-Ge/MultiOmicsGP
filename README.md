# Multi-Omics Genomic Prediction (MultiOmicsGP)

A Python-based software package for genomic prediction using multi-omics data integration. This tool combines genotype (SNP), transcriptome (RNA-seq), and phenotype data to improve prediction accuracy through biologically-informed feature weighting and selection.

## Features

- **Multi-omics Data Integration**: Seamlessly integrates genotype, transcriptome, and phenotype data
- **Quality Control**: Automated QC pipeline for transcriptome data (low-variance filtering and log2 transformation)
- **Genomic Annotation**: Maps SNPs to genes using physical genomic coordinates
- **Multiple Prediction Methods**:
  - **Standard**: Traditional genotype-only prediction
  - **Weighted**: Transcriptome-informed SNP weighting
  - **Feature Selection**: Top-gene-based SNP selection
- **Multiple ML Algorithms**: Supports Ridge Regression, Random Forest, Support Vector Regression (SVR), and Multi-Layer Perceptron (MLP)
- **Cross-Validation**: Supports k-fold cross-validation or single train-test split
- **Visualization**: Comparative plots showing prediction performance across methods (single split mode only)
- **Demo Mode**: Built-in mock data generation for testing

## Installation

### Requirements

- Python 3.7+
- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn

### Install Dependencies

```bash
pip install pandas numpy scikit-learn matplotlib seaborn
```

## Quick Start

### 1. Using Demo Data

Run the pipeline with automatically generated mock data:

```bash
python main.py --demo
```

This will:
- Generate mock cattle data (genotypes, transcriptomes, phenotypes, and coordinate maps)
- Run all three prediction methods
- Generate comparison plots

### 2. Using Your Own Data

Prepare your data files in CSV format:

**Required Files:**
- `geno.csv`: Genotype matrix (animals × SNPs), encoded as 0/1/2
- `trans.csv`: Transcriptome matrix (animals × genes), raw counts or normalized
- `pheno.csv`: Phenotype data (animals × traits)
- `snp_map.csv`: SNP coordinates with columns `[SNP_ID, Chromosome, Position]`
- `gene_map.csv`: Gene coordinates with columns `[Gene_ID, Chromosome, Start, End]`

**Example Command:**
```bash
python main.py --dir /path/to/data --trait YourTraitName --method all
```

## Usage

### Command-Line Help

To see all available options, run:

```bash
python main.py -h
```

This will display:

```
usage: main.py [-h] [--dir DIR] [--geno GENO] [--trans TRANS] [--pheno PHENO] [--map_snp MAP_SNP] [--map_gene MAP_GENE]
               [--trait TRAIT] [--method {standard,weighted,selected,all}] [--threshold THRESHOLD]
               [--ml {ridge,rf,svr,mlp}] [--cv CV] [--demo]

Multi-Omics Genomic Prediction Software v2.2

options:
  -h, --help            show this help message and exit
  --dir DIR             Directory containing input files
  --geno GENO           Genotype filename
  --trans TRANS         Transcriptome filename
  --pheno PHENO         Phenotype filename
  --map_snp MAP_SNP     SNP Position Map
  --map_gene MAP_GENE   Gene Position Map
  --trait TRAIT         Phenotype column to predict
  --method {standard,weighted,selected,all}
                        Analysis method to run
  --threshold THRESHOLD
                        Top % of genes to keep (for 'selected' method)
  --ml {ridge,rf,svr,mlp}
                        Machine Learning Algorithm
  --cv CV               Number of Cross-Validation folds (1 = Single Split)
  --demo                Generate demo data before running
```

### Command-Line Arguments

```
--dir DIR              Directory containing input files (default: test_data)
--geno FILE            Genotype filename (default: geno.csv)
--trans FILE           Transcriptome filename (default: trans.csv)
--pheno FILE           Phenotype filename (default: pheno.csv)
--map_snp FILE         SNP position map (default: snp_map.csv)
--map_gene FILE        Gene position map (default: gene_map.csv)
--trait COLUMN         Phenotype column to predict (default: CarcassWeight)
--method METHOD        Analysis method: standard, weighted, selected, or all (default: all)
--threshold PERCENT    Top % of genes to keep for 'selected' method (default: 0.10)
--ml ALGORITHM         ML algorithm: ridge, rf, svr, or mlp (default: ridge)
--cv FOLDS             Number of CV folds (1 = single split, default: 1)
--demo                 Generate demo data before running
```

### Example Workflows

**Run all methods:**
```bash
python main.py --dir test_data --method all
```

**Run only weighted method:**
```bash
python main.py --dir test_data --method weighted
```

**Run feature selection with custom threshold:**
```bash
python main.py --dir test_data --method selected --threshold 0.15
```

**Run with Random Forest and 5-fold cross-validation:**
```bash
python main.py --dir test_data --method all --ml rf --cv 5
```

**Run with MLP and single split (for visualization):**
```bash
python main.py --dir test_data --method all --ml mlp --cv 1
```

## Data Format

### Genotype File (`geno.csv`)
- Rows: Animal IDs
- Columns: SNP IDs
- Values: 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate)

### Transcriptome File (`trans.csv`)
- Rows: Animal IDs
- Columns: Gene IDs
- Values: Raw read counts or normalized expression values

### Phenotype File (`pheno.csv`)
- Rows: Animal IDs
- Columns: Trait names
- Values: Phenotypic measurements

### SNP Map (`snp_map.csv`)
```csv
SNP_ID,Chromosome,Position
SNP_0,Chr1,1000
SNP_1,Chr1,1100
...
```

### Gene Map (`gene_map.csv`)
```csv
Gene_ID,Chromosome,Start,End
Gene_0,Chr1,1000,2000
Gene_1,Chr1,3000,4000
...
```

## Methodology

### 1. Data Loading & Synchronization
- Loads genotype, transcriptome, and phenotype data
- Identifies common animals across all datasets
- Synchronizes data matrices to ensure consistent sample ordering

### 2. Quality Control
- **Transcriptome**: Removes low-variance genes (variance < 0.01) and applies log2(x+1) transformation
- Note: Genotype QC (MAF filtering) is available in the `QualityControl` class but not applied by default in the main pipeline

### 3. Genomic Annotation
- Maps SNPs to genes based on physical genomic coordinates
- Supports optional buffer regions (e.g., promoter regions)

### 4. Prediction Methods

#### Standard GBLUP
Traditional ridge regression using genotype data only:
```
y = Xβ + ε
```

#### Weighted GBLUP
SNPs are weighted by the variance of their associated genes in the transcriptome:
```
X_weighted = X × w
```
where `w` is derived from gene expression variance.

#### Feature Selection
Selects SNPs associated with top-variance genes (e.g., top 10% of genes):
```
X_selected = X[SNPs ∈ TopGenes]
```

### 5. Model Evaluation

#### Machine Learning Algorithms
- **Ridge Regression** (default): Ridge regression with alpha=10.0
- **Random Forest**: 100 estimators, max_depth=10
- **Support Vector Regression (SVR)**: RBF kernel, C=10.0, epsilon=0.1
- **Multi-Layer Perceptron (MLP)**: Two hidden layers (50, 25 neurons), max_iter=500

#### Validation Strategy
- **Single Split** (default, `--cv 1`): 70/30 train-test split for visualization
- **Cross-Validation** (`--cv > 1`): k-fold cross-validation for robust performance estimation
- Performance measured as Pearson correlation (r) between predicted and observed values

## Project Structure

```
MultiOmicsGP/
├── main.py                 # Main CLI entry point
├── main_pipeline.py        # Alternative pipeline script
├── data_loader.py          # Data loading and synchronization
├── quality_control.py      # QC pipelines
├── annotation.py           # SNP-to-gene mapping
├── models.py               # Prediction methods
├── visualization.py        # Plotting functions
├── test_data/             # Example data directory
│   ├── geno.csv
│   ├── trans.csv
│   ├── pheno.csv
│   ├── snp_map.csv
│   └── gene_map.csv
└── README.md
```

## Output

The pipeline generates:
- **Console output**: Prediction accuracy (correlation r) for each method
- **Visualization**: `prediction_results.png` comparing prediction methods (only when `--cv 1`)
- **Results**: Correlation coefficients printed to console

Example output (single split mode):
```
  -> Evaluating Standard (G-Only)...
     Single Split Accuracy (r): 0.7234
  -> Evaluating Weighted...
     Single Split Accuracy (r): 0.7891
  -> Evaluating Selected (Top 10%)...
     Single Split Accuracy (r): 0.7654
```

Example output (cross-validation mode):
```
  -> Evaluating Standard (G-Only)...
     Running 5-Fold CV... Done. Avg r = 0.7156
  -> Evaluating Weighted...
     Running 5-Fold CV... Done. Avg r = 0.7823
```

## License

[GPL-3]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

[higefee@gmail.com]

