# MultiOmicsGP User Manual

**Version 2.2**

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start Guide](#quick-start-guide)
4. [Data Preparation](#data-preparation)
5. [Running the Software](#running-the-software)
6. [Understanding the Output](#understanding-the-output)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)
9. [FAQ](#faq)
10. [Appendix](#appendix)

---

## Introduction

### What is MultiOmicsGP?

MultiOmicsGP is a Python-based software package for genomic prediction that integrates multiple types of omics data (genotype, transcriptome, and phenotype) to improve prediction accuracy. The software uses biologically-informed feature weighting and selection strategies to enhance traditional genomic prediction methods.

### Key Features

- **Multi-omics Integration**: Combines genotype (SNP), transcriptome (RNA-seq), and phenotype data
- **Multiple Prediction Methods**: Standard, weighted, and feature selection approaches
- **Flexible ML Algorithms**: Ridge Regression, Random Forest, SVR, and MLP
- **Cross-Validation Support**: Single split or k-fold cross-validation
- **Quality Control**: Automated transcriptome data processing
- **Genomic Annotation**: SNP-to-gene mapping based on physical coordinates

### When to Use MultiOmicsGP

- Livestock breeding (cattle, pigs, poultry)
- Plant breeding applications
- Any scenario where you have genotype, transcriptome, and phenotype data
- When you want to improve prediction accuracy using multi-omics information

---

## Installation

### System Requirements

- **Operating System**: Linux, macOS, or Windows
- **Python**: Version 3.7 or higher
- **RAM**: Minimum 4GB (8GB+ recommended for large datasets)
- **Disk Space**: ~500MB for software and dependencies

### Step-by-Step Installation

#### 1. Check Python Version

```bash
python --version
# Should show Python 3.7 or higher
```

#### 2. Install Dependencies

```bash
pip install pandas numpy scikit-learn matplotlib seaborn
```

Or install from a requirements file (if available):

```bash
pip install -r requirements.txt
```

#### 3. Verify Installation

```bash
python -c "import pandas, numpy, sklearn, matplotlib, seaborn; print('All packages installed successfully!')"
```

#### 4. Download/Clone the Software

If using Git:

```bash
git clone https://github.com/F-Ge/MultiOmicsGP.git
cd MultiOmicsGP
```

Or download and extract the ZIP file to your desired location.

#### 5. Test the Installation

Run the demo mode to verify everything works:

```bash
python main.py --demo
```

If successful, you should see output indicating data generation and analysis completion.

---

## Quick Start Guide

### First Run: Demo Mode

The fastest way to get started is using the built-in demo mode:

```bash
python main.py --demo
```

This command will:
1. Generate mock cattle data (100 animals, 500 SNPs, 50 genes)
2. Run all three prediction methods
3. Generate comparison plots
4. Display results in the console

**Expected Output:**
- Console messages showing data generation
- Data loading and processing messages
- Prediction accuracy results for each method
- A saved plot file: `prediction_results.png`

### Understanding Demo Results

After running `--demo`, you'll see:
- **Standard method**: Baseline prediction using only genotypes
- **Weighted method**: Improved prediction using transcriptome-weighted SNPs
- **Selected method**: Prediction using top-variance genes

The correlation coefficient (r) indicates prediction accuracy (higher is better, range: -1 to 1).

---

## Data Preparation

### Required Data Files

You need five CSV files to run MultiOmicsGP with your own data:

1. **Genotype file** (`geno.csv`)
2. **Transcriptome file** (`trans.csv`)
3. **Phenotype file** (`pheno.csv`)
4. **SNP map file** (`snp_map.csv`)
5. **Gene map file** (`gene_map.csv`)

### File Format Specifications

#### 1. Genotype File (`geno.csv`)

**Format**: CSV with animal IDs as row names and SNP IDs as column names

**Encoding**:
- `0` = Homozygous reference (AA)
- `1` = Heterozygous (Aa)
- `2` = Homozygous alternate (aa)

**Example**:
```csv
,SNP_0,SNP_1,SNP_2,SNP_3
Animal_1,0,1,2,0
Animal_2,1,1,0,2
Animal_3,2,0,1,1
```

**Requirements**:
- First column should be animal IDs (will be used as index)
- Missing values should be handled (imputed or removed) before analysis
- All values must be 0, 1, or 2

**Tips**:
- Use consistent animal ID naming across all files
- Ensure no duplicate animal IDs
- Check for missing data before running

#### 2. Transcriptome File (`trans.csv`)

**Format**: CSV with animal IDs as row names and gene IDs as column names

**Values**: Raw read counts or normalized expression values

**Example**:
```csv
,Gene_0,Gene_1,Gene_2,Gene_3
Animal_1,1250.5,342.1,89.3,1567.8
Animal_2,1180.2,398.5,95.1,1423.4
Animal_3,1320.7,315.9,102.6,1689.2
```

**Requirements**:
- First column should be animal IDs
- Values should be non-negative (counts or normalized expression)
- Missing values will cause issues - handle them beforehand

**Tips**:
- Can use raw counts or normalized values (TPM, FPKM, etc.)
- The software will apply log2(x+1) transformation automatically
- Not all animals need transcriptome data (software will synchronize)

#### 3. Phenotype File (`pheno.csv`)

**Format**: CSV with animal IDs as row names and trait names as column names

**Example**:
```csv
,CarcassWeight,MilkYield,FeedEfficiency
Animal_1,350.5,8200.3,2.45
Animal_2,385.2,9100.1,2.38
Animal_3,320.8,7800.5,2.52
```

**Requirements**:
- First column should be animal IDs
- Specify the trait column using `--trait` argument
- Missing values will be automatically removed for that animal

**Tips**:
- You can have multiple traits in one file
- Use descriptive trait names
- Ensure trait values are numeric

#### 4. SNP Map File (`snp_map.csv`)

**Format**: CSV with three columns: SNP_ID, Chromosome, Position

**Example**:
```csv
SNP_ID,Chromosome,Position
SNP_0,Chr1,1000
SNP_1,Chr1,1100
SNP_2,Chr1,1200
SNP_3,Chr2,5000
```

**Requirements**:
- Must have exactly these three columns
- SNP_ID must match column names in genotype file
- Position should be in base pairs
- Chromosome names should be consistent

**Tips**:
- Use standard chromosome naming (Chr1, Chr2, etc. or 1, 2, etc.)
- Ensure all SNPs in genotype file have entries in SNP map
- Positions should be numeric

#### 5. Gene Map File (`gene_map.csv`)

**Format**: CSV with four columns: Gene_ID, Chromosome, Start, End

**Example**:
```csv
Gene_ID,Chromosome,Start,End
Gene_0,Chr1,1000,2000
Gene_1,Chr1,3000,4000
Gene_2,Chr2,5000,6000
```

**Requirements**:
- Must have exactly these four columns
- Gene_ID must match column names in transcriptome file
- Start and End positions define gene boundaries
- Start < End for all genes

**Tips**:
- Gene boundaries should be accurate
- Overlapping genes are allowed
- Ensure all genes in transcriptome file have entries in gene map

### Data Synchronization

The software automatically:
- Finds animals present in ALL three data files (genotype, transcriptome, phenotype)
- Removes animals missing from any file
- Synchronizes row order across all datasets
- Reports how many common animals were found

**Important**: Make sure animal IDs are consistent across all files!

### Data Quality Checks

Before running analysis, check:

1. **Animal ID Consistency**: Same IDs in all files?
2. **Missing Data**: Are there many missing values?
3. **File Formats**: Are files properly formatted CSV?
4. **Coordinate Maps**: Do all SNPs/genes have map entries?
5. **Data Types**: Are numeric columns actually numeric?

---

## Running the Software

### Basic Command Structure

```bash
python main.py [OPTIONS]
```

### Common Workflows

#### Workflow 1: Quick Analysis with Default Settings

```bash
python main.py --dir test_data --trait CarcassWeight
```

This runs:
- All three methods (standard, weighted, selected)
- Ridge regression (default ML algorithm)
- Single train-test split (default validation)
- Generates visualization

#### Workflow 2: Compare Methods with Cross-Validation

```bash
python main.py --dir test_data --trait CarcassWeight --method all --cv 5
```

This runs:
- All three methods
- 5-fold cross-validation (more robust estimates)
- No visualization (CV mode doesn't generate plots)

#### Workflow 3: Test Different ML Algorithms

```bash
# Random Forest
python main.py --dir test_data --trait CarcassWeight --ml rf

# Support Vector Regression
python main.py --dir test_data --trait CarcassWeight --ml svr

# Multi-Layer Perceptron
python main.py --dir test_data --trait CarcassWeight --ml mlp
```

#### Workflow 4: Feature Selection with Custom Threshold

```bash
python main.py --dir test_data --trait CarcassWeight --method selected --threshold 0.15
```

This keeps top 15% of genes (instead of default 10%).

### Command-Line Options Explained

#### Data Input Options

| Option | Description | Default | Example |
|--------|-------------|---------|---------|
| `--dir` | Directory containing input files | `test_data` | `--dir /path/to/data` |
| `--geno` | Genotype filename | `geno.csv` | `--geno my_genotypes.csv` |
| `--trans` | Transcriptome filename | `trans.csv` | `--trans rna_seq.csv` |
| `--pheno` | Phenotype filename | `pheno.csv` | `--pheno traits.csv` |
| `--map_snp` | SNP position map | `snp_map.csv` | `--map_snp snp_coords.csv` |
| `--map_gene` | Gene position map | `gene_map.csv` | `--map_gene gene_coords.csv` |
| `--trait` | Phenotype column to predict | `CarcassWeight` | `--trait MilkYield` |

#### Analysis Options

| Option | Description | Default | Choices |
|--------|-------------|---------|---------|
| `--method` | Prediction method(s) to run | `all` | `standard`, `weighted`, `selected`, `all` |
| `--threshold` | Top % genes for selection | `0.10` | Any float 0.0-1.0 |
| `--ml` | Machine learning algorithm | `ridge` | `ridge`, `rf`, `svr`, `mlp` |
| `--cv` | Cross-validation folds | `1` | Integer ≥ 1 |

#### Utility Options

| Option | Description |
|--------|-------------|
| `--demo` | Generate demo data before running |
| `-h`, `--help` | Show help message |

### Method Selection Guide

#### Standard Method (`--method standard`)

**When to use:**
- Baseline comparison
- When transcriptome data quality is poor
- Quick analysis

**What it does:**
- Uses only genotype data
- Traditional GBLUP approach
- Fastest method

#### Weighted Method (`--method weighted`)

**When to use:**
- You have good quality transcriptome data
- Want to incorporate biological information
- Expect gene expression to be informative

**What it does:**
- Weights SNPs by their associated gene expression variance
- SNPs in highly variable genes get higher weights
- Maintains all SNPs (no filtering)

#### Selected Method (`--method selected`)

**When to use:**
- Want to reduce dimensionality
- Focus on most informative genes
- Computational efficiency is important

**What it does:**
- Selects SNPs in top-variance genes
- Filters out SNPs in low-variance genes
- Reduces feature space

#### All Methods (`--method all`)

**When to use:**
- Comparing all approaches
- First-time analysis
- Publication-ready results

**What it does:**
- Runs all three methods
- Provides comprehensive comparison

### ML Algorithm Selection Guide

#### Ridge Regression (`--ml ridge`)

**Best for:**
- Linear relationships
- Large datasets
- Baseline comparisons
- Fast computation

**Parameters:**
- Alpha (regularization) = 10.0

#### Random Forest (`--ml rf`)

**Best for:**
- Non-linear relationships
- Feature interactions
- Robust to outliers

**Parameters:**
- n_estimators = 100
- max_depth = 10

#### Support Vector Regression (`--ml svr`)

**Best for:**
- Non-linear patterns
- Small to medium datasets
- Complex relationships

**Parameters:**
- Kernel = RBF
- C = 10.0
- epsilon = 0.1

#### Multi-Layer Perceptron (`--ml mlp`)

**Best for:**
- Deep learning approach
- Complex non-linear patterns
- When other methods plateau

**Parameters:**
- Hidden layers: (50, 25)
- max_iter = 500

### Validation Strategy Guide

#### Single Split (`--cv 1`)

**Use when:**
- You want visualization plots
- Quick results needed
- Dataset is very large

**How it works:**
- 70% training, 30% testing
- One prediction per sample
- Generates scatter plots

#### Cross-Validation (`--cv 5` or higher)

**Use when:**
- Robust performance estimation needed
- Small to medium datasets
- Publication-quality results

**How it works:**
- k-fold cross-validation
- Average performance across folds
- More reliable estimates
- No visualization (too many plots)

**Recommended folds:**
- Small datasets (<100 samples): 5-fold
- Medium datasets (100-500): 5-10 fold
- Large datasets (>500): 10-fold

---

## Understanding the Output

### Console Output

#### Data Loading Phase

```
Loading data from 'test_data'...
Loading Genotypes from test_data/geno.csv...
Loaded Genotypes: (100, 500) (Animals x SNPs)
Loading Transcriptome from test_data/trans.csv...
Loaded Transcriptome: (80, 50) (Animals x Genes)
Loading Phenotypes from test_data/pheno.csv...
Loaded Phenotypes: (100, 1) (Animals)

--- Synchronizing Data ---
Found 80 common animals.
Data synchronized successfully.
```

**What to check:**
- Are the dimensions reasonable?
- How many common animals were found?
- Any warnings about missing data?

#### Quality Control Phase

```
--- Transcriptome QC & Normalization ---
Removed 2 low-variance genes.
Applying Log2(x + 1) transformation...
Final Transcriptome: 48 Genes
```

**What to check:**
- How many genes were filtered?
- Is the final gene count reasonable?

#### Annotation Phase

```
--- Loading Genomic Coordinates ---
Loaded 500 SNP positions and 50 Gene regions.
Mapping SNPs to Genes...
Successfully mapped 450 SNPs out of 500 (90.0%).
```

**What to check:**
- Mapping percentage (should be >80% typically)
- Are unmapped SNPs a concern?

#### Analysis Phase

**Single Split Mode:**
```
[Analysis] Running Method: ALL
[Model] Algorithm: RIDGE
  -> Evaluating Standard (G-Only)...
     Single Split Accuracy (r): 0.7234
  -> Evaluating Weighted...
     Single Split Accuracy (r): 0.7891
  -> Evaluating Selected (Top 10%)...
     Single Split Accuracy (r): 0.7654
```

**Cross-Validation Mode:**
```
[Analysis] Running Method: ALL
[Model] Algorithm: RIDGE
  -> Evaluating Standard (G-Only)...
     Running 5-Fold CV... Done. Avg r = 0.7156
  -> Evaluating Weighted...
     Running 5-Fold CV... Done. Avg r = 0.7823
  -> Evaluating Selected (Top 10%)...
     Running 5-Fold CV... Done. Avg r = 0.7589
```

### Interpreting Results

#### Correlation Coefficient (r)

- **Range**: -1 to +1
- **Interpretation**:
  - `r > 0.7`: Excellent prediction
  - `r = 0.5-0.7`: Good prediction
  - `r = 0.3-0.5`: Moderate prediction
  - `r < 0.3`: Poor prediction
  - `r < 0`: Negative correlation (model is worse than random)

#### Comparing Methods

- **Standard vs Weighted**: If weighted is higher, transcriptome information helps
- **Standard vs Selected**: If selected is higher, feature selection is beneficial
- **Weighted vs Selected**: Which biological strategy works better?

#### Visualization Output

When using `--cv 1`, a plot file `prediction_results.png` is generated showing:
- **Left panel**: Standard method predictions vs true values
- **Right panel**: Best alternative method (weighted or selected) vs true values
- **Correlation values**: Displayed on each plot
- **Trend lines**: Red lines showing prediction accuracy

**How to read the plots:**
- Points closer to the diagonal line = better predictions
- Tight clustering = consistent predictions
- Scattered points = high prediction error

---

## Advanced Usage

### Custom Data Directory Structure

Organize your data in a custom structure:

```
my_project/
├── data/
│   ├── genotypes/
│   │   └── cattle_geno.csv
│   ├── transcriptomes/
│   │   └── cattle_rna.csv
│   └── phenotypes/
│       └── cattle_traits.csv
├── maps/
│   ├── snp_coordinates.csv
│   └── gene_coordinates.csv
└── results/
```

Run with:
```bash
python main.py \
  --dir data \
  --geno genotypes/cattle_geno.csv \
  --trans transcriptomes/cattle_rna.csv \
  --pheno phenotypes/cattle_traits.csv \
  --map_snp ../maps/snp_coordinates.csv \
  --map_gene ../maps/gene_coordinates.csv \
  --trait CarcassWeight
```

### Batch Processing Multiple Traits

Create a script to process multiple traits:

```bash
#!/bin/bash
traits=("CarcassWeight" "MilkYield" "FeedEfficiency")

for trait in "${traits[@]}"; do
    echo "Processing $trait..."
    python main.py \
        --dir test_data \
        --trait "$trait" \
        --method all \
        --cv 5 \
        --ml ridge
done
```

### Comparing ML Algorithms

Test which algorithm works best for your data:

```bash
for ml in ridge rf svr mlp; do
    echo "Testing $ml..."
    python main.py \
        --dir test_data \
        --trait CarcassWeight \
        --method all \
        --ml "$ml" \
        --cv 5
done
```

### Feature Selection Sensitivity Analysis

Test different thresholds:

```bash
for threshold in 0.05 0.10 0.15 0.20 0.25; do
    echo "Testing threshold $threshold..."
    python main.py \
        --dir test_data \
        --trait CarcassWeight \
        --method selected \
        --threshold "$threshold" \
        --cv 5
done
```

### Integration with Other Tools

#### Export Results to CSV

Modify the code or use Python wrapper:

```python
import subprocess
import pandas as pd

# Run analysis
result = subprocess.run(
    ['python', 'main.py', '--dir', 'test_data', '--trait', 'CarcassWeight'],
    capture_output=True, text=True
)

# Parse and save results
# (You would need to parse the output or modify the code to save results)
```

---

## Troubleshooting

### Common Errors and Solutions

#### Error: "No common animals found!"

**Cause**: Animal IDs don't match across files

**Solution**:
1. Check animal ID format in all files
2. Ensure no extra spaces or special characters
3. Verify IDs are in the first column (index)
4. Use consistent naming (e.g., all uppercase or all lowercase)

**Check with**:
```python
import pandas as pd
geno = pd.read_csv('geno.csv', index_col=0)
trans = pd.read_csv('trans.csv', index_col=0)
pheno = pd.read_csv('pheno.csv', index_col=0)
print("Geno IDs:", geno.index.tolist()[:5])
print("Trans IDs:", trans.index.tolist()[:5])
print("Pheno IDs:", pheno.index.tolist()[:5])
```

#### Error: "Trait 'X' not found in file"

**Cause**: Trait name doesn't match column name in phenotype file

**Solution**:
1. Check exact column name (case-sensitive)
2. List available traits:
   ```python
   import pandas as pd
   pheno = pd.read_csv('pheno.csv', index_col=0)
   print("Available traits:", pheno.columns.tolist())
   ```
3. Use exact column name with `--trait`

#### Error: "Failed to load data"

**Cause**: File path issues or format problems

**Solution**:
1. Check file paths are correct
2. Verify files exist in specified directory
3. Ensure files are valid CSV format
4. Check file permissions

#### Warning: "No SNPs selected! Returning full matrix"

**Cause**: Feature selection threshold too strict or mapping issues

**Solution**:
1. Increase `--threshold` value
2. Check SNP-to-gene mapping (run annotation step separately)
3. Verify transcriptome data has variance

#### Poor Prediction Results (low r values)

**Possible causes and solutions**:

1. **Small sample size**
   - Solution: Need more animals (typically >100 for good results)

2. **Low heritability trait**
   - Solution: Some traits are inherently hard to predict

3. **Data quality issues**
   - Solution: Check for genotyping errors, missing data

4. **Wrong ML algorithm**
   - Solution: Try different algorithms (`--ml rf`, `--ml svr`)

5. **Trait not suitable for genomic prediction**
   - Solution: Some traits may not be predictable from genotypes

#### Memory Issues

**Symptoms**: Program crashes or runs very slowly

**Solutions**:
1. Reduce dataset size (fewer SNPs/genes)
2. Use feature selection method
3. Increase system RAM
4. Process in batches

#### Visualization Not Generated

**Cause**: Running in cross-validation mode (`--cv > 1`)

**Solution**: Use `--cv 1` for single split mode to generate plots

### Debugging Tips

#### 1. Check Data Dimensions

```python
import pandas as pd
geno = pd.read_csv('test_data/geno.csv', index_col=0)
trans = pd.read_csv('test_data/trans.csv', index_col=0)
pheno = pd.read_csv('test_data/pheno.csv', index_col=0)

print(f"Genotype: {geno.shape}")
print(f"Transcriptome: {trans.shape}")
print(f"Phenotype: {pheno.shape}")
```

#### 2. Verify Data Types

```python
print("Genotype dtypes:", geno.dtypes.value_counts())
print("Transcriptome dtypes:", trans.dtypes.value_counts())
print("Phenotype dtypes:", pheno.dtypes.value_counts())
```

#### 3. Check for Missing Data

```python
print("Missing in genotype:", geno.isnull().sum().sum())
print("Missing in transcriptome:", trans.isnull().sum().sum())
print("Missing in phenotype:", pheno.isnull().sum().sum())
```

#### 4. Test Individual Components

Run the pipeline step by step using `main_pipeline.py` or modify `main.py` to add debug prints.

---

## FAQ

### General Questions

**Q: Can I use this for human data?**
A: Yes, the software works with any species as long as you have the required data formats.

**Q: Do I need all three data types (genotype, transcriptome, phenotype)?**
A: Yes, all three are required. The software integrates them for improved prediction.

**Q: Can I use imputed genotypes?**
A: Yes, as long as they're encoded as 0/1/2 and in the correct format.

**Q: What if I don't have transcriptome data for all animals?**
A: The software will automatically use only animals with all three data types.

**Q: How long does analysis take?**
A: Depends on dataset size:
- Small (<100 animals, <1000 SNPs): Seconds
- Medium (100-500 animals, 1000-10000 SNPs): Minutes
- Large (>500 animals, >10000 SNPs): Hours

### Data Questions

**Q: What format should transcriptome data be in?**
A: Raw counts or normalized values (TPM, FPKM) are both fine. The software applies log2 transformation.

**Q: Can I use different chromosome naming conventions?**
A: Yes, but be consistent. Use the same format in both SNP and gene maps.

**Q: What if SNPs are not mapped to genes?**
A: Unmapped SNPs get default weight of 1.0 in weighted method and are excluded in selected method.

**Q: Can I use phased genotypes?**
A: The software expects unphased genotypes (0/1/2 encoding). Phased data would need conversion.

### Method Questions

**Q: Which method should I use?**
A: Start with `--method all` to compare all approaches, then choose the best for your data.

**Q: When should I use feature selection?**
A: When you have many SNPs and want to focus on the most informative ones, or for computational efficiency.

**Q: What's the difference between weighted and selected?**
A: Weighted keeps all SNPs but adjusts their importance. Selected filters to only top-gene SNPs.

**Q: Which ML algorithm is best?**
A: Depends on your data. Try all and compare. Ridge is often a good starting point.

### Technical Questions

**Q: Can I run this on a cluster?**
A: Yes, it's a command-line tool that can be submitted as a job.

**Q: How do I parallelize analysis?**
A: You can run different traits or methods in parallel using job arrays or Python multiprocessing.

**Q: Can I modify the code?**
A: Yes, it's open source. Modify as needed for your specific requirements.

**Q: How do I cite this software?**
A: See the README file for citation information.

### Results Questions

**Q: What's a good prediction accuracy?**
A: Depends on trait heritability. For highly heritable traits, r > 0.5 is good. For low heritability, r > 0.3 may be acceptable.

**Q: Why are my results different each time?**
A: If using random splits or cross-validation with shuffling, results will vary slightly. Use `random_state` for reproducibility.

**Q: Can I get prediction values for individual animals?**
A: Currently, the software reports correlation. You could modify the code to save individual predictions.

**Q: How do I interpret negative correlations?**
A: Negative correlations indicate the model is performing worse than random. Check data quality and method selection.

---

## Appendix

### A. File Format Examples

#### Complete Genotype File Example
```csv
,SNP_0,SNP_1,SNP_2,SNP_3,SNP_4
Cattle_1,0,1,2,0,1
Cattle_2,1,1,0,2,0
Cattle_3,2,0,1,1,2
Cattle_4,0,2,0,1,1
Cattle_5,1,1,2,0,2
```

#### Complete Transcriptome File Example
```csv
,Gene_0,Gene_1,Gene_2,Gene_3,Gene_4
Cattle_1,1250.5,342.1,89.3,1567.8,234.5
Cattle_2,1180.2,398.5,95.1,1423.4,267.8
Cattle_3,1320.7,315.9,102.6,1689.2,198.3
Cattle_4,1105.3,421.2,87.9,1512.6,289.1
Cattle_5,1289.6,356.7,94.5,1623.9,245.7
```

#### Complete Phenotype File Example
```csv
,CarcassWeight,MilkYield,FeedEfficiency
Cattle_1,350.5,8200.3,2.45
Cattle_2,385.2,9100.1,2.38
Cattle_3,320.8,7800.5,2.52
Cattle_4,365.7,8750.2,2.41
Cattle_5,340.2,8050.8,2.49
```

### B. Command Reference Card

```bash
# Basic usage
python main.py --dir DATA_DIR --trait TRAIT_NAME

# Full options
python main.py \
  --dir test_data \
  --geno geno.csv \
  --trans trans.csv \
  --pheno pheno.csv \
  --map_snp snp_map.csv \
  --map_gene gene_map.csv \
  --trait CarcassWeight \
  --method all \
  --threshold 0.10 \
  --ml ridge \
  --cv 1 \
  --demo

# Quick reference
python main.py -h  # Show all options
```

### C. Performance Benchmarks

Typical runtime on a standard laptop (8GB RAM, 4 cores):

| Dataset Size | Runtime (Single Split) | Runtime (5-Fold CV) |
|--------------|------------------------|---------------------|
| 100 animals, 1K SNPs | < 1 minute | 2-3 minutes |
| 500 animals, 10K SNPs | 5-10 minutes | 20-30 minutes |
| 1000 animals, 50K SNPs | 30-60 minutes | 2-3 hours |

### D. Glossary

- **GBLUP**: Genomic Best Linear Unbiased Prediction
- **MAF**: Minor Allele Frequency
- **SNP**: Single Nucleotide Polymorphism
- **QC**: Quality Control
- **CV**: Cross-Validation
- **ML**: Machine Learning
- **Ridge**: Ridge Regression
- **RF**: Random Forest
- **SVR**: Support Vector Regression
- **MLP**: Multi-Layer Perceptron

### E. Additional Resources

- **Documentation**: See README.md for overview
- **Code Repository**: Check GitHub for latest updates
- **Issues**: Report bugs on GitHub issues page
- **Contact**: higefee@gmail.com

---

**End of User Manual**

*Last Updated: Version 2.2*

