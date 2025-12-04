# 多组学基因组预测 (MultiOmicsGP)

一个基于 Python 的多组学数据整合基因组预测软件包。该工具结合基因型（SNP）、转录组（RNA-seq）和表型数据，通过生物学信息指导的特征加权和选择来提高预测准确性。

## 功能特点

- **多组学数据整合**：无缝整合基因型、转录组和表型数据
- **质量控制**：自动化的转录组数据质量控制流程（低方差过滤和 log2 转换）
- **基因组注释**：使用物理基因组坐标将 SNP 映射到基因
- **多种预测方法**：
  - **标准方法**：传统的仅使用基因型的预测方法
  - **加权方法**：基于转录组信息的 SNP 加权方法
  - **特征选择**：基于顶级基因的 SNP 选择方法
- **多种机器学习算法**：支持岭回归、随机森林、支持向量回归（SVR）和多层感知器（MLP）
- **交叉验证**：支持 k 折交叉验证或单次训练-测试分割
- **可视化**：比较不同方法预测性能的对比图（仅单次分割模式）
- **演示模式**：内置模拟数据生成功能，便于测试

## 安装

### 系统要求

- Python 3.7+
- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn

### 安装依赖

```bash
pip install pandas numpy scikit-learn matplotlib seaborn
```

## 快速开始

### 1. 使用演示数据

使用自动生成的模拟数据运行流程：

```bash
python main.py --demo
```

这将：
- 生成模拟牛群数据（基因型、转录组、表型和坐标图谱）
- 运行所有三种预测方法
- 生成对比图

### 2. 使用您自己的数据

准备 CSV 格式的数据文件：

**必需文件：**
- `geno.csv`：基因型矩阵（个体 × SNP），编码为 0/1/2
- `trans.csv`：转录组矩阵（个体 × 基因），原始计数或标准化值
- `pheno.csv`：表型数据（个体 × 性状）
- `snp_map.csv`：SNP 坐标，包含列 `[SNP_ID, Chromosome, Position]`
- `gene_map.csv`：基因坐标，包含列 `[Gene_ID, Chromosome, Start, End]`

**示例命令：**
```bash
python main.py --dir /path/to/data --trait YourTraitName --method all
```

## 使用方法

### 命令行帮助

要查看所有可用选项，请运行：

```bash
python main.py -h
```

这将显示：

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

### 命令行参数

```
--dir DIR              包含输入文件的目录（默认：test_data）
--geno FILE            基因型文件名（默认：geno.csv）
--trans FILE           转录组文件名（默认：trans.csv）
--pheno FILE           表型文件名（默认：pheno.csv）
--map_snp FILE         SNP 位置图谱（默认：snp_map.csv）
--map_gene FILE        基因位置图谱（默认：gene_map.csv）
--trait COLUMN         要预测的表型列（默认：CarcassWeight）
--method METHOD        分析方法：standard、weighted、selected 或 all（默认：all）
--threshold PERCENT    用于 'selected' 方法保留的顶级基因百分比（默认：0.10）
--ml ALGORITHM         机器学习算法：ridge、rf、svr 或 mlp（默认：ridge）
--cv FOLDS             交叉验证折数（1 = 单次分割，默认：1）
--demo                 运行前生成演示数据
```

### 示例工作流程

**运行所有方法：**
```bash
python main.py --dir test_data --method all
```

**仅运行加权方法：**
```bash
python main.py --dir test_data --method weighted
```

**使用自定义阈值运行特征选择：**
```bash
python main.py --dir test_data --method selected --threshold 0.15
```

**使用随机森林和 5 折交叉验证运行：**
```bash
python main.py --dir test_data --method all --ml rf --cv 5
```

**使用 MLP 和单次分割运行（用于可视化）：**
```bash
python main.py --dir test_data --method all --ml mlp --cv 1
```

## 数据格式

### 基因型文件 (`geno.csv`)
- 行：个体 ID
- 列：SNP ID
- 值：0（纯合参考型）、1（杂合型）、2（纯合替代型）

### 转录组文件 (`trans.csv`)
- 行：个体 ID
- 列：基因 ID
- 值：原始读段计数或标准化表达值

### 表型文件 (`pheno.csv`)
- 行：个体 ID
- 列：性状名称
- 值：表型测量值

### SNP 图谱 (`snp_map.csv`)
```csv
SNP_ID,Chromosome,Position
SNP_0,Chr1,1000
SNP_1,Chr1,1100
...
```

### 基因图谱 (`gene_map.csv`)
```csv
Gene_ID,Chromosome,Start,End
Gene_0,Chr1,1000,2000
Gene_1,Chr1,3000,4000
...
```

## 方法学

### 1. 数据加载与同步
- 加载基因型、转录组和表型数据
- 识别所有数据集中的共同个体
- 同步数据矩阵以确保样本顺序一致

### 2. 质量控制
- **转录组**：移除低方差基因（方差 < 0.01）并应用 log2(x+1) 转换
- 注意：基因型 QC（MAF 过滤）在 `QualityControl` 类中可用，但在主流程中默认不应用

### 3. 基因组注释
- 基于物理基因组坐标将 SNP 映射到基因
- 支持可选的缓冲区域（例如，启动子区域）

### 4. 预测方法

#### 标准 GBLUP
仅使用基因型数据的传统岭回归：
```
y = Xβ + ε
```

#### 加权 GBLUP
SNP 根据其关联基因在转录组中的方差进行加权：
```
X_weighted = X × w
```
其中 `w` 由基因表达方差推导得出。

#### 特征选择
选择与顶级方差基因相关的 SNP（例如，前 10% 的基因）：
```
X_selected = X[SNPs ∈ TopGenes]
```

### 5. 模型评估

#### 机器学习算法
- **岭回归**（默认）：alpha=10.0 的岭回归
- **随机森林**：100 个估计器，max_depth=10
- **支持向量回归（SVR）**：RBF 核，C=10.0，epsilon=0.1
- **多层感知器（MLP）**：两个隐藏层（50、25 个神经元），max_iter=500

#### 验证策略
- **单次分割**（默认，`--cv 1`）：70/30 训练-测试分割，用于可视化
- **交叉验证**（`--cv > 1`）：k 折交叉验证，用于稳健的性能估计
- 性能通过预测值与观测值之间的 Pearson 相关系数（r）来衡量

## 项目结构

```
MultiOmicsGP/
├── main.py                 # 主命令行入口
├── main_pipeline.py        # 替代流程脚本
├── data_loader.py          # 数据加载和同步
├── quality_control.py      # 质量控制流程
├── annotation.py           # SNP 到基因的映射
├── models.py               # 预测方法
├── visualization.py        # 绘图函数
├── test_data/             # 示例数据目录
│   ├── geno.csv
│   ├── trans.csv
│   ├── pheno.csv
│   ├── snp_map.csv
│   └── gene_map.csv
└── README.md
```

## 输出结果

流程生成：
- **控制台输出**：每种方法的预测准确性（相关系数 r）
- **可视化**：`prediction_results.png` 比较不同预测方法（仅当 `--cv 1` 时）
- **结果**：相关系数打印到控制台

示例输出（单次分割模式）：
```
  -> Evaluating Standard (G-Only)...
     Single Split Accuracy (r): 0.7234
  -> Evaluating Weighted...
     Single Split Accuracy (r): 0.7891
  -> Evaluating Selected (Top 10%)...
     Single Split Accuracy (r): 0.7654
```

示例输出（交叉验证模式）：
```
  -> Evaluating Standard (G-Only)...
     Running 5-Fold CV... Done. Avg r = 0.7156
  -> Evaluating Weighted...
     Running 5-Fold CV... Done. Avg r = 0.7823
```

## 许可证

[GPL-3]

## 贡献

欢迎贡献！请随时提交 Pull Request。

## 联系方式

[higefee@gmail.com]

