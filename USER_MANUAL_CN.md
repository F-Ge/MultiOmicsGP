# MultiOmicsGP 用户手册

**版本 2.2**

## 目录

1. [简介](#简介)
2. [安装](#安装)
3. [快速开始指南](#快速开始指南)
4. [数据准备](#数据准备)
5. [运行软件](#运行软件)
6. [理解输出结果](#理解输出结果)
7. [高级用法](#高级用法)
8. [故障排除](#故障排除)
9. [常见问题](#常见问题)
10. [附录](#附录)

---

## 简介

### 什么是 MultiOmicsGP？

MultiOmicsGP 是一个基于 Python 的基因组预测软件包，它整合多种类型的组学数据（基因型、转录组和表型）以提高预测准确性。该软件使用生物学信息指导的特征加权和选择策略来增强传统的基因组预测方法。

### 主要功能

- **多组学整合**：结合基因型（SNP）、转录组（RNA-seq）和表型数据
- **多种预测方法**：标准、加权和特征选择方法
- **灵活的机器学习算法**：岭回归、随机森林、SVR 和 MLP
- **交叉验证支持**：单次分割或 k 折交叉验证
- **质量控制**：自动化的转录组数据处理
- **基因组注释**：基于物理坐标的 SNP 到基因映射

### 何时使用 MultiOmicsGP

- 家畜育种（牛、猪、家禽）
- 植物育种应用
- 任何拥有基因型、转录组和表型数据的场景
- 当您想使用多组学信息提高预测准确性时

---

## 安装

### 系统要求

- **操作系统**：Linux、macOS 或 Windows
- **Python**：3.7 或更高版本
- **内存**：最少 4GB（大型数据集建议 8GB+）
- **磁盘空间**：约 500MB（用于软件和依赖项）

### 逐步安装

#### 1. 检查 Python 版本

```bash
python --version
# 应显示 Python 3.7 或更高版本
```

#### 2. 安装依赖项

```bash
pip install pandas numpy scikit-learn matplotlib seaborn
```

或从 requirements 文件安装（如果可用）：

```bash
pip install -r requirements.txt
```

#### 3. 验证安装

```bash
python -c "import pandas, numpy, sklearn, matplotlib, seaborn; print('所有包安装成功！')"
```

#### 4. 下载/克隆软件

如果使用 Git：

```bash
git clone https://github.com/F-Ge/MultiOmicsGP.git
cd MultiOmicsGP
```

或下载 ZIP 文件并解压到所需位置。

#### 5. 测试安装

运行演示模式以验证一切正常：

```bash
python main.py --demo
```

如果成功，您应该看到指示数据生成和分析完成的输出。

---

## 快速开始指南

### 首次运行：演示模式

最快开始的方法是使用内置的演示模式：

```bash
python main.py --demo
```

此命令将：
1. 生成模拟牛群数据（100 个个体，500 个 SNP，50 个基因）
2. 运行所有三种预测方法
3. 生成对比图
4. 在控制台显示结果

**预期输出：**
- 显示数据生成的控制台消息
- 数据加载和处理消息
- 每种方法的预测准确性结果
- 保存的图表文件：`prediction_results.png`

### 理解演示结果

运行 `--demo` 后，您将看到：
- **标准方法**：仅使用基因型的基线预测
- **加权方法**：使用转录组加权 SNP 的改进预测
- **选择方法**：使用顶级方差基因的预测

相关系数（r）表示预测准确性（越高越好，范围：-1 到 1）。

---

## 数据准备

### 所需数据文件

要使用自己的数据运行 MultiOmicsGP，您需要五个 CSV 文件：

1. **基因型文件** (`geno.csv`)
2. **转录组文件** (`trans.csv`)
3. **表型文件** (`pheno.csv`)
4. **SNP 图谱文件** (`snp_map.csv`)
5. **基因图谱文件** (`gene_map.csv`)

### 文件格式规范

#### 1. 基因型文件 (`geno.csv`)

**格式**：CSV 格式，个体 ID 作为行名，SNP ID 作为列名

**编码**：
- `0` = 纯合参考型（AA）
- `1` = 杂合型（Aa）
- `2` = 纯合替代型（aa）

**示例**：
```csv
,SNP_0,SNP_1,SNP_2,SNP_3
Animal_1,0,1,2,0
Animal_2,1,1,0,2
Animal_3,2,0,1,1
```

**要求**：
- 第一列应为个体 ID（将用作索引）
- 缺失值应在分析前处理（插补或删除）
- 所有值必须为 0、1 或 2

**提示**：
- 在所有文件中使用一致的个体 ID 命名
- 确保没有重复的个体 ID
- 运行前检查缺失数据

#### 2. 转录组文件 (`trans.csv`)

**格式**：CSV 格式，个体 ID 作为行名，基因 ID 作为列名

**值**：原始读段计数或标准化表达值

**示例**：
```csv
,Gene_0,Gene_1,Gene_2,Gene_3
Animal_1,1250.5,342.1,89.3,1567.8
Animal_2,1180.2,398.5,95.1,1423.4
Animal_3,1320.7,315.9,102.6,1689.2
```

**要求**：
- 第一列应为个体 ID
- 值应为非负数（计数或标准化表达）
- 缺失值会导致问题 - 请提前处理

**提示**：
- 可以使用原始计数或标准化值（TPM、FPKM 等）
- 软件将自动应用 log2(x+1) 转换
- 并非所有个体都需要转录组数据（软件将同步）

#### 3. 表型文件 (`pheno.csv`)

**格式**：CSV 格式，个体 ID 作为行名，性状名称作为列名

**示例**：
```csv
,CarcassWeight,MilkYield,FeedEfficiency
Animal_1,350.5,8200.3,2.45
Animal_2,385.2,9100.1,2.38
Animal_3,320.8,7800.5,2.52
```

**要求**：
- 第一列应为个体 ID
- 使用 `--trait` 参数指定性状列
- 缺失值将自动删除该个体

**提示**：
- 可以在一个文件中包含多个性状
- 使用描述性的性状名称
- 确保性状值为数值型

#### 4. SNP 图谱文件 (`snp_map.csv`)

**格式**：CSV 格式，三列：SNP_ID、Chromosome、Position

**示例**：
```csv
SNP_ID,Chromosome,Position
SNP_0,Chr1,1000
SNP_1,Chr1,1100
SNP_2,Chr1,1200
SNP_3,Chr2,5000
```

**要求**：
- 必须恰好有这三列
- SNP_ID 必须与基因型文件中的列名匹配
- Position 应为碱基对单位
- 染色体名称应一致

**提示**：
- 使用标准染色体命名（Chr1、Chr2 等或 1、2 等）
- 确保基因型文件中的所有 SNP 在 SNP 图谱中都有条目
- 位置应为数值型

#### 5. 基因图谱文件 (`gene_map.csv`)

**格式**：CSV 格式，四列：Gene_ID、Chromosome、Start、End

**示例**：
```csv
Gene_ID,Chromosome,Start,End
Gene_0,Chr1,1000,2000
Gene_1,Chr1,3000,4000
Gene_2,Chr2,5000,6000
```

**要求**：
- 必须恰好有这四列
- Gene_ID 必须与转录组文件中的列名匹配
- Start 和 End 位置定义基因边界
- 所有基因的 Start < End

**提示**：
- 基因边界应准确
- 允许重叠基因
- 确保转录组文件中的所有基因在基因图谱中都有条目

### 数据同步

软件自动：
- 查找所有三个数据文件（基因型、转录组、表型）中都存在的个体
- 删除任何文件中缺失的个体
- 同步所有数据集的行顺序
- 报告找到多少个共同个体

**重要**：确保所有文件中的个体 ID 一致！

### 数据质量检查

运行分析前，请检查：

1. **个体 ID 一致性**：所有文件中是否使用相同的 ID？
2. **缺失数据**：是否有大量缺失值？
3. **文件格式**：文件是否为正确格式的 CSV？
4. **坐标图谱**：所有 SNP/基因是否都有图谱条目？
5. **数据类型**：数值列是否确实是数值型？

---

## 运行软件

### 基本命令结构

```bash
python main.py [选项]
```

### 常见工作流程

#### 工作流程 1：使用默认设置的快速分析

```bash
python main.py --dir test_data --trait CarcassWeight
```

这将运行：
- 所有三种方法（标准、加权、选择）
- 岭回归（默认 ML 算法）
- 单次训练-测试分割（默认验证）
- 生成可视化

#### 工作流程 2：使用交叉验证比较方法

```bash
python main.py --dir test_data --trait CarcassWeight --method all --cv 5
```

这将运行：
- 所有三种方法
- 5 折交叉验证（更稳健的估计）
- 不生成可视化（CV 模式不生成图表）

#### 工作流程 3：测试不同的 ML 算法

```bash
# 随机森林
python main.py --dir test_data --trait CarcassWeight --ml rf

# 支持向量回归
python main.py --dir test_data --trait CarcassWeight --ml svr

# 多层感知器
python main.py --dir test_data --trait CarcassWeight --ml mlp
```

#### 工作流程 4：使用自定义阈值的特征选择

```bash
python main.py --dir test_data --trait CarcassWeight --method selected --threshold 0.15
```

这将保留前 15% 的基因（而不是默认的 10%）。

### 命令行选项说明

#### 数据输入选项

| 选项 | 说明 | 默认值 | 示例 |
|------|------|--------|------|
| `--dir` | 包含输入文件的目录 | `test_data` | `--dir /path/to/data` |
| `--geno` | 基因型文件名 | `geno.csv` | `--geno my_genotypes.csv` |
| `--trans` | 转录组文件名 | `trans.csv` | `--trans rna_seq.csv` |
| `--pheno` | 表型文件名 | `pheno.csv` | `--pheno traits.csv` |
| `--map_snp` | SNP 位置图谱 | `snp_map.csv` | `--map_snp snp_coords.csv` |
| `--map_gene` | 基因位置图谱 | `gene_map.csv` | `--map_gene gene_coords.csv` |
| `--trait` | 要预测的表型列 | `CarcassWeight` | `--trait MilkYield` |

#### 分析选项

| 选项 | 说明 | 默认值 | 选择 |
|------|------|--------|------|
| `--method` | 要运行的预测方法 | `all` | `standard`、`weighted`、`selected`、`all` |
| `--threshold` | 选择方法保留的顶级基因百分比 | `0.10` | 任何 0.0-1.0 的浮点数 |
| `--ml` | 机器学习算法 | `ridge` | `ridge`、`rf`、`svr`、`mlp` |
| `--cv` | 交叉验证折数 | `1` | 整数 ≥ 1 |

#### 实用选项

| 选项 | 说明 |
|------|------|
| `--demo` | 运行前生成演示数据 |
| `-h`, `--help` | 显示帮助消息 |

### 方法选择指南

#### 标准方法 (`--method standard`)

**何时使用：**
- 基线比较
- 转录组数据质量较差时
- 快速分析

**作用：**
- 仅使用基因型数据
- 传统 GBLUP 方法
- 最快的方法

#### 加权方法 (`--method weighted`)

**何时使用：**
- 您有高质量的转录组数据
- 想要整合生物学信息
- 期望基因表达具有信息性

**作用：**
- 根据关联基因的表达方差对 SNP 进行加权
- 高变异基因中的 SNP 获得更高权重
- 保留所有 SNP（不进行过滤）

#### 选择方法 (`--method selected`)

**何时使用：**
- 想要降低维度
- 专注于最具信息性的基因
- 计算效率很重要

**作用：**
- 选择顶级方差基因中的 SNP
- 过滤掉低方差基因中的 SNP
- 减少特征空间

#### 所有方法 (`--method all`)

**何时使用：**
- 比较所有方法
- 首次分析
- 准备发表的结果

**作用：**
- 运行所有三种方法
- 提供全面比较

### ML 算法选择指南

#### 岭回归 (`--ml ridge`)

**最适合：**
- 线性关系
- 大型数据集
- 基线比较
- 快速计算

**参数：**
- Alpha（正则化）= 10.0

#### 随机森林 (`--ml rf`)

**最适合：**
- 非线性关系
- 特征交互
- 对异常值稳健

**参数：**
- n_estimators = 100
- max_depth = 10

#### 支持向量回归 (`--ml svr`)

**最适合：**
- 非线性模式
- 中小型数据集
- 复杂关系

**参数：**
- Kernel = RBF
- C = 10.0
- epsilon = 0.1

#### 多层感知器 (`--ml mlp`)

**最适合：**
- 深度学习方法
- 复杂非线性模式
- 当其他方法达到平台期时

**参数：**
- 隐藏层：(50, 25)
- max_iter = 500

### 验证策略指南

#### 单次分割 (`--cv 1`)

**使用时机：**
- 想要可视化图表
- 需要快速结果
- 数据集非常大

**工作原理：**
- 70% 训练，30% 测试
- 每个样本一个预测
- 生成散点图

#### 交叉验证 (`--cv 5` 或更高)

**使用时机：**
- 需要稳健的性能估计
- 中小型数据集
- 发表质量的结果

**工作原理：**
- k 折交叉验证
- 跨折的平均性能
- 更可靠的估计
- 不生成可视化（图表太多）

**推荐折数：**
- 小型数据集（<100 样本）：5 折
- 中型数据集（100-500）：5-10 折
- 大型数据集（>500）：10 折

---

## 理解输出结果

### 控制台输出

#### 数据加载阶段

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

**检查要点：**
- 维度是否合理？
- 找到了多少个共同个体？
- 是否有关于缺失数据的警告？

#### 质量控制阶段

```
--- Transcriptome QC & Normalization ---
Removed 2 low-variance genes.
Applying Log2(x + 1) transformation...
Final Transcriptome: 48 Genes
```

**检查要点：**
- 过滤了多少个基因？
- 最终基因数量是否合理？

#### 注释阶段

```
--- Loading Genomic Coordinates ---
Loaded 500 SNP positions and 50 Gene regions.
Mapping SNPs to Genes...
Successfully mapped 450 SNPs out of 500 (90.0%).
```

**检查要点：**
- 映射百分比（通常应 >80%）
- 未映射的 SNP 是否值得关注？

#### 分析阶段

**单次分割模式：**
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

**交叉验证模式：**
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

### 解释结果

#### 相关系数 (r)

- **范围**：-1 到 +1
- **解释**：
  - `r > 0.7`：优秀预测
  - `r = 0.5-0.7`：良好预测
  - `r = 0.3-0.5`：中等预测
  - `r < 0.3`：较差预测
  - `r < 0`：负相关（模型比随机更差）

#### 比较方法

- **标准 vs 加权**：如果加权更高，转录组信息有帮助
- **标准 vs 选择**：如果选择更高，特征选择有益
- **加权 vs 选择**：哪种生物学策略效果更好？

#### 可视化输出

使用 `--cv 1` 时，会生成一个图表文件 `prediction_results.png`，显示：
- **左面板**：标准方法预测值 vs 真实值
- **右面板**：最佳替代方法（加权或选择）vs 真实值
- **相关系数值**：显示在每个图表上
- **趋势线**：红线显示预测准确性

**如何阅读图表：**
- 点越接近对角线 = 预测越好
- 紧密聚集 = 一致的预测
- 分散的点 = 高预测误差

---

## 高级用法

### 自定义数据目录结构

组织您的数据为自定义结构：

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

运行：
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

### 批量处理多个性状

创建脚本处理多个性状：

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

### 比较 ML 算法

测试哪种算法最适合您的数据：

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

### 特征选择敏感性分析

测试不同阈值：

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

### 与其他工具集成

#### 导出结果到 CSV

修改代码或使用 Python 包装器：

```python
import subprocess
import pandas as pd

# 运行分析
result = subprocess.run(
    ['python', 'main.py', '--dir', 'test_data', '--trait', 'CarcassWeight'],
    capture_output=True, text=True
)

# 解析并保存结果
# （您需要解析输出或修改代码以保存结果）
```

---

## 故障排除

### 常见错误和解决方案

#### 错误："No common animals found!"（未找到共同个体！）

**原因**：文件间的个体 ID 不匹配

**解决方案**：
1. 检查所有文件中的个体 ID 格式
2. 确保没有多余的空格或特殊字符
3. 验证 ID 在第一列（索引）
4. 使用一致的命名（例如，全部大写或全部小写）

**检查方法**：
```python
import pandas as pd
geno = pd.read_csv('geno.csv', index_col=0)
trans = pd.read_csv('trans.csv', index_col=0)
pheno = pd.read_csv('pheno.csv', index_col=0)
print("Geno IDs:", geno.index.tolist()[:5])
print("Trans IDs:", trans.index.tolist()[:5])
print("Pheno IDs:", pheno.index.tolist()[:5])
```

#### 错误："Trait 'X' not found in file"（文件中未找到性状 'X'）

**原因**：性状名称与表型文件中的列名不匹配

**解决方案**：
1. 检查确切的列名（区分大小写）
2. 列出可用性状：
   ```python
   import pandas as pd
   pheno = pd.read_csv('pheno.csv', index_col=0)
   print("Available traits:", pheno.columns.tolist())
   ```
3. 使用确切的列名和 `--trait`

#### 错误："Failed to load data"（加载数据失败）

**原因**：文件路径问题或格式问题

**解决方案**：
1. 检查文件路径是否正确
2. 验证指定目录中存在文件
3. 确保文件是有效的 CSV 格式
4. 检查文件权限

#### 警告："No SNPs selected! Returning full matrix"（未选择 SNP！返回完整矩阵）

**原因**：特征选择阈值过于严格或映射问题

**解决方案**：
1. 增加 `--threshold` 值
2. 检查 SNP 到基因的映射（单独运行注释步骤）
3. 验证转录组数据有方差

#### 预测结果较差（r 值低）

**可能原因和解决方案**：

1. **样本量小**
   - 解决方案：需要更多个体（通常 >100 才能获得良好结果）

2. **低遗传力性状**
   - 解决方案：某些性状本质上难以预测

3. **数据质量问题**
   - 解决方案：检查基因分型错误、缺失数据

4. **错误的 ML 算法**
   - 解决方案：尝试不同的算法（`--ml rf`、`--ml svr`）

5. **性状不适合基因组预测**
   - 解决方案：某些性状可能无法从基因型预测

#### 内存问题

**症状**：程序崩溃或运行非常缓慢

**解决方案**：
1. 减少数据集大小（更少的 SNP/基因）
2. 使用特征选择方法
3. 增加系统 RAM
4. 分批处理

#### 未生成可视化

**原因**：在交叉验证模式下运行（`--cv > 1`）

**解决方案**：使用 `--cv 1` 进行单次分割模式以生成图表

### 调试技巧

#### 1. 检查数据维度

```python
import pandas as pd
geno = pd.read_csv('test_data/geno.csv', index_col=0)
trans = pd.read_csv('test_data/trans.csv', index_col=0)
pheno = pd.read_csv('test_data/pheno.csv', index_col=0)

print(f"Genotype: {geno.shape}")
print(f"Transcriptome: {trans.shape}")
print(f"Phenotype: {pheno.shape}")
```

#### 2. 验证数据类型

```python
print("Genotype dtypes:", geno.dtypes.value_counts())
print("Transcriptome dtypes:", trans.dtypes.value_counts())
print("Phenotype dtypes:", pheno.dtypes.value_counts())
```

#### 3. 检查缺失数据

```python
print("Missing in genotype:", geno.isnull().sum().sum())
print("Missing in transcriptome:", trans.isnull().sum().sum())
print("Missing in phenotype:", pheno.isnull().sum().sum())
```

#### 4. 测试单个组件

使用 `main_pipeline.py` 逐步运行流程，或修改 `main.py` 添加调试打印。

---

## 常见问题

### 一般问题

**问：我可以将此用于人类数据吗？**
答：可以，只要您有所需的数据格式，该软件适用于任何物种。

**问：我需要所有三种数据类型（基因型、转录组、表型）吗？**
答：是的，所有三种都是必需的。软件整合它们以提高预测。

**问：我可以使用插补的基因型吗？**
答：可以，只要它们编码为 0/1/2 并且格式正确。

**问：如果我没有所有个体的转录组数据怎么办？**
答：软件将自动仅使用具有所有三种数据类型的个体。

**问：分析需要多长时间？**
答：取决于数据集大小：
- 小型（<100 个体，<1000 SNP）：秒级
- 中型（100-500 个体，1000-10000 SNP）：分钟级
- 大型（>500 个体，>10000 SNP）：小时级

### 数据问题

**问：转录组数据应该是什么格式？**
答：原始计数或标准化值（TPM、FPKM）都可以。软件应用 log2 转换。

**问：我可以使用不同的染色体命名约定吗？**
答：可以，但要保持一致。在 SNP 和基因图谱中使用相同的格式。

**问：如果 SNP 未映射到基因怎么办？**
答：未映射的 SNP 在加权方法中获得默认权重 1.0，在选择方法中被排除。

**问：我可以使用分型基因型吗？**
答：软件期望未分型基因型（0/1/2 编码）。分型数据需要转换。

### 方法问题

**问：我应该使用哪种方法？**
答：从 `--method all` 开始比较所有方法，然后为您的数据选择最佳方法。

**问：何时应该使用特征选择？**
答：当您有许多 SNP 并想专注于最具信息性的 SNP 时，或为了计算效率。

**问：加权和选择有什么区别？**
答：加权保留所有 SNP 但调整其重要性。选择过滤为仅顶级基因 SNP。

**问：哪种 ML 算法最好？**
答：取决于您的数据。尝试所有并比较。岭回归通常是一个好的起点。

### 技术问题

**问：我可以在集群上运行这个吗？**
答：可以，这是一个可以作为作业提交的命令行工具。

**问：如何并行化分析？**
答：您可以使用作业数组或 Python 多处理并行运行不同的性状或方法。

**问：我可以修改代码吗？**
答：可以，它是开源的。根据需要修改以满足您的特定要求。

**问：如何引用此软件？**
答：请参阅 README 文件了解引用信息。

### 结果问题

**问：什么是好的预测准确性？**
答：取决于性状遗传力。对于高遗传力性状，r > 0.5 是好的。对于低遗传力，r > 0.3 可能是可接受的。

**问：为什么我的结果每次都不同？**
答：如果使用随机分割或带打乱的交叉验证，结果会略有不同。使用 `random_state` 以获得可重复性。

**问：我可以获得单个个体的预测值吗？**
答：目前，软件报告相关性。您可以修改代码以保存单个预测。

**问：如何解释负相关？**
答：负相关表明模型表现比随机更差。检查数据质量和方法选择。

---

## 附录

### A. 文件格式示例

#### 完整基因型文件示例
```csv
,SNP_0,SNP_1,SNP_2,SNP_3,SNP_4
Cattle_1,0,1,2,0,1
Cattle_2,1,1,0,2,0
Cattle_3,2,0,1,1,2
Cattle_4,0,2,0,1,1
Cattle_5,1,1,2,0,2
```

#### 完整转录组文件示例
```csv
,Gene_0,Gene_1,Gene_2,Gene_3,Gene_4
Cattle_1,1250.5,342.1,89.3,1567.8,234.5
Cattle_2,1180.2,398.5,95.1,1423.4,267.8
Cattle_3,1320.7,315.9,102.6,1689.2,198.3
Cattle_4,1105.3,421.2,87.9,1512.6,289.1
Cattle_5,1289.6,356.7,94.5,1623.9,245.7
```

#### 完整表型文件示例
```csv
,CarcassWeight,MilkYield,FeedEfficiency
Cattle_1,350.5,8200.3,2.45
Cattle_2,385.2,9100.1,2.38
Cattle_3,320.8,7800.5,2.52
Cattle_4,365.7,8750.2,2.41
Cattle_5,340.2,8050.8,2.49
```

### B. 命令参考卡

```bash
# 基本用法
python main.py --dir DATA_DIR --trait TRAIT_NAME

# 完整选项
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

# 快速参考
python main.py -h  # 显示所有选项
```

### C. 性能基准

标准笔记本电脑上的典型运行时间（8GB RAM，4 核）：

| 数据集大小 | 运行时间（单次分割） | 运行时间（5 折 CV） |
|------------|---------------------|-------------------|
| 100 个体，1K SNP | < 1 分钟 | 2-3 分钟 |
| 500 个体，10K SNP | 5-10 分钟 | 20-30 分钟 |
| 1000 个体，50K SNP | 30-60 分钟 | 2-3 小时 |

### D. 术语表

- **GBLUP**：基因组最佳线性无偏预测（Genomic Best Linear Unbiased Prediction）
- **MAF**：次要等位基因频率（Minor Allele Frequency）
- **SNP**：单核苷酸多态性（Single Nucleotide Polymorphism）
- **QC**：质量控制（Quality Control）
- **CV**：交叉验证（Cross-Validation）
- **ML**：机器学习（Machine Learning）
- **Ridge**：岭回归（Ridge Regression）
- **RF**：随机森林（Random Forest）
- **SVR**：支持向量回归（Support Vector Regression）
- **MLP**：多层感知器（Multi-Layer Perceptron）

### E. 其他资源

- **文档**：查看 README.md 了解概述
- **代码仓库**：查看 GitHub 获取最新更新
- **问题报告**：在 GitHub issues 页面报告错误
- **联系方式**：higefee@gmail.com

---

**用户手册结束**

*最后更新：版本 2.2*

