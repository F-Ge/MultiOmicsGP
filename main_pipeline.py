import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge

from data_loader import DataLoader
from quality_control import QualityControl
from models import GenomicPredictor
from visualization import Visualizer
from annotation import GenomicAnnotator # <--- NEW IMPORT

# --- 1. SETUP & DATA GENERATION (With Coordinates) ---
def generate_mock_data():
    if True: 
        if not os.path.exists("test_data"):
            os.makedirs("test_data")
            
        animals = [f"Cattle_{i}" for i in range(1, 101)]
        n_snps = 500
        n_genes = 50
        
        # A. Create Genotypes & Phenotypes (Same as before)
        G = pd.DataFrame(np.random.randint(0, 3, size=(100, n_snps)), 
                         index=animals, columns=[f"SNP_{i}" for i in range(n_snps)])
        
        # Biology: First 5 genes are causal (SNPs 0-49)
        causal_snps = [f"SNP_{i}" for i in range(50)] 
        signal = np.sum(G[causal_snps].values, axis=1) * 2.0 
        noise = np.random.normal(0, 10, 100)
        P = pd.DataFrame(350 + signal + noise, index=animals, columns=["CarcassWeight"])
        
        # B. Create Transcriptome
        T_values = np.zeros((100, n_genes))
        T_values[:, 0:5] = np.random.normal(10, 5, size=(100, 5)) # Causal = High Var
        T_values[:, 5:] = np.random.normal(10, 0.5, size=(100, 45)) # Noise = Low Var
        T = pd.DataFrame(T_values[:80], index=animals[:80], columns=[f"Gene_{i}" for i in range(n_genes)])
        
        # --- C. GENERATE COORDINATE MAPS (NEW) ---
        # 1. Gene Map: We place genes along Chromosome 1
        # Each gene is 1000bp long, spaced 1000bp apart
        gene_rows = []
        for i in range(n_genes):
            start = i * 2000 + 1000 # 1000, 3000, 5000...
            end = start + 1000      # 2000, 4000, 6000...
            gene_rows.append([f"Gene_{i}", "Chr1", start, end])
        
        gene_map_df = pd.DataFrame(gene_rows, columns=["Gene_ID", "Chromosome", "Start", "End"])
        
        # 2. SNP Map: We place 10 SNPs inside each gene
        snp_rows = []
        snp_idx = 0
        for i in range(n_genes):
            # Get the start of the current gene
            gene_start = gene_rows[i][2]
            for j in range(10): # 10 SNPs per gene
                pos = gene_start + (j * 100) # Spread them out inside the gene
                snp_rows.append([f"SNP_{snp_idx}", "Chr1", pos])
                snp_idx += 1
                
        snp_map_df = pd.DataFrame(snp_rows, columns=["SNP_ID", "Chromosome", "Position"])

        # Save everything
        G.to_csv("test_data/geno.csv")
        T.to_csv("test_data/trans.csv")
        P.to_csv("test_data/pheno.csv")
        gene_map_df.to_csv("test_data/gene_map.csv", index=False)
        snp_map_df.to_csv("test_data/snp_map.csv", index=False)
        
        print("Biologically Realistic Cattle Data + Maps Generated.")

# --- 2. RUN THE PIPELINE ---
if __name__ == "__main__":
    generate_mock_data()
    
    # 1. Load Data
    loader = DataLoader("test_data")
    loader.load_genotype("geno.csv")
    loader.load_transcriptome("trans.csv")
    loader.load_phenotype("pheno.csv", target_column="CarcassWeight")
    X_clean, T_clean, y = loader.synchronize_data()
    
    # 2. QC (Optional: QC skips for now as data is clean)
    qc = QualityControl(X_clean, T_clean)
    T_clean = qc.process_transcriptome(log_transform=True)

    # 3. ANNOTATION (NEW STEP)
    # This replaces the "guesswork" mapping
    annotator = GenomicAnnotator()
    annotator.load_positions("test_data/snp_map.csv", "test_data/gene_map.csv")
    real_mapping = annotator.map_snps_to_genes(buffer_bp=0)

    # 4. Modeling
    print("\n--- Starting Genomic Prediction ---")
    predictor = GenomicPredictor()
    X_train, X_test, y_train, y_test = train_test_split(X_clean, y, test_size=0.3, random_state=42)
    y_test_flat = y_test.values.ravel()

    # Method 1: Standard
    model1 = Ridge(alpha=10.0)
    model1.fit(X_train, y_train)
    acc1 = np.corrcoef(y_test_flat, model1.predict(X_test).ravel())[0,1]

    # Method 2: Weighted (Pass the real mapping!)
    X_weighted = predictor.apply_bio_weights(X_clean, T_clean, real_mapping)
    Xw_train = X_weighted.loc[X_train.index]
    Xw_test = X_weighted.loc[X_test.index]
    
    model2 = Ridge(alpha=10.0)
    model2.fit(Xw_train, y_train)
    acc2 = np.corrcoef(y_test_flat, model2.predict(Xw_test).ravel())[0,1]

    # Method 3: Selected (Pass the real mapping!)
    X_selected = predictor.select_top_features(X_clean, T_clean, real_mapping, top_percent=0.10)
    Xs_train = X_selected.loc[X_train.index]
    Xs_test = X_selected.loc[X_test.index]
    
    model3 = Ridge(alpha=10.0)
    model3.fit(Xs_train, y_train)
    acc3 = np.corrcoef(y_test_flat, model3.predict(Xs_test).ravel())[0,1]

    print("\n================ RESULTS ================")
    print(f"Standard: r = {acc1:.4f}")
    print(f"Weighted: r = {acc2:.4f}")
    print(f"Selected: r = {acc3:.4f}")
    print("=========================================")
    
    # Visualization
    viz = Visualizer()
    viz.plot_comparison(y_test_flat, model1.predict(X_test).ravel(), model3.predict(Xs_test).ravel())