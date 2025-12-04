import argparse
import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold

# Import your modules
from data_loader import DataLoader
from quality_control import QualityControl
from models import GenomicPredictor
from annotation import GenomicAnnotator
from visualization import Visualizer

# --- SIMULATION FUNCTION ---
def generate_mock_data(output_dir):
    print(f"\n[DEMO MODE] Generating mock cattle data in '{output_dir}'...")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    animals = [f"Cattle_{i}" for i in range(1, 101)]
    n_snps = 500
    n_genes = 50
    
    # 1. Genotypes & Phenotypes
    G = pd.DataFrame(np.random.randint(0, 3, size=(100, n_snps)), 
                     index=animals, columns=[f"SNP_{i}" for i in range(n_snps)])
    
    # Biology: First 5 genes are causal
    causal_snps = [f"SNP_{i}" for i in range(50)] 
    signal = np.sum(G[causal_snps].values, axis=1) * 2.0 
    noise = np.random.normal(0, 10, 100)
    P = pd.DataFrame(350 + signal + noise, index=animals, columns=["CarcassWeight"])
    
    # 2. Transcriptome
    T_values = np.zeros((100, n_genes))
    T_values[:, 0:5] = np.random.normal(10, 5, size=(100, 5)) 
    T_values[:, 5:] = np.random.normal(10, 0.5, size=(100, 45)) 
    T = pd.DataFrame(T_values[:80], index=animals[:80], columns=[f"Gene_{i}" for i in range(n_genes)])
    
    # 3. Maps
    gene_rows = [[f"Gene_{i}", "Chr1", i*2000+1000, i*2000+2000] for i in range(n_genes)]
    gene_map_df = pd.DataFrame(gene_rows, columns=["Gene_ID", "Chromosome", "Start", "End"])
    
    snp_rows = []
    snp_idx = 0
    for i in range(n_genes):
        gene_start = gene_rows[i][2]
        for j in range(10): 
            snp_rows.append([f"SNP_{snp_idx}", "Chr1", gene_start + (j * 100)])
            snp_idx += 1
    snp_map_df = pd.DataFrame(snp_rows, columns=["SNP_ID", "Chromosome", "Position"])

    # Save
    G.to_csv(os.path.join(output_dir, "geno.csv"))
    T.to_csv(os.path.join(output_dir, "trans.csv"))
    P.to_csv(os.path.join(output_dir, "pheno.csv"))
    gene_map_df.to_csv(os.path.join(output_dir, "gene_map.csv"), index=False)
    snp_map_df.to_csv(os.path.join(output_dir, "snp_map.csv"), index=False)
    print("Mock data generated successfully.\n")

# --- MAIN CLI LOGIC ---
def main():
    parser = argparse.ArgumentParser(description="Multi-Omics Genomic Prediction Software v2.2")
    
    # 1. Data Arguments
    parser.add_argument("--dir", type=str, default="test_data", help="Directory containing input files")
    parser.add_argument("--geno", type=str, default="geno.csv", help="Genotype filename")
    parser.add_argument("--trans", type=str, default="trans.csv", help="Transcriptome filename")
    parser.add_argument("--pheno", type=str, default="pheno.csv", help="Phenotype filename")
    parser.add_argument("--map_snp", type=str, default="snp_map.csv", help="SNP Position Map")
    parser.add_argument("--map_gene", type=str, default="gene_map.csv", help="Gene Position Map")
    parser.add_argument("--trait", type=str, default="CarcassWeight", help="Phenotype column to predict")

    # 2. Analysis Settings
    parser.add_argument("--method", type=str, choices=['standard', 'weighted', 'selected', 'all'], 
                        default='all', help="Analysis method to run")
    # FIX: Changed '%' to '%%' in the help string below
    parser.add_argument("--threshold", type=float, default=0.10, help="Top %% of genes to keep (for 'selected' method)")
    
    # 3. ML Settings
    parser.add_argument("--ml", type=str, default="ridge", choices=['ridge', 'rf', 'svr', 'mlp'],
                        help="Machine Learning Algorithm")
    
    # 4. Validation Settings
    parser.add_argument("--cv", type=int, default=1, help="Number of Cross-Validation folds (1 = Single Split)")
    
    parser.add_argument("--demo", action="store_true", help="Generate demo data before running")
    
    args = parser.parse_args()

    # --- Step 0: Demo Mode ---
    if args.demo:
        generate_mock_data(args.dir)

    # --- Step 1: Load Data ---
    print(f"Loading data from '{args.dir}'...")
    try:
        loader = DataLoader(args.dir)
        loader.load_genotype(args.geno)
        loader.load_transcriptome(args.trans)
        loader.load_phenotype(args.pheno, target_column=args.trait)
        X_clean, T_clean, y = loader.synchronize_data()
    except Exception as e:
        print(f"[Error] Failed to load data: {e}")
        sys.exit(1)

    # --- Step 2: Quality Control ---
    qc = QualityControl(X_clean, T_clean)
    T_clean = qc.process_transcriptome(log_transform=True)

    # --- Step 3: Annotation ---
    print("\n[Annotation] Mapping SNPs to Genes...")
    annotator = GenomicAnnotator()
    snp_path = os.path.join(args.dir, args.map_snp)
    gene_path = os.path.join(args.dir, args.map_gene)
    
    annotator.load_positions(snp_path, gene_path)
    real_mapping = annotator.map_snps_to_genes()

    # --- Step 4: Modeling ---
    print(f"\n[Analysis] Running Method: {args.method.upper()}")
    print(f"[Model] Algorithm: {args.ml.upper()}")
    
    predictor = GenomicPredictor()
    results = {"True": []}

    # Helper function that handles EITHER CV or Single Split
    def evaluate_model(X_data, name):
        print(f"  -> Evaluating {name}...")
        y_vec = y.values.ravel()
        
        # A. CROSS-VALIDATION MODE
        if args.cv > 1:
            kf = KFold(n_splits=args.cv, shuffle=True, random_state=42)
            scores = []
            print(f"     Running {args.cv}-Fold CV...", end="")
            
            for train_idx, test_idx in kf.split(X_data):
                X_tr, X_te = X_data.iloc[train_idx], X_data.iloc[test_idx]
                y_tr, y_te = y_vec[train_idx], y_vec[test_idx]
                
                model = predictor.get_model(args.ml)
                model.fit(X_tr, y_tr)
                preds = model.predict(X_te)
                
                r = np.corrcoef(y_te, preds)[0,1]
                if np.isnan(r): r = 0.0
                scores.append(r)
                
            avg_r = np.mean(scores)
            print(f" Done. Avg r = {avg_r:.4f}")
            return None

        # B. SINGLE SPLIT MODE (For Plotting)
        else:
            X_tr, X_te, y_tr, y_te = train_test_split(X_data, y, test_size=0.3, random_state=42)
            
            model = predictor.get_model(args.ml)
            model.fit(X_tr, y_tr.values.ravel())
            preds = model.predict(X_te).ravel()
            
            acc = np.corrcoef(y_te.values.ravel(), preds)[0,1]
            if np.isnan(acc): acc = 0.0
            print(f"     Single Split Accuracy (r): {acc:.4f}")
            
            # Save for plotting
            if len(results["True"]) == 0:
                results["True"] = y_te.values.ravel()
            return preds

    # --- EXECUTE ---
    
    # 1. Standard
    if args.method in ['standard', 'all']:
        preds = evaluate_model(X_clean, "Standard (G-Only)")
        if preds is not None: results["Standard"] = preds

    # 2. Weighted
    if args.method in ['weighted', 'all']:
        X_weighted = predictor.apply_bio_weights(X_clean, T_clean, real_mapping)
        preds = evaluate_model(X_weighted, "Weighted")
        if preds is not None: results["Weighted"] = preds

    # 3. Selected
    if args.method in ['selected', 'all']:
        X_selected = predictor.select_top_features(X_clean, T_clean, real_mapping, top_percent=args.threshold)
        preds = evaluate_model(X_selected, f"Selected (Top {args.threshold:.0%})")
        if preds is not None: results["Selected"] = preds

    # --- Step 5: Visualization ---
    if args.cv == 1 and len(results) > 1:
        print("\n[Visualization] Generating plot...")
        viz = Visualizer()
        keys = list(results.keys())
        
        comp_method = "Standard"
        if "Selected" in keys: comp_method = "Selected"
        elif "Weighted" in keys: comp_method = "Weighted"
            
        if comp_method != "Standard":
            viz.plot_comparison(results["True"], results.get("Standard", results["True"]), results[comp_method])

if __name__ == "__main__":
    main()