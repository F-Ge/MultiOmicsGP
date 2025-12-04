import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor

class GenomicPredictor:
    def __init__(self):
        pass

    def get_model(self, model_name="ridge"):
        """
        Factory to return the correct ML model object.
        """
        model_name = model_name.lower()
        
        if model_name == "ridge":
            # Standard GBLUP
            return Ridge(alpha=10.0)
            
        elif model_name == "rf":
            # Random Forest (Captures non-linear interactions)
            return RandomForestRegressor(n_estimators=100, max_depth=10, random_state=42)
            
        elif model_name == "svr":
            # Support Vector Regression
            return SVR(kernel='rbf', C=10.0, epsilon=0.1)
            
        elif model_name == "mlp":
            # Multi-Layer Perceptron (Deep Learning)
            # Two hidden layers with 50 and 25 neurons
            return MLPRegressor(hidden_layer_sizes=(50, 25), max_iter=500, random_state=42)
            
        else:
            raise ValueError(f"Unknown model: {model_name}")

    def apply_bio_weights(self, X, T, mapping_dict):
        """
        Strategy B: Weighted GBLUP (Using Real Mapping)
        """
        print("\n--- Calculating Biological Weights (Real Map) ---")
        
        # 1. Gene Importance
        gene_weights = T.std(axis=0)
        gene_weights = gene_weights / gene_weights.mean()
        
        # 2. Create SNP Weight Vector
        snp_weights = []
        default_weight = 1.0 
        
        for snp in X.columns:
            if snp in mapping_dict:
                gene_id = mapping_dict[snp]
                if gene_id in gene_weights:
                    snp_weights.append(gene_weights[gene_id])
                else:
                    snp_weights.append(default_weight) 
            else:
                snp_weights.append(default_weight) 
                
        snp_weights = np.array(snp_weights)
        
        print(f"Applying weights: Min={snp_weights.min():.2f}, Max={snp_weights.max():.2f}")
        X_weighted = X * snp_weights
        return X_weighted

    def select_top_features(self, X, T, mapping_dict, top_percent=0.1):
        """
        Strategy C: Feature Selection (Using Real Mapping)
        """
        print(f"\n--- Strategy C: Feature Selection (Top {top_percent*100}%) ---")
        
        # 1. Rank Genes
        gene_variances = T.var(axis=0)
        sorted_genes = gene_variances.sort_values(ascending=False)
        num_keep = int(len(sorted_genes) * top_percent)
        if num_keep < 1: num_keep = 1
        top_genes = set(sorted_genes.index[:num_keep]) 
        
        print(f"Selected {len(top_genes)} 'Elite' Genes.")
        
        # 2. Filter SNPs
        snps_to_keep = []
        for snp in X.columns:
            if snp in mapping_dict:
                gene_id = mapping_dict[snp]
                if gene_id in top_genes:
                    snps_to_keep.append(snp)
        
        print(f"Filtering: Keeping {len(snps_to_keep)} SNPs out of {X.shape[1]}.")
        
        if len(snps_to_keep) == 0:
            print("WARNING: No SNPs selected! Returning full matrix.")
            return X
            
        return X[snps_to_keep]