import pandas as pd
import numpy as np
import os

class DataLoader:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.genotype = None      # Matrix X
        self.transcriptome = None # Matrix T
        self.phenotype = None     # Vector y
        self.common_ids = None    # List of animals present in ALL files

    def load_genotype(self, filename, format='csv'):
        path = os.path.join(self.data_dir, filename)
        print(f"Loading Genotypes from {path}...")
        self.genotype = pd.read_csv(path, index_col=0)
        print(f"Loaded Genotypes: {self.genotype.shape} (Animals x SNPs)")

    def load_transcriptome(self, filename):
        path = os.path.join(self.data_dir, filename)
        print(f"Loading Transcriptome from {path}...")
        self.transcriptome = pd.read_csv(path, index_col=0)
        print(f"Loaded Transcriptome: {self.transcriptome.shape} (Animals x Genes)")

    def load_phenotype(self, filename, target_column):
        path = os.path.join(self.data_dir, filename)
        print(f"Loading Phenotypes from {path}...")
        full_pheno = pd.read_csv(path, index_col=0)
        
        if target_column in full_pheno.columns:
            self.phenotype = full_pheno[[target_column]]
            self.phenotype = self.phenotype.dropna()
        else:
            raise ValueError(f"Trait '{target_column}' not found in file.")
            
        print(f"Loaded Phenotypes: {self.phenotype.shape} (Animals)")

    def synchronize_data(self):
        print("\n--- Synchronizing Data ---")
        g_ids = set(self.genotype.index)
        t_ids = set(self.transcriptome.index)
        p_ids = set(self.phenotype.index)
        
        # Find intersection
        self.common_ids = list(g_ids.intersection(t_ids).intersection(p_ids))
        
        if not self.common_ids:
            raise ValueError("No common animals found!")
            
        print(f"Found {len(self.common_ids)} common animals.")
        self.common_ids.sort()
        
        # Subset
        self.genotype = self.genotype.loc[self.common_ids]
        self.transcriptome = self.transcriptome.loc[self.common_ids]
        self.phenotype = self.phenotype.loc[self.common_ids]
        
        print("Data synchronized successfully.")
        return self.genotype, self.transcriptome, self.phenotype

# --- UPDATED TEST DATA GENERATOR ---
if __name__ == "__main__":
    # 1. Define IDs as Cattle
    animals = [f"Cattle_{i}" for i in range(1, 101)]
    
    # 2. Simulate Genotypes
    G = pd.DataFrame(np.random.randint(0, 3, size=(100, 500)), index=animals, columns=[f"SNP_{i}" for i in range(500)])
    
    # 3. Simulate Transcriptome (Only 80 cattle have RNA-seq)
    T_animals = animals[:80] 
    T = pd.DataFrame(np.random.rand(80, 50), index=T_animals, columns=[f"Gene_{i}" for i in range(50)])
    
    # 4. Simulate Phenotypes (Carcass Weight)
    # Simulate weights around 350kg with some variation
    weights = np.random.normal(350, 30, 100) 
    P = pd.DataFrame(weights, index=animals, columns=["CarcassWeight"])
    
    # Introduce missing data for Cattle #5
    P.iloc[5] = np.nan 
    
    # Save files
    os.makedirs("test_data", exist_ok=True)
    G.to_csv("test_data/geno.csv")
    T.to_csv("test_data/trans.csv")
    P.to_csv("test_data/pheno.csv")
    
    # --- RUN THE LOADER ---
    loader = DataLoader("test_data")
    loader.load_genotype("geno.csv")
    loader.load_transcriptome("trans.csv")
    
    # Update target column to CarcassWeight
    loader.load_phenotype("pheno.csv", target_column="CarcassWeight")
    
    X, T, y = loader.synchronize_data()
    
    print("\nFinal Dimensions:")
    print(f"Genotype: {X.shape}")
    print(f"Transcriptome: {T.shape}")
    print(f"Phenotype: {y.shape}")
    print("\nFirst 3 rows of Phenotype Data:")
    print(y.head(3))