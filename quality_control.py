import pandas as pd
import numpy as np

class QualityControl:
    def __init__(self, genotype_df, transcriptome_df):
        """
        Input: Synchronized DataFrames from the DataLoader.
        """
        self.genotype = genotype_df
        self.transcriptome = transcriptome_df

    def filter_genotypes(self, maf_threshold=0.05):
        """
        QC 1: Remove SNPs with low Minor Allele Frequency (MAF).
        Formula: MAF = min(p, q)
        """
        print(f"\n--- Genotype QC (MAF > {maf_threshold}) ---")
        original_count = self.genotype.shape[1]
        
        # Calculate Frequency of the coded allele (assuming 0, 1, 2 coding)
        # Mean dosage / 2 gives the frequency (p)
        freq = self.genotype.mean(axis=0) / 2
        
        # Calculate MAF (if p > 0.5, then q is the minor allele)
        maf = np.where(freq > 0.5, 1 - freq, freq)
        
        # Identify SNPs to keep
        # We also ensure MAF is not 0 (monomorphic)
        keep_snps = maf > maf_threshold
        
        self.genotype = self.genotype.loc[:, keep_snps]
        
        dropped = original_count - self.genotype.shape[1]
        print(f"Removed {dropped} SNPs. Remaining: {self.genotype.shape[1]}")
        
        return self.genotype

    def process_transcriptome(self, log_transform=True, min_variance=0.01):
        """
        QC 2: Process RNA-Seq Data.
        1. Remove genes with zero variance (constant expression).
        2. Log2 Transform: log2(counts + 1).
        """
        print("\n--- Transcriptome QC & Normalization ---")
        original_count = self.transcriptome.shape[1]

        # 1. Filter Low Variance Genes (Noise)
        variances = self.transcriptome.var(axis=0)
        keep_genes = variances > min_variance
        self.transcriptome = self.transcriptome.loc[:, keep_genes]
        
        dropped = original_count - self.transcriptome.shape[1]
        print(f"Removed {dropped} low-variance genes.")

        # 2. Log Transformation
        # We add +1 to avoid log(0) errors.
        if log_transform:
            print("Applying Log2(x + 1) transformation...")
            self.transcriptome = np.log2(self.transcriptome + 1)

        print(f"Final Transcriptome: {self.transcriptome.shape[1]} Genes")
        return self.transcriptome