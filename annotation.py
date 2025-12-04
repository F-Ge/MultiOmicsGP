import pandas as pd
import numpy as np

class GenomicAnnotator:
    def __init__(self):
        self.snp_map = None
        self.gene_map = None
        self.mapping = None

    def load_positions(self, snp_file, gene_file):
        """
        Load coordinate files.
        snp_file columns: [SNP_ID, Chromosome, Position]
        gene_file columns: [Gene_ID, Chromosome, Start, End]
        """
        print("\n--- Loading Genomic Coordinates ---")
        self.snp_map = pd.read_csv(snp_file)
        self.gene_map = pd.read_csv(gene_file)
        print(f"Loaded {len(self.snp_map)} SNP positions and {len(self.gene_map)} Gene regions.")

    def map_snps_to_genes(self, buffer_bp=0):
        """
        The Core Logic:
        Finds which Gene each SNP belongs to based on physical location.
        buffer_bp: Optional distance (e.g., 2000bp) to include promoter regions.
        """
        print("Mapping SNPs to Genes...")
        
        # Dictionary to store result: {SNP_ID: Gene_ID}
        snp_to_gene = {}
        
        # Optimization: Group by Chromosome to speed up search
        # (Comparing Chr1 SNPs only to Chr1 Genes)
        unique_chrs = self.snp_map['Chromosome'].unique()
        
        mapped_count = 0
        
        for chrom in unique_chrs:
            # 1. Get SNPs and Genes on this chromosome
            chr_snps = self.snp_map[self.snp_map['Chromosome'] == chrom]
            chr_genes = self.gene_map[self.gene_map['Chromosome'] == chrom]
            
            # 2. Iterate through genes on this chrom
            for _, gene in chr_genes.iterrows():
                gene_id = gene['Gene_ID']
                start = gene['Start'] - buffer_bp
                end = gene['End'] + buffer_bp
                
                # 3. Find SNPs inside this window
                # Logic: SNP > Start AND SNP < End
                mask = (chr_snps['Position'] >= start) & (chr_snps['Position'] <= end)
                found_snps = chr_snps[mask]['SNP_ID'].tolist()
                
                # 4. Save to dictionary
                for snp in found_snps:
                    snp_to_gene[snp] = gene_id
                    mapped_count += 1
                    
        self.mapping = snp_to_gene
        
        total_snps = len(self.snp_map)
        print(f"Successfully mapped {mapped_count} SNPs out of {total_snps} ({mapped_count/total_snps:.1%}).")
        
        return self.mapping