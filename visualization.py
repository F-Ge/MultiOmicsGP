import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class Visualizer:
    def __init__(self):
        # Set the aesthetic style of the plots
        sns.set_style("whitegrid")

    def plot_comparison(self, y_true, y_pred_base, y_pred_weighted, filename="prediction_results.png"):
        """
        Plots two scatter plots side-by-side:
        1. Standard GBLUP
        2. Weighted GBLUP
        """
        plt.figure(figsize=(14, 6))

        # --- Plot 1: Standard GBLUP ---
        plt.subplot(1, 2, 1)
        # Scatter points
        sns.scatterplot(x=y_true, y=y_pred_base, color='grey', alpha=0.6, s=80)
        
        # Add trend line (Red)
        m, b = np.polyfit(y_true, y_pred_base, 1)
        plt.plot(y_true, m*y_true + b, color='red', linewidth=2)
        
        # Labels
        plt.title("Method 1: Standard G-Only", fontsize=14, fontweight='bold')
        plt.xlabel("True Phenotype (kg)", fontsize=12)
        plt.ylabel("Predicted Phenotype (kg)", fontsize=12)
        
        # Calculate Correlation for label
        r_base = np.corrcoef(y_true, y_pred_base)[0,1]
        plt.text(0.05, 0.95, f'r = {r_base:.3f}', transform=plt.gca().transAxes, 
                 fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        # --- Plot 2: Weighted GBLUP ---
        plt.subplot(1, 2, 2)
        # Scatter points (Blue for better method)
        sns.scatterplot(x=y_true, y=y_pred_weighted, color='dodgerblue', alpha=0.6, s=80)
        
        # Add trend line (Red)
        m, b = np.polyfit(y_true, y_pred_weighted, 1)
        plt.plot(y_true, m*y_true + b, color='darkred', linewidth=2)
        
        # Labels
        plt.title("Method 2: Transcriptome-Weighted", fontsize=14, fontweight='bold')
        plt.xlabel("True Phenotype (kg)", fontsize=12)
        plt.ylabel("Predicted Phenotype (kg)", fontsize=12)
        
        # Calculate Correlation
        r_weight = np.corrcoef(y_true, y_pred_weighted)[0,1]
        plt.text(0.05, 0.95, f'r = {r_weight:.3f}', transform=plt.gca().transAxes, 
                 fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

        plt.tight_layout()
        plt.savefig(filename)
        print(f"\n[Success] Plot saved to {filename}")
        plt.show()