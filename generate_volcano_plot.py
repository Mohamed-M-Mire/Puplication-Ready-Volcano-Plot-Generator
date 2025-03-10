import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import text_annotate as ta

def generate_volcano_plot(DEA_res_noDubs: pd.DataFrame, lfcutoff: float, fdrcutoff: float, figureTitle: str, savefig_Path: str = None, maxGeneNameLen: int = 20):
    """
    Generate a volcano plot based on input data and thresholds.

    Parameters:
    DEA_res_noDubs (pd.DataFrame): DataFrame containing differential expression analysis results.
    lfcutoff (float): Log fold change cutoff value.
    fdrcutoff (float): False discovery rate cutoff value.
    figureTitle (str): Title for the generated plot.
    savefig_Path (str, optional): File path to save the plot. Defaults to None.
    maxGeneNameLen (int, optional): Maximum length of gene names to display. Defaults to 20.

    Returns:
    None

    Example usage:
    generate_volcano_plot(DEA_res_noDubs, lfcutoff=1.5, fdrcutoff=0.05, figureTitle="Differential Expression Volcano Plot", savefig_Path="./figures/volcano_plot.png", maxGeneNameLen=8)
    """
    
    # Process the input data to include necessary columns
    DEA_res_sig = DEA_res_noDubs.copy()
    DEA_res_sig["nlog10"] = -np.log10(DEA_res_sig["adj.P.Val"])

    # Categorize genes based on significance thresholds
    def categorize_genes(row):
        if row["adj.P.Val"] >= fdrcutoff and abs(row["logFC"]) <= lfcutoff:
            return 'Not Sig.'
        elif row["adj.P.Val"] >= fdrcutoff and abs(row["logFC"]) > lfcutoff:
            return 'Log2 FC'
        elif row["adj.P.Val"] < fdrcutoff and abs(row["logFC"]) <= lfcutoff:
            return 'FDR'
        elif row["adj.P.Val"] < fdrcutoff and abs(row["logFC"]) > lfcutoff:
            return 'Log2 FC & FDR'
        return 'Not Sig.'

    DEA_res_sig["color"] = DEA_res_sig.apply(categorize_genes, axis=1)

    # Highlight top genes
    top_genes = DEA_res_sig[DEA_res_sig["color"] == 'Log2 FC & FDR']
    top_genes = top_genes[top_genes.index.str.len() < maxGeneNameLen]
    top_genes = top_genes.sort_values(by="logFC", ascending=False).head(10).append(
                top_genes.sort_values(by="logFC", ascending=True).head(10))
    
    # Plotting
    plt.figure(figsize=(9, 6))
    sns.set_style("whitegrid")
    ax = sns.scatterplot(data=DEA_res_sig, x="logFC", y="nlog10", hue="color", hue_order=["Not Sig.", "Log2 FC", "FDR", "Log2 FC & FDR"], palette=["gray", "green", "blue", "red"], s=45, alpha=0.5, edgecolor="none")
    
    # Add threshold lines
    ax.axhline(-np.log10(fdrcutoff), ls="--", c="k", lw=2)
    ax.axvline(lfcutoff, ls="--", c="k", lw=2)
    ax.axvline(-lfcutoff, ls="--", c="k", lw=2)
    
    plt.legend(loc="center", bbox_to_anchor=(.33, 1.05), ncol=4, frameon=False, prop={"weight": "bold"})
    
    # Customize plot appearance
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.tick_params(width=2)
    
    plt.title(figureTitle, fontdict={"family": "sans-serif", "weight": "bold", "size": 18})
    plt.xlabel("$log_{2}$ Fold Change", fontdict={"family": "sans-serif", "weight": "bold", "size": 14})
    plt.ylabel("-$log_{10}$ FDR", fontdict={"family": "sans-serif", "weight": "bold", "size": 14})
    plt.xticks(size=12, weight="bold")
    plt.yticks(size=12, weight="bold")
    plt.tight_layout()
    
    # Annotate top genes
    ta.allocate_text(plt.gcf(), ax, top_genes["logFC"].values, top_genes["nlog10"].values, top_genes.index, draw_all=True, max_distance=0.8, min_distance=0.01, margin=0.01, linecolor="gray", avoid_label_lines_overlap=True, textsize=14, textcolor="black")
    
    # Save the figure
    if savefig_Path:
        plt.savefig(savefig_Path, dpi=600, bbox_inches="tight", facecolor='w')
    plt.show()
    plt.clf()
