# Puplication-Ready Volcano Plot Generator

This repository contains a function to generate volcano plots based on differential expression analysis (DEA) results. Volcano plots are useful for visualizing the relationship between fold change and statistical significance in gene expression studies.

## Function: `generate_volcano_plot`

### Description

The `generate_volcano_plot` function generates a volcano plot based on input transcriptomics data and specified thresholds. The function highlights significant genes based on log fold change (logFC) and false discovery rate (FDR) cutoffs.

### Parameters

- `DEA_res_noDubs` (pd.DataFrame): DataFrame containing differential expression analysis results.
- `lfcutoff` (float): Log fold change cutoff value.
- `fdrcutoff` (float): False discovery rate cutoff value.
- `figureTitle` (str): Title for the generated plot.
- `savefig_Path` (str, optional): File path to save the plot. Defaults to None.
- `maxGeneNameLen` (int, optional): Maximum length of gene names to display. Defaults to 20.

### Returns

- None

### Example Usage

```python
import pandas as pd

# Assuming DEA_res_noDubs is a DataFrame containing the necessary data
generate_volcano_plot(DEA_res_noDubs, lfcutoff=1.5, fdrcutoff=0.05, figureTitle="Differential Expression Volcano Plot", savefig_Path="./figures/volcano_plot.png", maxGeneNameLen=8)

### Notes

- The volcano plot visualizes the relationship between fold change (logFC) and statistical significance (FDR).
- Significant genes are highlighted based on the specified cutoffs.
- The color mapping function assigns gene names to categories: "Not Sig.", "Log2 FC", "FDR", or "Log2 FC & FDR".

## Dependencies

- pandas
- numpy
- matplotlib
- seaborn
- text_annotate
- 
## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Acknowledgments

Special thanks to the contributors of the libraries used in this project.

Inspired by the R library EnhancedVolcano.
