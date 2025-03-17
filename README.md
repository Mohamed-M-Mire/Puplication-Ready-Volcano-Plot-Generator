# Puplication-Ready Volcano Plot Generator

## Function: `generate_volcano_plot`

### Description

The `generate_volcano_plot` function generates a volcano plot based on input transcriptomics data and specified thresholds. The function highlights significant genes based on log fold change (logFC) and false discovery rate (FDR) cutoffs.

### Parameters

- `DEA_res_noDubs` (pd.DataFrame): DataFrame containing differential expression analysis results; must include 'logFC' and 'adj.P.Val' columns.
- `lfcutoff` (float): Log fold change cutoff value.
- `fdrcutoff` (float): False discovery rate cutoff value.
- `figureTitle` (str): Title for the generated plot.
- `topGeneCount` (int, optional): Number of top most significantly changed genes to subset the data by. Defaults to 15.
- `maxGeneNameLen` (int, optional): Maximum length of gene names to display. Defaults to 8.
- `savefig_Path` (str, optional): File path to save the plot. Defaults to None.


### Returns

- matplotlib figure object

### Example Usage

```python
import matplotlib.pyplot as plt
import textalloc as ta
import pandas as pd
import numpy as np

# Assuming DEA_res_noDubs is a DataFrame containing the necessary data
generate_volcano_plot(
    DEA_res_noDubs, 
    lfcutoff=1.5, 
    fdrcutoff=0.05, 
    figureTitle="Differential Expression Volcano Plot", 
    topGeneCount=15, 
    maxGeneNameLen=8, 
    savefig_Path="./figures/volcano_plot.png")
```
![alt text](https://github.com/Mohamed-M-Mire/Puplication-Ready-Volcano-Plot-Generator/blob/main/Figures/volcano_plot.png)

## Dependencies

- text_annotate
- matplotlib
- pandas
- numpy
  
## License

This project is licensed under the MIT License. See the LICENSE file for more details.

## Acknowledgments

Inspired by the R library EnhancedVolcano.
