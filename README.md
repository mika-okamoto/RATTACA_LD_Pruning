# Scripts for RATTACA LD Pruning Testing

### Main Notebook: plot_correlations.ipynb

Includes scripts to generate prediction performance data for various experiments and to visualize the data. 

Comparing phenotype prediction performance with different genome subsampling methods including:
- LD pruning parameters - $r^2$ and window size 
- Number of random SNPs
- Number of training rats
- Random vs LD Pruning
- LD clumping

And graphing:
- Prediction performance distributions
- Runtimes

/experimental: scripts for testing and new ideas (not public)
/old: in-between scripts and old ideas which were refined in other scripts (not public)
main folder: main scripts for plots in plot_correlations.ipynb
/pyrrBLUP: class to run rrBLUP in python, used for scikit-learn comparisons