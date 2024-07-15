# Scripts for RATTACA LD Pruning Testing

## [Associated with RATTACA: Genetic predictions in Heterogeneous Stock rats offer a new tool for genetic correlation and experimental design](https://doi.org/10.1101/2023.09.18.558279)
### Benjamin B. Johnson, Thiago M. Sanches, <ins>Mika H. Okamoto</ins>, Khai-Min Nguyen, Clara A. Ortez, Oksana Polesskaya, Abraham A. Palmer

### Main Notebook: plot_correlations.ipynb

Includes scripts written in R and Python to generate prediction performance data for various experiments and to visualize the data. 


Comparing phenotype prediction performance with different genome subsampling methods including:
- LD pruning parameters - $r^2$ and window size 
- Number of random SNPs
- Number of training rats
- Random vs LD Pruning
- LD clumping

And graphing:
- Prediction performance distributions
- Runtimes

main folder: main prediction pipeline and correlation plots<br>
full_pred_pipeline.r: code for general prediction pipeline<br>
plot_correlations.ipynb: correlation graphs from various experiments<br>
plot_runtimes_py.ipynb: runtime graphs from various experiments (to demonstrate cost of different methods)<br>
convex_hull.ipynb: code for testing various rat breeding algorithms for maximizing genetic diversity of offspring<br>
<br>
/experiments/code: pipelines for generating performance data for different genome subsampling methods<br>
/experiments/plots: notebooks of plots of different experiments of different genome subsampling methods<br>
/pyrrBLUP: class to run rrBLUP in python, used for scikit-learn comparisons<br>
/experimental: scripts for testing and new ideas (not public)<br>
/old: in-between scripts and old ideas which were refined in other scripts (not public)<br>
