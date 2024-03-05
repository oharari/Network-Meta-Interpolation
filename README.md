# Network-Meta-Interpolation
Code and data for our 2022 Research Synthesis Methods [paper](https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1608) "Network meta-interpolation: Effect modification adjustment in network meta-analysis using subgroup analyses". Go to [the tutorial](https://oharari.github.io/Network-Meta-Interpolation/) to learn about the method and how to run our code.

The *Code and Data* folder contains the following files:
1. **Example_AgD_ML_NMR.csv** – aggregate-level data at the arm level for running ML-NMR with dichotomous outcome and two effect modifiers ($x_1$ and $x_2$).
2. **Example_AgD_ML_NMR_3D.csv** – aggregate-level data at the arm level for running ML-NMR with dichotomous outcome and three effect modifiers ($x_1$, $x_2$ and $x_3$).
3. **Example_AgD_NMI.csv** – aggregate-level data at the study level for running NMI with dichotomous outcome and two effect modifiers ($x_1$ and $x_2$).
4. **Example_AgD_NMI_3D.csv** – aggregate-level data at the study level for running NMI with dichotomous outcome and three effect modifiers ($x_1$, $x_2$ and $x_3$).
5. **Example_AgD_NMR_NMA.csv** – aggregate-level data at the study level for running NMA/NMR with dichotomous outcome and two effect modifiers ($x_1$ and $x_2$).
6. **Example_AgD_NMR_NMA_3D.csv** – aggregate-level data at the study level for running NMA/NMR with dichotomous outcome and three effect modifiers ($x_1$, $x_2$ and $x_3$).
7. **Example_IPD.csv** – individual-level data with dichotomous outcome and two effect modifiers ($x_1$ and $x_2$). Used by both NMI and ML-NMR.
8. **Example_IPD_3D.csv** – individual-level data with dichotomous outcome and three effect modifiers ($x_1$, $x_2$ and $x_3$). Used by both NMI and ML-NMR.

In addition, included are the following R files – 
1. **NMI_Script.R** – an R script, to be run line-by-line, fitting NMA, NMR, ML-NMR and NMI to the simulated datasets included and comparing them to the ground-truth treatment effects that were used in their creation.  
2. **NMI_Functions.R** – the functions used in the abovementioned analysis. Sourced by NMI_Script.R.


$\text{\bf \underline{Instructions}:}$ keep the above files stored in the same directory, and run NMI_Script.R line-by-line for a full demo. 
