# RxTree

The superior therapeutic effects of synergistic combination therapy have been discovered in treating multiple types of cancers and it is of high clinical value to accurately identify the synergistic combinations via shared or similar underlying biological mechanisms. However, not all possible drug combinations can be realistically tested on patients in real clinical trials. In the manuscript, Probabilistic Learning of Treatment Trees in Cancer (Yao et al., 2022), we address this unmet translational research need by proposing a novel Bayesian probabilistic tree-based framework, referred to as treatment trees (Rx-tree), to investigate the hierarchical relationships between treatments. This repository contains the code for the manuscript and we further visualize the result through an R-shiny application in [https://yaots.shinyapps.io/RxTreePDX]

The model estimation is decoupled into (i) Euclidean parameters estimation by the approximated Bayesian computation (ABC) and (ii) tree parameters estimation by Metropolis-Hastings algorithm (MH) and we implemented four R scripts, names as 
* ABC_SyntheticData.R: generate the synthetic data for the ABC inference
* Inference.R: conduct the ABC and MH inference to draw posterior tree samples
* Tree_Summary.R: summarize the posterior tree samples through maximum a posteriori estimation (MAP) and integrated posterior co-clustering (iPCP)
* Visualization.R: visualize the iPCP and the MAP tree 

We demonstrate the code in the "RxTree_Reproduce.Rmd" with the data from Novartis Institutes for Biomedical Research PDX Encyclopedia (NIBR-PDXE; Gao et al., 2015 and Rashid et al., 2020) and replicate the Figures 5 and 6 in the Main paper.

