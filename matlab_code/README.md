This folder contains matlab code that computes the trajectory metrics considered in: 

"Fundamental dynamics of popularity-similarity trajectories in real networks", E. S. Papaefthymiou, C. Iordanou, and F. Papadopoulos, arXiv:2309.01675, September 2023.

The code also creates simulated counterparts of the real trajectories using the fBm model described in the above paper. A short description of each script is given below. 

- **popularity_sims.m**: reads the radial or expected degree trajectories and creates figures similar to Fig. 5 in the above paper.
    
- **similarity_sims.m**: reads the angular trajectories and creates figures similar to Fig. 5 in the above paper. It also creates angular histograms as in Fig. 2(f) in the above paper.  

- **popularity_predictions.m**: reads the radial or expected degree trajectories and performs predictions using simple heuristics, as in Figs. 38-42(a),(b) in the above paper.

- **similarity_predictions.m**: reads the angular trajectories and performs predictions using simple heuristics, as in Figs. 38-42(c) in the above paper.

- **mbm.m**: implements the fBm model described in the above paper.
 
Please see the comments inside each script for more details.
