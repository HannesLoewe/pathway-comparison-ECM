# pathway-comparison-ECM
MATLAB code for metabolic pathway comparison using Enzyme Cost Minimization

Description:
This repository contains the MATLAB files necessary to perform the Enzyme Cost Minimization algorithm that was developed by Noor et al. (2016) (https://doi.org/10.1371/journal.pcbi.1005167) and apply it to the comparison of different CO2 or C1 fixing pathways as presented in Löwe & Kremling (2021) (in review). Details on the purpose of this work can be found in the original publication.

Dependencies:
- Metabolic network toolbox by Wolfram Liebermeister (https://github.com/liebermeister/metabolic-network-toolbox)
- Enzyme cost minimization by Wolfram Liebermeister (https://github.com/liebermeister/enzyme-cost-minimization)
- All necessary dependencies of the two repositories

Installation/Start:
- Add the repository's folder to your MATLAB environment
- Run Compare_Pathways_Std_Plot2.m as a starting point
- modify pathway stoichiometry in PathwayStoichiometries.m
- Adjust CO2/HCO3-
- Adjust plotting section to your needs

Differences to the original ECM Matlab implementation:
- the ecm_enzyme_cost_minimization.m was modified: in case of occuring errors during ECM, the algorithm will not stop but rather terminate with empty return values. This was changed to be able to run multiple iterations without crashing. Crashes would occur when no feasible solution can be found.
