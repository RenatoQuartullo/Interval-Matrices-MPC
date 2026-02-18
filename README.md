# Interval Matrices Model Predictive Control

MATLAB code for implementing the Interval Matrices MPC technique.

## Required Toolboxes

- [MPT3](https://www.mpt3.org/) – for polyhedron operations  
- [CORA](https://it.mathworks.com/matlabcentral/fileexchange/68551-cora) – for interval matrices operations  
- [CVX](https://cvxr.com/cvx/) – for dealing with Linear Matrix Inequalities (LMI).  
  Not needed for replicating the paper results; all LMI results are already available in the data.

## Implemented Baselines

- [Polytopic Tube MPC](https://link.springer.com/book/10.1007/978-3-319-24853-0) – implemented following **Chapter 5** of:

  Kouvaritakis, Basil, and Mark Cannon. *Model Predictive Control*. Switzerland: Springer International Publishing, 2016, pp. 38.13–56.

## Other Baselines

- [Offline Tightening MPC](https://github.com/monimoyb/RMPCPy):
  Chen, S., Preciado, V. M., Morari, M., & Matni, N. (2024). *Robust model predictive control with polytopic model uncertainty through system level synthesis*. Automatica, 162, 111431.
- [System Level Synthesis MPC](https://github.com/ShaoruChen/Polytopic-SLSMPC): Bujarbaruah, M., Rosolia, U., Stürz, Y. R., Zhang, X., & Borrelli, F. (2022). *Robust MPC for LPV systems via a novel optimization-based constraint tightening*. Automatica, 143, 110459.


