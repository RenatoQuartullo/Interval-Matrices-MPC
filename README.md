# Interval Matrices Model Predictive Control

MATLAB code for implementing the Interval Matrices MPC technique proposed in this [paper](https://arxiv.org/abs/2602.17379):

Quartullo, R., Garulli, A., and Leomanni, M. (2026). *Robust Model Predictive Control for Linear Systems with Interval Matrix Model Uncertainty*. arXiv e-prints, arXiv: 2602.17379.


## Required Toolboxes

- [CORA](https://it.mathworks.com/matlabcentral/fileexchange/68551-cora) – for interval matrices operations

## Recommended Solver

- [Mosek](https://www.mosek.com/) - The runtime evaluation in the paper has been performed using Mosek. However, the code runs also without it. In the case the user does not have a Mosek licence, quadprog will refer to the MATLAB embedded function (without any modification)

## Implemented Baselines

- [Polytopic Tube MPC](https://link.springer.com/book/10.1007/978-3-319-24853-0) – implemented following **Chapter 5** of:

  Kouvaritakis, B., and Cannon, M. *Model Predictive Control*. Switzerland: Springer International Publishing, 2016, pp. 38.13–56.

## Other Baselines

- [Offline Tightening MPC](https://github.com/monimoyb/RMPCPy):
  Chen, S., Preciado, V. M., Morari, M., and Matni, N. (2024). *Robust model predictive control with polytopic model uncertainty through system level synthesis*. Automatica, 162, 111431.
- [System Level Synthesis MPC](https://github.com/ShaoruChen/Polytopic-SLSMPC): Bujarbaruah, M., Rosolia, U., Stürz, Y. R., Zhang, X., & Borrelli, F. (2022). *Robust MPC for LPV systems via a novel optimization-based constraint tightening*. Automatica, 143, 110459.


