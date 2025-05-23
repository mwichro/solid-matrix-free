# A matrix-free approach for finite-strain hyperelastic problems using geometric multigrid

This is the repository for the code associated with the paper 

M. Wichrowskia, M. Rezaee-Hajidehib, J. Korelcc, M. Kronbichlerd, S. Stupkiewicz
Matrix-Free Methods for Finite-Strain Elasticity: Automatic Code
Generation with No Performance Overhead

arXiv: https://arxiv.org/pdf/2505.15535

If you use this work, or find it useful in general, the citation of the aforementioned article would be much appreciated by the authors.

## Abstract
This study explores matrix-free tangent evaluations in finite-strain elasticity with the use of automatically-generated code for the quadrature-point level calculations. The code generation is done via automatic differentiation (AD) with AceGen. We compare hand-written and AD-generated codes under two computing strategies: on-the-fly evaluation and caching intermediate results. The comparison reveals that the AD-generated code achieves superior performance in matrix-free computations.
