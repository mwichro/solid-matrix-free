# Matrix-free FEM large-strain — kernel generation README

Short guide to generate and prepare Neo-Hookean evaluation kernels produced by AceGen for matrix-free finite‑element evaluation (large strain).

## Overview
1. Use the AceGen script `NeoHookean.wls` to generate C/C++ source for element/point evaluation kernels.  
2. The generated source needs three small modifications to be efficient and usable in a deal.II-based, vectorized, matrix‑free code:
    - replace AceGen "scratch" arrays with named scalar variables (to avoid large unused arrays and to enable compiler optimizations),
    - adapt function headers so arguments use deal.II Tensor types instead of raw C arrays,
    - make the function templated on `Number` so it can be instantiated for `double` or SIMD number types.

## Generate the source
Run AceGen (example):
```
# Open with Mathematica and run:
NeoHookean.wls   # produces NeoHookean.c / NeoHookean.cc (AceGen output)
```

## Post-process AceGen scratch data
Problem: AceGen emits a single large scratch array (e.g. `acegen_scratch[...]`) that prevents the compiler from optimizing and wastes cache.

Solution: use the provided helper script `ReplaceAceGenScratch.sh`. What it does:
- scans the generated source for all used occurrences of `acegen_scratch`,
- replaces them with uniquely named scalar variables `acegenscratch__00__`, `acegenscratch__01__`, ...
- inserts declarations for all these variables at the start of the function body, typed as `Number` so they are ready for templating.

Example usage:
```
./ReplaceAceGenScratch.sh NeoHookean.cc > NeoHookean_scratch_fixed.cc
# or inplace (if the script supports -i): ./ReplaceAceGenScratch.sh -i NeoHookean.cc
```

Notes:
- Verify the mapping (the script should print or log how scratch indices map to names).
- The resulting variables are plain scalars (not arrays) so unused ones do not allocate large contiguous memory.


## Templating by Number
Wrap the entire generated kernel function as a template:

Before:
```
void neo_hookean_eval(... /* uses double */) { ... }
```
After:
```
template <typename Number>
void neo_hookean_eval(... /* replace double arrays with Tensor<..,Number>& */) { ... }
```

Benefits:
- allows instantiation with `double` or with SIMD/vectorized number types (e.g. deal.II's `VectorizedArray<double>`),
- together with `ReplaceAceGenScratch.sh` (which declares scratch as `Number`) this makes the kernel SIMD‑capable.

Example instantiation / call:
```
using Number = double;
neo_hookean_eval<Number>(...);

using Number = dealii::VectorizedArray<double>;
neo_hookean_eval<Number>(...);   // for vectorized runs
```

If you prefer explicit instantiation to keep the object file small, add
```
template void neo_hookean_eval<double>(...); 
```
in a .cc file.


## Header / argument adjustments for deal.II
AceGen assumes raw C arrays (e.g. `double gradduIn[dim][dim]`). For easy use with deal.II replace such arguments by typed Tensor references:

- Replace
  ```
  double gradduIn[dim][dim]
  ```
  with
  ```
  Tensor<2, dim, Number> &gradduIn
  ```
- Similarly replace other multidimensional C arrays with `Tensor<rank, dim, Number> &` or `SymmetricTensor<>` where appropriate.

Reason: deal.II `Tensor` objects use the same operator syntax (`operator[][]`) as raw C arrays for read/write, so the body of the AceGen kernel typically needs no other changes after the header change.

Add required includes near the top of the file:
```
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>   // if used
```

## Build tips
- Ensure `DIM` or `dim` is provided via template or compile-time constant consistent with AceGen output.
- Compile with optimization flags (e.g. -O3 -march=native -ffast-math) for best SIMD performance.
- If using vectorized `Number`, enable relevant compile flags (e.g. AVX, AVX2).


## Quick checklist
- [ ] Run NeoHookean.wls → get .cc/.c file
- [ ] Run ReplaceAceGenScratch.sh → scalar scratch variables inserted
- [ ] Edit headers: raw C arrays → Tensor<..,Number> &
- [ ] Make function template<typename Number>
- [ ] Instantiate and test for SIMD Number type