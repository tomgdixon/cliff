# Changelog for `cliff`

## 0.1.0.0 (2026-08-14)
* Initial release of `cliff`, a zero-dependency high-performance Clifford Algebra library in pure Haskell.
* Statically-typed metric signatures `Signature (p :: Nat) (q :: Nat) (r :: Nat)` for 2D, 3D Euclidean, Spacetime Algebra, 3D PGA, CGA, and arbitrary 64D signatures.
* Dual multivector memory backends:
  * `DenseMV`: Flat unboxed arrays with in-place `ST` Cayley accumulation for dimensions $n \le 6$.
  * `SparseMV`: 64-bit basis blade bitmask maps for arbitrary high dimensions ($n \in [7, 64]$).
* Universal scalar field hierarchy: `Double`, `Float`, `Rational`, `Complex Double`, `GF2` ($\mathbb{F}_2$), `GFp` ($\mathbb{F}_p$), and `GF2m` ($\mathbb{F}_{2^m}$).
* Versors, rotors, motors, outermorphisms, and $2^n \times 2^n$ Gauss-Jordan general inverses.
* 3D Projective Geometric Algebra (PGA $Cl(3,0,1)$) with Poincaré complement duals, points, lines, planes, and motor kinematics.
* Analytic multivector functions: Numerically stable scaling-and-squaring exponential (`expMV`), trigonometric (`sinMV`, `cosMV`), hyperbolic (`sinhMV`, `coshMV`, `tanhMV`), and Horner polynomial evaluation (`hornerMV`).
* Multidimensional array buffers (`MVArray`), multivector GEMM (`matMulMV`), structure matrix isomorphisms, and gauge/rotor-invariant attention scoring.
* In-place `STUArray` sliding geometric convolutions (`convolve1D`, `convolve2D`, `convolveND`) with customizable boundary modes.
* Multidimensional spectral transforms:
  * Fast Clifford Fourier Transforms (FCFT) in $O(N \log N)$ radix-2 Cooley-Tukey butterfly.
  * Galois Field Number Theoretic Transforms (NTT) for prime fields $\mathbb{F}_p$.
  * Fast Walsh-Hadamard Transforms (FWHT) with exact dyadic (bitwise XOR) convolutions.
