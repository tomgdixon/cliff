# cliff: High-Performance Clifford Algebra (Geometric Algebra) in Pure Haskell

**`cliff`** is a foundational, high-performance Clifford Algebra (Geometric Algebra) library implemented in pure Haskell with **100% standard GHC boot libraries only** (`base`, `array`, `containers`, `deepseq`).

It is designed for:
1. **Geometric Calculus & Physics**: Exact spatial rotations, relativistic Spacetime Algebra $Cl(1,3,0)$, and 3D Projective Geometric Algebra $Cl(3,0,1)$ kinematics.
2. **Clifford Neural Networks & LLMs**: Multivector GEMM, structure matrix isomorphisms, gauge/rotor-invariant attention scoring, and self-adjoint gradient rules.
3. **Coding Theory & High Dimensions**: Arbitrary dimension algebras ($n$ up to 64 and beyond) over finite Galois fields ($\mathbb{F}_2, \mathbb{F}_p, \mathbb{F}_{2^m}$) using 64-bit basis blade bitmasks.
4. **Signal Processing & Transforms**: N-dimensional sliding geometric stencils, $O(N \log N)$ Fast Clifford Fourier Transforms (FCFT), Number Theoretic Transforms (NTT), and Fast Walsh-Hadamard Transforms (FWHT) with dyadic (XOR) convolutions.

---

## Key Features

- **Pure GHC Boot Footprint**: Zero external dependencies beyond standard GHC boot libraries (`base`, `array`, `containers`, `deepseq`). Zero Template Haskell, zero C FFI.
- **Arbitrary Metric Signatures**: Full type-level support for $Cl(p,q,r)$ from 2D Euclidean ($Cl(2,0,0)$) and Spacetime Algebra ($Cl(1,3,0)$) to 3D Projective Geometric Algebra ($Cl(3,0,1)$), Conformal Geometric Algebra ($Cl(4,1,0)$), and 64D algebras.
- **Dual Representation Strategy**:
  - `DenseMV`: Unboxed dense flat arrays with in-place `ST` Cayley accumulation ($n \le 6$).
  - `SparseMV`: Bitmask-indexed sparse multivectors for high dimensions ($n \in [7, 64]$).
- **Analytic Functions & Series Expansions**: Numerically stable scaling-and-squaring exponential (`expMV`), trigonometric (`sinMV`, `cosMV`), hyperbolic (`sinhMV`, `coshMV`, `tanhMV`), and Horner polynomial evaluations (`hornerMV`).
- **Galois Finite Fields**: First-class support for finite fields $\mathbb{F}_2$, $\mathbb{F}_p$, and $\mathbb{F}_{2^m}$ (AES polynomial $0x11B$) for error-correcting codes, cryptography, and combinatorics.
- **CliffordNet & Machine Learning Primitives**: Multivector GEMM (`matMulMV`), structure matrix isomorphisms ($M(ab) = M(a)M(b)$), rotor-invariant attention scoring, and self-adjoint gradient rules.
- **Multidimensional Spectral Transforms**: 1D/2D/ND geometric stencil convolutions, $O(N \log N)$ Fast Clifford Fourier Transforms (FCFT), Galois Number Theoretic Transforms (NTT), and Fast Walsh-Hadamard Transforms (FWHT) with exact dyadic (bitwise XOR) convolutions.

---

## Quickstart

```haskell
{-# LANGUAGE DataKinds #-}
import Clifford

-- 3D Euclidean GA (Cl(3,0,0))
type GA3 = DenseMV Cl300 Double

main :: IO ()
main = do
  -- Define vectors
  let v1 = 2 * e1 + 3 * e2 :: GA3
      v2 = 1 * e2 + 4 * e3 :: GA3

  -- Geometric product: v1 * v2 = v1 . v2 + v1 ^ v2
  let prod = v1 <*> v2
  putStrLn $ "v1 * v2 = " ++ show prod

  -- Rotor rotation of v1 by 90 degrees in the e12 plane
  let rotor = rotorExp (scaleMV (- pi / 4) e12)
      v1_rot = sandwich rotor v1
  putStrLn $ "Rotated v1 = " ++ show v1_rot

  -- 3D PGA translation
  let origin = pointPGA 0 0 0 :: DenseMV Cl301 Double
      motor = translatorPGA 1 2 3 :: DenseMV Cl301 Double
      pTranslated = sandwich motor origin
  putStrLn $ "Translated point: " ++ show pTranslated
```

---

## Benchmarks

Measured on a standard GHC 9.14 compiler using tasty-bench`(on a raspberry pi):

| Operation | Space / Scale | Latency |
|:---|:---|:---|
| **3D Multivector Product** | $Cl(3,0,0)$ Euclidean (Dense) | **261 ns** |
| **Spacetime Algebra Product** | $Cl(1,3,0)$ Minkowski (Dense) | **201 ns** |
| **Grade 8 ^ Grade 8 Wedge** | 32D Space (Sparse Bitmask) | **241 ns** |
| **Grade 16 ^ Grade 16 Wedge** | 64D Space (Sparse Bitmask) | **327 ns** |
| **Grade 32 * Grade 32 $\to$ Grade 64** | 64D Pseudoscalar Product | **468 ns** |
| **Grade 16 $\rfloor$ Grade 32 Contraction** | 64D Space (Sparse Bitmask) | **429 ns** |
| **64D Poincaré Dual $J(B_{32})$** | 64D Complement Dual | **370 ns** |
| **3D Horner Polynomial (Deg 4)** | $Cl(3,0,0)$ Horner Method | **893 ns** |
| **3D Multivector Sine (`sinMV`)** | $Cl(3,0,0)$ Scaling & Squaring | **2.69 μs** |
| **3D Multivector Cosine (`cosMV`)** | $Cl(3,0,0)$ Scaling & Squaring | **2.95 μs** |
| **3D Multivector Exp (`expMV`)** | $Cl(3,0,0)$ Scaling & Squaring | **6.04 μs** |
| **PGA 3D Screw Twist Exp** | $Cl(3,0,1)$ Rigid Body Twist | **5.55 μs** |
| **3D Gauss-Jordan Inverse** | $8 \times 8$ Augmented `ST` Matrix | **29.1 μs** |
| **Versor Inverse** | $V^{-1} = \widetilde{V} / \langle V \widetilde{V} \rangle_0$ | **~50 ns** |
| **16D Multivector Product (Real)** | $Cl(16,0,0)$ Double Sparse | **1.78 μs** |
| **32D Multivector Product (GF2)** | $Cl(32,0,0)$ over $\mathbb{F}_2$ ($2^{32}$ space) | **2.05 μs** |
| **64D Multivector Product (GF2)** | $Cl(64,0,0)$ over $\mathbb{F}_2$ ($2^{64}$ space) | **1.90 μs** |
| **Multivector GEMM** | $(8 \times 8) \times (8 \times 8)$ Matrix | **148 μs** |
| **2D Stencil Convolution** | $(16 \times 16)$ Field, $(3 \times 3)$ Kernel | **1.17 ms** |
| **2D Fast Clifford Fourier (FCFT)** | $(8 \times 8)$ Grid, $O(N \log N)$ | **407 μs** |
| **2D Number Theoretic Transform (NTT)**| $(8 \times 8)$ Grid over $\mathbb{F}_{17}$, $O(N \log N)$ | **482 μs** |
| **2D Direct Clifford Fourier (DCFT)** | $(8 \times 8)$ Grid, $O(N^2)$ Reference | **2.40 ms** |

---

## Building & Testing

```bash
# Build the library
cabal build

# Run all 42 tests
cabal test

# Run performance benchmarks
cabal bench
```

---

## Documentation & Tutorials

- [Beginner's Guide & Tutorial to Clifford Algebra](docs/TUTORIAL.md)
- [System & Mathematical Architecture](docs/ARCHITECTURE.md)
- [Complete API Reference](docs/API_REFERENCE.md)
- [Implementation Tasks & Roadmap](docs/TASKS.md)

---

## Interactive Examples

Run any of the self-contained tutorial examples directly with `cabal`:

```bash
# 1. Multivector Basics & Geometric Product
cabal run example-basics

# 2. 2D/3D Rotations & Rotors
cabal run example-rotors

# 3. 3D Projective Geometric Algebra (PGA) Kinematics
cabal run example-pga

# 4. Finite Galois Fields & High Dimensions (16D/32D)
cabal run example-galois

# 5. CliffordNet LLM Attention & Multivector GEMM
cabal run example-cliffordnet

# 6. 2D Geometric Stencils & Clifford Fourier Transforms
cabal run example-stencils
```

---

## License

BSD-3-Clause License. See [LICENSE](LICENSE) for details.
