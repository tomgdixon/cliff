{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : High-performance micro-benchmark suite for cliff
-- Note: AI helped develop this test bench
module Main (main) where

import Prelude hiding ((<*>))
import Test.Tasty.Bench
import Clifford

main :: IO ()
main = defaultMain
  [ bgroup "Clifford Micro-Operations"
      [ bench "3D Euclidean (Cl300) Multivector Product" $ nf (\(a, b) -> a <*> b)
          ( (2 * e1 + 3 * e2 + 5 * e12) :: DenseMV Cl300 Double
          , (1 * e2 + 4 * e3 + 7 * e123) :: DenseMV Cl300 Double
          )
      , bench "Spacetime Algebra (Cl130) Multivector Product" $ nf (\(a, b) -> a <*> b)
          ( (3 * basisMV 1 + 2 * basisMV 2) :: DenseMV Cl130 Double
          , (5 * basisMV 1 + 4 * basisMV 3) :: DenseMV Cl130 Double
          )
      , bench "Sparse 16D Algebra over GF(2) Product" $ nf (\(a, b) -> a <*> b)
          ( (basisMV 1 `addMV` basisMV 5 `addMV` basisMV 12) :: SparseMV Cl16_00 GF2
          , (basisMV 5 `addMV` basisMV 16) :: SparseMV Cl16_00 GF2
          )
      ]

  , bgroup "Higher-Order Grades & 64-Bit Sparse Operations"
      [ bench "Grade 8 ^ Grade 8 Wedge Product (32D)" $ nf (\(a, b) -> a /\ b)
          ( bladeMV (bladeFromIndices [1..8]) 1.0   :: SparseMV Cl32_00 Double
          , bladeMV (bladeFromIndices [9..16]) 1.0  :: SparseMV Cl32_00 Double
          )
      , bench "Grade 16 ^ Grade 16 Wedge Product (64D)" $ nf (\(a, b) -> a /\ b)
          ( bladeMV (bladeFromIndices [1..16]) 1.0  :: SparseMV Cl64_00 Double
          , bladeMV (bladeFromIndices [17..32]) 1.0 :: SparseMV Cl64_00 Double
          )
      , bench "Grade 32 * Grade 32 to Grade 64 Pseudoscalar (64D)" $ nf (\(a, b) -> a <*> b)
          ( bladeMV (bladeFromIndices [1..32]) 1.0  :: SparseMV Cl64_00 Double
          , bladeMV (bladeFromIndices [33..64]) 1.0 :: SparseMV Cl64_00 Double
          )
      , bench "Grade 16 . Grade 32 Left Contraction (64D)" $ nf (\(a, b) -> a <.> b)
          ( bladeMV (bladeFromIndices [1..16]) 1.0  :: SparseMV Cl64_00 Double
          , bladeMV (bladeFromIndices [1..32]) 1.0  :: SparseMV Cl64_00 Double
          )
      , bench "64D Poincaré Dual J(Grade 32 Blade)" $ nf (\a -> poincareDualMV a)
          ( bladeMV (bladeFromIndices [1..32]) 1.0  :: SparseMV Cl64_00 Double
          )
      , bench "16D Sparse Multivector Product (Double)" $ nf (\(a, b) -> a <*> b)
          ( fromBladeList
              [ (bladeFromIndices [1, 2, 3], 2.5)
              , (bladeFromIndices [4, 5, 6, 7], 1.8)
              , (bladeFromIndices [8, 9, 10, 11, 12], -3.1)
              ] :: SparseMV Cl16_00 Double
          , fromBladeList
              [ (bladeFromIndices [3, 4], 4.2)
              , (bladeFromIndices [6, 7, 8], -1.5)
              , (bladeFromIndices [12, 13, 14, 15, 16], 0.9)
              ] :: SparseMV Cl16_00 Double
          )
      , bench "32D Sparse Multivector Product (GF2, 2^32 space)" $ nf (\(a, b) -> a <*> b)
          ( fromBladeList
              [ (bladeFromIndices [1, 7, 15, 23], GF2 1)
              , (bladeFromIndices [4, 12, 20, 28], GF2 1)
              , (bladeFromIndices [2, 10, 18, 26, 31], GF2 1)
              ] :: SparseMV Cl32_00 GF2
          , fromBladeList
              [ (bladeFromIndices [1, 4, 9, 16, 25], GF2 1)
              , (bladeFromIndices [2, 8, 18, 30], GF2 1)
              , (bladeFromIndices [7, 14, 21, 28], GF2 1)
              ] :: SparseMV Cl32_00 GF2
          )
      , bench "64D Sparse Multivector Product (GF2, 2^64 space)" $ nf (\(a, b) -> a <*> b)
          ( fromBladeList
              [ (bladeFromIndices [1, 16, 32, 48, 64], GF2 1)
              , (bladeFromIndices [8, 24, 40, 56], GF2 1)
              , (bladeFromIndices [5, 10, 15, 20, 25, 30, 35, 40], GF2 1)
              ] :: SparseMV Cl64_00 GF2
          , fromBladeList
              [ (bladeFromIndices [2, 16, 30, 44, 58], GF2 1)
              , (bladeFromIndices [8, 16, 24, 32, 40, 48, 56, 64], GF2 1)
              ] :: SparseMV Cl64_00 GF2
          )
      ]

  , bgroup "Analytic Functions & Series Expansions"
      [ bench "3D Multivector Exponential (expMV)" $ nf (\m -> expMV m)
          ((scaleMV 0.5 e1 `addMV` scaleMV 0.8 e23 `subMV` scaleMV 0.4 (e123 :: DenseMV Cl300 Double)))
      , bench "3D Multivector Sine (sinMV)" $ nf (\m -> sinMV m)
          ((scaleMV 0.6 e1 `addMV` scaleMV 0.3 (e12 :: DenseMV Cl300 Double)))
      , bench "3D Multivector Cosine (cosMV)" $ nf (\m -> cosMV m)
          ((scaleMV 0.7 e12 `subMV` scaleMV 0.4 (e23 :: DenseMV Cl300 Double)))
      , bench "3D Multivector Hyperbolic Sine (sinhMV)" $ nf (\m -> sinhMV m)
          ((scaleMV 0.5 e1 `addMV` scaleMV 0.2 (e123 :: DenseMV Cl300 Double)))
      , bench "3D Multivector Horner Degree-4 Polynomial" $ nf (\m -> hornerMV [1.0, 2.0, -1.5, 0.5, 0.1] m)
          ((scaleMV 0.4 e1 `addMV` scaleMV 0.3 (e2 :: DenseMV Cl300 Double)))
      , bench "PGA 3D Screw Twist Exponential Map" $ nf (\twist -> expMV twist)
          ((scaleMV (- pi / 4) e12 `addMV` scaleMV 2.0 (basisMV 4 /\ basisMV 3)) :: DenseMV Cl301 Double)
      ]

  , bgroup "Tensor, Matrix & Transform Operations"
      [ bench "3D Gauss-Jordan General Inverse (8x8 matrix)" $ nf (\a -> generalInverse a)
          ((scalarMV 2 `addMV` scaleMV 3 e1 `addMV` scaleMV 4 e12 `addMV` scaleMV 1 (e123 :: DenseMV Cl300 Double)))
      , bench "2D Direct Clifford Fourier DCFT (8x8 grid, O(N^2))" $ nf (\(k, img) -> cftForward k img)
          ( defaultCFT2D 8 8
          , fromListND [8, 8] (replicate 64 ((1 * e1 + 2 * e2) :: DenseMV Cl200 Double))
          )
      , bench "2D Fast Clifford Fourier FCFT (8x8 grid, O(N log N))" $ nf (\img -> fcft2D (e12 :: DenseMV Cl200 Double) (e12 :: DenseMV Cl200 Double) img)
          ( fromListND [8, 8] (replicate 64 ((1 * e1 + 2 * e2) :: DenseMV Cl200 Double))
          )
      , bench "2D Number Theoretic Transform NTT (8x8 grid, GF17, O(N log N))" $ nf (\img -> ntt2D (GFp 2 :: GFp 17) (GFp 2 :: GFp 17) img)
          ( fromListND [8, 8] (replicate 64 ((scaleMV (GFp 1) e1 `addMV` scaleMV (GFp 2) e2) :: DenseMV Cl200 (GFp 17)))
          )
      , bench "2D Fast Walsh-Hadamard Transform FWHT (8x8 grid, Cl300)" $ nf (\img -> fwht2D img)
          ( fromListND [8, 8] (replicate 64 ((1 * e1 + 2 * e2) :: DenseMV Cl300 Double))
          )
      , bench "DenseMV 16-Blade FWHT (Cl130 STA)" $ nf (\v -> fwhtMV v)
          ((scalarMV 2 `addMV` scaleMV 3 e1 `addMV` scaleMV 4 e12 `addMV` scaleMV 1 (e123 :: DenseMV Cl130 Double)))
      , bench "DenseMV Dyadic (XOR) Convolution (Cl130 STA)" $ nf (\(a, b) -> dyadicConvolveMV a b)
          ( (scalarMV 2 `addMV` scaleMV 3 e1) :: DenseMV Cl130 Double
          , (scaleMV 4 e12 `addMV` scaleMV 1 e123) :: DenseMV Cl130 Double
          )
      , bench "1D Dyadic Array Convolution (N=16, Cl300)" $ nf (\(a, b) -> dyadicConvolve1D (<*>) a b)
          ( fromListND [16] (replicate 16 ((1 * e1 + 2 * e2) :: DenseMV Cl300 Double))
          , fromListND [16] (replicate 16 ((3 * e2 + 4 * e3) :: DenseMV Cl300 Double))
          )
      , bench "2D Stencil Convolution (16x16 with 3x3 kernel)" $ nf (\(k, img) -> convolve2D Periodic GeometricProd k img)
          ( fromListND [3, 3] (replicate 9 ((1 * e1 + 2 * e2) :: DenseMV Cl200 Double))
          , fromListND [16, 16] (replicate 256 ((3 * e1 + 4 * e12) :: DenseMV Cl200 Double))
          )
      , bench "Multivector GEMM (8x8 x 8x8)" $ nf (\(a, b) -> matMulMV a b)
          ( mkMatrix 8 8 (replicate 64 ((1 * e1 + 2 * e2) :: DenseMV Cl300 Double))
          , mkMatrix 8 8 (replicate 64 ((3 * e2 + 4 * e3) :: DenseMV Cl300 Double))
          )
      ]
  ]
