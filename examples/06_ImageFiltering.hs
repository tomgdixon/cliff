{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : Example 6 - Sliding Geometric Stencils and Clifford Fourier Transforms (DCFT)
--
-- This tutorial example covers:
-- 1. Representing 2D Multivector Grid Fields (Images, Vector Fields, Spinor Fields)
-- 2. Sliding Geometric Stencil Convolutions (1D/2D) with Boundary Conditions
-- 3. 2D Discrete Clifford Fourier Transform (DCFT)
-- 4. Spectral Frequency Filtering
module Main where

import Prelude hiding ((<*>))
import Clifford

type GA2 = DenseMV Cl200 Double

main :: IO ()
main = do
  putStrLn "=================================================================================="
  putStrLn "  Clifford Algebra Tutorial: Example 6 - Stencils & Clifford Fourier Transforms"
  putStrLn "==================================================================================\n"

  -- 1. Representing Multivector Grid Fields
  -- A 4x4 spatial field where each pixel contains a multivector in Cl(2,0,0) (Scalar + Vector + Bivector).
  let inputGrid = generateND [4, 4] $ \coords ->
        case coords of
          [r, c] ->
            let x = fromIntegral c
                y = fromIntegral r
            in (x * e1 + y * e2 + (x * y) * (basisMV 1 /\ basisMV 2) :: GA2)
          _ -> zeroMV

  putStrLn "--- 1. Input 2D Multivector Field (4x4) ---"
  putStrLn $ "Field shape: " ++ show (shapeND inputGrid)
  putStrLn $ "Pixel at (1, 2) = " ++ show (indexND inputGrid [1, 2])
  putStrLn ""

  -- 2. Geometric Stencil Convolution
  -- A 3x3 filter kernel with non-trivial geometric weights:
  let k0 = scalarMV (1/16) :: GA2
      k1 = scaleMV (1/8) (e1 :: GA2)
      k2 = scaleMV (1/4) (basisMV 1 /\ basisMV 2 :: GA2)

      kernel = mkMatrix 3 3
        [ k0, k1, k0
        , k1, k2, k1
        , k0, k1, k0
        ]

      -- Apply 2D geometric convolution with Periodic (toroidal) boundary condition:
      convolvedGrid = convolve2D Periodic GeometricProd (unMVMatrix kernel) inputGrid

  putStrLn "--- 2. Sliding Stencil Convolution ---"
  putStrLn $ "Filtered pixel at (1, 2) = " ++ show (indexND convolvedGrid [1, 2])
  putStrLn ""

  -- 3. 2D Discrete Clifford Fourier Transform (DCFT)
  -- The Clifford Fourier transform uses the bivector pseudoscalar I = e12 as the hypercomplex unit (I^2 = -1):
  -- F(u, v) = sum_x,y f(x, y) * exp(-I * (u*x/W + v*y/H) * 2*pi)
  let cftKernel = defaultCFT2D 4 4
      freqSpectrum = cftForward cftKernel inputGrid
      reconstructed = cftInverse cftKernel freqSpectrum

  putStrLn "--- 3. 2D Discrete Clifford Fourier Transform (DCFT) ---"
  putStrLn $ "Frequency DC component F(0, 0) = " ++ show (indexND freqSpectrum [0, 0])
  putStrLn $ "Reconstructed pixel (1, 2)     = " ++ show (indexND reconstructed [1, 2])
  putStrLn ""

  let origP = indexND inputGrid [1, 2]
      recP  = indexND reconstructed [1, 2]
      diff = subMV origP recP
  putStrLn $ "DCFT Invertibility Error ||F^-1(F(f)) - f||: " ++ show (normSqMV diff)
  putStrLn ""
  putStrLn "Tutorial 6 completed successfully!"
