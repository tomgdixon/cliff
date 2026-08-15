{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : Example 5 - Clifford Neural Networks (CliffordNet), Multivector GEMM & Rotor-Invariant Attention
--
-- This tutorial example covers:
-- 1. Multivector Matrices and General Matrix Multiplication (GEMM)
-- 2. Structure Matrix Homomorphism: M(a * b) = M(a) * M(b)
-- 3. Gauge / Rotor-Invariant Attention for Clifford Transformers & LLMs
-- 4. Self-Adjoint Backpropagation Gradient Rules
module Main where

import Prelude hiding ((<*>))
import Data.Array.IArray ((!))
import Clifford

type GA3 = DenseMV Cl300 Double

main :: IO ()
main = do
  putStrLn "================================================================================"
  putStrLn "  Clifford Algebra Tutorial: Example 5 - CliffordNet LLM Attention & GEMM"
  putStrLn "================================================================================\n"

  -- 1. Multivector Matrices & GEMM (Matrix Multiplication)
  -- Neural network layers operate on multivector feature channels.
  let v1 = 2 * e1 + 3 * e2 :: GA3
      v2 = 1 * e2 + 4 * e3 :: GA3
      v3 = scaleMV 2 e12    :: GA3
      v4 = scalarMV 5       :: GA3

      -- Construct 2x2 Multivector weight matrices:
      matA = mkMatrix 2 2 [v1, v2, v3, v4]
      matB = mkMatrix 2 2 [v4, v1, v2, v3]

      -- Multivector Matrix Multiply: C_ij = sum_k A_ik * B_kj
      matC = matMulMV matA matB

  putStrLn "--- 1. Multivector Matrix GEMM ---"
  putStrLn $ "Matrix A = " ++ show matA
  putStrLn $ "Matrix B = " ++ show matB
  putStrLn $ "Product C = A * B = " ++ show matC
  putStrLn $ "Entry C(0, 0) = " ++ show (matrixIndex matC 0 0)
  putStrLn ""

  -- 2. Structure Matrix Homomorphism
  -- Any multivector 'a' in Cl(3,0,0) corresponds to an 8x8 real matrix M(a).
  -- The map satisfies: M(a * b) == M(a) * M(b)
  let mA = toStructureMatrix v1
      mAb = toStructureMatrix (v1 <*> v2)

  putStrLn "--- 2. Structure Matrix Algebra Homomorphism ---"
  putStrLn $ "Multivector v1: " ++ show v1
  putStrLn $ "Structure matrix M(v1) entry (1,1): " ++ show (mA ! (1, 1))
  putStrLn $ "Structure matrix M(v1 * v2) entry (0,0): " ++ show (mAb ! (0, 0))
  putStrLn "This allows deploying Clifford models onto standard GPU GEMM accelerators without geometric overhead!"
  putStrLn ""

  -- 3. Rotor-Invariant Clifford Attention (CliffordNet LLM / Transformer)
  -- In physical systems and invariant ML models, attention weights between Query and Key tokens
  -- must not change when the coordinate frame is arbitrarily rotated by a rotor R.
  -- In Clifford attention: S(Q, K)_ij = sum_d <Q_id * ~K_jd>_0
  let q1 = 2 * e1 + 3 * e2  :: GA3
      q2 = 1 * e2 + 4 * e3  :: GA3
      k1 = 5 * e1 + 2 * e3  :: GA3
      k2 = 3 * e2 + 1 * e123:: GA3

      matQ = mkMatrix 2 2 [q1, q2, q2, q1]
      matK = mkMatrix 2 2 [k1, k2, k1, k2]

      -- Compute standard attention scores:
      scoresOrig = attentionScoresInvariant matQ matK

  putStrLn "--- 3. Rotor-Invariant Clifford Attention ---"
  putStrLn $ "Original Attention Score (0,0): " ++ show (scoresOrig ! (0, 0))
  putStrLn $ "Original Attention Score (0,1): " ++ show (scoresOrig ! (0, 1))

  -- Apply an arbitrary 3D rotation rotor R to all queries and keys:
  let rotor = normalizeRotor (scalarMV 3 `addMV` scaleMV 4 (e12 :: GA3))
      rotQ = mkMatrix 2 2 [sandwich rotor q1, sandwich rotor q2, sandwich rotor q2, sandwich rotor q1]
      rotK = mkMatrix 2 2 [sandwich rotor k1, sandwich rotor k2, sandwich rotor k1, sandwich rotor k2]
      scoresRot = attentionScoresInvariant rotQ rotK

  putStrLn $ "Rotated Attention Score  (0,0): " ++ show (scoresRot ! (0, 0))
  putStrLn $ "Rotated Attention Score  (0,1): " ++ show (scoresRot ! (0, 1))
  let maxDiff = maximum [abs (scoresOrig ! (i, j) - scoresRot ! (i, j)) | i <- [0 .. 1], j <- [0 .. 1]]
  putStrLn $ "Maximum difference under arbitrary rotation: " ++ show maxDiff
  putStrLn "Attention scores are mathematically 100% invariant under coordinate frame rotations!"
  putStrLn ""

  -- 4. Self-Adjoint Backpropagation Gradient Rules
  -- When computing gradients through multivector layers Z = X * Y:
  -- dL/dX = (dL/dZ) * ~Y
  -- dL/dY = ~X * (dL/dZ)
  let gradZ = 1 * e1 + 5 * e12 :: GA3
      gradX = adjointMulLeft gradZ v2
      gradY = adjointMulRight v1 gradZ

  putStrLn "--- 4. Self-Adjoint Backpropagation Gradients ---"
  putStrLn $ "Upstream Gradient dL/dZ = " ++ show gradZ
  putStrLn $ "Gradient dL/dX = dL/dZ * ~Y = " ++ show gradX
  putStrLn $ "Gradient dL/dY = ~X * dL/dZ = " ++ show gradY
  putStrLn ""
  putStrLn "Tutorial 5 completed successfully!"
