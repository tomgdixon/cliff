{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.MatrixSpec (matrixTests) where

import Prelude hiding ((<*>))
import Data.Array.IArray (Array, (!), array)
import Test.Tasty
import Test.Tasty.HUnit

import Clifford

matrixTests :: TestTree
matrixTests = testGroup "Multivector Matrices & LLM Primitives"
  [ testCase "Multivector Matrix Multiply (GEMM) associativity: (A * B) * C == A * (B * C)" $ do
      let v1 = (2 * e1 + 3 * e2) :: DenseMV Cl300 Double
          v2 = (1 * e2 + 4 * e3) :: DenseMV Cl300 Double
          v3 = scaleMV 2 e12 :: DenseMV Cl300 Double
          v4 = scalarMV 5 :: DenseMV Cl300 Double

          matA = mkMatrix 2 2 [v1, v2, v3, v4]
          matB = mkMatrix 2 2 [v4, v1, v2, v3]
          matC = mkMatrix 2 2 [v2, v3, v4, v1]

          lhs = matMulMV (matMulMV matA matB) matC
          rhs = matMulMV matA (matMulMV matB matC)

      lhs @?= rhs

  , testCase "Structure matrix algebra homomorphism: M(a * b) == M(a) . M(b)" $ do
      let a = (2 * e1 + 3 * e2) :: DenseMV Cl300 Double
          b = (1 * e2 + 4 * e3) :: DenseMV Cl300 Double
          ab = a <*> b

          mA = toStructureMatrix a
          mB = toStructureMatrix b
          mAb = toStructureMatrix ab

          -- Compute standard matrix product of mA and mB of size 8x8
          mMult :: Array (Int, Int) Double
          mMult = array ((0, 0), (7, 7))
            [ ((i, j), sum [mA ! (i, k) * mB ! (k, j) | k <- [0 .. 7]])
            | i <- [0 .. 7], j <- [0 .. 7]
            ]

      let isHomomorphic = all (\idx -> abs (mAb ! idx - mMult ! idx) < 1e-10) [((i, j)) | i <- [0 .. 7], j <- [0 .. 7]]
      isHomomorphic @?= True

  , testCase "Clifford Attention invariant scoring: S(R*Q, R*K) == S(Q, K)" $ do
      let rotor = normalizeRotor (scalarMV 3 `addMV` scaleMV 4 (e12 :: DenseMV Cl300 Double))
          q1 = 2 * e1 + 3 * e2 :: DenseMV Cl300 Double
          q2 = 1 * e2 + 4 * e3 :: DenseMV Cl300 Double
          k1 = 5 * e1 + 2 * e3 :: DenseMV Cl300 Double
          k2 = 3 * e2 + 1 * e123 :: DenseMV Cl300 Double

          matQ = mkMatrix 2 2 [q1, q2, q2, q1]
          matK = mkMatrix 2 2 [k1, k2, k1, k2]

          matRotQ = mkMatrix 2 2 [sandwich rotor q1, sandwich rotor q2, sandwich rotor q2, sandwich rotor q1]
          matRotK = mkMatrix 2 2 [sandwich rotor k1, sandwich rotor k2, sandwich rotor k1, sandwich rotor k2]

          scoresOrig = attentionScoresInvariant matQ matK
          scoresRot = attentionScoresInvariant matRotQ matRotK

      let diff = maximum [abs (scoresOrig ! (i, j) - scoresRot ! (i, j)) | i <- [0 .. 1], j <- [0 .. 1]]
      diff < 1e-10 @?= True

  , testCase "Self-adjoint gradient consistency: <X*Y, ~Z>_0 == <X, ~(Z*~Y)>_0" $ do
      let x = (2 * e1 + 3 * e2) :: DenseMV Cl300 Double
          y = (1 * e2 + 4 * e3) :: DenseMV Cl300 Double
          z = (5 * e1 + 2 * e2 + 7 * e12) :: DenseMV Cl300 Double

          lhs = (x <*> y) <|> reverseMV z
          rhs = x <|> reverseMV (adjointMulLeft z y)

      abs (lhs - rhs) < 1e-10 @?= True
  ]
