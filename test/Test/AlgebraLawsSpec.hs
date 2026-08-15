{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.AlgebraLawsSpec (algebraLawsTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit

import Clifford
import Test.Arbitrary ()

algebraLawsTests :: TestTree
algebraLawsTests = testGroup "Algebraic Laws & Axioms"
  [ testGroup "3D Euclidean GA (Cl(3,0,0)) Properties"
      [ testProperty "Geometric product associativity: (a * b) * c == a * (b * c)" $
          \(a :: DenseMV Cl300 Double) b c ->
            let lhs = (a <*> b) <*> c
                rhs = a <*> (b <*> c)
            in lhs ==~ rhs

      , testProperty "Left distributivity: a * (b + c) == (a * b) + (a * c)" $
          \(a :: DenseMV Cl300 Double) b c ->
            let lhs = a <*> (b `addMV` c)
                rhs = (a <*> b) `addMV` (a <*> c)
            in lhs ==~ rhs

      , testProperty "Right distributivity: (a + b) * c == (a * c) + (b * c)" $
          \(a :: DenseMV Cl300 Double) b c ->
            let lhs = (a `addMV` b) <*> c
                rhs = (a <*> c) `addMV` (b <*> c)
            in lhs ==~ rhs

      , testProperty "Reversion anti-automorphism: ~(a * b) == ~b * ~a" $
          \(a :: DenseMV Cl300 Double) b ->
            let lhs = reverseMV (a <*> b)
                rhs = reverseMV b <*> reverseMV a
            in lhs ==~ rhs

      , testProperty "Grade involution automorphism: ^(a * b) == ^a * ^b" $
          \(a :: DenseMV Cl300 Double) b ->
            let lhs = involMV (a <*> b)
                rhs = involMV a <*> involMV b
            in lhs ==~ rhs

      , testProperty "Double reversion identity: ~~a == a" $
          \(a :: DenseMV Cl300 Double) ->
            reverseMV (reverseMV a) ==~ a
      ]

  , testGroup "Spacetime Algebra (Cl(1,3,0)) Metric Verification"
      [ testCase "Positive time dimension: e1^2 == +1" $ do
          let e1_sta = basisMV 1 :: DenseMV Cl130 Double
          (e1_sta <*> e1_sta) @?= scalarMV 1

      , testCase "Negative spatial dimensions: e2^2 == -1, e3^2 == -1, e4^2 == -1" $ do
          let e2_sta = basisMV 2 :: DenseMV Cl130 Double
              e3_sta = basisMV 3 :: DenseMV Cl130 Double
              e4_sta = basisMV 4 :: DenseMV Cl130 Double
          (e2_sta <*> e2_sta) @?= scalarMV (-1)
          (e3_sta <*> e3_sta) @?= scalarMV (-1)
          (e4_sta <*> e4_sta) @?= scalarMV (-1)

      , testCase "Anticommutativity of orthogonal vectors: e1*e2 == -e2*e1" $ do
          let e1_sta = basisMV 1 :: DenseMV Cl130 Double
              e2_sta = basisMV 2 :: DenseMV Cl130 Double
          (e1_sta <*> e2_sta) @?= negateMV (e2_sta <*> e1_sta)
      ]

  , testGroup "Contraction & Wedge Laws"
      [ testCase "Wedge product of 1-vectors is antisymmetric" $ do
          let v1 = basisMV 1 :: DenseMV Cl300 Double
              v2 = basisMV 2 :: DenseMV Cl300 Double
          (v1 /\ v2) @?= negateMV (v2 /\ v1)

      , testCase "Wedge product of identical vectors is zero" $ do
          let v1 = basisMV 1 :: DenseMV Cl300 Double
          (v1 /\ v1) @?= (zeroMV :: DenseMV Cl300 Double)

      , testCase "Left contraction identity: e2 . (e1 ^ e2) == -e1" $ do
          let e1_val = basisMV 1 :: DenseMV Cl300 Double
              e2_val = basisMV 2 :: DenseMV Cl300 Double
              bivec = e1_val /\ e2_val
          (e2_val <.> bivec) @?= negateMV e1_val
      ]
  ]
