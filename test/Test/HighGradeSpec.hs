{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      : Test.HighGradeSpec
-- Description : Comprehensive test suite for High-Order Grades (Grades 1 to 64) and Sparse Bitmask Arithmetic
-- License     : BSD-3-Clause
--
-- This module tests:
-- 1. High-order basis blade construction, bitmask indexing, and grade evaluation up to Grade 64
-- 2. Wedge products of high-grade blades (disjoint additions, anticommutativity, nilpotency)
-- 3. Exact metric powers for every grade k in 1..64: B_k^2 == (-1)^(k*(k-1)/2)
-- 4. High-order contractions (Dorst/Lounesto Left & Right contractions)
-- 5. High-dimensional Poincaré duals in 16D, 32D, and 64D algebras
-- 6. QuickCheck property-based tests for SparseMV over Rational and GF(2) in 16D/32D
module Test.HighGradeSpec (highGradeTests) where

import Prelude hiding ((<*>))
import Data.Bits (testBit)
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Clifford

-- QuickCheck generators for High-Dimensional Sparse Multivectors
instance Arbitrary (SparseMV Cl16_00 Rational) where
  arbitrary = do
    numTerms <- (choose (1, 6) :: Gen Int)
    terms <- vectorOf numTerms $ do
      gradeVal <- (choose (0, 16) :: Gen Int)
      blades <- if gradeVal == 0
                  then return [scalarBlade]
                  else (vectorOf 3 (choose (1, 16) :: Gen Int) :: Gen [Int]) >>= return . pure . bladeFromIndices
      coeff <- (choose (-10, 10) :: Gen Integer)
      return (head blades, fromIntegral coeff)
    return (fromBladeList terms)

instance Arbitrary (SparseMV Cl32_00 GF2) where
  arbitrary = do
    numTerms <- (choose (1, 8) :: Gen Int)
    terms <- vectorOf numTerms $ do
      numIndices <- (choose (1, 8) :: Gen Int)
      indices <- (vectorOf numIndices (choose (1, 32) :: Gen Int) :: Gen [Int])
      return (bladeFromIndices indices, GF2 1)
    return (fromBladeList terms)

highGradeTests :: TestTree
highGradeTests = testGroup "High-Order Grades & 64-Bit Sparse Bitmask Suite"
  [ testGroup "High-Grade Basis Blade Construction (Grades 1 to 64)"
      [ testCase "Grade 8 blade has grade 8 and roundtrips indices" $ do
          let b8 = bladeFromIndices [1..8]
          bladeGrade b8 @?= 8
          bladeIndices b8 @?= [1..8]

      , testCase "Grade 16 blade has grade 16 and roundtrips indices" $ do
          let b16 = bladeFromIndices [1..16]
          bladeGrade b16 @?= 16
          bladeIndices b16 @?= [1..16]

      , testCase "Grade 32 blade has grade 32 and roundtrips indices" $ do
          let b32 = bladeFromIndices [1..32]
          bladeGrade b32 @?= 32
          bladeIndices b32 @?= [1..32]

      , testCase "Grade 48 blade has grade 48 and roundtrips indices" $ do
          let b48 = bladeFromIndices [1..48]
          bladeGrade b48 @?= 48
          bladeIndices b48 @?= [1..48]

      , testCase "Grade 64 pseudoscalar blade has grade 64" $ do
          let b64 = bladeFromIndices [1..64]
          bladeGrade b64 @?= 64
          bladeIndices b64 @?= [1..64]

      , testCase "64-bit boundary edge cases (bits 1, 32, 63, 64)" $ do
          let bEdges = bladeFromIndices [1, 32, 63, 64]
          bladeGrade bEdges @?= 4
          bladeIndices bEdges @?= [1, 32, 63, 64]
      ]

  , testGroup "High-Grade Wedge Product (Exterior Algebra)"
      [ testCase "Grade 8 ^ Grade 16 produces Grade 24 blade" $ do
          let b8  = bladeFromIndices [1..8]
              b16 = bladeFromIndices [9..24]
              mv8  = bladeMV b8 1  :: SparseMV Cl32_00 Double
              mv16 = bladeMV b16 1 :: SparseMV Cl32_00 Double
              res  = mv8 /\ mv16
          grades res @?= [24]
          getBlade (bladeFromIndices [1..24]) res @?= 1.0

      , testCase "Grade 16 ^ Grade 16 produces Grade 32 blade" $ do
          let b16_A = bladeFromIndices [1..16]
              b16_B = bladeFromIndices [17..32]
              mvA   = bladeMV b16_A 1 :: SparseMV Cl32_00 Double
              mvB   = bladeMV b16_B 1 :: SparseMV Cl32_00 Double
              res   = mvA /\ mvB
          grades res @?= [32]
          getBlade (bladeFromIndices [1..32]) res @?= 1.0

      , testCase "Grade 32 ^ Grade 32 produces Grade 64 pseudoscalar" $ do
          let b32_A = bladeFromIndices [1..32]
              b32_B = bladeFromIndices [33..64]
              mvA   = bladeMV b32_A 1 :: SparseMV Cl64_00 Double
              mvB   = bladeMV b32_B 1 :: SparseMV Cl64_00 Double
              res   = mvA /\ mvB
          grades res @?= [64]
          getBlade (bladeFromIndices [1..64]) res @?= 1.0

      , testCase "High-grade blade nilpotency: B_16 ^ B_16 == 0" $ do
          let b16 = bladeFromIndices [1..16]
              mv16 = bladeMV b16 1 :: SparseMV Cl32_00 Double
          (mv16 /\ mv16) @?= (zeroMV :: SparseMV Cl32_00 Double)

      , testCase "Overlapping high-grade blades wedge to zero: B[1..10] ^ B[10..20] == 0" $ do
          let bA = bladeFromIndices [1..10]
              bB = bladeFromIndices [10..20]
              mvA = bladeMV bA 1 :: SparseMV Cl32_00 Double
              mvB = bladeMV bB 1 :: SparseMV Cl32_00 Double
          (mvA /\ mvB) @?= (zeroMV :: SparseMV Cl32_00 Double)
      ]

  , testGroup "High-Grade Metric Squares: B_k^2 == (-1)^(k*(k-1)/2)"
      [ testCase ("Metric square of Grade " ++ show k ++ " blade") $ do
          let bk = bladeFromIndices [1..k]
              mvk = bladeMV bk 1 :: SparseMV Cl64_00 Double
              sq = mvk <*> mvk
              expectedSign = if testBit (k * (k - 1) `div` 2) 0 then -1.0 else 1.0
          sq @?= scalarMV expectedSign
      | k <- [1, 2, 3, 4, 7, 8, 11, 12, 15, 16, 23, 24, 31, 32, 47, 48, 63, 64]
      ]

  , testGroup "High-Grade Contractions (Left & Right)"
      [ testCase "Left contraction: Grade 8 . Grade 16 yields Grade 8" $ do
          let b8  = bladeFromIndices [1..8]
              b16 = bladeFromIndices [1..16]
              mv8  = bladeMV b8 1  :: SparseMV Cl32_00 Double
              mv16 = bladeMV b16 1 :: SparseMV Cl32_00 Double
              res  = mv8 <.> mv16
          grades res @?= [8]
          let expectedBlade = bladeFromIndices [9..16]
          getBlade expectedBlade res /= 0 @?= True

      , testCase "Left contraction of higher grade onto lower grade is strictly 0" $ do
          let b16 = bladeFromIndices [1..16]
              b8  = bladeFromIndices [1..8]
              mv16 = bladeMV b16 1 :: SparseMV Cl32_00 Double
              mv8  = bladeMV b8 1  :: SparseMV Cl32_00 Double
          (mv16 <.> mv8) @?= (zeroMV :: SparseMV Cl32_00 Double)
      ]

  , testGroup "High-Grade Poincaré Duals in 16D, 32D, 64D"
      [ testCase "Dual of 16D scalar is 16D pseudoscalar" $ do
          let s = scalarMV 1 :: SparseMV Cl16_00 Double
              dual_s = poincareDualMV s
          grades dual_s @?= [16]
          getBlade (bladeFromIndices [1..16]) dual_s @?= 1.0

      , testCase "Dual of 32D scalar is 32D pseudoscalar" $ do
          let s = scalarMV 1 :: SparseMV Cl32_00 Double
              dual_s = poincareDualMV s
          grades dual_s @?= [32]
          getBlade (bladeFromIndices [1..32]) dual_s @?= 1.0

      , testCase "Dual of 64D scalar is 64D pseudoscalar" $ do
          let s = scalarMV 1 :: SparseMV Cl64_00 Double
              dual_s = poincareDualMV s
          grades dual_s @?= [64]
          getBlade (bladeFromIndices [1..64]) dual_s @?= 1.0

      , testCase "Double Poincaré dual identity J(J(A)) == (-1)^(k*(n-k)) A in 16D" $ do
          let b7 = bladeFromIndices [1, 3, 5, 7, 9, 11, 13]
              mv7 = bladeMV b7 1 :: SparseMV Cl16_00 Double
              dual1 = poincareDualMV mv7
              dual2 = poincareDualMV dual1
              k = 7 :: Int
              n = 16 :: Int
              expectedSign = if testBit (k * (n - k)) 0 then -1.0 else 1.0
          dual2 @?= scaleMV expectedSign mv7
      ]

  , testGroup "High-Dimensional QuickCheck Properties (16D & 32D Sparse)"
      [ testProperty "16D SparseMV Associativity over Rational: (a * b) * c == a * (b * c)" $
          \(a :: SparseMV Cl16_00 Rational) b c ->
            (a <*> b) <*> c == a <*> (b <*> c)

      , testProperty "16D SparseMV Distributivity over Rational: a * (b + c) == a*b + a*c" $
          \(a :: SparseMV Cl16_00 Rational) b c ->
            a <*> (b `addMV` c) == (a <*> b) `addMV` (a <*> c)

      , testProperty "16D SparseMV Reversion Anti-Automorphism: ~(a * b) == ~b * ~a" $
          \(a :: SparseMV Cl16_00 Rational) b ->
            reverseMV (a <*> b) == reverseMV b <*> reverseMV a

      , testProperty "32D SparseMV Associativity over GF(2): (a * b) * c == a * (b * c)" $
          \(a :: SparseMV Cl32_00 GF2) b c ->
            (a <*> b) <*> c == a <*> (b <*> c)

      , testProperty "32D SparseMV Distributivity over GF(2): a * (b + c) == a*b + a*c" $
          \(a :: SparseMV Cl32_00 GF2) b c ->
            a <*> (b `addMV` c) == (a <*> b) `addMV` (a <*> c)
      ]
  ]
