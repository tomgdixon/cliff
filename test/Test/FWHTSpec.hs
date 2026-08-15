{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Test.FWHTSpec
-- Description : Comprehensive test suite for Fast Walsh-Hadamard Transforms (FWHT)
-- License     : BSD-3-Clause
module Test.FWHTSpec (fwhtTests) where

import Prelude hiding ((<*>))
import Data.Bits (xor)
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Clifford

fwhtTests :: TestTree
fwhtTests = testGroup "Fast Walsh-Hadamard Transform (FWHT) & Dyadic Convolution Suite"
  [ testGroup "Dense Multivector (Intra-MV) FWHT"
      [ testCase "Exact 4-element Sylvester Hadamard basis check on Cl(2,0,0)" $ do
          -- v = 1*1 + 2*e1 + 3*e2 + 4*e12
          let v = fromBladeList (zip (allBlades 2) [1, 2, 3, 4]) :: DenseMV Cl200 Double
              fw = fwhtMV v
              coeffs = [getBlade b fw | b <- allBlades 2]
          -- Expected: H4 * [1, 2, 3, 4] = [10, -2, -4, 0]
          coeffs @?= [10.0, -2.0, -4.0, 0.0]

      , testCase "DenseMV FWHT Invertibility Roundtrip on Cl(2,0,0) (N=4)" $ do
          let v = fromBladeList (zip (allBlades 2) [3.5, -2.1, 7.8, 1.4]) :: DenseMV Cl200 Double
              roundtrip = ifwhtMV (fwhtMV v)
          roundtrip ==~ v @?= True

      , testCase "DenseMV FWHT Invertibility Roundtrip on Cl(3,0,0) (N=8)" $ do
          let v = fromBladeList (zip (allBlades 3) [1.0, -2.5, 3.2, 0.8, -4.1, 5.6, -1.9, 6.3]) :: DenseMV Cl300 Double
              roundtrip = ifwhtMV (fwhtMV v)
          roundtrip ==~ v @?= True

      , testCase "DenseMV FWHT Invertibility Roundtrip on Cl(1,3,0) STA (N=16)" $ do
          let vals = [fromIntegral i * 1.1 - 8.0 | i <- [1..16]]
              v = fromBladeList (zip (allBlades 4) vals) :: DenseMV Cl130 Double
              roundtrip = ifwhtMV (fwhtMV v)
          roundtrip ==~ v @?= True

      , testCase "DenseMV Dyadic (XOR) Convolution Theorem vs Direct Quadratic Sum" $ do
          let aVals = [2.0, -1.5, 3.0, 4.5]
              bVals = [1.0, 2.5, -0.5, 3.0]
              a = fromBladeList (zip (allBlades 2) aVals) :: DenseMV Cl200 Double
              b = fromBladeList (zip (allBlades 2) bVals) :: DenseMV Cl200 Double
              fastConv = dyadicConvolveMV a b
              -- Direct quadratic sum: c_k = sum_{i xor j == k} a_i * b_j
              directCoeff k = sum [ (aVals !! i) * (bVals !! j) | i <- [0..3], j <- [0..3], i `xor` j == k ]
              expected = fromBladeList (zip (allBlades 2) [directCoeff k | k <- [0..3]]) :: DenseMV Cl200 Double
          fastConv ==~ expected @?= True
      ]

  , testGroup "Multivector Array (Inter-MV) FWHT"
      [ testCase "1D Array FWHT Invertibility (N=4, Cl300)" $ do
          let mvs = [ scalarMV (fromIntegral i) :: DenseMV Cl300 Double | i <- [1..4] ]
              arr = fromListND [4] mvs
              roundtrip = ifwht1D (fwht1D arr)
          and [ indexND roundtrip [i] ==~ indexND arr [i] | i <- [0..3] ] @?= True

      , testCase "1D Array FWHT Invertibility (N=16, Cl300)" $ do
          let mvs = [ scalarMV (sin (fromIntegral i)) `addMV` bladeMV (basisBlade 1) (cos (fromIntegral i))
                    | i <- [1..16] :: [Int]
                    ] :: [DenseMV Cl300 Double]
              arr = fromListND [16] mvs
              roundtrip = ifwht1D (fwht1D arr)
          and [ indexND roundtrip [i] ==~ indexND arr [i] | i <- [0..15] ] @?= True

      , testCase "2D Grid FWHT Invertibility (4x4, Cl200)" $ do
          let grid = [ scalarMV (fromIntegral (r * 4 + c)) :: DenseMV Cl200 Double
                     | r <- [0..3], c <- [0..3]
                     ]
              arr = fromListND [4, 4] grid
              roundtrip = ifwht2D (fwht2D arr)
          and [ indexND roundtrip [r, c] ==~ indexND arr [r, c] | r <- [0..3], c <- [0..3] ] @?= True

      , testCase "2D Grid FWHT Invertibility (8x8, Cl300)" $ do
          let grid = [ scalarMV (fromIntegral (r * 8 + c)) `addMV` bladeMV (basisBlade 2) (fromIntegral (r - c))
                     | r <- [0..7], c <- [0..7]
                     ] :: [DenseMV Cl300 Double]
              arr = fromListND [8, 8] grid
              roundtrip = ifwht2D (fwht2D arr)
          and [ indexND roundtrip [r, c] ==~ indexND arr [r, c] | r <- [0..7], c <- [0..7] ] @?= True

      , testCase "1D Array Dyadic Convolution Theorem vs Direct Convolution" $ do
          let aList = [ scalarMV 2 `addMV` bladeMV (basisBlade 1) 3
                      , scalarMV (-1) `addMV` bladeMV (basisBlade 2) 4
                      , scalarMV 5 `addMV` bladeMV (basisBlade 1) (-2)
                      , scalarMV 0 `addMV` bladeMV (basisBlade 3) 1
                      ] :: [DenseMV Cl300 Double]
              bList = [ scalarMV 1 `addMV` bladeMV (basisBlade 2) (-3)
                      , scalarMV 4 `addMV` bladeMV (basisBlade 1) 2
                      , scalarMV (-2) `addMV` bladeMV (basisBlade 3) 5
                      , scalarMV 3 `addMV` bladeMV (basisBlade 1) 1
                      ] :: [DenseMV Cl300 Double]
              arrA = fromListND [4] aList
              arrB = fromListND [4] bList
              fastConv = dyadicConvolve1D (<*>) arrA arrB
              -- Direct calculation: C(k) = sum_{i xor j == k} A(i) <*> B(j)
              directAt k = foldr addMV zeroMV [ (aList !! i) <*> (bList !! j) | i <- [0..3], j <- [0..3], i `xor` j == k ]
          and [ indexND fastConv [k] ==~ directAt k | k <- [0..3] ] @?= True
      ]

  , testGroup "Galois Field Exact FWHT (F_p)"
      [ testCase "DenseMV FWHT Invertibility in F_17 (Cl(2,0,0))" $ do
          let v = fromBladeList (zip (allBlades 2) [GFp 3, GFp 7, GFp 12, GFp 1]) :: DenseMV Cl200 (GFp 17)
              roundtrip = ifwhtMV (fwhtMV v)
          roundtrip @?= v

      , testCase "1D Array FWHT Invertibility in F_257 (N=8, Cl300)" $ do
          let mvs = [ scalarMV (GFp (fromIntegral i * 13 `mod` 257)) :: DenseMV Cl300 (GFp 257) | i <- [0..7] ]
              arr = fromListND [8] mvs
              roundtrip = ifwht1D (fwht1D arr)
          and [ indexND roundtrip [i] == indexND arr [i] | i <- [0..7] ] @?= True
      ]
  ]
