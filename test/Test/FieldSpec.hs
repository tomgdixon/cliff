{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.FieldSpec (fieldTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Clifford

fieldTests :: TestTree
fieldTests = testGroup "Finite Fields & Higher-Order Algebras"
  [ testGroup "Binary Field GF(2)"
      [ testCase "Addition is XOR (1 + 1 == 0, 1 + 0 == 1)" $ do
          (GF2 1 + GF2 1) @?= GF2 0
          (GF2 1 + GF2 0) @?= GF2 1

      , testCase "Multiplication is AND (1 * 1 == 1, 1 * 0 == 0)" $ do
          (GF2 1 * GF2 1) @?= GF2 1
          (GF2 1 * GF2 0) @?= GF2 0

      , testCase "Characteristic 2 negation: -1 == 1" $ do
          negate (GF2 1) @?= GF2 1

      , testCase "Metric scale injection: fromMetricScale (-1) == GF2 1" $ do
          (fromMetricScale (-1) :: GF2) @?= GF2 1
          (fromMetricScale 0 :: GF2) @?= GF2 0
          (fromMetricScale 1 :: GF2) @?= GF2 1
      ]

  , testGroup "Prime Modular Field GF(p)"
      [ testCase "Modular inverse in GF(7): 3 * 3^(-1) == 1 (mod 7)" $ do
          let x = GFp 3 :: GFp 7
              xInv = recip x
          (x * xInv) @?= GFp 1

      , testCase "Modular inverse in GF(13): 5 * 5^(-1) == 1 (mod 13)" $ do
          let x = GFp 5 :: GFp 13
              xInv = recip x
          (x * xInv) @?= GFp 1
      ]

  , testGroup "Galois Extension Field GF(2^8) (AES Field)"
      [ testCase "Multiplicative inverse in GF(2^8) with poly 0x11B" $ do
          let a = GF2m 0x53 :: GF2m 8 0x11B
              aInv = recip a
          (a * aInv) @?= GF2m 1
      ]

  , testGroup "Higher-Order Clifford Algebras (n=16, n=32) over GF(2)"
      [ testCase "Geometric product associativity in 16D Clifford algebra over GF(2)" $ do
          let v1 = basisMV 1 `addMV` basisMV 5 `addMV` basisMV 12 :: SparseMV Cl16_00 GF2
              v2 = basisMV 5 `addMV` basisMV 16 :: SparseMV Cl16_00 GF2
              v3 = basisMV 1 `addMV` basisMV 16 :: SparseMV Cl16_00 GF2
              lhs = (v1 <*> v2) <*> v3
              rhs = v1 <*> (v2 <*> v3)
          lhs @?= rhs

      , testCase "Blade XOR reduction in 32D Clifford algebra over GF(2)" $ do
          let b1 = basisMV 3 :: SparseMV Cl32_00 GF2
              b2 = basisMV 3 :: SparseMV Cl32_00 GF2
          (b1 <*> b2) @?= scalarMV (GF2 1)
      ]
  ]
