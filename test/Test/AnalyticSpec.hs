{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Test.AnalyticSpec
-- Description : Comprehensive test suite for Multivector Analytic Functions and Taylor Series
-- License     : BSD-3-Clause
module Test.AnalyticSpec (analyticTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit

import Clifford

type GA3 = DenseMV Cl300 Double
type STA = DenseMV Cl130 Double
type PGA3 = DenseMV Cl301 Double

analyticTests :: TestTree
analyticTests = testGroup "Multivector Analytic Functions & Taylor Series Suite"
  [ testGroup "Multivector Exponential Map (expMV)"
      [ testCase "exp(0) == 1 in 3D Euclidean GA" $ do
          let z = zeroMV :: GA3
              ez = expMV z
          ez ==~ scalarMV 1 @?= True

      , testCase "exp(0) == 1 in Spacetime Algebra" $ do
          let z = zeroMV :: STA
              ez = expMV z
          ez ==~ scalarMV 1 @?= True

      , testCase "exp(B) matches exact Euler closed-form rotor on pure bivectors" $ do
          let theta = pi / 3 :: Double
              bivec = scaleMV (- theta / 2) (e12 :: GA3)
              taylorRotor = expMV bivec
              eulerRotor = rotorExp bivec
          taylorRotor ==~ eulerRotor @?= True

      , testCase "exp(A + B) == exp(A) * exp(B) for commuting multivectors" $ do
          let a = scaleMV 0.2 (scalarMV 1) `addMV` scaleMV 0.3 (e12 :: GA3)
              b = scaleMV 0.5 (scalarMV 1) `subMV` scaleMV 0.1 (e12 :: GA3)
              expSum = expMV (a `addMV` b)
              expProd = expMV a <*> expMV b
          expSum ==~ expProd @?= True

      , testCase "exp(M) * exp(-M) == 1 (Inverse exponential law)" $ do
          let m = scaleMV 0.5 e1 `addMV` (scaleMV 0.8 e23 `subMV` scaleMV 0.4 (e123 :: GA3))
              posExp = expMV m
              negExp = expMV (negateMV m)
              identity = posExp <*> negExp
          identity ==~ scalarMV 1 @?= True

      , testCase "Nilpotent exponential in PGA: exp(e0) == 1 + e0" $ do
          let e0_pga = e0 :: PGA3
              expE0 = expMV e0_pga
              expected = scalarMV 1 `addMV` e0_pga
          expE0 ==~ expected @?= True
      ]

  , testGroup "Trigonometric Functions (sinMV, cosMV)"
      [ testCase "sin(0) == 0 and cos(0) == 1" $ do
          let z = zeroMV :: GA3
          sinMV z ==~ zeroMV @?= True
          cosMV z ==~ scalarMV 1 @?= True

      , testCase "Odd and even parity: sin(-M) == -sin(M), cos(-M) == cos(M)" $ do
          let m = scaleMV 0.6 e1 `addMV` scaleMV 0.3 (e12 :: GA3)
              negM = negateMV m
          sinMV negM ==~ negateMV (sinMV m) @?= True
          cosMV negM ==~ cosMV m @?= True

      , testCase "Pythagorean identity: cos^2(M) + sin^2(M) == 1" $ do
          let m = scaleMV 0.7 e12 `subMV` scaleMV 0.4 (e23 :: GA3)
              sVal = sinMV m
              cVal = cosMV m
              pythag = (cVal <*> cVal) `addMV` (sVal <*> sVal)
          pythag ==~ scalarMV 1 @?= True
      ]

  , testGroup "Hyperbolic Functions (sinhMV, coshMV, tanhMV)"
      [ testCase "sinh(0) == 0 and cosh(0) == 1" $ do
          let z = zeroMV :: GA3
          sinhMV z ==~ zeroMV @?= True
          coshMV z ==~ scalarMV 1 @?= True

      , testCase "Hyperbolic identity: cosh^2(M) - sinh^2(M) == 1" $ do
          let m = scaleMV 0.5 e1 `addMV` scaleMV 0.2 (e123 :: GA3)
              sh = sinhMV m
              ch = coshMV m
              hypId = (ch <*> ch) `subMV` (sh <*> sh)
          hypId ==~ scalarMV 1 @?= True

      , testCase "tanh(0) == Just 0" $ do
          let z = zeroMV :: GA3
          case tanhMV z of
            Just th -> th ==~ zeroMV @?= True
            Nothing -> assertFailure "tanh(0) returned Nothing"
      ]

  , testGroup "Horner Polynomial Evaluation (hornerMV)"
      [ testCase "hornerMV [1, 2, 3] e1 == 4 + 2*e1 (since e1^2 = 1)" $ do
          let v = e1 :: GA3
              pVal = hornerMV [1, 2, 3] v
              expected = (scalarMV 4 `addMV` scaleMV 2 e1) :: GA3
          pVal ==~ expected @?= True

      , testCase "hornerMV [5, 0, -2] e12 == 7 (since e12^2 = -1)" $ do
          let b = e12 :: GA3
              pVal = hornerMV [5, 0, -2] b
              expected = scalarMV 7 :: GA3
          pVal ==~ expected @?= True
      ]

  , testGroup "PGA 3D Continuous Screw Twist Integration"
      [ testCase "exp(omega * e12 + v * e03) integrates rigid body screw motion" $ do
          let theta = pi / 2
              dz = 4.0
              -- Twist: rotation in e12 plane + translation along Z (e03)
              twist = (scaleMV (- theta / 2) e12 `addMV` scaleMV (dz / 2) (basisMV 4 /\ basisMV 3)) :: PGA3
              screwMotor = expMV twist
              origin = pointPGA 0 0 0 :: PGA3
              transformed = sandwich screwMotor origin
              expectedPoint = pointPGA 0 0 4 :: PGA3
          transformed ==~ expectedPoint @?= True
      ]

  , testGroup "Finite Field Frobenius Endomorphism"
      [ testCase "Frobenius map Phi_2(M) = M^2 in 16D GF(2) algebra" $ do
          let a = basisMV 1 `addMV` basisMV 5 :: SparseMV Cl16_00 GF2
              frobA = frobeniusMV 2 a
              expected = a <*> a
          frobA @?= expected
      ]
  ]
