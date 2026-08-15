{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.PGASpec (pgaTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit

import Clifford

pgaTests :: TestTree
pgaTests = testGroup "3D Projective Geometric Algebra (PGA) & Conformal GA (CGA)"
  [ testGroup "3D PGA Cl(3,0,1) Kinematics & Incidence"
      [ testCase "Degenerate null basis: e0^2 == 0" $ do
          let e0_pga = e0 :: DenseMV Cl301 Double
          (e0_pga <*> e0_pga) @?= scalarMV 0

      , testCase "Poincaré dual on degenerate blade e0 produces -e123 without error" $ do
          let e0_pga = e0 :: DenseMV Cl301 Double
              dual_e0 = poincareDualMV e0_pga
              expected = e123 :: DenseMV Cl301 Double
          dual_e0 @?= negateMV expected

      , testCase "PGA translation motor translates origin point (0,0,0) to (1,2,3)" $ do
          let origin = pointPGA 0 0 0 :: DenseMV Cl301 Double
              tMotor = translatorPGA 1 2 3 :: DenseMV Cl301 Double
              translatedPoint = sandwich tMotor origin
              expectedPoint = pointPGA 1 2 3 :: DenseMV Cl301 Double
          translatedPoint ==~ expectedPoint @?= True

      , testCase "PGA rotation motor rotates point (1,0,0) by 90 degrees around Z axis" $ do
          let pInit = pointPGA 1 0 0 :: DenseMV Cl301 Double
              zAxisBivec = e12 :: DenseMV Cl301 Double
              rMotor = rotatorPGA (pi / 2) zAxisBivec
              pRot = sandwich rMotor pInit
              pExpected = pointPGA 0 1 0 :: DenseMV Cl301 Double
          pRot ==~ pExpected @?= True

      , testCase "PGA combined rigid body motor (M = T * R)" $ do
          let pInit = pointPGA 1 0 0 :: DenseMV Cl301 Double
              tMotor = translatorPGA 0 0 5 :: DenseMV Cl301 Double
              rRotor = rotatorPGA (pi / 2) (e12 :: DenseMV Cl301 Double)
              m = motorPGA tMotor rRotor
              pFinal = sandwich m pInit
              pExpected = pointPGA 0 1 5 :: DenseMV Cl301 Double
          pFinal ==~ pExpected @?= True

      , testCase "PGA 3D Line join (P1 \\/ P2) and Plane meet (pi1 /\\ pi2)" $ do
          let p1 = pointPGA 0 0 0 :: DenseMV Cl301 Double
              p2 = pointPGA 0 0 1 :: DenseMV Cl301 Double
              lineZ = linePGA p1 p2
          grades lineZ @?= [2]
          let planeX = planePGA 1 0 0 0 :: DenseMV Cl301 Double
              planeY = planePGA 0 1 0 0 :: DenseMV Cl301 Double
              intersectLine = intersectPlanesPGA planeX planeY
          grades intersectLine @?= [2]

      , testCase "LaTeX export formatting" $ do
          let v = (2 * e1 - 3 * e12) :: DenseMV Cl300 Double
              latexStr = toLatex v
          latexStr @?= "2.0 e_{1} - 3.0 e_{12}"
      ]

  , testGroup "Conformal Geometric Algebra CGA Cl(4,1,0)"
      [ testCase "CGA null bases: e_inf^2 == 0 and e_0^2 == 0" $ do
          let eInf = eInfCGA :: DenseMV Cl410 Double
              e0_cga = e0CGA :: DenseMV Cl410 Double
          (eInf <*> eInf) ==~ scalarMV 0 @?= True
          (e0_cga <*> e0_cga) ==~ scalarMV 0 @?= True

      , testCase "CGA inner product: e_inf . e_0 == -1" $ do
          let eInf = eInfCGA :: DenseMV Cl410 Double
              e0_cga = e0CGA :: DenseMV Cl410 Double
              dot = eInf <|> e0_cga
          abs (dot - (-1.0)) < 1e-12 @?= True

      , testCase "CGA point Euclidean distance: d^2(P1, P2) == -2 (P1 . P2)" $ do
          let p1 = pointCGA 1 2 3 :: DenseMV Cl410 Double
              p2 = pointCGA 4 6 3 :: DenseMV Cl410 Double
              sqDist = sqDistanceCGA p1 p2
          -- (4-1)^2 + (6-2)^2 + (3-3)^2 = 9 + 16 = 25
          abs (sqDist - 25.0) < 1e-12 @?= True

      , testCase "CGA sphere null boundary: P on sphere S(c, r) satisfies P . S == 0" $ do
          let center = pointCGA 0 0 0 :: DenseMV Cl410 Double
              sphere = sphereCGA center 5.0 :: DenseMV Cl410 Double
              ptOnSphere = pointCGA 3 4 0 :: DenseMV Cl410 Double  -- 3^2 + 4^2 = 25 = 5^2
              dotProd = ptOnSphere <|> sphere
          abs dotProd < 1e-12 @?= True
      ]
  ]
