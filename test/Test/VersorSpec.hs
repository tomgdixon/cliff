{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.VersorSpec (versorTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit

import Clifford

versorTests :: TestTree
versorTests = testGroup "Versors, Rotors, Projections & Inverses"
  [ testCase "Bivector exponential rotor rotates e1 to e2 by 90 degrees" $ do
      let theta = pi / 2
          bivec = scaleMV (- theta / 2) (e12 :: DenseMV Cl300 Double)
          rotor = rotorExp bivec
          v1 = e1 :: DenseMV Cl300 Double
          rotated = sandwich rotor v1
          expected = e2 :: DenseMV Cl300 Double
      rotated ==~ expected @?= True

  , testCase "Rotor norm preservation: ||R v ~R|| == ||v||" $ do
      let rotor = normalizeRotor (scalarMV 3 `addMV` scaleMV 4 (e12 :: DenseMV Cl300 Double))
          v = (scaleMV 2 e1 `addMV` scaleMV 5 e2 `addMV` scaleMV 7 e3) :: DenseMV Cl300 Double
          vRot = sandwich rotor v
          normOrig = normSqMV v
          normRot = normSqMV vRot
      abs (normOrig - normRot) < 1e-10 @?= True

  , testCase "Versor inverse satisfies V * V^(-1) == 1" $ do
      let rotor = (scalarMV 0.8 `addMV` scaleMV 0.6 (e12 :: DenseMV Cl300 Double))
      case versorInverse rotor of
        Nothing -> assertFailure "Versor inverse calculation failed"
        Just rInv -> do
          let prod = rotor <*> rInv
          prod ==~ scalarMV 1 @?= True

  , testCase "General Gauss-Jordan inverse satisfies A * A^(-1) == 1" $ do
      let a = (scalarMV 2 `addMV` scaleMV 3 e1 `addMV` scaleMV 4 e12 `addMV` scaleMV 1 e123) :: DenseMV Cl300 Double
      case generalInverse a of
        Nothing -> assertFailure "General inverse calculation failed"
        Just aInv -> do
          let prod = a <*> aInv
          prod ==~ scalarMV 1 @?= True

  , testCase "Subspace Projection & Rejection onto 2D plane: v == Proj_B(v) + Rej_B(v)" $ do
      let v = (scaleMV 3 e1 `addMV` scaleMV 4 e2 `addMV` scaleMV 5 e3) :: DenseMV Cl300 Double
          bPlane = e12 :: DenseMV Cl300 Double
      case (projectMV v bPlane, rejectMV v bPlane) of
        (Just projV, Just rejV) -> do
          let expectedProj = (scaleMV 3 e1 `addMV` scaleMV 4 e2) :: DenseMV Cl300 Double
              expectedRej = (scaleMV 5 e3) :: DenseMV Cl300 Double
          projV ==~ expectedProj @?= True
          rejV ==~ expectedRej @?= True
          let vReconstructed = projV `addMV` rejV
          vReconstructed ==~ v @?= True
        _ -> assertFailure "projectMV or rejectMV returned Nothing on invertible 2-blade"

  , testCase "Hyperplane Reflection: reflect(v, n) reflects across normal n" $ do
      let v = (scaleMV 3 e1 `addMV` scaleMV 4 e2) :: DenseMV Cl300 Double
          normal = e1 :: DenseMV Cl300 Double  -- Normal to YZ plane
      case reflectMV v normal of
        Just vReflected -> do
          let expected = (scaleMV (-3) e1 `addMV` scaleMV 4 e2) :: DenseMV Cl300 Double
          vReflected ==~ expected @?= True
        Nothing -> assertFailure "reflectMV failed on unit vector normal"

  , testCase "Even and Odd Subalgebra Decomposition: a == evenPart(a) + oddPart(a)" $ do
      let a = (scalarMV 2 `addMV` scaleMV 3 e1 `addMV` scaleMV 4 e12 `addMV` scaleMV 5 e123) :: DenseMV Cl300 Double
          aEven = evenPart a
          aOdd = oddPart a
          expectedEven = (scalarMV 2 `addMV` scaleMV 4 e12) :: DenseMV Cl300 Double
          expectedOdd = (scaleMV 3 e1 `addMV` scaleMV 5 e123) :: DenseMV Cl300 Double
      aEven ==~ expectedEven @?= True
      aOdd ==~ expectedOdd @?= True
      (aEven `addMV` aOdd) ==~ a @?= True

  , testCase "Hestenes Inner Product (hestenesDot) symmetric contraction" $ do
      let v1 = scaleMV 2 e1 `addMV` scaleMV 3 e2 :: DenseMV Cl300 Double
          v2 = scaleMV 4 e1 `addMV` scaleMV 5 e2 :: DenseMV Cl300 Double
          dot12 = hestenesDot v1 v2
          expectedScalar = scalarMV (2*4 + 3*5) :: DenseMV Cl300 Double
      dot12 ==~ expectedScalar @?= True
  ]
