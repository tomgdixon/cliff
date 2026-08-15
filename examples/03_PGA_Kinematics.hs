{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : Example 3 - 3D Projective Geometric Algebra (PGA Cl(3,0,1)) and Rigid Body Kinematics
--
-- This tutorial example covers:
-- 1. Understanding 3D PGA: 3 Euclidean dimensions + 1 Null / Degenerate dimension (e0^2 = 0)
-- 2. Geometric Primitives: Points (Trivectors) and Planes (Vectors)
-- 3. Translation Motors: T = 1 + 0.5 * (dx e01 + dy e02 + dz e03)
-- 4. Rotation Motors: R = cos(theta/2) - sin(theta/2) B
-- 5. Rigid Body Screw Motion (SE(3) Kinematics)
module Main where

import Prelude hiding ((<*>))
import Clifford

-- PGA Cl(3,0,1): 3 positive metric dimensions, 0 negative, 1 degenerate null dimension (e0^2 = 0)
type PGA3 = DenseMV Cl301 Double

main :: IO ()
main = do
  putStrLn "========================================================================="
  putStrLn "  Clifford Algebra Tutorial: Example 3 - 3D PGA & Rigid Body Kinematics"
  putStrLn "=========================================================================\n"

  -- 1. Concept: In 3D PGA, space has 4 basis vectors: e1, e2, e3 (Euclidean) and e0 (Null plane at infinity)
  -- The null dimension satisfies e0 * e0 = 0.
  -- This allows translations and rotations to be unified into a single algebraic group: Motors (SE(3))!

  -- 2. Constructing Geometric Points and Planes
  -- In PGA:
  -- - Origin Point (0,0,0) is represented by the trivector: e123
  -- - A Point (x,y,z) is represented by: e123 + x e023 + y e031 + z e012
  -- - A Plane (ax + by + cz + d = 0) is represented by: a e1 + b e2 + c e3 + d e0
  let p0 = pointPGA 0 0 0 :: PGA3
      p1 = pointPGA 1 2 3 :: PGA3
      groundPlane = planePGA 0 0 1 0 :: PGA3 -- z = 0 plane

  putStrLn "--- 1. Points and Planes ---"
  putStrLn $ "Origin point P0     = " ++ show p0
  putStrLn $ "Point P1 (1, 2, 3)  = " ++ show p1
  putStrLn $ "Ground plane (z=0)  = " ++ show groundPlane
  putStrLn ""

  -- 3. Translation Motors (Pure Translation)
  -- A translation motor translates any geometric entity (points, lines, planes) by displacement (dx, dy, dz):
  -- T = 1 + 0.5 * (dx e01 + dy e02 + dz e03)
  let tMotor = translatorPGA 5 10 (-2) :: PGA3
      p0_translated = sandwich tMotor p0

  putStrLn "--- 2. Translation Motor ---"
  putStrLn $ "Translation Motor T(5, 10, -2) = " ++ show tMotor
  putStrLn $ "Translating Origin P0 -> " ++ show p0_translated
  putStrLn $ "Expected point (5, 10, -2)    = " ++ show (pointPGA 5 10 (-2) :: PGA3)
  putStrLn $ "Match: " ++ show (p0_translated == pointPGA 5 10 (-2))
  putStrLn ""

  -- 4. Rotation Motors (Pure Rotation around Line/Axis)
  -- Rotate by 90 degrees around the Z axis:
  let rMotor = rotatorPGA (pi / 2) (e12 :: PGA3)
      p1_rotated = sandwich rMotor p1

  putStrLn "--- 3. Rotation Motor ---"
  putStrLn $ "Rotation Motor R(90° around Z) = " ++ show rMotor
  putStrLn $ "Point P1 before rotation       = " ++ show p1 ++ "  (x=1, y=2, z=3)"
  putStrLn $ "Point P1 after 90° Z rotation  = " ++ show p1_rotated ++ "  (x=-2, y=1, z=3)"
  putStrLn ""

  -- 5. Rigid Body Screw Motion (Translation + Rotation)
  -- Composing a translation and a rotation creates a general Rigid Body Motor: M = T * R
  let screwMotor = tMotor <*> rMotor
      p0_screw = sandwich screwMotor p0

  putStrLn "--- 4. Full Rigid Body Screw Motion (SE(3)) ---"
  putStrLn $ "Combined Screw Motor M = T * R = " ++ show screwMotor
  putStrLn $ "Transformed Origin P0 -> " ++ show p0_screw
  putStrLn ""
  putStrLn "PGA allows transforming points, lines, planes, and camera rays with ONE identical formula: M * X * ~M"
  putStrLn ""
  putStrLn "Tutorial 3 completed successfully!"
