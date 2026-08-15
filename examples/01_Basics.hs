{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : Example 1 - Basics of Multivectors and Geometric Algebra
--
-- This tutorial example covers:
-- 1. Constructing multivectors in 3D Euclidean space Cl(3,0,0)
-- 2. The Geometric Product: ab = a.b + a^b
-- 3. Wedge (Exterior) Product and Bivectors
-- 4. Involutions: Reversion (~a), Grade Involution, and Poincaré Duals
module Main where

import Prelude hiding ((<*>))
import Clifford

-- We work in 3D Euclidean Geometric Algebra Cl(3,0,0) with Double precision.
type GA3 = DenseMV Cl300 Double

main :: IO ()
main = do
  putStrLn "============================================================"
  putStrLn "  Clifford Algebra Tutorial: Example 1 - Multivector Basics"
  putStrLn "============================================================\n"

  -- 1. Constructing Basis Vectors
  -- In 3D Euclidean space, we have 3 orthogonal unit vectors: e1, e2, e3
  -- satisfying e1^2 = 1, e2^2 = 1, e3^2 = 1 and e_i * e_j = - e_j * e_i (for i /= j)
  let u = 2 * e1 + 3 * e2          :: GA3
      v = 1 * e2 + 4 * e3          :: GA3

  putStrLn "--- 1. Vector Definitions ---"
  putStrLn $ "Vector u = " ++ show u
  putStrLn $ "Vector v = " ++ show v
  putStrLn ""

  -- 2. The Fundamental Geometric Product (<*>)
  -- In Clifford Algebra, multiplying two vectors combines their dot product
  -- (symmetric/scalar part) and their wedge product (antisymmetric/bivector part):
  -- u * v = (u . v) + (u ^ v)
  let uv = u <*> v

  putStrLn "--- 2. Geometric Product (u * v) ---"
  putStrLn $ "u <*> v = " ++ show uv
  putStrLn "Notice the result contains both a scalar part and bivector parts (e12, e13, e23)!"
  putStrLn ""

  -- 3. Exterior (Wedge) Product (/\) and Dot Product (<.>)
  -- The wedge product produces an oriented area segment (Bivector):
  let u_wedge_v = u /\ v
      u_dot_v   = u <.> v
      u_scalar  = u <|> v

  putStrLn "--- 3. Wedge vs Dot vs Scalar Products ---"
  putStrLn $ "Wedge product (u /\\ v)  = " ++ show u_wedge_v ++ "  (Oriented Area)"
  putStrLn $ "Dot product   (u <.> v) = " ++ show u_dot_v   ++ "  (Symmetric Inner Product)"
  putStrLn $ "Scalar product(u <|> v) = " ++ show u_scalar  ++ "  (Scalar Value)"
  putStrLn ""

  -- Verify the fundamental identity: u * v == (u . v) + (u ^ v)
  let reconstructed = u_dot_v `addMV` u_wedge_v
  putStrLn $ "Is u*v == (u.v) + (u^v)? " ++ show (uv == reconstructed)
  putStrLn ""

  -- 4. Grade Projections
  -- A general multivector can contain terms of grade 0 (scalars), 1 (vectors),
  -- 2 (bivectors), and 3 (trivectors / pseudoscalars).
  putStrLn "--- 4. Grade Decomposition of (u * v) ---"
  putStrLn $ "Active grades: " ++ show (grades uv)
  putStrLn $ "Grade 0 (Scalar):   " ++ show (gradeProj 0 uv)
  putStrLn $ "Grade 1 (Vector):   " ++ show (gradeProj 1 uv)
  putStrLn $ "Grade 2 (Bivector): " ++ show (gradeProj 2 uv)
  putStrLn $ "Grade 3 (Volume):   " ++ show (gradeProj 3 uv)
  putStrLn ""

  -- 5. Involutions: Reversion and Poincaré Dual
  -- Reversion (~A) reverses the order of vector factors: ~(e1 * e2) = e2 * e1 = -e12
  let rev_uv = reverseMV uv
      dual_u = poincareDualMV u

  putStrLn "--- 5. Algebraic Involutions ---"
  putStrLn $ "Reversion ~(u * v)      = " ++ show rev_uv
  putStrLn $ "Poincaré Dual J(u)      = " ++ show dual_u ++ " (Vectors map to orthogonal bivectors)"
  putStrLn ""
  putStrLn "Tutorial 1 completed successfully!"
