{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}

-- |
-- Module      : Clifford
-- Description : Top-level re-export of the Clifford Algebra (Geometric Algebra) library
-- License     : BSD-3-Clause
--
-- 'Clifford' is a high-performance, general-purpose Geometric Algebra library
-- in pure Haskell.
--
-- = Quickstart Example
--
-- @
-- import Clifford
--
-- -- 3D Euclidean Geometric Algebra with Double scalars
-- type GA3 = DenseMV Cl300 Double
--
-- v1, v2 :: GA3
-- v1 = 2 * e1 + 3 * e2
-- v2 = 1 * e2 + 4 * e3
--
-- -- Geometric product: v1 * v2 = v1 . v2 + v1 ^ v2
-- prod :: GA3
-- prod = v1 <*> v2
--
-- -- Rotor rotation of vector v1 by angle theta in the e12 plane
-- rotor :: GA3
-- rotor = rotorExp (scaleMV (pi / 4) e12)
--
-- rotatedV1 :: GA3
-- rotatedV1 = sandwich rotor v1
-- @
module Clifford
  ( -- * Re-exported Modules
    module Clifford.Signature
  , module Clifford.Blade
  , module Clifford.Field
  , module Clifford.Class
  , module Clifford.Dense
  , module Clifford.Sparse
  , module Clifford.Versor
  , module Clifford.Array.Tensor
  , module Clifford.Array.Matrix
  , module Clifford.Array.Convolution
  , module Clifford.Array.Fourier
  , module Clifford.Analytic
    -- * Basis Vector Shorthands
  , e1
  , e2
  , e3
  , e4
  , e5
  , e0
    -- * Common Blade Shorthands
  , e12
  , e23
  , e31
  , e123
    -- * 3D Projective Geometric Algebra (PGA Cl(3,0,1)) Helpers
  , pointPGA
  , planePGA
  , linePGA
  , intersectPlanesPGA
  , intersectPlaneLinePGA
  , translatorPGA
  , rotatorPGA
  , motorPGA
    -- * Conformal Geometric Algebra (CGA Cl(4,1,0)) Helpers
  , eInfCGA
  , e0CGA
  , pointCGA
  , sphereCGA
  , sqDistanceCGA
  ) where

import Prelude hiding ((<*>))
import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class
import Clifford.Dense
import Clifford.Sparse
import Clifford.Versor
import Clifford.Array.Tensor
import Clifford.Array.Matrix
import Clifford.Array.Convolution
import Clifford.Array.Fourier
import Clifford.Analytic

--------------------------------------------------------------------------------
-- Basis Shorthands
--------------------------------------------------------------------------------

-- | Basis vector \(e_1\).
e1 :: (CliffordVectorSpace sig k mv) => mv
e1 = basisMV 1
{-# INLINE e1 #-}

-- | Basis vector \(e_2\).
e2 :: (CliffordVectorSpace sig k mv) => mv
e2 = basisMV 2
{-# INLINE e2 #-}

-- | Basis vector \(e_3\).
e3 :: (CliffordVectorSpace sig k mv) => mv
e3 = basisMV 3
{-# INLINE e3 #-}

-- | Basis vector \(e_4\).
e4 :: (CliffordVectorSpace sig k mv) => mv
e4 = basisMV 4
{-# INLINE e4 #-}

-- | Basis vector \(e_5\).
e5 :: (CliffordVectorSpace sig k mv) => mv
e5 = basisMV 5
{-# INLINE e5 #-}

-- | Degenerate / Origin basis vector \(e_0\) (mapped to \(e_4\) in PGA \(Cl(3,0,1)\)).
e0 :: (CliffordVectorSpace sig k mv) => mv
e0 = basisMV 4
{-# INLINE e0 #-}

-- | Bivector \(e_{12} = e_1 \wedge e_2\).
e12 :: (CliffordAlgebra sig k mv) => mv
e12 = basisMV 1 /\ basisMV 2
{-# INLINE e12 #-}

-- | Bivector \(e_{23} = e_2 \wedge e_3\).
e23 :: (CliffordAlgebra sig k mv) => mv
e23 = basisMV 2 /\ basisMV 3
{-# INLINE e23 #-}

-- | Bivector \(e_{31} = e_3 \wedge e_1\).
e31 :: (CliffordAlgebra sig k mv) => mv
e31 = basisMV 3 /\ basisMV 1
{-# INLINE e31 #-}

-- | Trivector / Pseudoscalar \(e_{123} = e_1 \wedge e_2 \wedge e_3\).
e123 :: (CliffordAlgebra sig k mv) => mv
e123 = basisMV 1 /\ basisMV 2 /\ basisMV 3
{-# INLINE e123 #-}

--------------------------------------------------------------------------------
-- 3D PGA Kinematics Helpers (Cl(3,0,1))
--------------------------------------------------------------------------------

-- | Construct a 3D point in PGA \(Cl(3,0,1)\):
-- \(P = e_{123} + x e_{023} + y e_{031} + z e_{012}\).
pointPGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => k -> k -> k -> mv
pointPGA x y z =
  let e023 = basisMV 4 /\ basisMV 2 /\ basisMV 3
      e031 = basisMV 4 /\ basisMV 3 /\ basisMV 1
      e012 = basisMV 4 /\ basisMV 1 /\ basisMV 2
  in e123 `addMV` (scaleMV x e023 `addMV` (scaleMV y e031 `addMV` scaleMV z e012))

-- | Construct a plane in PGA \(Cl(3,0,1)\):
-- \(\pi = a e_1 + b e_2 + c e_3 + d e_0\).
planePGA :: forall k mv. (CliffordVectorSpace Cl301 k mv) => k -> k -> k -> k -> mv
planePGA a b c d =
  scaleMV a (basisMV 1) `addMV` (scaleMV b (basisMV 2) `addMV` (scaleMV c (basisMV 3) `addMV` scaleMV d (basisMV 4)))

-- | Construct a 3D line in PGA \(Cl(3,0,1)\) connecting two points:
-- \(L = P_1 \vee P_2\).
linePGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
linePGA p1 p2 = p1 \/ p2

-- | Compute the intersection line of two planes in PGA \(Cl(3,0,1)\):
-- \(L = \pi_1 \wedge \pi_2\).
intersectPlanesPGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
intersectPlanesPGA pi1 pi2 = pi1 /\ pi2

-- | Compute the intersection point of a plane and a line in PGA \(Cl(3,0,1)\):
-- \(P = \pi \wedge L\).
intersectPlaneLinePGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
intersectPlaneLinePGA piPlane line = piPlane /\ line

-- | Construct a pure translation motor in PGA \(Cl(3,0,1)\):
-- \(T = 1 + \frac{1}{2} (dx e_{01} + dy e_{02} + dz e_{03})\).
translatorPGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => k -> k -> k -> mv
translatorPGA dx dy dz =
  let e01 = basisMV 4 /\ basisMV 1
      e02 = basisMV 4 /\ basisMV 2
      e03 = basisMV 4 /\ basisMV 3
      half = fromRational (1/2)
      tBivec = scaleMV (half * dx) e01 `addMV` (scaleMV (half * dy) e02 `addMV` scaleMV (half * dz) e03)
  in scalarMV 1 `addMV` tBivec

-- | Construct a pure rotation rotor in PGA \(Cl(3,0,1)\) around coordinate axes:
-- \(R = \cos(\theta/2) - \sin(\theta/2) B\).
rotatorPGA :: forall k mv. (Floating k, CliffordAlgebra Cl301 k mv) => k -> mv -> mv
rotatorPGA theta axisBivec =
  let halfTheta = theta / 2
      cosTerm = scalarMV (cos halfTheta)
      sinTerm = scaleMV (sin halfTheta) axisBivec
  in cosTerm `subMV` sinTerm

-- | Compose a translation motor and rotation rotor into a rigid body motor:
-- \(M = T \star R\).
motorPGA :: forall k mv. (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
motorPGA tMotor rRotor = tMotor <*> rRotor

--------------------------------------------------------------------------------
-- Conformal Geometric Algebra (CGA Cl(4,1,0)) Helpers
--------------------------------------------------------------------------------

-- | Conformal point at infinity \(e_\infty = e_4 + e_5\) (\(e_\infty^2 = 0\)).
eInfCGA :: (CliffordAlgebra Cl410 k mv) => mv
eInfCGA = basisMV 4 `addMV` basisMV 5

-- | Conformal point at origin \(e_0 = \frac{1}{2} (e_5 - e_4)\) (\(e_0^2 = 0, e_\infty \cdot e_0 = -1\)).
e0CGA :: (CliffordAlgebra Cl410 k mv) => mv
e0CGA = scaleMV (fromRational (1/2)) (basisMV 5 `subMV` basisMV 4)

-- | Construct a conformal point in CGA \(Cl(4,1,0)\):
-- \(P(\mathbf{x}) = \mathbf{x} + \frac{1}{2} \|\mathbf{x}\|^2 e_\infty + e_0\).
pointCGA :: forall k mv. (Floating k, CliffordAlgebra Cl410 k mv) => k -> k -> k -> mv
pointCGA x y z =
  let xVec = scaleMV x (basisMV 1) `addMV` (scaleMV y (basisMV 2) `addMV` scaleMV z (basisMV 3))
      sqDist = x*x + y*y + z*z
      infTerm = scaleMV (fromRational (1/2) * sqDist) eInfCGA
  in xVec `addMV` infTerm `addMV` e0CGA

-- | Construct a conformal sphere with given center point and radius in CGA \(Cl(4,1,0)\):
-- \(S(\mathbf{c}, r) = P(\mathbf{c}) - \frac{1}{2} r^2 e_\infty\).
sphereCGA :: forall k mv. (Floating k, CliffordAlgebra Cl410 k mv) => mv -> k -> mv
sphereCGA centerPt radius =
  let rSqTerm = scaleMV (fromRational (1/2) * radius * radius) eInfCGA
  in centerPt `subMV` rSqTerm

-- | Compute Euclidean squared distance between two conformal points:
-- \(d^2(P_1, P_2) = -2 (P_1 \cdot P_2)\).
sqDistanceCGA :: forall k mv. (CliffordAlgebra Cl410 k mv) => mv -> mv -> k
sqDistanceCGA p1 p2 = fromMetricScale (-2) * (p1 <|> p2)
