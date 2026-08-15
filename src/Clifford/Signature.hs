{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Clifford.Signature
-- Description : Type-level representations of Clifford algebra metric signatures Cl(p,q,r)
-- License     : BSD-3-Clause
--
-- This module defines the type-level metric signature @Signature p q r@ where:
--
-- * @p@ is the number of positive basis vectors (\(e_i^2 = +1\))
-- * @q@ is the number of negative basis vectors (\(e_i^2 = -1\))
-- * @r@ is the number of degenerate / null basis vectors (\(e_i^2 = 0\))
module Clifford.Signature
  ( -- * Type-Level Signature
    Signature(..)
  , KnownSignature(..)
    -- * Type Families
  , Dim
  , TotalBlades
  , IsEuclidean
  , IsDegenerate
    -- * Common Signature Aliases
  , Cl200
  , Cl300
  , Cl130
  , Cl310
  , Cl301
  , Cl410
  , Cl8_00
  , Cl16_00
  , Cl24_00
  , Cl32_00
  , Cl48_00
  , Cl64_00
    -- * Reflection Functions
  , signatureVal
  , metricSign
  , metricSignOf
  ) where

import Data.Kind (Type)
import GHC.TypeLits
import Data.Proxy

-- | Metric signature \(Cl(p,q,r)\) parameterized at the type level.
data Signature (p :: Nat) (q :: Nat) (r :: Nat) = Signature

-- | Total dimension \(n = p + q + r\).
type family Dim (sig :: Type) :: Nat where
  Dim (Signature p q r) = p + q + r

-- | Total number of basis blades \(2^n\).
type family TotalBlades (sig :: Type) :: Nat where
  TotalBlades sig = 2 ^ Dim sig

-- | Check if signature is purely Euclidean (\(q = 0, r = 0\)).
type family IsEuclidean (sig :: Type) :: Bool where
  IsEuclidean (Signature _ 0 0) = 'True
  IsEuclidean _                 = 'False

-- | Check if signature has degenerate/null dimensions (\(r > 0\)).
type family IsDegenerate (sig :: Type) :: Bool where
  IsDegenerate (Signature _ _ 0) = 'False
  IsDegenerate _                 = 'True

-- | 2D Euclidean (Isomorphic to Complex numbers in even subalgebra)
type Cl200 = Signature 2 0 0

-- | 3D Euclidean Geometric Algebra (Quaternions as even subalgebra)
type Cl300 = Signature 3 0 0

-- | Spacetime Algebra (STA, Minkowski metric +---)
type Cl130 = Signature 1 3 0

-- | Spacetime Algebra (Minkowski metric -+++)
type Cl310 = Signature 3 1 0

-- | 3D Projective Geometric Algebra (PGA, 3D Euclidean kinematics)
type Cl301 = Signature 3 0 1

-- | Conformal Geometric Algebra (CGA)
type Cl410 = Signature 4 1 0

-- | 8D Euclidean Clifford algebra
type Cl8_00  = Signature 8 0 0

-- | 16D Euclidean Clifford algebra for high-dimensional and coding experiments
type Cl16_00 = Signature 16 0 0

-- | 24D Euclidean Clifford algebra (Golay / Leech lattice space)
type Cl24_00 = Signature 24 0 0

-- | 32D Euclidean Clifford algebra
type Cl32_00 = Signature 32 0 0

-- | 48D Euclidean Clifford algebra
type Cl48_00 = Signature 48 0 0

-- | 64D Euclidean Clifford algebra for 64-bit sparse bitmask limits
type Cl64_00 = Signature 64 0 0

-- | Class for reflecting type-level signature values to runtime integers.
class (KnownNat p, KnownNat q, KnownNat r) => KnownSignature (p :: Nat) (q :: Nat) (r :: Nat) where
  signatureP   :: proxy (Signature p q r) -> Int
  signatureQ   :: proxy (Signature p q r) -> Int
  signatureR   :: proxy (Signature p q r) -> Int
  signatureDim :: proxy (Signature p q r) -> Int

instance (KnownNat p, KnownNat q, KnownNat r) => KnownSignature p q r where
  signatureP _ = fromIntegral (natVal (Proxy :: Proxy p))
  signatureQ _ = fromIntegral (natVal (Proxy :: Proxy q))
  signatureR _ = fromIntegral (natVal (Proxy :: Proxy r))
  signatureDim _ = fromIntegral (natVal (Proxy :: Proxy p) + natVal (Proxy :: Proxy q) + natVal (Proxy :: Proxy r))
  {-# INLINE signatureP #-}
  {-# INLINE signatureQ #-}
  {-# INLINE signatureR #-}
  {-# INLINE signatureDim #-}

-- | Extract runtime @(p, q, r)@ tuple from a proxy.
signatureVal :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> (Int, Int, Int)
signatureVal p = (signatureP p, signatureQ p, signatureR p)
{-# INLINE signatureVal #-}

-- | Evaluate metric signature quadratic form for basis vector \(e_i\) (1-based index):
--
-- * Returns @+1@ for \(1 \le i \le p\)
-- * Returns @-1@ for \(p < i \le p+q\)
-- * Returns @0@ for \(p+q < i \le p+q+r\)
metricSign :: forall p q r. KnownSignature p q r => Int -> Int
metricSign i
  | i <= 0     = error $ "Clifford.Signature.metricSign: index must be positive, got " ++ show i
  | i <= p     = 1
  | i <= p+q   = -1
  | i <= p+q+r = 0
  | otherwise  = error $ "Clifford.Signature.metricSign: index " ++ show i ++ " exceeds dimension " ++ show (p+q+r)
  where
    (p, q, r) = signatureVal (Proxy :: Proxy (Signature p q r))
{-# INLINE metricSign #-}

-- | Evaluate metric quadratic form given an explicit proxy.
metricSignOf :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Int -> Int
metricSignOf p i
  | i <= 0     = error $ "Clifford.Signature.metricSignOf: index must be positive, got " ++ show i
  | i <= pDim  = 1
  | i <= pDim+qDim = -1
  | i <= pDim+qDim+rDim = 0
  | otherwise  = error $ "Clifford.Signature.metricSignOf: index " ++ show i ++ " exceeds dimension " ++ show (pDim+qDim+rDim)
  where
    (pDim, qDim, rDim) = signatureVal p
{-# INLINE metricSignOf #-}
