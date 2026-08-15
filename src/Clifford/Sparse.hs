{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Clifford.Sparse
-- Description : Map-backed sparse multivectors for high-dimensional Clifford algebras (n in 7..64)
-- License     : BSD-3-Clause
--
-- This module implements 'SparseMV', storing only non-zero basis blade coefficients in a balanced
-- 'Map Blade k'. It is designed for high-dimensional spaces where \(2^n\) is too large to represent
-- densely (e.g. \(n=16, 32, 64\)).
module Clifford.Sparse
  ( -- * Sparse Multivector Type
    SparseMV(..)
    -- * Conversions with DenseMV
  , toSparse
  , toDense
    -- * Specialized Sparse Operations
  , sparseMul
  , sparseWedge
  , sparseDotLeft
  , sparseDotRight
  , sparseScalarProd
  ) where

import Prelude hiding ((<*>))
import Control.DeepSeq (NFData(..))
import Data.Bits (testBit)
import Data.List (intercalate)
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Proxy (Proxy(..))

import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class
import Clifford.Dense

-- | Sparse multivector parameterized by signature @sig@ and scalar field @k@.
newtype SparseMV sig k = SparseMV { unSparseMV :: Map Blade k }

instance NFData k => NFData (SparseMV sig k) where
  rnf (SparseMV m) = rnf m

instance (Field k, Eq k) => Eq (SparseMV (Signature p q r) k) where
  (SparseMV m1) == (SparseMV m2) =
    let clean1 = Map.filter (not . isZero) m1
        clean2 = Map.filter (not . isZero) m2
    in clean1 == clean2

instance (Field k, Show k) => Show (SparseMV (Signature p q r) k) where
  show (SparseMV m) =
    let clean = Map.filter (not . isZero) m
        showTerm (b, v)
          | b == scalarBlade = show v
          | otherwise        = show v ++ "*" ++ show b
    in if Map.null clean
         then "0"
         else intercalate " + " (map showTerm (Map.toList clean))

--------------------------------------------------------------------------------
-- Typeclass Instances
--------------------------------------------------------------------------------

instance Field k => CliffordVectorSpace (Signature p q r) k (SparseMV (Signature p q r) k) where
  zeroMV = SparseMV Map.empty
  {-# INLINE zeroMV #-}

  scalarMV v
    | isZero v  = SparseMV Map.empty
    | otherwise = SparseMV (Map.singleton scalarBlade v)
  {-# INLINE scalarMV #-}

  bladeMV b v
    | isZero v  = SparseMV Map.empty
    | otherwise = SparseMV (Map.singleton b v)
  {-# INLINE bladeMV #-}

  gradeProj k (SparseMV m) =
    SparseMV (Map.filterWithKey (\b _ -> bladeGrade b == k) m)
  {-# INLINE gradeProj #-}

  grades (SparseMV m) =
    let clean = Map.filter (not . isZero) m
        active = [bladeGrade b | b <- Map.keys clean]
    in foldl (\acc g -> if g `elem` acc then acc else acc ++ [g]) [] active
  {-# INLINE grades #-}

  getBlade b (SparseMV m) = Map.findWithDefault (fromMetricScale 0) b m
  {-# INLINE getBlade #-}

  fromBladeList pairs =
    let clean = filter (not . isZero . snd) pairs
        m = Map.fromListWith (+) clean
    in SparseMV (Map.filter (not . isZero) m)
  {-# INLINE fromBladeList #-}

  toBladeList (SparseMV m) =
    Map.toList (Map.filter (not . isZero) m)
  {-# INLINE toBladeList #-}

  addMV (SparseMV m1) (SparseMV m2) =
    SparseMV (Map.filter (not . isZero) (Map.unionWith (+) m1 m2))
  {-# INLINE addMV #-}

  subMV (SparseMV m1) (SparseMV m2) =
    SparseMV (Map.filter (not . isZero) (Map.unionWith (+) m1 (Map.map negate m2)))
  {-# INLINE subMV #-}

  scaleMV s (SparseMV m)
    | isZero s  = SparseMV Map.empty
    | otherwise = SparseMV (Map.filter (not . isZero) (Map.map (s *) m))
  {-# INLINE scaleMV #-}

instance (KnownSignature p q r, Field k) => CliffordAlgebra (Signature p q r) k (SparseMV (Signature p q r) k) where
  (<*>) = sparseMul
  {-# INLINE (<*>) #-}

  (/\) = sparseWedge
  {-# INLINE (/\) #-}

  (<.>) = sparseDotLeft
  {-# INLINE (<.>) #-}

  (>.>) = sparseDotRight
  {-# INLINE (>.>) #-}

  (<|>) = sparseScalarProd
  {-# INLINE (<|>) #-}

  reverseMV (SparseMV m) =
    let revSign g = if testBit (g * (g - 1) `div` 2) 0 then -1 else 1
        m' = Map.mapWithKey (\b v -> fromMetricScale (revSign (bladeGrade b)) * v) m
    in SparseMV (Map.filter (not . isZero) m')
  {-# INLINE reverseMV #-}

  involMV (SparseMV m) =
    let m' = Map.mapWithKey (\b v -> fromMetricScale (if testBit (bladeGrade b) 0 then -1 else 1) * v) m
    in SparseMV (Map.filter (not . isZero) m')
  {-# INLINE involMV #-}

  poincareDualMV (SparseMV m) =
    let dim = signatureDim (Proxy :: Proxy (Signature p q r))
        pairs = [ let (sgn, dualBlade) = poincareDual dim b
                  in (dualBlade, fromMetricScale sgn * v)
                | (b, v) <- Map.toList m
                ]
    in SparseMV (Map.filter (not . isZero) (Map.fromListWith (+) pairs))
  {-# INLINE poincareDualMV #-}

instance (KnownSignature p q r, Field k) => Num (SparseMV (Signature p q r) k) where
  (+) = addMV
  (-) = subMV
  (*) = (<*>)
  negate = negateMV
  abs x = x
  signum x = x
  fromInteger n = scalarMV (fromInteger n)
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}
  {-# INLINE negate #-}
  {-# INLINE fromInteger #-}

--------------------------------------------------------------------------------
-- Inter-conversion with DenseMV
--------------------------------------------------------------------------------

-- | Convert a dense multivector to sparse representation.
toSparse :: (CliffordVectorSpace sig k (DenseMV sig k), CliffordVectorSpace sig k (SparseMV sig k))
         => DenseMV sig k -> SparseMV sig k
toSparse dense = fromBladeList (toBladeList dense)
{-# INLINE toSparse #-}

-- | Convert a sparse multivector to dense representation.
toDense :: (CliffordVectorSpace sig k (DenseMV sig k), CliffordVectorSpace sig k (SparseMV sig k))
        => SparseMV sig k -> DenseMV sig k
toDense sparse = fromBladeList (toBladeList sparse)
{-# INLINE toDense #-}

--------------------------------------------------------------------------------
-- Explicit Sparse Products
--------------------------------------------------------------------------------

-- | Sparse geometric product: \(A \star B\).
sparseMul :: forall p q r k. (KnownSignature p q r, Field k)
          => SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k
sparseMul (SparseMV m1) (SparseMV m2) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      terms = [ (rBlade, a * b * fromMetricScale s)
              | (b1, a) <- Map.toList m1
              , (b2, b) <- Map.toList m2
              , let (s, rBlade) = bladeMul sigProxy b1 b2
              , s /= 0
              ]
      m = Map.fromListWith (+) terms
  in SparseMV (Map.filter (not . isZero) m)
{-# INLINE sparseMul #-}

-- | Sparse exterior / wedge product: \(A \wedge B\).
sparseWedge :: forall p q r k. Field k
            => SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k
sparseWedge (SparseMV m1) (SparseMV m2) =
  let terms = [ (rBlade, a * b * fromMetricScale s)
              | (b1, a) <- Map.toList m1
              , (b2, b) <- Map.toList m2
              , Just (s, rBlade) <- [bladeWedge b1 b2]
              ]
      m = Map.fromListWith (+) terms
  in SparseMV (Map.filter (not . isZero) m)
{-# INLINE sparseWedge #-}

-- | Sparse left contraction.
sparseDotLeft :: forall p q r k. (KnownSignature p q r, Field k)
              => SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k
sparseDotLeft (SparseMV m1) (SparseMV m2) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      terms = [ (rBlade, a * b * fromMetricScale s)
              | (b1, a) <- Map.toList m1
              , (b2, b) <- Map.toList m2
              , Just (s, rBlade) <- [bladeDotLeft sigProxy b1 b2]
              ]
      m = Map.fromListWith (+) terms
  in SparseMV (Map.filter (not . isZero) m)
{-# INLINE sparseDotLeft #-}

-- | Sparse right contraction.
sparseDotRight :: forall p q r k. (KnownSignature p q r, Field k)
               => SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k
sparseDotRight (SparseMV m1) (SparseMV m2) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      terms = [ (rBlade, a * b * fromMetricScale s)
              | (b1, a) <- Map.toList m1
              , (b2, b) <- Map.toList m2
              , Just (s, rBlade) <- [bladeDotRight sigProxy b1 b2]
              ]
      m = Map.fromListWith (+) terms
  in SparseMV (Map.filter (not . isZero) m)
{-# INLINE sparseDotRight #-}

-- | Sparse scalar inner product.
sparseScalarProd :: forall p q r k. (KnownSignature p q r, Field k)
                 => SparseMV (Signature p q r) k -> SparseMV (Signature p q r) k -> k
sparseScalarProd (SparseMV m1) (SparseMV m2) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      shared = Map.intersectionWith (*) m1 m2
      terms = [ v * fromMetricScale (bladeScalarProd sigProxy b b)
              | (b, v) <- Map.toList shared
              ]
  in foldr (+) (fromMetricScale 0) terms
{-# INLINE sparseScalarProd #-}
