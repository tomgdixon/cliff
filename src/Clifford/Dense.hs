{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Clifford.Dense
-- Description : Statically unrolled dense multivectors with unboxed array storage for low dimensions
-- License     : BSD-3-Clause
--
-- This module implements 'DenseMV', a dense multivector holding exactly \(2^n\) unboxed scalar
-- coefficients. It is optimized for dimensions \(n \le 6\) (up to 64 blades), providing zero-allocation
-- unrolled geometric products and cache-friendly operations.
module Clifford.Dense
  ( -- * Dense Multivector Type
    DenseMV(..)
    -- * Constructors & Extraction
  , mkDense
  , denseCoeffs
  , denseIndex
  , denseBladeAt
    -- * Specialized Fast Products
  , denseMul
  , denseWedge
  , denseDotLeft
  , denseDotRight
  , denseScalarProd
  ) where

import Prelude hiding ((<*>))
import Control.DeepSeq (NFData(..))
import Control.Monad (forM_, unless)
import Control.Monad.ST (runST)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray (IArray, (!), array, listArray, assocs)
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (STUArray, runSTUArray)
import Data.Array.Unboxed (UArray)
import Data.Bits (shiftL, testBit)
import Data.List (intercalate)
import Data.Proxy (Proxy(..))
import GHC.TypeLits

import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class

-- | Statically dimensioned dense multivector holding exactly \(2^n\) unboxed scalar coefficients.
data DenseMV sig k = DenseMV !(UArray Int k)

instance NFData (DenseMV sig k) where
  rnf (DenseMV arr) = arr `seq` ()

instance (KnownSignature p q r, UnboxedField k, Eq k) => Eq (DenseMV (Signature p q r) k) where
  (DenseMV arr1) == (DenseMV arr2) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
    in all (\i -> unsafeAt arr1 i == unsafeAt arr2 i) [0 .. tot - 1]

instance (KnownSignature p q r, UnboxedField k, Show k) => Show (DenseMV (Signature p q r) k) where
  show mv =
    let terms = toBladeList mv
        showTerm (b, v)
          | b == scalarBlade = show v
          | otherwise        = show v ++ "*" ++ show b
    in if null terms
         then "0"
         else intercalate " + " (map showTerm terms)

-- | Construct a dense multivector from a list of blade-coefficient pairs.
mkDense :: forall p q r k. (KnownSignature p q r, UnboxedField k) => [(Blade, k)] -> DenseMV (Signature p q r) k
mkDense blades =
  let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      updated = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ blades $ \(b, v) -> do
          let idx = denseIndex b
          if idx < tot
            then writeArray marr idx v
            else return ()
        return marr
  in DenseMV updated
{-# INLINE mkDense #-}

-- | Extract the underlying flat unboxed array of coefficients.
denseCoeffs :: DenseMV sig k -> UArray Int k
denseCoeffs (DenseMV arr) = arr
{-# INLINE denseCoeffs #-}

-- | Convert a basis blade to its slot index in the dense buffer.
denseIndex :: Blade -> Int
denseIndex (Blade w) = fromIntegral w
{-# INLINE denseIndex #-}

-- | Convert a slot index in the dense buffer to its basis blade.
denseBladeAt :: Int -> Blade
denseBladeAt i = Blade (fromIntegral i)
{-# INLINE denseBladeAt #-}

--------------------------------------------------------------------------------
-- Typeclass Instances
--------------------------------------------------------------------------------

instance (KnownSignature p q r, UnboxedField k) => CliffordVectorSpace (Signature p q r) k (DenseMV (Signature p q r) k) where
  zeroMV =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
    in DenseMV (array (0, tot - 1) [(i, fromMetricScale 0) | i <- [0 .. tot - 1]])
  {-# INLINE zeroMV #-}

  scalarMV v =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
    in DenseMV (array (0, tot - 1) ((0, v) : [(i, fromMetricScale 0) | i <- [1 .. tot - 1]]))
  {-# INLINE scalarMV #-}

  bladeMV b v =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        idx = denseIndex b
    in DenseMV (array (0, tot - 1) [(i, if i == idx then v else fromMetricScale 0) | i <- [0 .. tot - 1]])
  {-# INLINE bladeMV #-}

  gradeProj k (DenseMV arr) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        updated = runSTUArray $ do
          marr <- newArray (0, tot - 1) (fromMetricScale 0)
          forM_ [0 .. tot - 1] $ \i -> do
            let b = denseBladeAt i
            if bladeGrade b == k
              then writeArray marr i (unsafeAt arr i)
              else writeArray marr i (fromMetricScale 0)
          return marr
    in DenseMV updated
  {-# INLINE gradeProj #-}

  grades (DenseMV arr) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        active = [bladeGrade (denseBladeAt i) | i <- [0 .. tot - 1], not (isZero (unsafeAt arr i))]
    in foldl' (\acc g -> if g `elem` acc then acc else acc ++ [g]) [] active
    where
      foldl' f z []     = z
      foldl' f z (x:xs) = foldl' f (f z x) xs
  {-# INLINE grades #-}

  getBlade b (DenseMV arr) =
    let idx = denseIndex b
        tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
    in if idx < tot then unsafeAt arr idx else fromMetricScale 0
  {-# INLINE getBlade #-}

  fromBladeList = mkDense
  {-# INLINE fromBladeList #-}

  toBladeList (DenseMV arr) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
    in [(denseBladeAt i, unsafeAt arr i) | i <- [0 .. tot - 1], not (isZero (unsafeAt arr i))]
  {-# INLINE toBladeList #-}

  addMV (DenseMV a) (DenseMV b) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        res = listArray (0, tot - 1) [unsafeAt a i + unsafeAt b i | i <- [0 .. tot - 1]]
    in DenseMV res
  {-# INLINE addMV #-}

  subMV (DenseMV a) (DenseMV b) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        res = listArray (0, tot - 1) [unsafeAt a i - unsafeAt b i | i <- [0 .. tot - 1]]
    in DenseMV res
  {-# INLINE subMV #-}

  scaleMV s (DenseMV a) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        res = listArray (0, tot - 1) [s * unsafeAt a i | i <- [0 .. tot - 1]]
    in DenseMV res
  {-# INLINE scaleMV #-}

instance (KnownSignature p q r, UnboxedField k) => CliffordAlgebra (Signature p q r) k (DenseMV (Signature p q r) k) where
  (<*>) = denseMul
  {-# INLINE (<*>) #-}

  (/\) = denseWedge
  {-# INLINE (/\) #-}

  (<.>) = denseDotLeft
  {-# INLINE (<.>) #-}

  (>.>) = denseDotRight
  {-# INLINE (>.>) #-}

  (<|>) = denseScalarProd
  {-# INLINE (<|>) #-}

  reverseMV (DenseMV arr) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        revSign g = if testBit (g * (g - 1) `div` 2) 0 then -1 else 1
        res = listArray (0, tot - 1)
          [ let g = bladeGrade (denseBladeAt i)
            in fromMetricScale (revSign g) * unsafeAt arr i
          | i <- [0 .. tot - 1]
          ]
    in DenseMV res
  {-# INLINE reverseMV #-}

  involMV (DenseMV arr) =
    let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
        res = listArray (0, tot - 1)
          [ let g = bladeGrade (denseBladeAt i)
            in fromMetricScale (if testBit g 0 then -1 else 1) * unsafeAt arr i
          | i <- [0 .. tot - 1]
          ]
    in DenseMV res
  {-# INLINE involMV #-}

  poincareDualMV (DenseMV arr) =
    let dim = signatureDim (Proxy :: Proxy (Signature p q r))
        tot = 1 `shiftL` dim
        res = runSTUArray $ do
          marr <- newArray (0, tot - 1) (fromMetricScale 0)
          forM_ [0 .. tot - 1] $ \i -> do
            let val = unsafeAt arr i
            unless (isZero val) $ do
              let (sgn, dualBlade) = poincareDual dim (denseBladeAt i)
                  dualIdx = denseIndex dualBlade
              writeArray marr dualIdx (fromMetricScale sgn * val)
          return marr
    in DenseMV res
  {-# INLINE poincareDualMV #-}

instance (KnownSignature p q r, UnboxedField k) => Num (DenseMV (Signature p q r) k) where
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
-- Explicit Product Implementations
--------------------------------------------------------------------------------

-- | High-performance geometric product for dense multivectors via in-place ST array.
denseMul :: forall p q r k. (KnownSignature p q r, UnboxedField k)
         => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
denseMul (DenseMV arrA) (DenseMV arrB) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      tot = 1 `shiftL` signatureDim sigProxy
      res = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \i -> do
          let ai = unsafeAt arrA i
          unless (isZero ai) $ do
            let bi = denseBladeAt i
            forM_ [0 .. tot - 1] $ \j -> do
              let bj = unsafeAt arrB j
              unless (isZero bj) $ do
                let bjBlade = denseBladeAt j
                    (s, rBlade) = bladeMul sigProxy bi bjBlade
                unless (s == 0) $ do
                  let rIdx = denseIndex rBlade
                      term = ai * bj * fromMetricScale s
                  curr <- readArray marr rIdx
                  writeArray marr rIdx (curr + term)
        return marr
  in DenseMV res
{-# INLINE denseMul #-}

-- | Exterior / Wedge product for dense multivectors.
denseWedge :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
denseWedge (DenseMV arrA) (DenseMV arrB) =
  let tot = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      res = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \i -> do
          let ai = unsafeAt arrA i
          unless (isZero ai) $ do
            let bi = denseBladeAt i
            forM_ [0 .. tot - 1] $ \j -> do
              let bj = unsafeAt arrB j
              unless (isZero bj) $ do
                let bjBlade = denseBladeAt j
                case bladeWedge bi bjBlade of
                  Nothing -> return ()
                  Just (s, rBlade) -> do
                    let rIdx = denseIndex rBlade
                        term = ai * bj * fromMetricScale s
                    curr <- readArray marr rIdx
                    writeArray marr rIdx (curr + term)
        return marr
  in DenseMV res
{-# INLINE denseWedge #-}

-- | Left contraction for dense multivectors.
denseDotLeft :: forall p q r k. (KnownSignature p q r, UnboxedField k)
             => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
denseDotLeft (DenseMV arrA) (DenseMV arrB) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      tot = 1 `shiftL` signatureDim sigProxy
      res = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \i -> do
          let ai = unsafeAt arrA i
          unless (isZero ai) $ do
            let bi = denseBladeAt i
            forM_ [0 .. tot - 1] $ \j -> do
              let bj = unsafeAt arrB j
              unless (isZero bj) $ do
                let bjBlade = denseBladeAt j
                case bladeDotLeft sigProxy bi bjBlade of
                  Nothing -> return ()
                  Just (s, rBlade) -> do
                    let rIdx = denseIndex rBlade
                        term = ai * bj * fromMetricScale s
                    curr <- readArray marr rIdx
                    writeArray marr rIdx (curr + term)
        return marr
  in DenseMV res
{-# INLINE denseDotLeft #-}

-- | Right contraction for dense multivectors.
denseDotRight :: forall p q r k. (KnownSignature p q r, UnboxedField k)
              => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
denseDotRight (DenseMV arrA) (DenseMV arrB) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      tot = 1 `shiftL` signatureDim sigProxy
      res = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \i -> do
          let ai = unsafeAt arrA i
          unless (isZero ai) $ do
            let bi = denseBladeAt i
            forM_ [0 .. tot - 1] $ \j -> do
              let bj = unsafeAt arrB j
              unless (isZero bj) $ do
                let bjBlade = denseBladeAt j
                case bladeDotRight sigProxy bi bjBlade of
                  Nothing -> return ()
                  Just (s, rBlade) -> do
                    let rIdx = denseIndex rBlade
                        term = ai * bj * fromMetricScale s
                    curr <- readArray marr rIdx
                    writeArray marr rIdx (curr + term)
        return marr
  in DenseMV res
{-# INLINE denseDotRight #-}

-- | Scalar inner product for dense multivectors.
denseScalarProd :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> k
denseScalarProd (DenseMV arrA) (DenseMV arrB) =
  let sigProxy = Proxy :: Proxy (Signature p q r)
      tot = 1 `shiftL` signatureDim sigProxy
      terms = [ let ai = unsafeAt arrA i
                    bi = unsafeAt arrB i
                in if isZero ai || isZero bi
                     then fromMetricScale 0
                     else ai * bi * fromMetricScale (bladeScalarProd sigProxy (denseBladeAt i) (denseBladeAt i))
              | i <- [0 .. tot - 1]
              ]
  in foldr (+) (fromMetricScale 0) terms
{-# INLINE denseScalarProd #-}
