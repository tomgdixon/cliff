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
-- Module      : Clifford.Array.Tensor
-- Description : Contiguous flat unboxed array buffers for N-dimensional geometric multivector fields
-- License     : BSD-3-Clause
--
-- This module implements 'MVArray', a cache-aligned flat unboxed contiguous buffer
-- holding an \(N\)-dimensional grid of Clifford multivectors.
--
-- Layout: Each cell holds \(2^n\) contiguous scalar blade coefficients in canonical bitmask order.
module Clifford.Array.Tensor
  ( -- * Multivector Array Type
    MVArray(..)
    -- * Array Constructors
  , generateND
  , fromListND
  , zerosND
  , constND
    -- * Indexing & Slicing
  , indexND
  , coordsToFlat
  , flatToCoords
  , computeStrides
  , totalElements
  , sliceGrade
    -- * Element-wise Mappings & Reductions
  , mapMV
  , zipWithMV
  , dotProductArray
  ) where

import Control.DeepSeq (NFData(..))
import Control.Monad (forM_)
import Control.Monad.ST (runST)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray (IArray, (!), array, listArray)
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (STUArray, runSTUArray)
import Data.Array.Unboxed (UArray)
import Data.Bits (shiftL)
import Data.Proxy (Proxy(..))
import GHC.TypeLits

import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class
import Clifford.Dense

-- | Contiguous flat buffer storing an ND grid of multivectors.
data MVArray sig k = MVArray
  { shapeND      :: ![Int]         -- ^ Dimensions @[d_0, d_1, ..., d_{m-1}]@
  , mvaStrides   :: ![Int]         -- ^ Coordinate strides
  , mvaBuffer    :: !(UArray Int k)-- ^ Contiguous flat unboxed buffer
  }

instance NFData (MVArray sig k) where
  rnf (MVArray s st buf) = s `seq` st `seq` buf `seq` ()

instance (KnownSignature p q r, UnboxedField k, Eq k) => Eq (MVArray (Signature p q r) k) where
  (MVArray s1 _ b1) == (MVArray s2 _ b2) =
    s1 == s2 && b1 == b2

instance (KnownSignature p q r, UnboxedField k, Show k) => Show (MVArray (Signature p q r) k) where
  show (MVArray s _ _) = "MVArray { shape = " ++ show s ++ " }"

-- | Compute standard row-major strides for a given shape.
computeStrides :: [Int] -> [Int]
computeStrides [] = []
computeStrides dims =
  let initStrides = reverse (scanl (*) 1 (reverse (tail dims)))
  in initStrides
{-# INLINE computeStrides #-}

-- | Compute total number of multivector cells.
totalElements :: MVArray sig k -> Int
totalElements (MVArray s _ _) = product s
{-# INLINE totalElements #-}

-- | Map ND coordinate list to a 1D linear cell index.
coordsToFlat :: [Int] -> [Int] -> Int
coordsToFlat strides coords = sum (zipWith (*) strides coords)
{-# INLINE coordsToFlat #-}

-- | Map a 1D linear cell index to ND coordinates.
flatToCoords :: [Int] -> Int -> [Int]
flatToCoords strides idx = go strides idx
  where
    go [] _ = []
    go (s:ss) curr =
      let (q, r) = curr `quotRem` s
      in q : go ss r
{-# INLINE flatToCoords #-}

-- | Generate an ND multivector array by applying a generator function to each coordinate.
generateND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => [Int] -> ([Int] -> DenseMV (Signature p q r) k) -> MVArray (Signature p q r) k
generateND dims genFn =
  let strides = computeStrides dims
      totCells = product dims
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      totScalars = totCells * bladesPerCell

      buf = runSTUArray $ do
        marr <- newArray (0, totScalars - 1) (fromMetricScale 0)
        forM_ [0 .. totCells - 1] $ \cellIdx -> do
          let coords = flatToCoords strides cellIdx
              DenseMV cellArr = genFn coords
              cellOffset = cellIdx * bladesPerCell
          forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
            writeArray marr (cellOffset + bIdx) (unsafeAt cellArr bIdx)
        return marr
  in MVArray dims strides buf
{-# INLINE generateND #-}

-- | Construct an ND multivector array from a flat list of multivectors.
fromListND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => [Int] -> [DenseMV (Signature p q r) k] -> MVArray (Signature p q r) k
fromListND dims mvs =
  let strides = computeStrides dims
      totCells = product dims
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      totScalars = totCells * bladesPerCell

      buf = runSTUArray $ do
        marr <- newArray (0, totScalars - 1) (fromMetricScale 0)
        forM_ (zip [0 .. totCells - 1] mvs) $ \(cellIdx, DenseMV cellArr) -> do
          let cellOffset = cellIdx * bladesPerCell
          forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
            writeArray marr (cellOffset + bIdx) (unsafeAt cellArr bIdx)
        return marr
  in MVArray dims strides buf
{-# INLINE fromListND #-}

-- | Construct an ND array initialized to zero multivectors.
zerosND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => [Int] -> MVArray (Signature p q r) k
zerosND dims =
  let strides = computeStrides dims
      totCells = product dims
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      totScalars = totCells * bladesPerCell
      buf = array (0, totScalars - 1) [(i, fromMetricScale 0) | i <- [0 .. totScalars - 1]]
  in MVArray dims strides buf
{-# INLINE zerosND #-}

-- | Construct an ND array with a constant multivector in every cell.
constND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => [Int] -> DenseMV (Signature p q r) k -> MVArray (Signature p q r) k
constND dims val = generateND dims (const val)
{-# INLINE constND #-}

-- | Index a multivector at given ND coordinates.
indexND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => MVArray (Signature p q r) k -> [Int] -> DenseMV (Signature p q r) k
indexND (MVArray dims strides buf) coords =
  let cellIdx = coordsToFlat strides coords
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      cellOffset = cellIdx * bladesPerCell
      coeffs = [unsafeAt buf (cellOffset + i) | i <- [0 .. bladesPerCell - 1]]
  in DenseMV (listArray (0, bladesPerCell - 1) coeffs)
{-# INLINE indexND #-}

-- | Extract a scalar sub-tensor for a specific blade grade @g@.
sliceGrade :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => Int -> MVArray (Signature p q r) k -> MVArray (Signature p q r) k
sliceGrade g arr = mapMV (gradeProj g) arr
{-# INLINE sliceGrade #-}

-- | Pointwise mapping of a multivector function across all cells.
mapMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
      => (DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k)
      -> MVArray (Signature p q r) k
      -> MVArray (Signature p q r) k
mapMV f input =
  let dims = shapeND input
      strides = mvaStrides input
      totCells = totalElements input
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      buf = mvaBuffer input

      outBuf = runSTUArray $ do
        marr <- newArray (0, totCells * bladesPerCell - 1) (fromMetricScale 0)
        forM_ [0 .. totCells - 1] $ \cellIdx -> do
          let cellOffset = cellIdx * bladesPerCell
              inMV = DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt buf (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])
              DenseMV outCellArr = f inMV
          forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
            writeArray marr (cellOffset + bIdx) (unsafeAt outCellArr bIdx)
        return marr
  in MVArray dims strides outBuf
{-# INLINE mapMV #-}

-- | Pointwise binary operation across two conforming multivector arrays.
zipWithMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
          => (DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k)
          -> MVArray (Signature p q r) k
          -> MVArray (Signature p q r) k
          -> MVArray (Signature p q r) k
zipWithMV f arrA arrB
  | shapeND arrA /= shapeND arrB = error "Clifford.Array.Tensor.zipWithMV: shape mismatch"
  | otherwise =
      let dims = shapeND arrA
          strides = mvaStrides arrA
          totCells = totalElements arrA
          bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          bufA = mvaBuffer arrA
          bufB = mvaBuffer arrB

          outBuf = runSTUArray $ do
            marr <- newArray (0, totCells * bladesPerCell - 1) (fromMetricScale 0)
            forM_ [0 .. totCells - 1] $ \cellIdx -> do
              let cellOffset = cellIdx * bladesPerCell
                  inA = DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufA (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])
                  inB = DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufB (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])
                  DenseMV outCellArr = f inA inB
              forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
                writeArray marr (cellOffset + bIdx) (unsafeAt outCellArr bIdx)
            return marr
      in MVArray dims strides outBuf
{-# INLINE zipWithMV #-}

-- | Scalar dot product across two multivector array fields: \(\sum_{\mathbf{x}} \langle A(\mathbf{x}) B(\mathbf{x}) \rangle_0\).
dotProductArray :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                => MVArray (Signature p q r) k -> MVArray (Signature p q r) k -> k
dotProductArray arrA arrB
  | shapeND arrA /= shapeND arrB = error "Clifford.Array.Tensor.dotProductArray: shape mismatch"
  | otherwise =
      let totCells = totalElements arrA
          strides = mvaStrides arrA
          bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          bufA = mvaBuffer arrA
          bufB = mvaBuffer arrB

          terms = [ let cellOffset = cellIdx * bladesPerCell
                        inA = DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufA (cellOffset + i) | i <- [0 .. bladesPerCell - 1]]) :: DenseMV (Signature p q r) k
                        inB = DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufB (cellOffset + i) | i <- [0 .. bladesPerCell - 1]]) :: DenseMV (Signature p q r) k
                    in inA <|> inB
                  | cellIdx <- [0 .. totCells - 1]
                  ]
      in foldr (+) (fromMetricScale 0) terms
{-# INLINE dotProductArray #-}
