{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Array.Convolution
-- Description : Sliding geometric stencils, boundary padding modes, and algebraic convolution
-- License     : BSD-3-Clause
--
-- This module implements sliding stencils for 1D, 2D, and N-dimensional multivector fields:
--
-- * In-place ST array accelerated kernel sliding
-- * Configurable boundary conditions ('ZeroPad', 'Periodic', 'Replicate', 'Mirror')
-- * Configurable algebraic convolution operations (Geometric, Wedge, Left Contraction, Scalar Product)
module Clifford.Array.Convolution
  ( -- * Boundary & Product Configurations
    BoundaryMode(..)
  , ProductMode(..)
    -- * Sliding Stencils
  , convolve1D
  , convolve2D
  , convolveND
    -- * Index Mapping Helper
  , mapBoundaryIndex
  ) where

import Prelude hiding ((<*>))
import Control.Monad (forM_)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray (listArray)
import Data.Array.MArray (newArray, writeArray)
import Data.Array.ST (runSTUArray)
import Data.Bits (shiftL)
import Data.Proxy (Proxy(..))

import Clifford.Signature
import Clifford.Field
import Clifford.Class
import Clifford.Dense
import Clifford.Array.Tensor

-- | Boundary padding strategy when kernel extends outside the array boundary.
data BoundaryMode
  = ZeroPad   -- ^ Out-of-bounds cells evaluate to \(0\)
  | Periodic  -- ^ Toroidal wrap-around / circular padding
  | Replicate -- ^ Clamp to edge coordinate
  | Mirror    -- ^ Reflect coordinates across edge
  deriving (Eq, Show)

-- | Algebraic product used to combine kernel weight and input field.
data ProductMode
  = GeometricProd    -- ^ Geometric product: \(W \star X\)
  | WedgeProd        -- ^ Wedge / Exterior product: \(W \wedge X\)
  | LeftContractProd -- ^ Left contraction: \(W \rfloor X\)
  deriving (Eq, Show)

-- | Map a 1D coordinate with respect to boundary condition.
mapBoundaryIndex :: BoundaryMode -> Int -> Int -> Maybe Int
mapBoundaryIndex mode dim idx
  | idx >= 0 && idx < dim = Just idx
  | otherwise = case mode of
      ZeroPad -> Nothing
      Periodic -> Just ((idx `mod` dim + dim) `mod` dim)
      Replicate -> Just (max 0 (min (dim - 1) idx))
      Mirror ->
        let period = 2 * dim - 2
            wrapped = (idx `mod` period + period) `mod` period
        in Just (if wrapped < dim then wrapped else period - wrapped)
{-# INLINE mapBoundaryIndex #-}

-- | 1D Clifford Convolution:
--
-- \(Y(x) = \sum_{k} W(k) \star X(x + k - \text{center})\)
convolve1D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => BoundaryMode
           -> ProductMode
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
convolve1D bMode pMode kernel input = convolveND bMode pMode kernel input
{-# INLINE convolve1D #-}

-- | 2D Clifford Convolution for spatial multivector fields / images:
--
-- \(Y(r, c) = \sum_{kr, kc} W(kr, kc) \star X(r + kr - cr, c + kc - cc)\)
convolve2D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => BoundaryMode
           -> ProductMode
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
convolve2D bMode pMode kernel input = convolveND bMode pMode kernel input
{-# INLINE convolve2D #-}

-- | General N-Dimensional Clifford Convolution across conforming grids.
convolveND :: forall p q r k. (KnownSignature p q r, UnboxedField k)
           => BoundaryMode
           -> ProductMode
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
convolveND bMode pMode kernel input =
  let inDims = shapeND input
      kDims = shapeND kernel
      inStrides = mvaStrides input
      kStrides = mvaStrides kernel
      totInCells = totalElements input
      totKCells = totalElements kernel
      bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
      bufIn = mvaBuffer input
      bufK = mvaBuffer kernel

      kCenter = [d `div` 2 | d <- kDims]

      op :: DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
      op = case pMode of
        GeometricProd    -> (<*>)
        WedgeProd        -> (/\)
        LeftContractProd -> (<.>)

      getInMV :: Int -> DenseMV (Signature p q r) k
      getInMV cellIdx =
        let cellOffset = cellIdx * bladesPerCell
        in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufIn (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

      getKMV :: Int -> DenseMV (Signature p q r) k
      getKMV cellIdx =
        let cellOffset = cellIdx * bladesPerCell
        in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufK (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

      outBuf = runSTUArray $ do
        marr <- newArray (0, totInCells * bladesPerCell - 1) (fromMetricScale 0)
        forM_ [0 .. totInCells - 1] $ \outIdx -> do
          let outCoords = flatToCoords inStrides outIdx
              accumK kIdx acc =
                let kCoords = flatToCoords kStrides kIdx
                    sampleCoords = zipWith3 (\o k c -> o + k - c) outCoords kCoords kCenter
                    validCoords = sequence (zipWith (mapBoundaryIndex bMode) inDims sampleCoords)
                in case validCoords of
                     Nothing -> acc
                     Just validSample ->
                       let inCellIdx = coordsToFlat inStrides validSample
                           wVal = getKMV kIdx
                           xVal = getInMV inCellIdx
                       in acc `addMV` (wVal `op` xVal)
              resMV = foldr accumK (zeroMV :: DenseMV (Signature p q r) k) [0 .. totKCells - 1]
              DenseMV resArr = resMV
              outCellOffset = outIdx * bladesPerCell
          forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
            writeArray marr (outCellOffset + bIdx) (unsafeAt resArr bIdx)
        return marr
  in MVArray inDims inStrides outBuf
{-# INLINE convolveND #-}
