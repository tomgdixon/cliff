{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Versor
-- Description : Rotors, Motors, Outermorphisms, Exponentials, and General Inverses
-- License     : BSD-3-Clause
--
-- This module implements versor and rotor transformations in Geometric Algebra:
--
-- * Sandwich transformations (\(R v \widetilde{R}\))
-- * Versor exponential map (\(\exp(B) = \cos \|B\| + \frac{B}{\|B\|} \sin \|B\|\))
-- * Rotor logarithms, interpolation, and rotor construction between vectors
-- * General multivector inversion via \(2^n \times 2^n\) Gauss-Jordan elimination on structure matrices
module Clifford.Versor
  ( -- * Versor Actions & Outermorphism
    sandwich
  , outermorphism
  , reflectMV
    -- * Subspace Projection & Rejection
  , projectMV
  , rejectMV
    -- * Rotor Construction & Exponential Map
  , rotorNorm
  , normalizeRotor
  , rotorExp
  , rotorLog
  , rotorBetween
    -- * Inversion
  , versorInverse
  , generalInverse
  ) where

import Prelude hiding ((<*>))
import Control.Monad (forM, forM_, unless)
import Control.Monad.ST (ST, runST)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray ((!))
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (STArray, runSTArray)
import Data.Bits (shiftL)
import Data.Proxy (Proxy(..))

import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class
import Clifford.Dense

-- | Apply sandwich transformation \(V x \widetilde{V}\).
sandwich :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
sandwich v x = v <*> x <*> reverseMV v
{-# INLINE sandwich #-}

-- | Reflect a vector / multivector \(V\) across a hyperplane with normal 1-vector \(N\):
--
-- \[ V' = - N \star V \star N^{-1} \]
reflectMV :: CliffordAlgebra sig k mv => mv -> mv -> Maybe mv
reflectMV v n = case versorInverse n of
  Just invN -> Just (negateMV (n <*> v <*> invN))
  Nothing   -> Nothing
{-# INLINE reflectMV #-}

-- | Project a multivector \(A\) onto an invertible subspace blade \(B\):
--
-- \[ \operatorname{Proj}_B(A) = (A \rfloor B) \star B^{-1} \]
projectMV :: CliffordAlgebra sig k mv => mv -> mv -> Maybe mv
projectMV a b = case versorInverse b of
  Just invB -> Just ((a <.> b) <*> invB)
  Nothing   -> Nothing
{-# INLINE projectMV #-}

-- | Reject a multivector \(A\) from an invertible subspace blade \(B\) (perpendicular component):
--
-- \[ \operatorname{Rej}_B(A) = (A \wedge B) \star B^{-1} \]
--
-- Satisfies the orthogonal decomposition identity: \(A = \operatorname{Proj}_B(A) + \operatorname{Rej}_B(A)\).
rejectMV :: CliffordAlgebra sig k mv => mv -> mv -> Maybe mv
rejectMV a b = case versorInverse b of
  Just invB -> Just ((a /\ b) <*> invB)
  Nothing   -> Nothing
{-# INLINE rejectMV #-}

-- | Apply outermorphism induced by linear map \(f\) on arbitrary blade/multivector \(A\).
outermorphism :: (CliffordAlgebra sig k mv) => (mv -> mv) -> mv -> mv
outermorphism f a =
  let blades = toBladeList a
      transformBlade (b, coeff) =
        let idxs = bladeIndices b
        in if null idxs
             then scalarMV coeff
             else
               let mappedVectors = map (\i -> f (basisMV i)) idxs
                   wedgeAll []     = scalarMV 1
                   wedgeAll (v:vs) = foldl (/\) v vs
               in scaleMV coeff (wedgeAll mappedVectors)
  in foldr addMV zeroMV (map transformBlade blades)
{-# INLINE outermorphism #-}

-- | Rotor magnitude norm \(\|R\| = \sqrt{|\langle R \widetilde{R} \rangle_0|}\).
rotorNorm :: (CliffordAlgebra sig k mv, Floating k) => mv -> k
rotorNorm r = sqrt (abs (normSqMV r))
{-# INLINE rotorNorm #-}

-- | Normalize a rotor so that \(R \widetilde{R} = 1\).
normalizeRotor :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv
normalizeRotor r =
  let n = rotorNorm r
  in if isZero n then r else scaleMV (1 / n) r
{-# INLINE normalizeRotor #-}

-- | Closed-form exponential of a simple 2-blade / bivector \(B\):
--
-- \(\exp(B) = \cos(\theta) + \frac{B}{\theta} \sin(\theta)\) where \(\theta = \sqrt{|\langle B \widetilde{B} \rangle_0|}\).
rotorExp :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv
rotorExp b =
  let thetaSq = normSqMV b
      theta = sqrt (abs thetaSq)
  in if isZero theta
       then scalarMV 1 `addMV` b
       else
         let cosTerm = scalarMV (cos theta)
             sinTerm = scaleMV (sin theta / theta) b
         in cosTerm `addMV` sinTerm
{-# INLINE rotorExp #-}

-- | Logarithm of a rotor: returns the generator bivector \(B\) such that \(\exp(B) = R\).
rotorLog :: (CliffordAlgebra sig k mv, RealFloat k) => mv -> mv
rotorLog r =
  let s = getBlade scalarBlade r
      b = gradeProj 2 r
      bNorm = rotorNorm b
  in if isZero bNorm
       then zeroMV
       else
         let theta = atan2 bNorm s
         in scaleMV (theta / bNorm) b
{-# INLINE rotorLog #-}

-- | Construct the shortest-arc rotor taking unit vector \(u\) to unit vector \(v\):
--
-- \(R = \frac{1 + v u}{\|1 + v u\|}\).
rotorBetween :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv -> mv
rotorBetween u v =
  let r = scalarMV 1 `addMV` (v <*> u)
  in normalizeRotor r
{-# INLINE rotorBetween #-}

-- | Multiplicative inverse of a versor \(V\):
--
-- \(V^{-1} = \frac{\widetilde{V}}{\langle V \widetilde{V} \rangle_0}\).
versorInverse :: (CliffordAlgebra sig k mv) => mv -> Maybe mv
versorInverse v =
  let denom = normSqMV v
  in case fieldInverse denom of
       Nothing -> Nothing
       Just dInv -> Just (scaleMV dInv (reverseMV v))
{-# INLINE versorInverse #-}

-- | General multivector inversion via Gauss-Jordan elimination on its \(2^n \times 2^n\) Cayley representation matrix.
generalInverse :: forall p q r k. (KnownSignature p q r, UnboxedField k)
               => DenseMV (Signature p q r) k -> Maybe (DenseMV (Signature p q r) k)
generalInverse a =
  let dim = signatureDim (Proxy :: Proxy (Signature p q r))
      tot = 1 `shiftL` dim

      -- Build augmented matrix [M | I] of size tot x (2*tot)
      augMat = runSTArray $ do
        arr <- newArray ((0, 0), (tot - 1, 2 * tot - 1)) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \j -> do
          let ej = bladeMV (denseBladeAt j) (fromMetricScale 1) :: DenseMV (Signature p q r) k
              colMV = a <*> ej
              DenseMV colArr = colMV
          forM_ [0 .. tot - 1] $ \i -> do
            writeArray arr (i, j) (unsafeAt colArr i)
        forM_ [0 .. tot - 1] $ \i -> do
          writeArray arr (i, tot + i) (fromMetricScale 1)
        return arr

      -- Run in-place Gauss-Jordan elimination in ST
      result = runST $ do
        marr <- newArray ((0, 0), (tot - 1, 2 * tot - 1)) (fromMetricScale 0) :: ST s (STArray s (Int, Int) k)
        forM_ [0 .. tot - 1] $ \i ->
          forM_ [0 .. 2 * tot - 1] $ \j ->
            writeArray marr (i, j) (augMat ! (i, j))

        let eliminate col = if col >= tot then return True else do
              -- Find pivot
              pivotRow <- findPivot col col
              case pivotRow of
                Nothing -> return False
                Just pRow -> do
                  swapRows pRow col
                  pivotVal <- readArray marr (col, col)
                  case fieldInverse pivotVal of
                    Nothing -> return False
                    Just pInv -> do
                      -- Scale pivot row
                      forM_ [0 .. 2 * tot - 1] $ \j -> do
                        val <- readArray marr (col, j)
                        writeArray marr (col, j) (val * pInv)
                      -- Eliminate other rows
                      forM_ [0 .. tot - 1] $ \rIdx -> do
                        unless (rIdx == col) $ do
                          factor <- readArray marr (rIdx, col)
                          unless (isZero factor) $ do
                            forM_ [0 .. 2 * tot - 1] $ \j -> do
                              vP <- readArray marr (col, j)
                              vR <- readArray marr (rIdx, j)
                              writeArray marr (rIdx, j) (vR - factor * vP)
                      eliminate (col + 1)

            findPivot r c
              | r >= tot = return Nothing
              | otherwise = do
                  v <- readArray marr (r, c)
                  if isZero v then findPivot (r + 1) c else return (Just r)

            swapRows r1 r2
              | r1 == r2 = return ()
              | otherwise = do
                  forM_ [0 .. 2 * tot - 1] $ \j -> do
                    v1 <- readArray marr (r1, j)
                    v2 <- readArray marr (r2, j)
                    writeArray marr (r1, j) v2
                    writeArray marr (r2, j) v1

        success <- eliminate 0
        if not success
          then return Nothing
          else do
            -- Extract first column of inverse (since a * aInv = 1 => col 0 corresponds to basis blade 1)
            invCoeffs <- forM [0 .. tot - 1] $ \i -> readArray marr (i, tot)
            return (Just (mkDense (zip (allBlades dim) invCoeffs)))
  in result
