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
-- Module      : Clifford.Array.Matrix
-- Description : Multivector Matrix Multiply (GEMM), Structure Matrices, and Rotor-Invariant Attention
-- License     : BSD-3-Clause
--
-- This module implements matrix operations over Clifford multivectors:
--
-- * Multivector GEMM (\(C = A \star B\))
-- * Multivector Matrix-Vector multiplication
-- * Structure Matrix algebra homomorphism (\(M(a \star b) = M(a) \cdot M(b)\))
-- * Rotor-invariant Clifford attention scoring
-- * Self-adjoint Clifford backpropagation gradient rules
module Clifford.Array.Matrix
  ( -- * Multivector Matrix Type
    MVMatrix(..)
  , mkMatrix
  , matrixIndex
    -- * Multivector GEMM
  , matMulMV
  , matVecMulMV
    -- * Clifford Structure Matrix Homomorphism
  , toStructureMatrix
    -- * CliffordNet LLM Attention & Gradients
  , attentionScoresInvariant
  , adjointMulLeft
  , adjointMulRight
  ) where

import Prelude hiding ((<*>))
import Control.DeepSeq (NFData(..))
import Control.Monad (forM_)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray (Array, IArray, (!), array, listArray)
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (runSTUArray)
import Data.Bits (shiftL)
import Data.Proxy (Proxy(..))

import Clifford.Signature
import Clifford.Blade (allBlades)
import Clifford.Field
import Clifford.Class
import Clifford.Dense
import Clifford.Array.Tensor

-- | Multivector Matrix with dimension \((M \times N)\) stored as a contiguous 2D 'MVArray'.
newtype MVMatrix sig k = MVMatrix { unMVMatrix :: MVArray sig k }

instance NFData (MVMatrix sig k) where
  rnf (MVMatrix arr) = rnf arr

instance (KnownSignature p q r, UnboxedField k, Eq k) => Eq (MVMatrix (Signature p q r) k) where
  (MVMatrix a) == (MVMatrix b) = a == b

instance Show (MVMatrix (Signature p q r) k) where
  show (MVMatrix arr) =
    case shapeND arr of
      [rows, cols] -> "MVMatrix (" ++ show rows ++ "x" ++ show cols ++ ")"
      _            -> "MVMatrix (invalid shape)"

-- | Construct an \((M \times N)\) multivector matrix from a row-major list of multivectors.
mkMatrix :: forall p q r k. (KnownSignature p q r, UnboxedField k)
         => Int -> Int -> [DenseMV (Signature p q r) k] -> MVMatrix (Signature p q r) k
mkMatrix rows cols mvs
  | length mvs /= rows * cols = error "Clifford.Array.Matrix.mkMatrix: multivector count does not match rows * cols"
  | otherwise = MVMatrix (fromListND [rows, cols] mvs)
{-# INLINE mkMatrix #-}

-- | Index into a multivector matrix at \((row, col)\).
matrixIndex :: forall p q r k. (KnownSignature p q r, UnboxedField k)
            => MVMatrix (Signature p q r) k -> Int -> Int -> DenseMV (Signature p q r) k
matrixIndex (MVMatrix arr) r c = indexND arr [r, c]
{-# INLINE matrixIndex #-}

-- | Multivector Matrix-Matrix Multiplication (GEMM): \(C_{ij} = \sum_k A_{ik} \star B_{kj}\).
matMulMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
         => MVMatrix (Signature p q r) k -> MVMatrix (Signature p q r) k -> MVMatrix (Signature p q r) k
matMulMV (MVMatrix arrA) (MVMatrix arrB) =
  case (shapeND arrA, shapeND arrB) of
    ([m, kA], [kB, n])
      | kA /= kB  -> error $ "Clifford.Array.Matrix.matMulMV: inner dimension mismatch: " ++ show kA ++ " vs " ++ show kB
      | otherwise ->
          let k = kA
              bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
              bufA = mvaBuffer arrA
              bufB = mvaBuffer arrB
              stridesA = mvaStrides arrA
              stridesB = mvaStrides arrB
              outDims = [m, n]
              outStrides = computeStrides outDims

              getMVA :: Int -> Int -> DenseMV (Signature p q r) k
              getMVA rIdx cIdx =
                let cellOffset = coordsToFlat stridesA [rIdx, cIdx] * bladesPerCell
                in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufA (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

              getMVB :: Int -> Int -> DenseMV (Signature p q r) k
              getMVB rIdx cIdx =
                let cellOffset = coordsToFlat stridesB [rIdx, cIdx] * bladesPerCell
                in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt bufB (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

              outBuf = runSTUArray $ do
                marr <- newArray (0, m * n * bladesPerCell - 1) (fromMetricScale 0)
                forM_ [0 .. m - 1] $ \i ->
                  forM_ [0 .. n - 1] $ \j -> do
                    let sumTerms pIdx acc =
                          let a_ip = getMVA i pIdx
                              b_pj = getMVB pIdx j
                          in acc `addMV` (a_ip <*> b_pj)
                        resMV = foldr sumTerms (zeroMV :: DenseMV (Signature p q r) k) [0 .. k - 1]
                        DenseMV resArr = resMV
                        outCellOffset = coordsToFlat outStrides [i, j] * bladesPerCell
                    forM_ [0 .. bladesPerCell - 1] $ \bIdx ->
                      writeArray marr (outCellOffset + bIdx) (unsafeAt resArr bIdx)
                return marr
          in MVMatrix (MVArray outDims outStrides outBuf)
    _ -> error "Clifford.Array.Matrix.matMulMV: inputs must be 2D matrices"
{-# INLINE matMulMV #-}

-- | Multivector Matrix-Vector Multiplication: \(y_i = \sum_k A_{ik} \star x_k\).
matVecMulMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
            => MVMatrix (Signature p q r) k -> [DenseMV (Signature p q r) k] -> [DenseMV (Signature p q r) k]
matVecMulMV mat vec =
  let [m, n] = shapeND (unMVMatrix mat)
  in if length vec /= n
       then error "Clifford.Array.Matrix.matVecMulMV: vector dimension mismatch"
       else [ foldr (\j acc -> acc `addMV` (matrixIndex mat i j <*> (vec !! j))) (zeroMV :: DenseMV (Signature p q r) k) [0 .. n - 1]
            | i <- [0 .. m - 1]
            ]
{-# INLINE matVecMulMV #-}

-- | Construct the \(2^n \times 2^n\) real/field structure matrix \(M(a)\) of a multivector \(a\).
--
-- Satisfies the fundamental algebraic isomorphism: \(M(a \star b) = M(a) \cdot M(b)\).
toStructureMatrix :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                  => DenseMV (Signature p q r) k -> Array (Int, Int) k
toStructureMatrix a =
  let dim = signatureDim (Proxy :: Proxy (Signature p q r))
      tot = 1 `shiftL` dim
      blades = allBlades dim
      indices = [((i, j), coeff i j) | i <- [0 .. tot - 1], j <- [0 .. tot - 1]]
      coeff i j =
        let ej = bladeMV (blades !! j) (fromMetricScale 1) :: DenseMV (Signature p q r) k
            colMV = a <*> ej
        in getBlade (blades !! i) colMV
  in array ((0, 0), (tot - 1, tot - 1)) indices
{-# INLINE toStructureMatrix #-}

-- | Rotor-invariant Clifford attention scoring:
--
-- \(S(Q, K)_{ij} = \langle Q_i \widetilde{K_j} \rangle_0\).
--
-- Satisfies exact invariance under arbitrary rotor group transformations \(Q \mapsto R Q \widetilde{R}\) and \(K \mapsto R K \widetilde{R}\).
attentionScoresInvariant :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                         => MVMatrix (Signature p q r) k
                         -> MVMatrix (Signature p q r) k
                         -> Array (Int, Int) k
attentionScoresInvariant (MVMatrix qArr) (MVMatrix kArr) =
  case (shapeND qArr, shapeND kArr) of
    ([numQ, dimQ], [numK, dimK])
      | dimQ /= dimK -> error "Clifford.Array.Matrix.attentionScoresInvariant: query and key feature dimension mismatch"
      | otherwise ->
          let indices = [ ((i, j), score i j) | i <- [0 .. numQ - 1], j <- [0 .. numK - 1] ]
              score i j =
                let pairDots = [ let qi = indexND qArr [i, d]
                                     kj = indexND kArr [j, d]
                                 in (qi :: DenseMV (Signature p q r) k) <|> (kj :: DenseMV (Signature p q r) k)
                               | d <- [0 .. dimQ - 1]
                               ]
                in foldr (+) (fromMetricScale 0) pairDots
          in array ((0, 0), (numQ - 1, numK - 1)) indices
    _ -> error "Clifford.Array.Matrix.attentionScoresInvariant: inputs must be 2D matrices"
{-# INLINE attentionScoresInvariant #-}

-- | Adjoint gradient rule for left multiplication \(Z = X \star Y\):
--
-- \(\frac{\partial L}{\partial X} = \frac{\partial L}{\partial Z} \star \widetilde{Y}\).
adjointMulLeft :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
adjointMulLeft gradZ y = gradZ <*> reverseMV y
{-# INLINE adjointMulLeft #-}

-- | Adjoint gradient rule for right multiplication \(Z = X \star Y\):
--
-- \(\frac{\partial L}{\partial Y} = \widetilde{X} \star \frac{\partial L}{\partial Z}\).
adjointMulRight :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
adjointMulRight x gradZ = reverseMV x <*> gradZ
{-# INLINE adjointMulRight #-}
