{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Array.Fourier
-- Description : Fast Clifford Fourier Transforms (FCFT), Number Theoretic Transforms (NTT), and Fast Walsh-Hadamard Transforms (FWHT)
-- License     : BSD-3-Clause
--
-- This module implements multidimensional Fourier and Walsh-Hadamard spectral transforms
-- across continuous, rational, and finite fields:
--
-- * __Continuous Fast Clifford Fourier Transform (FCFT)__: Radix-2 Cooley-Tukey butterfly
--   using commuting bivector/pseudoscalar roots of unity (\(\exp(-I \frac{2\pi k}{N})\)).
-- * __Galois Field Number Theoretic Transform (NTT)__: Exact \(O(N \log N)\) modular butterfly
--   over prime fields \(\mathbb{F}_p\) using primitive roots of unity with zero precision loss.
-- * __Fast Walsh-Hadamard Transform (FWHT)__: Additive butterfly over the Boolean hypercube
--   character group, enabling fast multivector Walsh spectrum analysis and dyadic (XOR) convolutions.
module Clifford.Array.Fourier
  ( -- * Fast Clifford Fourier Transforms (O(N log N))
    fcft1D
  , ifcft1D
  , fcft2D
  , ifcft2D
    -- * Galois Field Number Theoretic Transforms (F_p)
  , ntt1D
  , intt1D
  , ntt2D
  , intt2D
    -- * Fast Walsh-Hadamard Transforms (FWHT & Dyadic Convolutions)
  , fwhtMV
  , ifwhtMV
  , dyadicConvolveMV
  , fwht1D
  , ifwht1D
  , fwht2D
  , ifwht2D
  , dyadicConvolve1D
    -- * Direct Discrete Transforms & Kernels (Reference)
  , CFTKernel(..)
  , mkCFTKernel
  , defaultCFT2D
  , defaultCFT3D
  , cftForward
  , cftInverse
  , spectralFilter
  ) where

import Prelude hiding ((<*>))
import Control.DeepSeq (NFData(..))
import Control.Monad (forM_)
import Control.Monad.ST (ST)
import Data.Array.Base (unsafeAt)
import Data.Array.IArray (listArray)
import Data.Array.MArray (newArray, readArray, writeArray)
import Data.Array.ST (STUArray, runSTUArray)
import Data.Bits (shiftL, shiftR, testBit, setBit, (.&.))
import Data.List (zipWith4)
import Data.Proxy (Proxy(..))
import GHC.TypeLits

import Clifford.Signature
import Clifford.Blade
import Clifford.Field
import Clifford.Class
import Clifford.Dense
import Clifford.Array.Tensor

--------------------------------------------------------------------------------
-- Bit Reversal & Power of 2 Helpers
--------------------------------------------------------------------------------

isPowerOfTwo :: Int -> Bool
isPowerOfTwo x = x > 0 && (x .&. (x - 1)) == 0

log2Int :: Int -> Int
log2Int n
  | n <= 1    = 0
  | otherwise = 1 + log2Int (n `shiftR` 1)

bitReverse :: Int -> Int -> Int
bitReverse x bits = go 0 0
  where
    go !acc !idx
      | idx >= bits = acc
      | testBit x idx = go (setBit acc (bits - 1 - idx)) (idx + 1)
      | otherwise     = go acc (idx + 1)

--------------------------------------------------------------------------------
-- 1D & 2D Fast Clifford Fourier Transform (Radix-2 Cooley-Tukey FCFT)
--------------------------------------------------------------------------------

-- | 1D Forward Fast Clifford Fourier Transform (FCFT) in \(O(N \log N)\) operations
-- using Radix-2 Cooley-Tukey decimation-in-time on an invertible unit bivector / pseudoscalar \(I\) (\(I^2 = -1\)).
fcft1D :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
       => DenseMV (Signature p q r) k  -- ^ Unit bivector / pseudoscalar \(I\) (\(I^2 = -1\))
       -> MVArray (Signature p q r) k
       -> MVArray (Signature p q r) k
fcft1D iGen input = runRadix2FCFT (-1) iGen input

-- | 1D Inverse Fast Clifford Fourier Transform (IFCFT) in \(O(N \log N)\) operations.
ifcft1D :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
        => DenseMV (Signature p q r) k
        -> MVArray (Signature p q r) k
        -> MVArray (Signature p q r) k
ifcft1D iGen input =
  let fwd = runRadix2FCFT 1 iGen input
      n = totalElements input
      scaleFactor = 1 / fromIntegral n
  in mapMV (scaleMV scaleFactor) fwd

runRadix2FCFT :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
              => k -> DenseMV (Signature p q r) k -> MVArray (Signature p q r) k -> MVArray (Signature p q r) k
runRadix2FCFT signMult iGen input
  | not (isPowerOfTwo n) = error "Clifford.Array.Fourier.fcft1D: array length must be a power of 2"
  | otherwise =
      let bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          numBits = log2Int n
          buf = mvaBuffer input

          calcTwiddle len kIdx =
            let angle = signMult * 2 * pi * fromIntegral kIdx / fromIntegral len
                cosTerm = scalarMV (cos angle)
                sinTerm = scaleMV (sin angle) iGen
            in cosTerm `addMV` sinTerm

          outBuf = runSTUArray $ do
            (marr :: STUArray s Int k) <- newArray (0, n * bladesPerCell - 1) (fromMetricScale 0)
            -- 1. Bit-reversal permutation copy
            forM_ [0 .. n - 1] $ \i -> do
              let revI = bitReverse i numBits
                  srcOffset = i * bladesPerCell
                  dstOffset = revI * bladesPerCell
              forM_ [0 .. bladesPerCell - 1] $ \b ->
                writeArray marr (dstOffset + b) (unsafeAt buf (srcOffset + b))

            let readMV :: STUArray s Int k -> Int -> ST s (DenseMV (Signature p q r) k)
                readMV arr cellIdx = do
                  let offset = cellIdx * bladesPerCell
                  coeffs <- mapM (\b -> readArray arr (offset + b)) [0 .. bladesPerCell - 1]
                  return (DenseMV (listArray (0, bladesPerCell - 1) coeffs))

                writeMV :: STUArray s Int k -> Int -> DenseMV (Signature p q r) k -> ST s ()
                writeMV arr cellIdx (DenseMV dArr) = do
                  let offset = cellIdx * bladesPerCell
                  forM_ [0 .. bladesPerCell - 1] $ \b ->
                    writeArray arr (offset + b) (unsafeAt dArr b)

            -- 2. Cooley-Tukey Butterfly stages
            let loopStages s
                  | s > numBits = return ()
                  | otherwise = do
                      let len = 1 `shiftL` s
                          half = len `shiftR` 1
                          step = len
                          loopBlocks i
                            | i >= n = return ()
                            | otherwise = do
                                forM_ [0 .. half - 1] $ \j -> do
                                  let uIdx = i + j
                                      vIdx = i + j + half
                                      w = calcTwiddle len j
                                  uMV <- readMV marr uIdx
                                  vMV <- readMV marr vIdx
                                  let vw = vMV <*> w
                                      newU = uMV `addMV` vw
                                      newV = uMV `subMV` vw
                                  writeMV marr uIdx newU
                                  writeMV marr vIdx newV
                                loopBlocks (i + step)
                      loopBlocks 0
                      loopStages (s + 1)
            loopStages 1
            return marr
      in MVArray (shapeND input) (mvaStrides input) outBuf
  where
    n = totalElements input

-- | 2D Forward Fast Clifford Fourier Transform (FCFT) in \(O(H W \log(HW))\) operations
-- via row-column separable 1D Cooley-Tukey butterflies.
fcft2D :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
       => DenseMV (Signature p q r) k  -- ^ Row bivector kernel \(I_1\)
       -> DenseMV (Signature p q r) k  -- ^ Column bivector kernel \(I_2\)
       -> MVArray (Signature p q r) k
       -> MVArray (Signature p q r) k
fcft2D iRow iCol input = case shapeND input of
  [h, w]
    | not (isPowerOfTwo h && isPowerOfTwo w) -> error "Clifford.Array.Fourier.fcft2D: dimensions must be powers of 2"
    | otherwise ->
        let rowTransformed = [ fcft1D iRow (fromListND [w] [ indexND input [r, c] | c <- [0 .. w - 1] ]) | r <- [0 .. h - 1] ]
            intermediateGrid = fromListND [h, w] [ indexND (rowTransformed !! r) [c] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
            colTransformed = [ fcft1D iCol (fromListND [h] [ indexND intermediateGrid [r, c] | r <- [0 .. h - 1] ]) | c <- [0 .. w - 1] ]
        in fromListND [h, w] [ indexND (colTransformed !! c) [r] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
  _ -> error "Clifford.Array.Fourier.fcft2D: input must be 2-dimensional"

-- | 2D Inverse Fast Clifford Fourier Transform (IFCFT).
ifcft2D :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
        => DenseMV (Signature p q r) k
        -> DenseMV (Signature p q r) k
        -> MVArray (Signature p q r) k
        -> MVArray (Signature p q r) k
ifcft2D iRow iCol input = case shapeND input of
  [h, w] ->
    let fwd = fcft2D (negateMV iRow) (negateMV iCol) input
        scaleFactor = 1 / fromIntegral (h * w)
    in mapMV (scaleMV scaleFactor) fwd
  _ -> error "Clifford.Array.Fourier.ifcft2D: input must be 2-dimensional"

--------------------------------------------------------------------------------
-- Number Theoretic Transform (NTT) for Prime Fields GF(p)
--------------------------------------------------------------------------------

-- | 1D Forward Number Theoretic Transform (NTT) over prime modular field \(\mathbb{F}_p\) in \(O(N \log N)\)
-- exact integer arithmetic using a primitive \(N\)-th root of unity \(\omega \in \mathbb{F}_p\).
ntt1D :: forall p q r (prime :: Nat). (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime))
      => GFp prime  -- ^ Primitive \(N\)-th root of unity \(\omega\) satisfying \(\omega^N \equiv 1 \pmod p\)
      -> MVArray (Signature p q r) (GFp prime)
      -> MVArray (Signature p q r) (GFp prime)
ntt1D omega input
  | not (isPowerOfTwo n) = error "Clifford.Array.Fourier.ntt1D: array length must be a power of 2"
  | otherwise =
      let bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          numBits = log2Int n
          buf = mvaBuffer input

          getOmegaLen len =
            let pVal = natVal (Proxy :: Proxy prime)
                expVal = fromIntegral (n `div` len)
                baseVal = fromIntegral (unGFp omega)
                wVal = modPow baseVal expVal pVal
            in GFp (fromIntegral wVal)

          outBuf = runSTUArray $ do
            (marr :: STUArray s Int (GFp prime)) <- newArray (0, n * bladesPerCell - 1) (GFp 0)
            -- 1. Bit-reversal permutation copy
            forM_ [0 .. n - 1] $ \i -> do
              let revI = bitReverse i numBits
                  srcOffset = i * bladesPerCell
                  dstOffset = revI * bladesPerCell
              forM_ [0 .. bladesPerCell - 1] $ \b ->
                writeArray marr (dstOffset + b) (unsafeAt buf (srcOffset + b))

            let readMV :: STUArray s Int (GFp prime) -> Int -> ST s (DenseMV (Signature p q r) (GFp prime))
                readMV arr cellIdx = do
                  let offset = cellIdx * bladesPerCell
                  coeffs <- mapM (\b -> readArray arr (offset + b)) [0 .. bladesPerCell - 1]
                  return (DenseMV (listArray (0, bladesPerCell - 1) coeffs))

                writeMV :: STUArray s Int (GFp prime) -> Int -> DenseMV (Signature p q r) (GFp prime) -> ST s ()
                writeMV arr cellIdx (DenseMV dArr) = do
                  let offset = cellIdx * bladesPerCell
                  forM_ [0 .. bladesPerCell - 1] $ \b ->
                    writeArray arr (offset + b) (unsafeAt dArr b)

            -- 2. Modular Butterfly stages
            let loopStages s
                  | s > numBits = return ()
                  | otherwise = do
                      let len = 1 `shiftL` s
                          half = len `shiftR` 1
                          step = len
                          wLen = getOmegaLen len
                          loopBlocks i
                            | i >= n = return ()
                            | otherwise = do
                                let loopJ j !wAcc
                                      | j >= half = return ()
                                      | otherwise = do
                                          let uIdx = i + j
                                              vIdx = i + j + half
                                          uMV <- readMV marr uIdx
                                          vMV <- readMV marr vIdx
                                          let vw = scaleMV wAcc vMV
                                              newU = uMV `addMV` vw
                                              newV = uMV `subMV` vw
                                          writeMV marr uIdx newU
                                          writeMV marr vIdx newV
                                          let nextW = wAcc * wLen
                                          loopJ (j + 1) nextW
                                loopJ 0 (GFp 1)
                                loopBlocks (i + step)
                      loopBlocks 0
                      loopStages (s + 1)
            loopStages 1
            return marr
      in MVArray (shapeND input) (mvaStrides input) outBuf
  where
    n = totalElements input

-- | 1D Inverse Number Theoretic Transform (INTT) over \(\mathbb{F}_p\).
intt1D :: forall p q r (prime :: Nat). (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime))
       => GFp prime  -- ^ Primitive \(N\)-th root of unity \(\omega\)
       -> MVArray (Signature p q r) (GFp prime)
       -> MVArray (Signature p q r) (GFp prime)
intt1D omega input =
  case fieldInverse omega of
    Nothing -> error "Clifford.Array.Fourier.intt1D: omega is not invertible"
    Just omegaInv ->
      let fwd = ntt1D omegaInv input
          n = totalElements input
          nMod = fromIntegral n :: GFp prime
      in case fieldInverse nMod of
           Just invN -> mapMV (scaleMV invN) fwd
           Nothing   -> error "Clifford.Array.Fourier.intt1D: length N is not invertible modulo p"

-- | 2D Forward Number Theoretic Transform (NTT) over \(\mathbb{F}_p\).
ntt2D :: forall p q r (prime :: Nat). (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime))
      => GFp prime  -- ^ Row primitive root of unity \(\omega_1\)
      -> GFp prime  -- ^ Column primitive root of unity \(\omega_2\)
      -> MVArray (Signature p q r) (GFp prime)
      -> MVArray (Signature p q r) (GFp prime)
ntt2D wRow wCol input = case shapeND input of
  [h, w]
    | not (isPowerOfTwo h && isPowerOfTwo w) -> error "Clifford.Array.Fourier.ntt2D: dimensions must be powers of 2"
    | otherwise ->
        let rowTransformed = [ ntt1D wRow (fromListND [w] [ indexND input [r, c] | c <- [0 .. w - 1] ]) | r <- [0 .. h - 1] ]
            intermediateGrid = fromListND [h, w] [ indexND (rowTransformed !! r) [c] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
            colTransformed = [ ntt1D wCol (fromListND [h] [ indexND intermediateGrid [r, c] | r <- [0 .. h - 1] ]) | c <- [0 .. w - 1] ]
        in fromListND [h, w] [ indexND (colTransformed !! c) [r] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
  _ -> error "Clifford.Array.Fourier.ntt2D: input must be 2-dimensional"

-- | 2D Inverse Number Theoretic Transform (INTT) over \(\mathbb{F}_p\).
intt2D :: forall p q r (prime :: Nat). (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime))
       => GFp prime  -- ^ Row primitive root of unity \(\omega_1\)
       -> GFp prime  -- ^ Column primitive root of unity \(\omega_2\)
       -> MVArray (Signature p q r) (GFp prime)
       -> MVArray (Signature p q r) (GFp prime)
intt2D wRow wCol input =
  case (fieldInverse wRow, fieldInverse wCol) of
    (Just invRow, Just invCol) ->
      let fwd = ntt2D invRow invCol input
      in case shapeND input of
           [h, w] ->
             let nTotal = fromIntegral (h * w) :: GFp prime
             in case fieldInverse nTotal of
                  Just invN -> mapMV (scaleMV invN) fwd
                  Nothing   -> error "Clifford.Array.Fourier.intt2D: grid size is not invertible modulo p"
           _ -> error "Clifford.Array.Fourier.intt2D: input must be 2-dimensional"
    _ -> error "Clifford.Array.Fourier.intt2D: root of unity is not invertible"

--------------------------------------------------------------------------------
-- Fast Walsh-Hadamard Transforms (FWHT & Dyadic Convolutions)
--------------------------------------------------------------------------------

-- | Fast Walsh-Hadamard Transform of the \(2^n\) blade coefficients of a dense multivector.
-- Operates in \(O(n 2^n)\) additions and subtractions via in-place 'ST' array mutation
-- over the Boolean hypercube character lattice.
fwhtMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
       => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
fwhtMV (DenseMV inArr) =
  let dim = signatureDim (Proxy :: Proxy (Signature p q r))
      tot = 1 `shiftL` dim
      res = runSTUArray $ do
        marr <- newArray (0, tot - 1) (fromMetricScale 0)
        forM_ [0 .. tot - 1] $ \i -> writeArray marr i (unsafeAt inArr i)
        let loopStages s
              | s > dim = return ()
              | otherwise = do
                  let len = 1 `shiftL` s
                      half = len `shiftR` 1
                      step = len
                      loopBlocks i
                        | i >= tot = return ()
                        | otherwise = do
                            forM_ [0 .. half - 1] $ \j -> do
                              let uIdx = i + j
                                  vIdx = i + j + half
                              u <- readArray marr uIdx
                              v <- readArray marr vIdx
                              writeArray marr uIdx (u + v)
                              writeArray marr vIdx (u - v)
                            loopBlocks (i + step)
                  loopBlocks 0
                  loopStages (s + 1)
        loopStages 1
        return marr
  in DenseMV res
{-# INLINE fwhtMV #-}

-- | Inverse Fast Walsh-Hadamard Transform for dense multivector blade coefficients.
-- Satisfies @ifwhtMV (fwhtMV v) == v@ across all fields where \(2\) is invertible.
ifwhtMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
ifwhtMV v =
  let dim = signatureDim (Proxy :: Proxy (Signature p q r))
      tot = 1 `shiftL` dim
      scaleFactor = fromRational (1 / fromIntegral tot) :: k
      DenseMV arr = fwhtMV v
      scaled = listArray (0, tot - 1) [scaleFactor * unsafeAt arr i | i <- [0 .. tot - 1]]
  in DenseMV scaled
{-# INLINE ifwhtMV #-}

-- | Fast Dyadic (bitwise XOR) Convolution of two dense multivectors:
--
-- \((A \star_{\text{xor}} B)_k = \sum_{i \oplus j = k} A_i B_j\)
--
-- Evaluated in \(O(n 2^n)\) operations via the Walsh-Hadamard Convolution Theorem:
-- \(\text{FWHT}(A \star_{\text{xor}} B) = \text{FWHT}(A) \odot \text{FWHT}(B)\).
dyadicConvolveMV :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                 => DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k
dyadicConvolveMV a b =
  let DenseMV fA = fwhtMV a
      DenseMV fB = fwhtMV b
      dim = signatureDim (Proxy :: Proxy (Signature p q r))
      tot = 1 `shiftL` dim
      fProd = listArray (0, tot - 1) [unsafeAt fA i * unsafeAt fB i | i <- [0 .. tot - 1]]
  in ifwhtMV (DenseMV fProd)
{-# INLINE dyadicConvolveMV #-}

-- | 1D Fast Walsh-Hadamard Transform across an array of \(N = 2^m\) multivectors.
-- Transforms spatial multivector cells in \(O(N \log N)\) operations using additive butterflies.
fwht1D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
       => MVArray (Signature p q r) k -> MVArray (Signature p q r) k
fwht1D input
  | not (isPowerOfTwo n) = error "Clifford.Array.Fourier.fwht1D: array length must be a power of 2"
  | otherwise =
      let bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          numBits = log2Int n
          buf = mvaBuffer input

          outBuf = runSTUArray $ do
            marr <- newArray (0, n * bladesPerCell - 1) (fromMetricScale 0)
            forM_ [0 .. n * bladesPerCell - 1] $ \i -> writeArray marr i (unsafeAt buf i)

            let readMV :: STUArray s Int k -> Int -> ST s (DenseMV (Signature p q r) k)
                readMV arr cellIdx = do
                  let offset = cellIdx * bladesPerCell
                  coeffs <- mapM (\b -> readArray arr (offset + b)) [0 .. bladesPerCell - 1]
                  return (DenseMV (listArray (0, bladesPerCell - 1) coeffs))

                writeMV :: STUArray s Int k -> Int -> DenseMV (Signature p q r) k -> ST s ()
                writeMV arr cellIdx (DenseMV dArr) = do
                  let offset = cellIdx * bladesPerCell
                  forM_ [0 .. bladesPerCell - 1] $ \b ->
                    writeArray arr (offset + b) (unsafeAt dArr b)

            let loopStages s
                  | s > numBits = return ()
                  | otherwise = do
                      let len = 1 `shiftL` s
                          half = len `shiftR` 1
                          step = len
                          loopBlocks i
                            | i >= n = return ()
                            | otherwise = do
                                forM_ [0 .. half - 1] $ \j -> do
                                  let uIdx = i + j
                                      vIdx = i + j + half
                                  uMV <- readMV marr uIdx
                                  vMV <- readMV marr vIdx
                                  let newU = uMV `addMV` vMV
                                      newV = uMV `subMV` vMV
                                  writeMV marr uIdx newU
                                  writeMV marr vIdx newV
                                loopBlocks (i + step)
                      loopBlocks 0
                      loopStages (s + 1)
            loopStages 1
            return marr
      in MVArray (shapeND input) (mvaStrides input) outBuf
  where
    n = totalElements input
{-# INLINE fwht1D #-}

-- | 1D Inverse Fast Walsh-Hadamard Transform for multivector arrays.
-- Satisfies @ifwht1D (fwht1D arr) == arr@.
ifwht1D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => MVArray (Signature p q r) k -> MVArray (Signature p q r) k
ifwht1D input =
  let fwd = fwht1D input
      n = totalElements input
      scaleFactor = fromRational (1 / fromIntegral n) :: k
  in mapMV (scaleMV scaleFactor) fwd
{-# INLINE ifwht1D #-}

-- | 2D Fast Walsh-Hadamard Transform across an \(H \times W\) multivector grid where \(H, W\) are powers of 2.
fwht2D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
       => MVArray (Signature p q r) k -> MVArray (Signature p q r) k
fwht2D input = case shapeND input of
  [h, w]
    | not (isPowerOfTwo h && isPowerOfTwo w) -> error "Clifford.Array.Fourier.fwht2D: dimensions must be powers of 2"
    | otherwise ->
        let rowTransformed = [ fwht1D (fromListND [w] [ indexND input [r, c] | c <- [0 .. w - 1] ]) | r <- [0 .. h - 1] ]
            intermediateGrid = fromListND [h, w] [ indexND (rowTransformed !! r) [c] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
            colTransformed = [ fwht1D (fromListND [h] [ indexND intermediateGrid [r, c] | r <- [0 .. h - 1] ]) | c <- [0 .. w - 1] ]
        in fromListND [h, w] [ indexND (colTransformed !! c) [r] | r <- [0 .. h - 1], c <- [0 .. w - 1] ]
  _ -> error "Clifford.Array.Fourier.fwht2D: input must be 2-dimensional"
{-# INLINE fwht2D #-}

-- | 2D Inverse Fast Walsh-Hadamard Transform across an \(H \times W\) multivector grid.
ifwht2D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
        => MVArray (Signature p q r) k -> MVArray (Signature p q r) k
ifwht2D input =
  let fwd = fwht2D input
  in case shapeND input of
       [h, w] ->
         let scaleFactor = fromRational (1 / fromIntegral (h * w)) :: k
         in mapMV (scaleMV scaleFactor) fwd
       _ -> error "Clifford.Array.Fourier.ifwht2D: input must be 2-dimensional"
{-# INLINE ifwht2D #-}

-- | Fast Dyadic (XOR) Convolution of two 1D multivector array fields via FWHT.
--
-- \((f \star_{\text{xor}} g)(k) = \sum_{i \oplus j = k} f(i) \star g(j)\)
dyadicConvolve1D :: forall p q r k. (KnownSignature p q r, UnboxedField k)
                 => (DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k -> DenseMV (Signature p q r) k)
                 -- ^ Multivector product operation (e.g. geometric @(<*>)@ or wedge @(/\\)@)
                 -> MVArray (Signature p q r) k
                 -> MVArray (Signature p q r) k
                 -> MVArray (Signature p q r) k
dyadicConvolve1D prodOp a b
  | shapeND a /= shapeND b = error "Clifford.Array.Fourier.dyadicConvolve1D: shape mismatch"
  | otherwise =
      let fA = fwht1D a
          fB = fwht1D b
          fProd = zipWithMV prodOp fA fB
      in ifwht1D fProd
{-# INLINE dyadicConvolve1D #-}

--------------------------------------------------------------------------------
-- Direct Discrete Clifford Fourier Transform (DCFT Reference)
--------------------------------------------------------------------------------

-- | Multidimensional Clifford Fourier Transform kernel specified by bivectors/pseudoscalars.
data CFTKernel sig k = CFTKernel
  { cftBivectors  :: ![DenseMV sig k]   -- ^ Bivectors / Pseudoscalars satisfying \(I_j^2 = -1\)
  , cftDimensions :: ![Int]             -- ^ Grid dimensions @[N_0, N_1, ...]@
  }

instance NFData (CFTKernel sig k) where
  rnf (CFTKernel bv d) = bv `seq` d `seq` ()

-- | Construct and validate a CFT kernel.
mkCFTKernel :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
            => [DenseMV (Signature p q r) k] -> [Int] -> Maybe (CFTKernel (Signature p q r) k)
mkCFTKernel bivecs dims
  | length bivecs /= length dims = Nothing
  | not (all isValidI bivecs)    = Nothing
  | otherwise                    = Just (CFTKernel bivecs dims)
  where
    isValidI iVal =
      let iSq = iVal <*> iVal
          expected = scalarMV (-1)
      in iSq ==~ expected

-- | Construct default 2D Euclidean CFT kernel for \(Cl(2,0,0)\) with \(I = e_{12}\).
defaultCFT2D :: forall k. (Floating k, UnboxedField k)
             => Int -> Int -> CFTKernel Cl200 k
defaultCFT2D h w =
  let e12_blade = bladeMV (bladeFromIndices [1, 2]) (fromMetricScale 1)
  in CFTKernel [e12_blade, e12_blade] [h, w]

-- | Construct default 3D Euclidean CFT kernel for \(Cl(3,0,0)\) with central pseudoscalar \(I = e_{123}\).
defaultCFT3D :: forall k. (Floating k, UnboxedField k)
             => Int -> Int -> Int -> CFTKernel Cl300 k
defaultCFT3D d h w =
  let e123_blade = bladeMV (bladeFromIndices [1, 2, 3]) (fromMetricScale 1)
  in CFTKernel [e123_blade, e123_blade, e123_blade] [d, h, w]

-- | Forward Discrete Clifford Fourier Transform (DCFT):
--
-- \(F(\mathbf{u}) = \sum_{\mathbf{x}} f(\mathbf{x}) \prod_{j} \exp\left(-I_j \, 2\pi \frac{x_j u_j}{N_j}\right)\)
cftForward :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
           => CFTKernel (Signature p q r) k
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
cftForward (CFTKernel bivecs dims) input
  | shapeND input /= dims = error "Clifford.Array.Fourier.cftForward: array dimensions do not match kernel"
  | otherwise =
      let strides = mvaStrides input
          totCells = totalElements input
          bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          buf = mvaBuffer input

          getInMV cellIdx =
            let cellOffset = cellIdx * bladesPerCell
            in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt buf (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

          calcPhase signMult uCoords xCoords =
            let phaseTerms = zipWith4 (\iVal dim uVal xVal ->
                  let angle = signMult * 2 * pi * fromIntegral (uVal * xVal) / fromIntegral dim
                      cosTerm = scalarMV (cos angle)
                      sinTerm = scaleMV (sin angle) iVal
                  in cosTerm `addMV` sinTerm) bivecs dims uCoords xCoords
            in foldl (<*>) (scalarMV 1) phaseTerms

          outBuf = runSTUArray $ do
            marr <- newArray (0, totCells * bladesPerCell - 1) (fromMetricScale 0)
            forM_ [0 .. totCells - 1] $ \uIdx -> do
              let uCoords = flatToCoords strides uIdx
                  accumX xIdx acc =
                    let xCoords = flatToCoords strides xIdx
                        f_x = getInMV xIdx
                        phase = calcPhase (-1) uCoords xCoords
                    in acc `addMV` (f_x <*> phase)
                  resMV = foldr accumX zeroMV [0 .. totCells - 1]
                  DenseMV resArr = resMV
                  cellOffset = uIdx * bladesPerCell
              forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
                writeArray marr (cellOffset + bIdx) (unsafeAt resArr bIdx)
            return marr
      in MVArray dims strides outBuf
{-# INLINE cftForward #-}

-- | Inverse Discrete Clifford Fourier Transform:
--
-- \(f(\mathbf{x}) = \frac{1}{\prod N_j} \sum_{\mathbf{u}} F(\mathbf{u}) \prod_{j} \exp\left(+I_j \, 2\pi \frac{x_j u_j}{N_j}\right)\)
cftInverse :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
           => CFTKernel (Signature p q r) k
           -> MVArray (Signature p q r) k
           -> MVArray (Signature p q r) k
cftInverse (CFTKernel bivecs dims) input
  | shapeND input /= dims = error "Clifford.Array.Fourier.cftInverse: array dimensions do not match kernel"
  | otherwise =
      let strides = mvaStrides input
          totCells = totalElements input
          bladesPerCell = 1 `shiftL` signatureDim (Proxy :: Proxy (Signature p q r))
          buf = mvaBuffer input
          normFactor = fromRational (1 / fromIntegral totCells)

          getInMV cellIdx =
            let cellOffset = cellIdx * bladesPerCell
            in DenseMV (listArray (0, bladesPerCell - 1) [unsafeAt buf (cellOffset + i) | i <- [0 .. bladesPerCell - 1]])

          calcPhase signMult xCoords uCoords =
            let phaseTerms = zipWith4 (\iVal dim xVal uVal ->
                  let angle = signMult * 2 * pi * fromIntegral (xVal * uVal) / fromIntegral dim
                      cosTerm = scalarMV (cos angle)
                      sinTerm = scaleMV (sin angle) iVal
                  in cosTerm `addMV` sinTerm) bivecs dims xCoords uCoords
            in foldl (<*>) (scalarMV 1) phaseTerms

          outBuf = runSTUArray $ do
            marr <- newArray (0, totCells * bladesPerCell - 1) (fromMetricScale 0)
            forM_ [0 .. totCells - 1] $ \xIdx -> do
              let xCoords = flatToCoords strides xIdx
                  accumU uIdx acc =
                    let uCoords = flatToCoords strides uIdx
                        f_u = getInMV uIdx
                        phase = calcPhase 1 xCoords uCoords
                    in acc `addMV` (f_u <*> phase)
                  resMV = foldr accumU zeroMV [0 .. totCells - 1]
                  scaledMV = scaleMV normFactor resMV
                  DenseMV resArr = scaledMV
                  cellOffset = xIdx * bladesPerCell
              forM_ [0 .. bladesPerCell - 1] $ \bIdx -> do
                writeArray marr (cellOffset + bIdx) (unsafeAt resArr bIdx)
            return marr
      in MVArray dims strides outBuf
{-# INLINE cftInverse #-}

-- | Frequency-domain Clifford spectral filtering:
--
-- \(f_{\text{filtered}} = \text{CFT}^{-1}(H \cdot \text{CFT}(f))\)
spectralFilter :: forall p q r k. (KnownSignature p q r, Floating k, UnboxedField k)
               => CFTKernel (Signature p q r) k
               -> (MVArray (Signature p q r) k -> MVArray (Signature p q r) k)
               -> MVArray (Signature p q r) k
               -> MVArray (Signature p q r) k
spectralFilter kernel transferFn input =
  let freqDomain = cftForward kernel input
      filteredFreq = transferFn freqDomain
  in cftInverse kernel filteredFreq
