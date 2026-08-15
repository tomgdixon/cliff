{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Blade
-- Description : Basis blade representations, bitmask combinatorics, and canonical sign evaluation
-- License     : BSD-3-Clause
--
-- This module implements basis blades using 64-bit bitmasks, supporting spaces
-- up to \(n = 64\) dimensions (\(2^{64}\) basis blades).
--
-- In this representation, basis vector \(e_i\) corresponds to bit index \((i-1)\).
-- For example:
--
-- * \(1\) (scalar) \(\rightarrow\) @Blade 0@
-- * \(e_1\) \(\rightarrow\) @Blade 1@ (@0b001@)
-- * \(e_2\) \(\rightarrow\) @Blade 2@ (@0b010@)
-- * \(e_{12} = e_1 \wedge e_2\) \(\rightarrow\) @Blade 3@ (@0b011@)
-- * \(e_{123} = e_1 \wedge e_2 \wedge e_3\) \(\rightarrow\) @Blade 7@ (@0b111@)
module Clifford.Blade
  ( -- * Basis Blade Representation
    Blade(..)
  , scalarBlade
    -- * Blade Constructors
  , basisBlade
  , pseudoscalarBlade
  , bladeFromIndices
  , bladeIndices
  , bladeGrade
    -- * Blade Arithmetic & Sign Rules
  , canonicalSwapParity
  , metricFactor
  , bladeMul
  , bladeWedge
  , bladeDotLeft
  , bladeDotRight
  , bladeScalarProd
  , poincareDual
    -- * Blade Enumerations
  , allBlades
  , bladesOfGrade
  ) where

import Control.DeepSeq (NFData(..))
import Data.Bits
import Data.Word (Word64)
import Clifford.Signature

-- | Basis blade represented as a 64-bit integer bitmask.
newtype Blade = Blade { unBlade :: Word64 }
  deriving (Eq, Ord, NFData)

instance Show Blade where
  show (Blade 0) = "1"
  show b = "e" ++ concatMap show (bladeIndices b)

-- | The scalar unit blade \(1\) (empty bitmask).
scalarBlade :: Blade
scalarBlade = Blade 0
{-# INLINE scalarBlade #-}

-- | Construct a 1-vector basis blade \(e_i\) (1-indexed, \(1 \le i \le 64\)).
basisBlade :: Int -> Blade
basisBlade i
  | i <= 0 || i > 64 = error $ "Clifford.Blade.basisBlade: index out of range [1..64]: " ++ show i
  | otherwise        = Blade (1 `shiftL` (i - 1))
{-# INLINE basisBlade #-}

-- | Construct the unit pseudoscalar \(I = e_1 e_2 \dots e_n\) for an \(n\)-dimensional algebra.
pseudoscalarBlade :: Int -> Blade
pseudoscalarBlade n
  | n < 0 || n > 64 = error $ "Clifford.Blade.pseudoscalarBlade: dimension out of range [0..64]: " ++ show n
  | n == 64         = Blade maxBound
  | otherwise       = Blade ((1 `shiftL` n) - 1)
{-# INLINE pseudoscalarBlade #-}

-- | Construct a basis blade from an ordered list of 1-based basis indices:
--
-- >>> bladeFromIndices [1, 2] == Blade 3
-- True
bladeFromIndices :: [Int] -> Blade
bladeFromIndices = foldr (\i (Blade acc) -> Blade (acc .|. (1 `shiftL` (i - 1)))) (Blade 0)
{-# INLINE bladeFromIndices #-}

-- | Deconstruct a basis blade into its sorted 1-based basis vector indices.
bladeIndices :: Blade -> [Int]
bladeIndices (Blade w) = [i + 1 | i <- [0 .. 63], testBit w i]
{-# INLINE bladeIndices #-}

-- | Grade \(k\) of a basis blade, equal to its Hamming weight (population count).
bladeGrade :: Blade -> Int
bladeGrade (Blade w) = popCount w
{-# INLINE bladeGrade #-}

-- | Parity sign factor \((-1)^{\text{swaps}}\) required to reorder basis vectors into canonical order.
canonicalSwapParity :: Blade -> Blade -> Int
canonicalSwapParity (Blade w1) (Blade w2) =
  let swaps = evalSwaps (w1 `shiftR` 1) w2 0
  in if testBit swaps 0 then -1 else 1
  where
    evalSwaps 0 !_ !acc = acc
    evalSwaps !s1 !s2 !acc =
      let count = popCount (s2 .&. s1)
      in evalSwaps (s1 `shiftR` 1) s2 (acc + count)
{-# INLINE canonicalSwapParity #-}

-- | Metric contraction factor for shared basis vectors between @b1@ and @b2@.
-- Returns \(\prod_{i \in b_1 \cap b_2} e_i^2 \in \{-1, 0, +1\}\).
metricFactor :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Blade -> Blade -> Int
metricFactor p (Blade w1) (Blade w2) =
  let shared = w1 .&. w2
  in if shared == 0
       then 1
       else evalShared shared 1 1
  where
    evalShared 0 !_ !acc = acc
    evalShared !s !idx !acc =
      if testBit s 0
        then
          let sgn = metricSignOf p idx
          in if sgn == 0
               then 0
               else evalShared (s `shiftR` 1) (idx + 1) (acc * sgn)
        else evalShared (s `shiftR` 1) (idx + 1) acc
{-# INLINE metricFactor #-}

-- | Full geometric product of two basis blades:
--
-- @b1 * b2 = (signFactor, resultingBlade)@ where @signFactor@ \(\in \{-1, 0, +1\}\).
bladeMul :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Blade -> Blade -> (Int, Blade)
bladeMul p b1 b2 =
  let !mFactor = metricFactor p b1 b2
  in if mFactor == 0
       then (0, scalarBlade)
       else
         let !pSign = canonicalSwapParity b1 b2
             !resBlade = Blade (unBlade b1 `xor` unBlade b2)
         in (pSign * mFactor, resBlade)
{-# INLINE bladeMul #-}

-- | Exterior / Wedge product of two basis blades:
-- Returns @Just (sign, b1 ^ b2)@ if blades are linearly independent (disjoint basis sets),
-- or @Nothing@ if they share any basis vectors.
bladeWedge :: Blade -> Blade -> Maybe (Int, Blade)
bladeWedge b1 b2
  | unBlade b1 .&. unBlade b2 /= 0 = Nothing
  | otherwise =
      let !pSign = canonicalSwapParity b1 b2
      in Just (pSign, Blade (unBlade b1 .|. unBlade b2))
{-# INLINE bladeWedge #-}

-- | Left contraction (Dorst/Lounesto) of two basis blades:
-- \(A \rfloor B = \langle A B \rangle_{\text{grade}(B) - \text{grade}(A)}\).
--
-- Non-zero only if all basis vectors of @b1@ are contained in @b2@.
bladeDotLeft :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Blade -> Blade -> Maybe (Int, Blade)
bladeDotLeft p b1 b2
  | (unBlade b1 .&. unBlade b2) /= unBlade b1 = Nothing
  | otherwise =
      let (s, res) = bladeMul p b1 b2
      in if s == 0 then Nothing else Just (s, res)
{-# INLINE bladeDotLeft #-}

-- | Right contraction of two basis blades:
-- \(A \lfloor B = \langle A B \rangle_{\text{grade}(A) - \text{grade}(B)}\).
bladeDotRight :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Blade -> Blade -> Maybe (Int, Blade)
bladeDotRight p b1 b2
  | (unBlade b1 .&. unBlade b2) /= unBlade b2 = Nothing
  | otherwise =
      let (s, res) = bladeMul p b1 b2
      in if s == 0 then Nothing else Just (s, res)
{-# INLINE bladeDotRight #-}

-- | Scalar inner product of two basis blades:
-- \(\langle b_1 b_2 \rangle_0\). Non-zero only if @b1 == b2@.
bladeScalarProd :: forall p q r proxy. KnownSignature p q r => proxy (Signature p q r) -> Blade -> Blade -> Int
bladeScalarProd p b1 b2
  | b1 /= b2  = 0
  | otherwise = fst (bladeMul p b1 b1)
{-# INLINE bladeScalarProd #-}

-- | Poincaré complement dual \(J(A)\) for dimension \(n\):
-- Grade-reversing complement bijection defined for all metric signatures (including degenerate PGA \(r > 0\)):
--
-- \(J(e_K) = \text{sign}(e_K \wedge e_{\bar{K}} = I) \, e_{\bar{K}}\)
poincareDual :: Int -> Blade -> (Int, Blade)
poincareDual n (Blade w) =
  let fullMask = if n >= 64 then maxBound else (1 `shiftL` n) - 1
      compMask = fullMask `xor` w
      compBlade = Blade compMask
      pSign = canonicalSwapParity (Blade w) compBlade
  in (pSign, compBlade)
{-# INLINE poincareDual #-}

-- | Generate all \(2^n\) basis blades in canonical bitmask order.
allBlades :: Int -> [Blade]
allBlades n
  | n <= 0    = [Blade 0]
  | n > 64    = error "Clifford.Blade.allBlades: n cannot exceed 64"
  | n == 64   = map Blade [0 .. maxBound]
  | otherwise = [Blade w | w <- [0 .. (1 `shiftL` n) - 1]]

-- | Generate all basis blades of a specific grade \(k\) in dimension \(n\).
bladesOfGrade :: Int -> Int -> [Blade]
bladesOfGrade n k = filter (\b -> bladeGrade b == k) (allBlades n)
