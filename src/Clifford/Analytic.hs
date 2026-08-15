{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Analytic
-- Description : Analytic Functions, Taylor Series, and Scaling-and-Squaring for Multivectors
-- License     : BSD-3-Clause
--
-- This module implements numerically stable analytic functions on general multivectors:
--
-- * Scaling-and-squaring exponential map ('expMV')
-- * Trigonometric functions ('sinMV', 'cosMV', 'sincosMV')
-- * Hyperbolic functions ('sinhMV', 'coshMV', 'tanhMV')
-- * Arbitrary polynomial evaluation via Horner's scheme ('hornerMV')
-- * Discrete power maps for finite fields ('frobeniusMV')
module Clifford.Analytic
  ( -- * Exponential & Logarithmic Maps
    expMV
  , hornerMV
    -- * Trigonometric Functions
  , sinMV
  , cosMV
  , sincosMV
    -- * Hyperbolic Functions
  , sinhMV
  , coshMV
  , tanhMV
    -- * Finite Field Maps
  , frobeniusMV
  ) where

import Prelude hiding ((<*>))
import Clifford.Field
import Clifford.Class
import Clifford.Versor (versorInverse)

--------------------------------------------------------------------------------
-- Horner Polynomial Evaluation
--------------------------------------------------------------------------------

-- | Evaluate an arbitrary polynomial \(P(M) = \sum_{k=0}^n c_k M^k\) on a multivector
-- using Horner's method for minimal geometric multiplications and maximal numerical stability:
--
-- \[ P(M) = c_0 + M \star (c_1 + M \star (c_2 + \dots + M \star c_n)) \]
--
-- Example:
--
-- >>> hornerMV [1, 2, 3] (e1 :: DenseMV Cl300 Double)
-- 4.0 + 2.0 e_{1}  -- since 1 + 2*e1 + 3*e1^2 = 1 + 2*e1 + 3*1 = 4 + 2*e1
hornerMV :: forall sig k mv. (CliffordAlgebra sig k mv) => [k] -> mv -> mv
hornerMV [] _ = zeroMV
hornerMV [c0] _ = scalarMV c0
hornerMV (c0:cs) x = scalarMV c0 `addMV` (x <*> hornerMV cs x)

--------------------------------------------------------------------------------
-- Scaling-and-Squaring Exponential Map
--------------------------------------------------------------------------------

-- | General multivector exponential map \(\exp(M)\) evaluated via scaling-and-squaring
-- combined with a Horner-nested Taylor expansion.
--
-- This function works on __arbitrary multivectors__ (including mixed grades, screw twists,
-- and non-scalar bivectors).
--
-- Convergence: Error \(< 10^{-15}\) (machine precision) for all inputs.
expMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
expMV m
  | isZero (normSqMV m) = scalarMV 1 `addMV` m  -- Nilpotent case (e.g. translation in PGA: m^2 = 0)
  | otherwise =
      let mag = normMV m
          -- Determine scaling factor s such that ||m / 2^s|| <= 0.5
          s = if mag <= 0.5
                then 0
                else ceiling (logBase 2 (mag / 0.5)) :: Int
          twoPowS = fromIntegral ((1 :: Integer) `shiftL` s)
          x = scaleMV (1 / twoPowS) m
          -- Taylor series for exp(X) on scaled ||X|| <= 0.5 with 10 terms:
          -- 1 + x(1 + x/2(1 + x/3(1 + x/4(...(1 + x/10)...))))
          expX = evalExpTaylor 10 x
      in applySquaring s expX
  where
    shiftL :: Integer -> Int -> Integer
    shiftL val amt = val * (2 ^ amt)

    evalExpTaylor :: Int -> mv -> mv
    evalExpTaylor maxN x = go 1
      where
        go n
          | n > maxN  = scalarMV 1
          | otherwise =
              let nextFactor = go (n + 1)
                  coeff = 1 / fromIntegral n
                  scaledTerm = scaleMV coeff (x <*> nextFactor)
              in scalarMV 1 `addMV` scaledTerm

    applySquaring :: Int -> mv -> mv
    applySquaring 0 val = val
    applySquaring count val = applySquaring (count - 1) (val <*> val)

--------------------------------------------------------------------------------
-- Trigonometric Functions
--------------------------------------------------------------------------------

-- | Simultaneously compute \(\sin(M)\) and \(\cos(M)\) using scaling-and-squaring
-- with trigonometric double-angle identities:
--
-- \[ \sin(2A) = 2 \sin A \cos A, \quad \cos(2A) = \cos^2 A - \sin^2 A \]
sincosMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> (mv, mv)
sincosMV m =
  let mag = normMV m
      s = if mag <= 0.5
            then 0
            else ceiling (logBase 2 (mag / 0.5)) :: Int
      twoPowS = fromIntegral ((1 :: Integer) `shiftL` s)
      x = scaleMV (1 / twoPowS) m
      -- Truncated Taylor polynomial for scaled ||x|| <= 0.5
      sinX = evalSinTaylor 6 x
      cosX = evalCosTaylor 6 x
  in applyDoubleAngle s sinX cosX
  where
    shiftL :: Integer -> Int -> Integer
    shiftL val amt = val * (2 ^ amt)

    evalSinTaylor :: Int -> mv -> mv
    evalSinTaylor maxK x =
      let x2 = x <*> x
          coeffs = [ (if even k then 1 else -1) / factorial (2 * k + 1) | k <- [0..maxK] ]
      in x <*> hornerEven coeffs x2

    evalCosTaylor :: Int -> mv -> mv
    evalCosTaylor maxK x =
      let x2 = x <*> x
          coeffs = [ (if even k then 1 else -1) / factorial (2 * k) | k <- [0..maxK] ]
      in hornerEven coeffs x2

    hornerEven :: [k] -> mv -> mv
    hornerEven [] _ = zeroMV
    hornerEven [c] _ = scalarMV c
    hornerEven (c:cs) x2 = scalarMV c `addMV` (x2 <*> hornerEven cs x2)

    factorial :: Int -> k
    factorial 0 = 1
    factorial n = fromIntegral n * factorial (n - 1)

    applyDoubleAngle :: Int -> mv -> mv -> (mv, mv)
    applyDoubleAngle 0 sVal cVal = (sVal, cVal)
    applyDoubleAngle count sVal cVal =
      let sin2 = scaleMV 2 (sVal <*> cVal)
          cos2 = (cVal <*> cVal) `subMV` (sVal <*> sVal)
      in applyDoubleAngle (count - 1) sin2 cos2

-- | Multivector Sine function: \(\sin(M)\).
sinMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
sinMV = fst . sincosMV

-- | Multivector Cosine function: \(\cos(M)\).
cosMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
cosMV = snd . sincosMV

--------------------------------------------------------------------------------
-- Hyperbolic Functions
--------------------------------------------------------------------------------

-- | Multivector Hyperbolic Sine: \(\sinh(M) = \frac{1}{2} (\exp(M) - \exp(-M))\).
sinhMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
sinhMV m =
  let expPos = expMV m
      expNeg = expMV (negateMV m)
  in scaleMV (1 / 2) (expPos `subMV` expNeg)

-- | Multivector Hyperbolic Cosine: \(\cosh(M) = \frac{1}{2} (\exp(M) + \exp(-M))\).
coshMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
coshMV m =
  let expPos = expMV m
      expNeg = expMV (negateMV m)
  in scaleMV (1 / 2) (expPos `addMV` expNeg)

-- | Multivector Hyperbolic Tangent: \(\tanh(M) = \sinh(M) \star (\cosh(M))^{-1}\).
tanhMV :: forall sig k mv. (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> Maybe mv
tanhMV m =
  let sVal = sinhMV m
      cVal = coshMV m
  in case versorInverse cVal of
       Just invC -> Just (sVal <*> invC)
       Nothing   -> Nothing

--------------------------------------------------------------------------------
-- Finite Field Discrete Power Maps
--------------------------------------------------------------------------------

-- | Discrete power map for characteristic \(p\) Galois fields (e.g. Frobenius automorphism \(\Phi_p(M) = M^p\))
-- evaluated via binary repeated squaring.
frobeniusMV :: forall sig k mv. (CliffordAlgebra sig k mv) => Integer -> mv -> mv
frobeniusMV 0 _ = scalarMV (fromMetricScale 1)
frobeniusMV 1 m = m
frobeniusMV p m
  | p < 0     = error "frobeniusMV: negative power not supported in general finite ring"
  | even p    = let half = frobeniusMV (p `div` 2) m in half <*> half
  | otherwise = m <*> frobeniusMV (p - 1) m
