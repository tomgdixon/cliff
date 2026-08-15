{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Clifford.Field
-- Description : Unified scalar field abstractions for continuous, rational, and finite Galois fields
-- License     : BSD-3-Clause
--
-- This module defines the 'Field' typeclass, generalizing numeric scalar arithmetic
-- across standard continuous numbers ('Double', 'Float', 'Complex Double'), exact rationals ('Rational'),
-- binary fields ('GF2'), modular prime fields ('GFp'), and Galois extension fields ('GF2m').
module Clifford.Field
  ( -- * Scalar Field Typeclass
    Field(..)
  , UnboxedField
    -- * Finite & Galois Fields
  , GF2(..)
  , GFp(..)
  , GF2m(..)
    -- * Galois Field Arithmetic Utilities
  , extGCD
  , modPow
  , gf2mMult
  , gf2mInv
  ) where

import Control.DeepSeq (NFData(..))
import Control.Monad.ST (ST)
import Data.Array.Base (IArray(..), MArray(..), UArray(..), STUArray(..))
import Data.Bits
import Data.Complex (Complex(..), conjugate)
import Data.Proxy (Proxy(..))
import Data.Ratio (Ratio, (%), numerator, denominator, Rational)
import Data.Word (Word8, Word64)
import GHC.TypeLits
import Unsafe.Coerce (unsafeCoerce)

-- | Generalized scalar field abstraction.
class (Eq k, Fractional k) => Field k where
  -- | Inject integer metric scaling factor (\(-1, 0, +1\)) into the scalar field.
  fromMetricScale :: Int -> k
  fromMetricScale n = fromIntegral n
  {-# INLINE fromMetricScale #-}

  -- | Multiplicative inverse (returns @Nothing@ for zero).
  fieldInverse :: k -> Maybe k
  fieldInverse x
    | isZero x  = Nothing
    | otherwise = Just (recip x)
  {-# INLINE fieldInverse #-}

  -- | Check if scalar is zero.
  isZero :: k -> Bool
  isZero x = x == 0
  {-# INLINE isZero #-}

  -- | Field involution / complex conjugate (identity for real and finite fields).
  fieldConj :: k -> k
  fieldConj x = x
  {-# INLINE fieldConj #-}

  -- | Squared field norm \(x \cdot \overline{x}\).
  fieldNormSq :: k -> k
  fieldNormSq x = x * fieldConj x
  {-# INLINE fieldNormSq #-}

-- | Class for fields that can be stored in unboxed arrays.
class (Field k, IArray UArray k, forall s. MArray (STUArray s) k (ST s)) => UnboxedField k
instance (Field k, IArray UArray k, forall s. MArray (STUArray s) k (ST s)) => UnboxedField k

-- | Instance for 'Double'.
instance Field Double where
  isZero x = abs x < 1e-9
  {-# INLINE isZero #-}

-- | Instance for 'Float'.
instance Field Float where
  isZero x = abs x < 1e-5
  {-# INLINE isZero #-}

-- | Instance for exact 'Rational'.
instance Field Rational

-- | Instance for 'Complex Double'.
instance Field (Complex Double) where
  fieldConj = conjugate
  isZero (r :+ i) = (r*r + i*i) < 1e-18
  {-# INLINE fieldConj #-}
  {-# INLINE isZero #-}

--------------------------------------------------------------------------------
-- Binary Galois Field GF(2)
--------------------------------------------------------------------------------

-- | Binary Field \(\mathbb{F}_2 = \{0, 1\}\).
--
-- Addition is bitwise XOR (\(\oplus\)), multiplication is bitwise AND (\(\land\)).
-- In \(\mathbb{F}_2\), \(1 + 1 = 0\) and \(-1 \equiv +1\).
newtype GF2 = GF2 { unGF2 :: Word8 }
  deriving (Eq, Ord, NFData)

instance Show GF2 where
  show (GF2 0) = "0"
  show (GF2 1) = "1"
  show (GF2 w) = show (w .&. 1)

instance Num GF2 where
  (GF2 a) + (GF2 b) = GF2 ((a `xor` b) .&. 1)
  (GF2 a) - (GF2 b) = GF2 ((a `xor` b) .&. 1)
  (GF2 a) * (GF2 b) = GF2 ((a .&. b) .&. 1)
  negate x          = x
  abs x             = x
  signum x          = x
  fromInteger n     = GF2 (fromIntegral (n .&. 1))
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}
  {-# INLINE negate #-}
  {-# INLINE fromInteger #-}

instance Fractional GF2 where
  (GF2 a) / (GF2 b)
    | b .&. 1 == 0 = error "GF2: division by zero"
    | otherwise    = GF2 (a .&. 1)
  recip (GF2 0) = error "GF2: division by zero"
  recip (GF2 1) = GF2 1
  recip (GF2 w) = recip (GF2 (w .&. 1))
  fromRational r = fromInteger (numerator r) / fromInteger (denominator r)
  {-# INLINE (/) #-}
  {-# INLINE recip #-}

instance Field GF2 where
  fromMetricScale s = if s == 0 then GF2 0 else GF2 1
  fieldInverse (GF2 0) = Nothing
  fieldInverse (GF2 _) = Just (GF2 1)
  isZero (GF2 w) = (w .&. 1) == 0
  {-# INLINE fromMetricScale #-}
  {-# INLINE fieldInverse #-}
  {-# INLINE isZero #-}

-- Unboxed UArray instances for GF2
instance IArray UArray GF2 where
  {-# INLINE bounds #-}
  bounds (UArray l u _ _) = (l, u)
  {-# INLINE numElements #-}
  numElements (UArray _ _ n _) = n
  {-# INLINE unsafeArray #-}
  unsafeArray lu ies = unsafeCoerce ((unsafeArray (unsafeCoerce lu :: (Int, Int)) (unsafeCoerce ies :: [(Int, Word8)])) :: UArray Int Word8)
  {-# INLINE unsafeAt #-}
  unsafeAt arr i = GF2 (unsafeAt (unsafeCoerce arr :: UArray Int Word8) i)
  {-# INLINE unsafeReplace #-}
  unsafeReplace arr ies = unsafeCoerce (unsafeReplace (unsafeCoerce arr :: UArray Int Word8) (unsafeCoerce ies))
  {-# INLINE unsafeAccum #-}
  unsafeAccum f arr ies = unsafeCoerce (unsafeAccum (\a b -> unGF2 (f (GF2 a) (unsafeCoerce b))) (unsafeCoerce arr :: UArray Int Word8) ies)
  {-# INLINE unsafeAccumArray #-}
  unsafeAccumArray f initVal lu ies = unsafeCoerce ((unsafeAccumArray (\a b -> unGF2 (f (GF2 a) (unsafeCoerce b))) (unGF2 initVal) (unsafeCoerce lu :: (Int, Int)) ies) :: UArray Int Word8)

instance MArray (STUArray s) GF2 (ST s) where
  {-# INLINE getBounds #-}
  getBounds (STUArray l u _ _) = return (l, u)
  {-# INLINE getNumElements #-}
  getNumElements (STUArray _ _ n _) = return n
  {-# INLINE newArray #-}
  newArray (l, u) (GF2 initVal) = do
    STUArray _ _ n mba <- (newArray (unsafeCoerce (l, u) :: (Int, Int)) initVal :: ST s (STUArray s Int Word8))
    return (STUArray l u n mba)
  {-# INLINE newArray_ #-}
  newArray_ (l, u) = do
    STUArray _ _ n mba <- (newArray_ (unsafeCoerce (l, u) :: (Int, Int)) :: ST s (STUArray s Int Word8))
    return (STUArray l u n mba)
  {-# INLINE unsafeRead #-}
  unsafeRead (STUArray _ _ _ mba) i = do
    w <- unsafeRead (STUArray 0 0 0 mba :: STUArray s Int Word8) i
    return (GF2 w)
  {-# INLINE unsafeWrite #-}
  unsafeWrite (STUArray _ _ _ mba) i (GF2 w) =
    unsafeWrite (STUArray 0 0 0 mba :: STUArray s Int Word8) i w

--------------------------------------------------------------------------------
-- Prime Modular Field GF(p)
--------------------------------------------------------------------------------

-- | Prime Galois Field \(\mathbb{F}_p\) for type-level prime @p@.
newtype GFp (p :: Nat) = GFp { unGFp :: Word64 }
  deriving (Eq, Ord, NFData)

instance KnownNat p => Show (GFp p) where
  show (GFp a) = show a ++ " (mod " ++ show (natVal (Proxy :: Proxy p)) ++ ")"

instance KnownNat p => Num (GFp p) where
  (GFp a) + (GFp b) = GFp ((a + b) `mod` pVal) where pVal = fromIntegral (natVal (Proxy :: Proxy p))
  (GFp a) - (GFp b) = GFp ((a - b + pVal) `mod` pVal) where pVal = fromIntegral (natVal (Proxy :: Proxy p))
  (GFp a) * (GFp b) = GFp ((a * b) `mod` pVal) where pVal = fromIntegral (natVal (Proxy :: Proxy p))
  negate (GFp 0)    = GFp 0
  negate (GFp a)    = GFp (fromIntegral (natVal (Proxy :: Proxy p)) - (a `mod` fromIntegral (natVal (Proxy :: Proxy p))))
  abs x             = x
  signum x          = x
  fromInteger n     = GFp (fromIntegral ((n `mod` pVal + pVal) `mod` pVal)) where pVal = natVal (Proxy :: Proxy p)
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}
  {-# INLINE negate #-}
  {-# INLINE fromInteger #-}

instance KnownNat p => Fractional (GFp p) where
  recip (GFp 0) = error "GFp: division by zero"
  recip (GFp a) = GFp (fromIntegral (modInverse (fromIntegral a) (natVal (Proxy :: Proxy p))))
  (GFp a) / (GFp b)
    | b == 0    = error "GFp: division by zero"
    | otherwise = (GFp a) * recip (GFp b)
  fromRational r = fromInteger (numerator r) / fromInteger (denominator r)
  {-# INLINE recip #-}
  {-# INLINE (/) #-}

instance KnownNat p => Field (GFp p) where
  fieldInverse (GFp 0) = Nothing
  fieldInverse x       = Just (recip x)
  isZero (GFp a)       = (a `mod` fromIntegral (natVal (Proxy :: Proxy p))) == 0
  {-# INLINE fieldInverse #-}
  {-# INLINE isZero #-}

instance KnownNat p => IArray UArray (GFp p) where
  {-# INLINE bounds #-}
  bounds (UArray l u _ _) = (l, u)
  {-# INLINE numElements #-}
  numElements (UArray _ _ n _) = n
  {-# INLINE unsafeArray #-}
  unsafeArray lu ies = unsafeCoerce ((unsafeArray (unsafeCoerce lu :: (Int, Int)) (unsafeCoerce ies :: [(Int, Word64)])) :: UArray Int Word64)
  {-# INLINE unsafeAt #-}
  unsafeAt arr i = GFp (unsafeAt (unsafeCoerce arr :: UArray Int Word64) i)
  {-# INLINE unsafeReplace #-}
  unsafeReplace arr ies = unsafeCoerce (unsafeReplace (unsafeCoerce arr :: UArray Int Word64) (unsafeCoerce ies))
  {-# INLINE unsafeAccum #-}
  unsafeAccum f arr ies = unsafeCoerce (unsafeAccum (\a b -> unGFp (f (GFp a) (unsafeCoerce b))) (unsafeCoerce arr :: UArray Int Word64) ies)
  {-# INLINE unsafeAccumArray #-}
  unsafeAccumArray f initVal lu ies = unsafeCoerce ((unsafeAccumArray (\a b -> unGFp (f (GFp a) (unsafeCoerce b))) (unGFp initVal) (unsafeCoerce lu :: (Int, Int)) ies) :: UArray Int Word64)

instance KnownNat p => MArray (STUArray s) (GFp p) (ST s) where
  {-# INLINE getBounds #-}
  getBounds (STUArray l u _ _) = return (l, u)
  {-# INLINE getNumElements #-}
  getNumElements (STUArray _ _ n _) = return n
  {-# INLINE newArray #-}
  newArray (l, u) (GFp initVal) = do
    STUArray _ _ n mba <- (newArray (unsafeCoerce (l, u) :: (Int, Int)) initVal :: ST s (STUArray s Int Word64))
    return (STUArray l u n mba)
  {-# INLINE newArray_ #-}
  newArray_ (l, u) = do
    STUArray _ _ n mba <- (newArray_ (unsafeCoerce (l, u) :: (Int, Int)) :: ST s (STUArray s Int Word64))
    return (STUArray l u n mba)
  {-# INLINE unsafeRead #-}
  unsafeRead (STUArray _ _ _ mba) i = do
    w <- unsafeRead (STUArray 0 0 0 mba :: STUArray s Int Word64) i
    return (GFp w)
  {-# INLINE unsafeWrite #-}
  unsafeWrite (STUArray _ _ _ mba) i (GFp w) =
    unsafeWrite (STUArray 0 0 0 mba :: STUArray s Int Word64) i w

--------------------------------------------------------------------------------
-- Galois Extension Field GF(2^m)
--------------------------------------------------------------------------------

-- | Galois Extension Field \(\mathbb{F}_{2^m}\) parameterized by:
--
-- * @m@: degree of extension (\(1 \le m \le 64\))
-- * @poly@: irreducible reduction polynomial represented as an integer bitmask
--
-- For example, @GF2m 8 0x11B@ represents \(\text{GF}(2^8) = \text{GF}(256)\)
-- with irreducible polynomial \(x^8 + x^4 + x^3 + x + 1\).
newtype GF2m (m :: Nat) (poly :: Nat) = GF2m { unGF2m :: Word64 }
  deriving (Eq, Ord, NFData)

instance KnownNat m => Show (GF2m m poly) where
  show (GF2m w) = "GF(2^" ++ show (natVal (Proxy :: Proxy m)) ++ ":" ++ show w ++ ")"

instance (KnownNat m, KnownNat poly) => Num (GF2m m poly) where
  (GF2m a) + (GF2m b) = GF2m (a `xor` b)
  (GF2m a) - (GF2m b) = GF2m (a `xor` b)
  (GF2m a) * (GF2m b) = GF2m (gf2mMult a b (fromIntegral (natVal (Proxy :: Proxy m))) (fromIntegral (natVal (Proxy :: Proxy poly))))
  negate x            = x
  abs x               = x
  signum x            = x
  fromInteger n       = GF2m (fromIntegral n .&. mask)
    where mask = if m >= 64 then maxBound else (1 `shiftL` m) - 1
          m = fromIntegral (natVal (Proxy :: Proxy m))
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}

instance (KnownNat m, KnownNat poly) => Fractional (GF2m m poly) where
  recip (GF2m 0) = error "GF2m: division by zero"
  recip (GF2m a) = GF2m (gf2mInv a (fromIntegral (natVal (Proxy :: Proxy m))) (fromIntegral (natVal (Proxy :: Proxy poly))))
  (GF2m a) / (GF2m b)
    | b == 0    = error "GF2m: division by zero"
    | otherwise = (GF2m a) * recip (GF2m b)
  fromRational r = fromInteger (numerator r) / fromInteger (denominator r)
  {-# INLINE recip #-}
  {-# INLINE (/) #-}

instance (KnownNat m, KnownNat poly) => Field (GF2m m poly) where
  fromMetricScale s = if s == 0 then GF2m 0 else GF2m 1
  fieldInverse (GF2m 0) = Nothing
  fieldInverse x        = Just (recip x)
  isZero (GF2m w)       = w == 0
  {-# INLINE fromMetricScale #-}
  {-# INLINE fieldInverse #-}
  {-# INLINE isZero #-}

instance (KnownNat m, KnownNat poly) => IArray UArray (GF2m m poly) where
  {-# INLINE bounds #-}
  bounds (UArray l u _ _) = (l, u)
  {-# INLINE numElements #-}
  numElements (UArray _ _ n _) = n
  {-# INLINE unsafeArray #-}
  unsafeArray lu ies = unsafeCoerce ((unsafeArray (unsafeCoerce lu :: (Int, Int)) (unsafeCoerce ies :: [(Int, Word64)])) :: UArray Int Word64)
  {-# INLINE unsafeAt #-}
  unsafeAt arr i = GF2m (unsafeAt (unsafeCoerce arr :: UArray Int Word64) i)
  {-# INLINE unsafeReplace #-}
  unsafeReplace arr ies = unsafeCoerce (unsafeReplace (unsafeCoerce arr :: UArray Int Word64) (unsafeCoerce ies))
  {-# INLINE unsafeAccum #-}
  unsafeAccum f arr ies = unsafeCoerce (unsafeAccum (\a b -> unGF2m (f (GF2m a) (unsafeCoerce b))) (unsafeCoerce arr :: UArray Int Word64) ies)
  {-# INLINE unsafeAccumArray #-}
  unsafeAccumArray f initVal lu ies = unsafeCoerce ((unsafeAccumArray (\a b -> unGF2m (f (GF2m a) (unsafeCoerce b))) (unGF2m initVal) (unsafeCoerce lu :: (Int, Int)) ies) :: UArray Int Word64)

instance (KnownNat m, KnownNat poly) => MArray (STUArray s) (GF2m m poly) (ST s) where
  {-# INLINE getBounds #-}
  getBounds (STUArray l u _ _) = return (l, u)
  {-# INLINE getNumElements #-}
  getNumElements (STUArray _ _ n _) = return n
  {-# INLINE newArray #-}
  newArray (l, u) (GF2m initVal) = do
    STUArray _ _ n mba <- (newArray (unsafeCoerce (l, u) :: (Int, Int)) initVal :: ST s (STUArray s Int Word64))
    return (STUArray l u n mba)
  {-# INLINE newArray_ #-}
  newArray_ (l, u) = do
    STUArray _ _ n mba <- (newArray_ (unsafeCoerce (l, u) :: (Int, Int)) :: ST s (STUArray s Int Word64))
    return (STUArray l u n mba)
  {-# INLINE unsafeRead #-}
  unsafeRead (STUArray _ _ _ mba) i = do
    w <- unsafeRead (STUArray 0 0 0 mba :: STUArray s Int Word64) i
    return (GF2m w)
  {-# INLINE unsafeWrite #-}
  unsafeWrite (STUArray _ _ _ mba) i (GF2m w) =
    unsafeWrite (STUArray 0 0 0 mba :: STUArray s Int Word64) i w

--------------------------------------------------------------------------------
-- Mathematical Utilities
--------------------------------------------------------------------------------

-- | Extended Euclidean Algorithm: computes @(g, x, y)@ such that @a*x + b*y = g = gcd(a, b)@.
extGCD :: Integer -> Integer -> (Integer, Integer, Integer)
extGCD a 0 = (a, 1, 0)
extGCD a b =
  let (q, r) = a `quotRem` b
      (g, s, t) = extGCD b r
  in (g, t, s - q * t)

-- | Compute modular inverse of @a@ modulo @m@.
modInverse :: Integer -> Integer -> Integer
modInverse a m =
  let (g, x, _) = extGCD a m
  in if g /= 1
       then error $ "Clifford.Field.modInverse: " ++ show a ++ " has no inverse modulo " ++ show m
       else (x `mod` m + m) `mod` m

-- | Modular exponentiation: computes \((b^e) \pmod m\).
modPow :: Integer -> Integer -> Integer -> Integer
modPow _ 0 _ = 1
modPow b e m =
  let half = modPow b (e `shiftR` 1) m
      sq = (half * half) `mod` m
  in if testBit e 0
       then (sq * (b `mod` m)) `mod` m
       else sq

-- | Carryless multiplication modulo irreducible polynomial in \(\mathbb{F}_{2^m}\).
gf2mMult :: Word64 -> Word64 -> Int -> Word64 -> Word64
gf2mMult !a !b !m !poly = go b 0 a
  where
    !highBit = 1 `shiftL` (m - 1)
    !mask = if m >= 64 then maxBound else (1 `shiftL` m) - 1
    go 0 !res !_ = res
    go !currB !res !currA =
      let !newRes = if testBit currB 0 then res `xor` currA else res
          !carry = currA .&. highBit
          !shiftedA = (currA `shiftL` 1) .&. mask
          !reducedA = if carry /= 0 then shiftedA `xor` (poly .&. mask) else shiftedA
      in go (currB `shiftR` 1) newRes reducedA
{-# INLINE gf2mMult #-}

-- | Inversion in \(\mathbb{F}_{2^m}\) via Fermat's Little Theorem: \(a^{-1} = a^{2^m - 2}\).
gf2mInv :: Word64 -> Int -> Word64 -> Word64
gf2mInv 0 _ _ = error "Clifford.Field.gf2mInv: division by zero"
gf2mInv a m poly =
  let expVal = (1 `shiftL` m) - 2 :: Integer
      mult x y = gf2mMult x y m poly
      pow _ 0 = 1
      pow base e =
        let half = pow base (e `shiftR` 1)
            sq = mult half half
        in if testBit e 0 then mult sq base else sq
  in pow a expVal
