{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Clifford.Class
-- Description : Foundational algebraic typeclasses, products, involutions, and operator fixities
-- License     : BSD-3-Clause
--
-- This module defines the core typeclasses 'CliffordVectorSpace' and 'CliffordAlgebra',
-- establishing a unified interface for all multivector representations (dense unrolled and sparse).
module Clifford.Class
  ( -- * Core Algebraic Typeclasses
    CliffordVectorSpace(..)
  , CliffordAlgebra(..)
    -- * Operator Precedence & Aliases
  , (==~)
  , contractLeft
  , contractRight
  , scalarProd
  , normMV
  , normalizeMV
  , toLatex
  , evenPart
  , oddPart
  , hestenesDot
  ) where

import Prelude hiding ((<*>))
import Clifford.Signature
import Clifford.Blade
import Clifford.Field

-- Fixities for geometric algebra operators
infixl 7 <*>, /\, \/, <.>, >.>, <|>
infixl 6 `addMV`, `subMV`
infix 4 ==~

-- | Vector space over multivectors with signature @sig@ and scalar field @k@.
class (Field k) => CliffordVectorSpace sig k mv | mv -> sig k where
  -- | The zero multivector (\(0\)).
  zeroMV :: mv

  -- | Lift a scalar into grade-0 multivector (\(\alpha \cdot 1\)).
  scalarMV :: k -> mv

  -- | Construct 1-vector basis element \(e_i\) (1-indexed).
  basisMV :: Int -> mv
  basisMV i = bladeMV (basisBlade i) (fromMetricScale 1)

  -- | Construct a single-blade multivector \(\alpha \cdot e_K\).
  bladeMV :: Blade -> k -> mv

  -- | Grade projection \(\langle A \rangle_k\), retaining only components of grade @k@.
  gradeProj :: Int -> mv -> mv

  -- | List of active grade indices present in the multivector.
  grades :: mv -> [Int]

  -- | Lookup the scalar coefficient for a specific basis blade.
  getBlade :: Blade -> mv -> k

  -- | Construct a multivector from a list of blade-coefficient pairs.
  fromBladeList :: [(Blade, k)] -> mv

  -- | Convert a multivector to a list of non-zero blade-coefficient pairs.
  toBladeList :: mv -> [(Blade, k)]

  -- | Multivector addition (\(A + B\)).
  addMV :: mv -> mv -> mv

  -- | Multivector subtraction (\(A - B\)).
  subMV :: mv -> mv -> mv

  -- | Scalar multiplication (\(\alpha \cdot A\)).
  scaleMV :: k -> mv -> mv

  -- | Multivector negation (\(-A\)).
  negateMV :: mv -> mv
  negateMV = scaleMV (fromMetricScale (-1))

-- | Clifford Algebra structure providing products, involutions, and dualities.
class (CliffordVectorSpace sig k mv) => CliffordAlgebra sig k mv where
  -- | Geometric product: \(A \star B\).
  (<*>) :: mv -> mv -> mv

  -- | Exterior / Wedge product: \(A \wedge B\).
  (/\) :: mv -> mv -> mv

  -- | Regressive / Meet product: \(A \vee B = J^{-1}(J(A) \wedge J(B))\).
  (\/) :: mv -> mv -> mv
  a \/ b = poincareDualMV (poincareDualMV a /\ poincareDualMV b)

  -- | Left contraction (Dorst/Lounesto): \(A \rfloor B = \sum_{r,s} \langle \langle A \rangle_r \langle B \rangle_s \rangle_{s-r}\).
  (<.>) :: mv -> mv -> mv

  -- | Right contraction: \(A \lfloor B = \sum_{r,s} \langle \langle A \rangle_r \langle B \rangle_s \rangle_{r-s}\).
  (>.>) :: mv -> mv -> mv

  -- | Scalar inner product: \(\langle A B \rangle_0 \in k\).
  (<|>) :: mv -> mv -> k

  -- | Reversion anti-automorphism: \(\widetilde{A} = \sum_k (-1)^{k(k-1)/2} \langle A \rangle_k\).
  --
  -- Satisfies \(\widetilde{A B} = \widetilde{B} \widetilde{A}\).
  reverseMV :: mv -> mv

  -- | Grade involution automorphism: \(\widehat{A} = \sum_k (-1)^k \langle A \rangle_k\).
  --
  -- Satisfies \(\widehat{A B} = \widehat{A} \widehat{B}\).
  involMV :: mv -> mv

  -- | Clifford conjugation: \(\overline{A} = \widetilde{\widehat{A}} = \sum_k (-1)^{k(k+1)/2} \langle A \rangle_k\).
  conjMV :: mv -> mv
  conjMV = reverseMV . involMV

  -- | Poincaré complement dual: \(J(A)\).
  --
  -- Grade-reversing complement bijection well-defined for all signatures including PGA (\(r > 0\)).
  poincareDualMV :: mv -> mv

  -- | Squared magnitude norm: \(\langle A \widetilde{A} \rangle_0\).
  normSqMV :: mv -> k
  normSqMV a = a <|> reverseMV a

-- | Commutator: \([A, B] = \frac{1}{2} (A B - B A)\).
  commutator :: mv -> mv -> mv
  commutator a b = scaleMV (fromRational (1/2)) (subMV (a <*> b) (b <*> a))

  -- | Anti-commutator: \(\{A, B\} = \frac{1}{2} (A B + B A)\).
  antiCommutator :: mv -> mv -> mv
  antiCommutator a b = scaleMV (fromRational (1/2)) (addMV (a <*> b) (b <*> a))

-- | Multivector Frobenius magnitude norm: \(\|A\| = \sqrt{\langle A \widetilde{A} \rangle_0}\).
normMV :: (Floating k, CliffordAlgebra sig k mv) => mv -> k
normMV a = sqrt (abs (normSqMV a))
{-# INLINE normMV #-}

-- | Normalize a non-zero multivector to unit norm.
normalizeMV :: (Floating k, CliffordAlgebra sig k mv) => mv -> mv
normalizeMV a =
  let n = normMV a
  in if n == 0 then a else scaleMV (1 / n) a
{-# INLINE normalizeMV #-}

-- | Format any multivector as a publication-ready LaTeX mathematical string.
--
-- Example:
--
-- >>> toLatex (2 * e1 - 3 * e12 :: DenseMV Cl300 Double)
-- "2.0 e_{1} - 3.0 e_{12}"
toLatex :: forall sig k mv. (CliffordVectorSpace sig k mv, Show k, Eq k, Num k) => mv -> String
toLatex mv =
  let terms = toBladeList mv
      nonZero = filter (\(_, c) -> c /= 0) terms
  in if null nonZero
       then "0"
       else unwords (zipWith formatTerm (True : repeat False) nonZero)
  where
    formatTerm isFirst (b, c) =
      let bStr = if unBlade b == 0
                   then ""
                   else "e_{" ++ concatMap show (bladeIndices b) ++ "}"
          cStr = show (abs c)
          valStr = if null bStr then cStr else (if abs c == 1 then "" else cStr ++ " ") ++ bStr
          signStr
            | isFirst   = if show c /= "" && head (show c) == '-' then "-" else ""
            | show c /= "" && head (show c) == '-' = "- "
            | otherwise = "+ "
      in signStr ++ valStr

-- | Extract the even-grade subalgebra component (\(Cl^+\)): \(\sum_{k \text{ even}} \langle A \rangle_k\).
evenPart :: (CliffordVectorSpace sig k mv) => mv -> mv
evenPart a = foldr addMV zeroMV [ gradeProj g a | g <- grades a, even g ]
{-# INLINE evenPart #-}

-- | Extract the odd-grade subspace component (\(Cl^-\)): \(\sum_{k \text{ odd}} \langle A \rangle_k\).
oddPart :: (CliffordVectorSpace sig k mv) => mv -> mv
oddPart a = foldr addMV zeroMV [ gradeProj g a | g <- grades a, odd g ]
{-# INLINE oddPart #-}

-- | Hestenes symmetric inner product ("fat-dot" \(\bullet\)):
--
-- \[ \langle A \rangle_r \bullet \langle B \rangle_s = \begin{cases} \langle A B \rangle_{|r - s|} & \text{if } r, s > 0 \\ 0 & \text{otherwise} \end{cases} \]
hestenesDot :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
hestenesDot a b =
  foldr addMV zeroMV
    [ gradeProj (abs (r - s)) (gradeProj r a <*> gradeProj s b)
    | r <- filter (> 0) (grades a)
    , s <- filter (> 0) (grades b)
    ]

-- | Named alias for left contraction (\(A \rfloor B\)).
contractLeft :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
contractLeft = (<.>)
{-# INLINE contractLeft #-}

-- | Named alias for right contraction (\(A \lfloor B\)).
contractRight :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
contractRight = (>.>)
{-# INLINE contractRight #-}

-- | Named alias for scalar inner product (\(\langle A B \rangle_0\)).
scalarProd :: (CliffordAlgebra sig k mv) => mv -> mv -> k
scalarProd = (<|>)
{-# INLINE scalarProd #-}

-- | Approximate equality for multivectors.
(==~) :: forall sig k mv. (CliffordVectorSpace sig k mv) => mv -> mv -> Bool
a ==~ b =
  let diff = subMV a b
      blades = toBladeList diff
  in all (isZero . snd) blades
{-# INLINE (==~) #-}
