# `cliff` API Reference

## Module Index

1. [`Clifford.Signature`](#1-cliffordsignature) — Compile-time metric signatures and type families.
2. [`Clifford.Blade`](#2-cliffordblade) — 64-bit basis blade bitmask combinatorics and signs.
3. [`Clifford.Field`](#3-cliffordfield) — Scalar field abstractions (Real, Complex, Rational, $\mathbb{F}_2, \mathbb{F}_p, \mathbb{F}_{2^m}$).
4. [`Clifford.Class`](#4-cliffordclass) — Core algebraic typeclasses and infix operators.
5. [`Clifford.Dense`](#5-clifforddense) — Unboxed dense multivectors with Cayley unrolling.
6. [`Clifford.Sparse`](#6-cliffordsparse) — Bitmask-indexed sparse multivectors.
7. [`Clifford.Versor`](#7-cliffordversor) — Rotors, motors, exponentials, and inverses.
8. [`Clifford.Array.Tensor`](#8-cliffordarraytensor) — $N$-dimensional multivector buffers.
9. [`Clifford.Array.Matrix`](#9-cliffordarraymatrix) — Multivector GEMM, structure matrices, and attention.
10. [`Clifford.Array.Convolution`](#10-cliffordarrayconvolution) — Sliding geometric stencils in `ST`.
11. [`Clifford.Array.Fourier`](#11-cliffordarrayfourier) — Discrete & Fast Clifford Fourier Transforms.
12. [`Clifford.Analytic`](#12-cliffordanalytic) — Analytic functions, Taylor series, and polynomial evaluation.
13. [`Clifford`](#13-clifford) — Top-level re-export and basis shorthands.

---

## 1. `Clifford.Signature`

```haskell
-- | Metric signature Cl(p,q,r)
data Signature (p :: Nat) (q :: Nat) (r :: Nat)

-- | Total dimension n = p + q + r
type family Dim (sig :: Signature p q r) :: Nat

-- | Total blade count 2^n
type family TotalBlades (sig :: Signature p q r) :: Nat

-- | Runtime reflection class
class (KnownNat p, KnownNat q, KnownNat r) => KnownSignature (p :: Nat) (q :: Nat) (r :: Nat) where
  signatureP   :: Proxy ('Signature p q r) -> Int
  signatureQ   :: Proxy ('Signature p q r) -> Int
  signatureR   :: Proxy ('Signature p q r) -> Int
  signatureDim :: Proxy ('Signature p q r) -> Int

-- | Extract (p, q, r) tuple
signatureVal :: forall p q r proxy. KnownSignature p q r => proxy ('Signature p q r) -> (Int, Int, Int)

-- | Metric quadratic form for basis vector e_i (1-based index)
metricSign :: forall p q r. KnownSignature p q r => Int -> Int

-- | Common Signatures
type Cl200 = 'Signature 2 0 0  -- 2D Euclidean
type Cl300 = 'Signature 3 0 0  -- 3D Euclidean GA
type Cl130 = 'Signature 1 3 0  -- Spacetime Algebra (STA, +---)
type Cl310 = 'Signature 3 1 0  -- Spacetime Algebra (-+++)
type Cl301 = 'Signature 3 0 1  -- 3D Projective Geometric Algebra (PGA)
type Cl410 = 'Signature 4 1 0  -- Conformal Geometric Algebra (CGA)
```

---

## 2. `Clifford.Blade`

```haskell
-- | 64-bit Basis Blade Bitmask
newtype Blade = Blade { unBlade :: Word64 } deriving (Eq, Ord, Show)

scalarBlade        :: Blade
basisBlade         :: Int -> Blade
pseudoscalarBlade  :: Int -> Blade
bladeFromIndices   :: [Int] -> Blade
bladeIndices       :: Blade -> [Int]
bladeGrade         :: Blade -> Int

canonicalSwapParity :: Blade -> Blade -> Int
metricFactor        :: forall p q r. KnownSignature p q r => Blade -> Blade -> Int
bladeMul            :: forall p q r. KnownSignature p q r => Blade -> Blade -> (Int, Blade)
bladeWedge          :: Blade -> Blade -> Maybe Blade
bladeDot            :: forall p q r. KnownSignature p q r => Blade -> Blade -> Maybe (Int, Blade)
poincareDual        :: Int -> Blade -> (Int, Blade)
allBlades           :: Int -> [Blade]
bladesOfGrade       :: Int -> Int -> [Blade]
```

---

## 3. `Clifford.Field`

```haskell
class (Eq k, Fractional k) => Field k where
  fromMetricScale :: Int -> k
  fieldInverse    :: k -> Maybe k
  isZero          :: k -> Bool
  fieldConj       :: k -> k

-- | Instances: Double, Float, Rational, Complex Double

-- | Binary Field F_2 = {0, 1}
newtype GF2 = GF2 { unGF2 :: Word8 } deriving (Eq, Ord, Show)

-- | Prime Modular Field F_p
newtype GFp (p :: Nat) = GFp { unGFp :: Integer } deriving (Eq, Ord, Show)

-- | Galois Extension Field F_(2^m) with irreducible polynomial
newtype GF2m (m :: Nat) (poly :: Nat) = GF2m { unGF2m :: Word64 } deriving (Eq, Ord, Show)

extGCD  :: Integer -> Integer -> (Integer, Integer, Integer)
modPow  :: Integer -> Integer -> Integer -> Integer
gf2mInv :: Word64 -> Int -> Word64 -> Word64
```

---

## 4. `Clifford.Class`

```haskell
-- Operator Precedences
infixl 7 <*>, /\, \/, <.>, >.>, <|>
infixl 6 `addMV`, `subMV`

class (Field k, KnownSignature p q r) => CliffordVectorSpace (sig :: Signature p q r) k mv | mv -> sig k where
  zeroMV        :: mv
  scalarMV      :: k -> mv
  basisMV       :: Int -> mv
  bladeMV       :: Blade -> k -> mv
  gradeProj     :: Int -> mv -> mv
  grades        :: mv -> [Int]
  getBlade      :: Blade -> mv -> k
  fromBladeList :: [(Blade, k)] -> mv
  toBladeList   :: mv -> [(Blade, k)]
  addMV         :: mv -> mv -> mv
  subMV         :: mv -> mv -> mv
  scaleMV       :: k -> mv -> mv
  negateMV      :: mv -> mv

class (CliffordVectorSpace sig k mv) => CliffordAlgebra sig k mv where
  (<*>)          :: mv -> mv -> mv       -- Geometric product
  (/\)           :: mv -> mv -> mv       -- Exterior wedge product
  (\/)           :: mv -> mv -> mv       -- Regressive meet product
  (<.>)          :: mv -> mv -> mv       -- Left contraction (Dorst)
  (>.>)          :: mv -> mv -> mv       -- Right contraction
  (<|>)          :: mv -> mv -> k        -- Scalar inner product
  reverseMV      :: mv -> mv             -- Reversion ~A
  involMV        :: mv -> mv             -- Grade involution ^A
  conjMV         :: mv -> mv             -- Clifford conjugation -A~
  poincareDualMV :: mv -> mv             -- Poincaré complement dual J(A)
  normSqMV       :: mv -> k              -- <A ~A>_0
  commutator     :: mv -> mv -> mv       -- 1/2 (AB - BA)
  antiCommutator :: mv -> mv -> mv       -- 1/2 (AB + BA)

-- | Subalgebra Projections & Involutions
evenPart    :: (CliffordVectorSpace sig k mv) => mv -> mv  -- Even subalgebra Cl^+
oddPart     :: (CliffordVectorSpace sig k mv) => mv -> mv  -- Odd subspace Cl^-
hestenesDot :: (CliffordAlgebra sig k mv) => mv -> mv -> mv -- Symmetric inner product
normMV      :: (Floating k, CliffordAlgebra sig k mv) => mv -> k
normalizeMV :: (Floating k, CliffordAlgebra sig k mv) => mv -> mv
toLatex     :: (CliffordVectorSpace sig k mv, Show k, Eq k, Num k) => mv -> String
```

---

## 5. `Clifford.Dense`

```haskell
data DenseMV (sig :: Signature p q r) k = DenseMV !(UArray Int k) deriving (Eq, Show)

mkDense       :: (KnownSignature p q r, Field k, IArray UArray k) => [(Blade, k)] -> DenseMV sig k
denseCoeffs   :: DenseMV sig k -> UArray Int k
denseIndex    :: Blade -> Int
denseBladeAt  :: Int -> Blade
```

---

## 6. `Clifford.Sparse`

```haskell
newtype SparseMV (sig :: Signature p q r) k = SparseMV { unSparse :: Map.Map Blade k } deriving (Eq, Show)

mkSparse          :: (Field k) => [(Blade, k)] -> SparseMV sig k
sparseLookup      :: (Field k) => Blade -> SparseMV sig k -> k
sparseInsert      :: (Field k) => Blade -> k -> SparseMV sig k -> SparseMV sig k
sparseFilterZeros :: (Field k) => SparseMV sig k -> SparseMV sig k
toDense           :: (KnownSignature p q r, Field k, IArray UArray k) => SparseMV sig k -> DenseMV sig k
toSparse          :: (KnownSignature p q r, Field k, IArray UArray k) => DenseMV sig k -> SparseMV sig k
```

---

## 7. `Clifford.Versor`

```haskell
-- | Versor Sandwiches & Hyperplane Reflections
sandwich       :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
reflectMV      :: (CliffordAlgebra sig k mv, Eq k, Fractional k) => mv -> mv -> Maybe mv

-- | Subspace Projections & Rejections
projectMV      :: (CliffordAlgebra sig k mv, Eq k, Fractional k) => mv -> mv -> Maybe mv
rejectMV       :: (CliffordAlgebra sig k mv, Eq k, Fractional k) => mv -> mv -> Maybe mv

-- | Rotor Construction & Exponentials
rotorNorm      :: (CliffordAlgebra sig k mv, Floating k) => mv -> k
normalizeRotor :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv
rotorExp       :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv
rotorLog       :: (CliffordAlgebra sig k mv, Floating k) => mv -> Maybe mv
rotorBetween   :: (CliffordAlgebra sig k mv, Floating k) => mv -> mv -> Maybe mv
outermorphism  :: (CliffordAlgebra sig k mv) => (mv -> mv) -> mv -> mv
versorInverse  :: (Field k, CliffordAlgebra sig k mv) => mv -> Maybe mv
generalInverse :: (Field k, KnownSignature p q r, IArray UArray k) => DenseMV sig k -> Maybe (DenseMV sig k)
```

---

## 8. `Clifford.Array.Tensor`

```haskell
data MVArray (sig :: Signature p q r) k = MVArray
  { mvaShape   :: ![Int]
  , mvaStrides :: ![Int]
  , mvaBuffer  :: !(UArray Int k)
  } deriving (Eq, Show)

fromListND      :: (KnownSignature p q r, Field k, IArray UArray k) => [Int] -> [DenseMV sig k] -> MVArray sig k
generateND      :: (KnownSignature p q r, Field k, IArray UArray k) => [Int] -> ([Int] -> DenseMV sig k) -> MVArray sig k
indexND         :: (KnownSignature p q r, Field k, IArray UArray k) => MVArray sig k -> [Int] -> DenseMV sig k
mapMV           :: (KnownSignature p q r, Field k, IArray UArray k) => (DenseMV sig k -> DenseMV sig k) -> MVArray sig k -> MVArray sig k
zipWithMV       :: (KnownSignature p q r, Field k, IArray UArray k) => (DenseMV sig k -> DenseMV sig k -> DenseMV sig k) -> MVArray sig k -> MVArray sig k -> MVArray sig k
sliceGrade      :: (KnownSignature p q r, Field k, IArray UArray k) => Int -> MVArray sig k -> MVArray sig k
dotProductArray :: (KnownSignature p q r, Field k, IArray UArray k) => MVArray sig k -> MVArray sig k -> k
```

---

## 9. `Clifford.Array.Matrix`

```haskell
data MVMatrix (sig :: Signature p q r) k = MVMatrix
  { mvmRows   :: !Int
  , mvmCols   :: !Int
  , mvmBuffer :: !(UArray Int k)
  } deriving (Eq, Show)

matMulMV                 :: (KnownSignature p q r, Field k, IArray UArray k) => MVMatrix sig k -> MVMatrix sig k -> MVMatrix sig k
matVecMulMV              :: (KnownSignature p q r, Field k, IArray UArray k) => MVMatrix sig k -> [DenseMV sig k] -> [DenseMV sig k]
toStructureMatrix        :: (KnownSignature p q r, Field k, IArray UArray k) => DenseMV sig k -> UArray (Int, Int) k
attentionScoresInvariant :: (KnownSignature p q r, Field k, IArray UArray k) => MVMatrix sig k -> MVMatrix sig k -> UArray (Int, Int) k
adjointMulLeft           :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
adjointMulRight          :: (CliffordAlgebra sig k mv) => mv -> mv -> mv
```

---

## 10. `Clifford.Array.Convolution`

```haskell
data BoundaryMode = ZeroPad | Periodic | Replicate | Mirror deriving (Eq, Show)
data ConvProduct  = GeometricProd | WedgeProd | LeftContract | RightContract deriving (Eq, Show)

convolve1D :: (KnownSignature p q r, Field k, IArray UArray k) => BoundaryMode -> ConvProduct -> MVArray sig k -> MVArray sig k -> MVArray sig k
convolve2D :: (KnownSignature p q r, Field k, IArray UArray k) => BoundaryMode -> ConvProduct -> MVArray sig k -> MVArray sig k -> MVArray sig k
convolveND :: (KnownSignature p q r, Field k, IArray UArray k) => BoundaryMode -> ConvProduct -> MVArray sig k -> MVArray sig k -> MVArray sig k
```

---

## 11. `Clifford.Array.Fourier`

```haskell
-- | Fast Clifford Fourier Transforms (O(N log N) Radix-2 Cooley-Tukey)
fcft1D  :: (KnownSignature p q r, Floating k, UnboxedField k) => DenseMV sig k -> MVArray sig k -> MVArray sig k
ifcft1D :: (KnownSignature p q r, Floating k, UnboxedField k) => DenseMV sig k -> MVArray sig k -> MVArray sig k
fcft2D  :: (KnownSignature p q r, Floating k, UnboxedField k) => DenseMV sig k -> DenseMV sig k -> MVArray sig k -> MVArray sig k
ifcft2D :: (KnownSignature p q r, Floating k, UnboxedField k) => DenseMV sig k -> DenseMV sig k -> MVArray sig k -> MVArray sig k

-- | Galois Field Number Theoretic Transforms (F_p exact modular arithmetic)
ntt1D   :: (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime)) => GFp prime -> MVArray sig (GFp prime) -> MVArray sig (GFp prime)
intt1D  :: (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime)) => GFp prime -> MVArray sig (GFp prime) -> MVArray sig (GFp prime)
ntt2D   :: (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime)) => GFp prime -> GFp prime -> MVArray sig (GFp prime) -> MVArray sig (GFp prime)
intt2D  :: (KnownSignature p q r, KnownNat prime, UnboxedField (GFp prime)) => GFp prime -> GFp prime -> MVArray sig (GFp prime) -> MVArray sig (GFp prime)

-- | Fast Walsh-Hadamard Transforms (FWHT & Dyadic Convolutions)
fwhtMV           :: (KnownSignature p q r, UnboxedField k) => DenseMV sig k -> DenseMV sig k
ifwhtMV          :: (KnownSignature p q r, UnboxedField k) => DenseMV sig k -> DenseMV sig k
dyadicConvolveMV :: (KnownSignature p q r, UnboxedField k) => DenseMV sig k -> DenseMV sig k -> DenseMV sig k
fwht1D           :: (KnownSignature p q r, UnboxedField k) => MVArray sig k -> MVArray sig k
ifwht1D          :: (KnownSignature p q r, UnboxedField k) => MVArray sig k -> MVArray sig k
fwht2D           :: (KnownSignature p q r, UnboxedField k) => MVArray sig k -> MVArray sig k
ifwht2D          :: (KnownSignature p q r, UnboxedField k) => MVArray sig k -> MVArray sig k
dyadicConvolve1D :: (KnownSignature p q r, UnboxedField k) => (DenseMV sig k -> DenseMV sig k -> DenseMV sig k) -> MVArray sig k -> MVArray sig k -> MVArray sig k

-- | Direct Reference Discrete Transforms (Arbitrary Dims & Kernels)
data CFTKernel (sig :: Signature p q r) k = CFTKernel
  { cftBivectors  :: ![DenseMV sig k]
  , cftDimensions :: ![Int]
  } deriving (Eq, Show)

mkCFTKernel    :: (KnownSignature p q r, Floating k, UnboxedField k) => [DenseMV sig k] -> [Int] -> Maybe (CFTKernel sig k)
cftForward     :: (KnownSignature p q r, Floating k, UnboxedField k) => CFTKernel sig k -> MVArray sig k -> MVArray sig k
cftInverse     :: (KnownSignature p q r, Floating k, UnboxedField k) => CFTKernel sig k -> MVArray sig k -> MVArray sig k
spectralFilter :: (KnownSignature p q r, Floating k, UnboxedField k) => CFTKernel sig k -> (MVArray sig k -> MVArray sig k) -> MVArray sig k -> MVArray sig k
```

---

---

## 12. `Clifford.Analytic`

```haskell
-- | Scaling-and-Squaring General Multivector Exponential
expMV     :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv

-- | Trigonometric Functions
sinMV     :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
cosMV     :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
sincosMV  :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> (mv, mv)

-- | Hyperbolic Functions
sinhMV    :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
coshMV    :: (RealFrac k, Floating k, CliffordAlgebra sig k mv) => mv -> mv
tanhMV    :: (RealFrac k, Floating k, Eq k, CliffordAlgebra sig k mv) => mv -> Maybe mv

-- | Horner Polynomial Evaluation: P(M) = sum c_k M^k
hornerMV  :: (CliffordAlgebra sig k mv) => [k] -> mv -> mv

-- | Finite Field Discrete Power Map (Frobenius Automorphism Phi_p(M) = M^p)
frobeniusMV :: (CliffordAlgebra sig k mv) => Integer -> mv -> mv
```

---

## 13. `Clifford`

```haskell
-- Standard Basis Helpers
e1, e2, e3, e4, e5, e0 :: (CliffordVectorSpace sig k mv) => mv
e12, e23, e31, e123   :: (CliffordAlgebra sig k mv) => mv

-- 3D PGA Kinematics (Cl(3,0,1))
pointPGA              :: (CliffordAlgebra Cl301 k mv) => k -> k -> k -> mv
planePGA              :: (CliffordVectorSpace Cl301 k mv) => k -> k -> k -> k -> mv
linePGA               :: (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
intersectPlanesPGA    :: (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
intersectPlaneLinePGA :: (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv
translatorPGA         :: (CliffordAlgebra Cl301 k mv) => k -> k -> k -> mv
rotatorPGA            :: (Floating k, CliffordAlgebra Cl301 k mv) => k -> mv -> mv
motorPGA              :: (CliffordAlgebra Cl301 k mv) => mv -> mv -> mv

-- Conformal Geometric Algebra (CGA Cl(4,1,0))
eInfCGA, e0CGA        :: (CliffordAlgebra Cl410 k mv) => mv
pointCGA              :: (Floating k, CliffordAlgebra Cl410 k mv) => k -> k -> k -> mv
sphereCGA             :: (Floating k, CliffordAlgebra Cl410 k mv) => mv -> k -> mv
sqDistanceCGA         :: (CliffordAlgebra Cl410 k mv) => mv -> mv -> k
```
