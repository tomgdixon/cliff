{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.StencilSpec (stencilTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit

import Clifford

stencilTests :: TestTree
stencilTests = testGroup "Sliding Stencils & Fast Clifford Fourier Transforms"
  [ testCase "2D Convolution with centered delta impulse reproduces input image" $ do
      let v1 = (scaleMV 2 e1 `addMV` scaleMV 3 e2) :: DenseMV Cl200 Double
          v2 = (scaleMV 1 e2) :: DenseMV Cl200 Double
          v3 = (scaleMV 4 e12) :: DenseMV Cl200 Double
          v4 = scalarMV 5 :: DenseMV Cl200 Double

          img = fromListND [2, 2] [v1, v2, v3, v4]

          -- Centered 3x3 delta kernel with 1 at center (1,1)
          deltaKernel = fromListND [3, 3]
            [ zeroMV, zeroMV, zeroMV
            , zeroMV, scalarMV 1, zeroMV
            , zeroMV, zeroMV, zeroMV
            ]

          convolved = convolve2D Periodic GeometricProd deltaKernel img

      convolved @?= img

  , testCase "2D Discrete Clifford Fourier Transform (DCFT) invertibility: CFT^-1(CFT(f)) == f" $ do
      let v1 = (scaleMV 2 e1 `addMV` scaleMV 3 e2) :: DenseMV Cl200 Double
          v2 = (scaleMV 1 e2) :: DenseMV Cl200 Double
          v3 = (scaleMV 4 e12) :: DenseMV Cl200 Double
          v4 = scalarMV 5 :: DenseMV Cl200 Double

          grid = fromListND [2, 2] [v1, v2, v3, v4]
          kernel = defaultCFT2D 2 2

          freqDomain = cftForward kernel grid
          reconstructed = cftInverse kernel freqDomain

      let diff = maximum [ normSqMV (indexND grid [r, c] `subMV` indexND reconstructed [r, c])
                         | r <- [0 .. 1], c <- [0 .. 1]
                         ]
      diff < 1e-10 @?= True

  , testCase "1D Fast Clifford Fourier Transform (FCFT): ifcft1D(fcft1D(f)) == f (N=4)" $ do
      let v0 = (scaleMV 1 e1 `addMV` scaleMV 2 e2) :: DenseMV Cl200 Double
          v1 = (scaleMV 3 e1 `subMV` scaleMV 1 e12) :: DenseMV Cl200 Double
          v2 = (scaleMV 4 e2 `addMV` scalarMV 2) :: DenseMV Cl200 Double
          v3 = (scaleMV (-2) e1 `addMV` scaleMV 5 e12) :: DenseMV Cl200 Double
          arr = fromListND [4] [v0, v1, v2, v3]
          e12_blade = e12 :: DenseMV Cl200 Double

          freq = fcft1D e12_blade arr
          reconstructed = ifcft1D e12_blade freq

      let diff = maximum [ normSqMV (indexND arr [i] `subMV` indexND reconstructed [i])
                         | i <- [0 .. 3]
                         ]
      diff < 1e-10 @?= True

  , testCase "2D Fast Clifford Fourier Transform (FCFT): ifcft2D(fcft2D(f)) == f (4x4 Grid)" $ do
      let cell r c = (scaleMV (fromIntegral (r + 1)) e1 `addMV` scaleMV (fromIntegral (c * 2)) e12) :: DenseMV Cl200 Double
          gridList = [ cell r c | r <- [0 .. 3], c <- [0 .. 3] ]
          grid = fromListND [4, 4] gridList
          e12_blade = e12 :: DenseMV Cl200 Double

          freq = fcft2D e12_blade e12_blade grid
          reconstructed = ifcft2D e12_blade e12_blade freq

      let diff = maximum [ normSqMV (indexND grid [r, c] `subMV` indexND reconstructed [r, c])
                         | r <- [0 .. 3], c <- [0 .. 3]
                         ]
      diff < 1e-10 @?= True

  , testCase "1D Number Theoretic Transform (NTT) over GF(17): intt1D(ntt1D(f)) == f (N=8, omega=2)" $ do
      -- In F_17, 2 is a primitive 8-th root of unity (2^8 = 256 = 15*17 + 1 = 1 mod 17)
      let omega = GFp 2 :: GFp 17
          cell i = (scaleMV (GFp (fromIntegral (i + 1))) e1
                    `addMV` scaleMV (GFp (fromIntegral (3 * i + 2))) e2) :: DenseMV Cl200 (GFp 17)
          arr = fromListND [8] [ cell i | i <- [0 .. 7] ]
          fwd = ntt1D omega arr
          reconstructed = intt1D omega fwd

      reconstructed @?= arr

  , testCase "2D Number Theoretic Transform (NTT) over GF(17): intt2D(ntt2D(f)) == f (4x4, omega=4)" $ do
      -- In F_17, 4 is a primitive 4-th root of unity (4^4 = 256 = 1 mod 17)
      let omega = GFp 4 :: GFp 17
          cell r c = (scaleMV (GFp (fromIntegral (r + 1))) e1
                      `addMV` scaleMV (GFp (fromIntegral (2 * c + 3))) e2) :: DenseMV Cl200 (GFp 17)
          grid = fromListND [4, 4] [ cell r c | r <- [0 .. 3], c <- [0 .. 3] ]
          fwd = ntt2D omega omega grid
          reconstructed = intt2D omega omega fwd

      reconstructed @?= grid
  ]
