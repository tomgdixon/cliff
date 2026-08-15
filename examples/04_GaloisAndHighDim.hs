{-# LANGUAGE DataKinds #-}

-- |
-- Module      : Main
-- Description : Example 4 - Finite Galois Fields and High-Dimensional Clifford Algebras
--
-- This tutorial example covers:
-- 1. Characteristic 2 Binary Arithmetic over GF(2)
-- 2. Prime Modular Fields GF(p)
-- 3. Galois Extension Fields GF(2^8) (AES Rijndael Field)
-- 4. High-Dimensional 16D and 32D Clifford Algebras using SparseMV
module Main where

import Prelude hiding ((<*>))
import Clifford

main :: IO ()
main = do
  putStrLn "=========================================================================="
  putStrLn "  Clifford Algebra Tutorial: Example 4 - Galois Fields & High Dimensions"
  putStrLn "==========================================================================\n"

  -- 1. Binary Galois Field GF(2) = {0, 1}
  -- In GF(2), addition is XOR (1 + 1 = 0) and multiplication is AND (1 * 1 = 1).
  let g1 = GF2 1
      g0 = GF2 0

  putStrLn "--- 1. Binary Galois Field GF(2) ---"
  putStrLn $ "1 + 1 (in GF(2)) = " ++ show (g1 + g1) ++ "  (Addition is XOR)"
  putStrLn $ "1 + 0 (in GF(2)) = " ++ show (g1 + g0)
  putStrLn $ "1 * 1 (in GF(2)) = " ++ show (g1 * g1) ++ "  (Multiplication is AND)"
  putStrLn $ "-1    (in GF(2)) = " ++ show (negate g1) ++ "  (In char 2, -1 == 1)"
  putStrLn ""

  -- 2. 16-Dimensional Clifford Algebra over GF(2)
  -- Total basis blades in 16D = 2^16 = 65,536 blades.
  -- Using SparseMV, only non-zero active blades are stored.
  let v1_16 = basisMV 1 :: SparseMV Cl16_00 GF2
      v5_16 = basisMV 5 :: SparseMV Cl16_00 GF2
      v12_16 = basisMV 12 :: SparseMV Cl16_00 GF2
      v16_16 = basisMV 16 :: SparseMV Cl16_00 GF2

      mvA = v1_16 `addMV` v5_16 `addMV` v12_16
      mvB = v5_16 `addMV` v16_16

      prod16D = mvA <*> mvB

  putStrLn "--- 2. High-Dimensional 16D Clifford Algebra over GF(2) ---"
  putStrLn $ "Multivector A (in 16D) = " ++ show mvA
  putStrLn $ "Multivector B (in 16D) = " ++ show mvB
  putStrLn $ "Product A * B (in 16D) = " ++ show prod16D
  putStrLn "Notice how shared basis vector e5 contracted via e5^2 = 1 in GF(2)!"
  putStrLn ""

  -- 3. 32-Dimensional Clifford Algebra over GF(2)
  -- Total basis blades in 32D = 2^32 = 4,294,967,296 blades!
  let v7_32 = basisMV 7 :: SparseMV Cl32_00 GF2
      v29_32 = basisMV 29 :: SparseMV Cl32_00 GF2
      prod32D = (v7_32 `addMV` v29_32) <*> (v7_32 `addMV` v29_32)

  putStrLn "--- 3. 32D Clifford Algebra (2^32 = 4.29 Billion Blades) ---"
  putStrLn $ "(e7 + e29)^2 = " ++ show prod32D ++ " (Scalar 0 because e7^2=1, e29^2=1 => 1+1=0 in GF(2))"
  putStrLn ""

  -- 4. Galois Extension Field GF(2^8) (AES Cryptographic Field)
  -- GF(256) with irreducible reduction polynomial x^8 + x^4 + x^3 + x + 1 (0x11B)
  let a_gf8 = GF2m 0x53 :: GF2m 8 0x11B
      inv_a = recip a_gf8

  putStrLn "--- 4. Galois Extension Field GF(2^8) (AES Field) ---"
  putStrLn $ "Element a       = " ++ show a_gf8 ++ " (Byte 0x53)"
  putStrLn $ "Inverse a^(-1)  = " ++ show inv_a ++ " (Computed via Fermat exponentiation)"
  putStrLn $ "Verification a * a^(-1) = " ++ show (a_gf8 * inv_a) ++ " (Should be exactly 1)"
  putStrLn ""
  putStrLn "Tutorial 4 completed successfully!"
