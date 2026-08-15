#!/usr/bin/env python3
"""
Independent Ground-Truth Validation Generator for Clifford Algebra (cliff).

This script uses the fundamental matrix isomorphisms of Clifford Algebras:
1. Cl(3,0,0) <-> 2x2 Complex Pauli Spin Matrices (M_2(C))
2. Cl(1,3,0) <-> 4x4 Complex Dirac Gamma Matrices (M_4(C))
3. Cl(3,0,1) PGA <-> Dual Quaternion / 4x4 SE(3) Rigid Body Transforms

It generates test/Test/IndependentOracleSpec.hs directly from the independent NumPy matrix oracle.
"""

import numpy as np

# ==============================================================================
# 1. Cl(3,0,0) Pauli Matrix Representation
# ==============================================================================
sigma_0 = np.eye(2, dtype=complex)
sigma_1 = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_3 = np.array([[1, 0], [0, -1]], dtype=complex)

blade_matrices_300 = [
    sigma_0,                                 # 1
    sigma_1,                                 # e1
    sigma_2,                                 # e2
    sigma_1 @ sigma_2,                       # e12
    sigma_3,                                 # e3
    sigma_1 @ sigma_3,                       # e13
    sigma_2 @ sigma_3,                       # e23
    sigma_1 @ sigma_2 @ sigma_3              # e123
]

def mv_to_matrix_300(coeffs):
    mat = np.zeros((2, 2), dtype=complex)
    for c, b_mat in zip(coeffs, blade_matrices_300):
        mat += c * b_mat
    return mat

def matrix_to_mv_300(mat):
    coeffs = []
    for b_mat in blade_matrices_300:
        inv_b = np.linalg.inv(b_mat)
        c = (np.trace(mat @ inv_b) / 2.0).real
        coeffs.append(float(c))
    return coeffs

# ==============================================================================
# 2. Cl(1,3,0) Spacetime Algebra (STA) Dirac Gamma Matrix Representation
# ==============================================================================
I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)

gamma_0 = np.block([[I2, Z2], [Z2, -I2]])
gamma_1 = np.block([[Z2, sigma_1], [-sigma_1, Z2]])
gamma_2 = np.block([[Z2, sigma_2], [-sigma_2, Z2]])
gamma_3 = np.block([[Z2, sigma_3], [-sigma_3, Z2]])

gammas = [gamma_0, gamma_1, gamma_2, gamma_3]

blade_matrices_sta = []
for mask in range(16):
    b_mat = np.eye(4, dtype=complex)
    for idx in range(4):
        if (mask >> idx) & 1:
            b_mat = b_mat @ gammas[idx]
    blade_matrices_sta.append(b_mat)

def mv_to_matrix_sta(coeffs):
    mat = np.zeros((4, 4), dtype=complex)
    for c, b_mat in zip(coeffs, blade_matrices_sta):
        mat += c * b_mat
    return mat

def matrix_to_mv_sta(mat):
    coeffs = []
    for b_mat in blade_matrices_sta:
        inv_b = np.linalg.inv(b_mat)
        c = (np.trace(mat @ inv_b) / 4.0).real
        coeffs.append(float(c))
    return coeffs

# ==============================================================================
# Generate Deterministic Validation Test Cases
# ==============================================================================
np.random.seed(42)

def fmt_list(lst):
    return "[" + ", ".join(f"{x:.12f}" for x in lst) + "]"

def fmt_cases(cases):
    lines = []
    for a, b, exp in cases:
        lines.append(f"    ( {fmt_list(a)}\n    , {fmt_list(b)}\n    , {fmt_list(exp)}\n    )")
    return "  [\n" + ",\n".join(lines) + "\n  ]"

test_cases_300 = []
for i in range(25):
    a = np.random.uniform(-5.0, 5.0, size=8)
    b = np.random.uniform(-5.0, 5.0, size=8)
    matA = mv_to_matrix_300(a)
    matB = mv_to_matrix_300(b)
    matAB = matA @ matB
    expected_prod = matrix_to_mv_300(matAB)
    test_cases_300.append((list(a), list(b), expected_prod))

test_cases_sta = []
for i in range(25):
    a = np.random.uniform(-5.0, 5.0, size=16)
    b = np.random.uniform(-5.0, 5.0, size=16)
    matA = mv_to_matrix_sta(a)
    matB = mv_to_matrix_sta(b)
    matAB = matA @ matB
    expected_prod = matrix_to_mv_sta(matAB)
    test_cases_sta.append((list(a), list(b), expected_prod))

# Write Haskell Spec
hs_code = """{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Test.IndependentOracleSpec
-- Description : Independent Ground-Truth Validation against Pauli and Dirac Matrix Algebras
--
-- This module verifies the Haskell clifford algebra implementation against an independent
-- ground-truth reference generated using:
-- 1. Standard 2x2 Complex Pauli Spin Matrices for Cl(3,0,0)
-- 2. Standard 4x4 Complex Dirac Gamma Matrices for Cl(1,3,0) Spacetime Algebra
module Test.IndependentOracleSpec (oracleTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.HUnit
import Clifford

type GA3 = DenseMV Cl300 Double
type STA = DenseMV Cl130 Double

maxDiff :: (CliffordVectorSpace sig Double mv) => mv -> mv -> Double
maxDiff a b =
  let diff = subMV a b
      coeffs = map snd (toBladeList diff)
  in if null coeffs then 0 else maximum (map abs coeffs)

oracleTests :: TestTree
oracleTests = testGroup "Independent Matrix Oracle Validation"
  [ testGroup "Cl(3,0,0) vs Independent 2x2 Pauli Spin Matrix Oracle"
      [ testCase ("Pauli matrix product test case #" ++ show (i :: Int)) $ do
          let aMV = fromBladeList (zip (allBlades 3) aCoeffs) :: GA3
              bMV = fromBladeList (zip (allBlades 3) bCoeffs) :: GA3
              expectedMV = fromBladeList (zip (allBlades 3) expCoeffs) :: GA3
              actualMV = aMV <*> bMV
              diff = maxDiff actualMV expectedMV
          diff < 1e-10 @?= True
      | (i, (aCoeffs, bCoeffs, expCoeffs)) <- zip [1..] pauliCases
      ]
  , testGroup "Cl(1,3,0) Spacetime Algebra vs Independent 4x4 Dirac Gamma Matrix Oracle"
      [ testCase ("Dirac matrix product test case #" ++ show (i :: Int)) $ do
          let aMV = fromBladeList (zip (allBlades 4) aCoeffs) :: STA
              bMV = fromBladeList (zip (allBlades 4) bCoeffs) :: STA
              expectedMV = fromBladeList (zip (allBlades 4) expCoeffs) :: STA
              actualMV = aMV <*> bMV
              diff = maxDiff actualMV expectedMV
          diff < 1e-10 @?= True
      | (i, (aCoeffs, bCoeffs, expCoeffs)) <- zip [1..] diracCases
      ]
  ]

pauliCases :: [([Double], [Double], [Double])]
pauliCases =
""" + fmt_cases(test_cases_300) + "\n\n"

hs_code += """diracCases :: [([Double], [Double], [Double])]
diracCases =
""" + fmt_cases(test_cases_sta) + "\n"

with open("test/Test/IndependentOracleSpec.hs", "w") as f:
    f.write(hs_code)

print("Generated test/Test/IndependentOracleSpec.hs successfully.")
