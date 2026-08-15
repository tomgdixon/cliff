{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module Test.Arbitrary () where

import Test.Tasty.QuickCheck
import Clifford

instance Arbitrary (DenseMV Cl300 Double) where
  arbitrary = do
    c0 <- choose (-10, 10)
    c1 <- choose (-10, 10)
    c2 <- choose (-10, 10)
    c3 <- choose (-10, 10)
    c4 <- choose (-10, 10)
    c5 <- choose (-10, 10)
    c6 <- choose (-10, 10)
    c7 <- choose (-10, 10)
    return $ fromBladeList
      [ (Blade 0, c0), (Blade 1, c1), (Blade 2, c2), (Blade 4, c3)
      , (Blade 3, c4), (Blade 5, c5), (Blade 6, c6), (Blade 7, c7)
      ]

instance Arbitrary (DenseMV Cl130 Double) where
  arbitrary = do
    coeffs <- vectorOf 16 (choose (-5, 5))
    return $ fromBladeList (zip (allBlades 4) coeffs)

instance Arbitrary (SparseMV Cl300 Double) where
  arbitrary = do
    dense <- (arbitrary :: Gen (DenseMV Cl300 Double))
    return (toSparse dense)
