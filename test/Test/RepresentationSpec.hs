{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Test.RepresentationSpec (representationTests) where

import Prelude hiding ((<*>))
import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit

import Clifford
import Test.Arbitrary ()

representationTests :: TestTree
representationTests = testGroup "Representation Equivalence (DenseMV vs SparseMV)"
  [ testProperty "Roundtrip: toDense (toSparse d) == d" $
      \(d :: DenseMV Cl300 Double) ->
        toDense (toSparse d) == d

  , testProperty "Roundtrip: toSparse (toDense s) == s" $
      \(s :: SparseMV Cl300 Double) ->
        toSparse (toDense s) == s

  , testProperty "Geometric product equivalence: toSparse (a * b) == toSparse a * toSparse b" $
      \(a :: DenseMV Cl300 Double) b ->
        toSparse (a <*> b) == (toSparse a <*> toSparse b)

  , testProperty "Wedge product equivalence: toSparse (a ^ b) == toSparse a ^ toSparse b" $
      \(a :: DenseMV Cl300 Double) b ->
        toSparse (a /\ b) == (toSparse a /\ toSparse b)

  , testProperty "Left contraction equivalence: toSparse (a . b) == toSparse a . toSparse b" $
      \(a :: DenseMV Cl300 Double) b ->
        toSparse (a <.> b) == (toSparse a <.> toSparse b)

  , testProperty "Reversion equivalence: toSparse (~a) == ~(toSparse a)" $
      \(a :: DenseMV Cl300 Double) ->
        toSparse (reverseMV a) == reverseMV (toSparse a)

  , testProperty "Poincaré dual equivalence: toSparse (J(a)) == J(toSparse a)" $
      \(a :: DenseMV Cl300 Double) ->
        toSparse (poincareDualMV a) == poincareDualMV (toSparse a)
  ]
