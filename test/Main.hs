module Main (main) where

import Test.Tasty

import Test.AlgebraLawsSpec (algebraLawsTests)
import Test.RepresentationSpec (representationTests)
import Test.FieldSpec (fieldTests)
import Test.VersorSpec (versorTests)
import Test.PGASpec (pgaTests)
import Test.MatrixSpec (matrixTests)
import Test.StencilSpec (stencilTests)
import Test.IndependentOracleSpec (oracleTests)
import Test.HighGradeSpec (highGradeTests)
import Test.AnalyticSpec (analyticTests)
import Test.FWHTSpec (fwhtTests)

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "cliff Library Test Suite"
  [ algebraLawsTests
  , representationTests
  , fieldTests
  , versorTests
  , pgaTests
  , matrixTests
  , stencilTests
  , oracleTests
  , highGradeTests
  , analyticTests
  , fwhtTests
  ]
