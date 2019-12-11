



test_that('inweb can get interactors',{
  expect_equal(sum(interactors('CACNA1C', F)$significant), 53)
})