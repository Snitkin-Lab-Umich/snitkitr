test_that("ConvertTimeToHours returns expected times", {
  time <- "10:00:00"
  out <- ConvertTimeToHours(time)
  expect_equal(10, out)
  
  time <- "10:30:00"
  out <- ConvertTimeToHours(time)
  expect_equal(10.5, out)
  
  time <- "1-10:00:00"
  out <- ConvertTimeToHours(time)
  expect_equal(34, out)
})

test_that("ConvertTimeToHours gives error with incorrect input", {
  time <- 10
  expect_error(ConvertTimeToHours(time))
  
  time <- "10:30"
  expect_error(ConvertTimeToHours(time))
})