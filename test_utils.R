#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

source("utilities.R")

context("mdy to days")

test_that("incorrect format is not accepted", {
    dates <- c("md10", "10d", "d.1", "who knows", "mdy23")
    sapply(dates, function(x) { 
        expect_error(mdy_to_day(x))
    })
})

ave_month <- 30.5
ave_year <- 366

test_that("months processed", {
    dates <- c("m1", "m10", "m100")
    res <- ave_month * c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("years processed", {
    dates <- c("y1", "y10", "y100")
    res <- ave_year * c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("days processed", {
    dates <- c("d1", "d10", "d100")
    res <- c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("suffix text is removed", {
    dates <- c("d0.tdx2", "d0.tdx2.bag1", "d0.tdx2.bag2", "d0.untdx")
    res <- c(0, 0, 0, 0)
    expect_equal(mdy_to_day(dates), res)
})

test_that("suffix post is removed", {
    dates <- c("m2post", "m2.5post", "y2post")
    res <- c(ave_month*2, ave_month*2.5, ave_year*2)
    expect_equal(mdy_to_day(dates), res)
})

test_that("dot is accepted", {
    dates <- c("d1.", "d10.", "d100.")
    res <- c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("number dot number is accepted", {
    dates <- c("m1.1", "m10.1", "m100.1")
    res <- ave_month * c(1.1, 10.1, 100.1)
    expect_equal(mdy_to_day(dates), res)
})

test_that("can process mix of mdy", {
    dates <- c("d1", "d10", "d100")
    res <- c(1, 10, 100)
    dates <- c(dates, c("m1.1", "m10.1", "m100.1"))
    res <- c(res, ave_month * c(1.1, 10.1, 100.1))
    dates <- c(dates, c("y1", "y10", "y100"))
    res <- c(res, ave_year * c(1, 10, 100))
    expect_equal(mdy_to_day(dates), res)
})

test_that("can process input where mdy is not ordered as 'dmy'",{
  dates <- c("y1", "m4", "d10", "m3", "y2", "d99")
  res <- c(ave_year, 4*ave_month, 10.0, 3*ave_month, 2*ave_year, 99.0)
  expect_equal(mdy_to_day(dates), res)
})
