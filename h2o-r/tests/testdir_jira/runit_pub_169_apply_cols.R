######################################################################
# Test for PUB-169
# apply should work across columns
######################################################################

# setwd("/Users/tomk/0xdata/ws/h2o/R/tests/testdir_jira")

setwd(normalizePath(dirname(R.utils::commandArgs(asValues=TRUE)$"f")))
options(echo=TRUE)
source('../h2o-runit.R')

test.colapply <- function(conn) {
  Log.info('Uploading cebbinom.csv to H2O...')
  cebH2O <- h2o.importFile(conn, normalizePath(locate('smalldata/jira/cebbinom.csv')))

  Log.info('Print head of dataset')
  Log.info(head(cebH2O))

  Log.info("Slice a subset of columns")
  cebstdH2O = cebH2O[,8:23]

  Log.info("Slice again and apply scale over columns")
  scaledH2O = apply(cebstdH2O[,2:16], 2, scale)
  Log.info(head(scaledH2O))

  cebR <- read.csv(locate("smalldata/jira/cebbinom.csv"))
  cebstdR <- cebR[,8:23]
  scaledR <- apply(cebstdR[,2:16], 2, scale)

  Log.info("Comparing results to R")
  scaledH2O.df = as.data.frame(scaledH2O)
  expect_equal(scaledH2O.df, (as.data.frame(scaledR)))

  testEnd()
}

doTest("PUB-169 Test: Apply scale over columns", test.colapply)