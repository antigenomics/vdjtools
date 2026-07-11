# inext_oracle.R — regenerate the iNEXT (v3.0.2) oracle constants used by
# tests/python/test_inext.py. Requires the iNEXT package. Point .libPaths at a
# library that has iNEXT 3.0.2 installed, then run:
#   Rscript tests/assets/inext_oracle.R
# The Python test does NOT call R; it reads the committed spider_*.txt vectors
# and the hardcoded constants this script prints.
#
# 2026-07-11
suppressMessages(library(iNEXT))
data(spider)
options(digits = 10)

for (nm in c("Girdled", "Logged")) {
  x <- spider[[nm]]
  cat("==== ", nm, " ====\n", sep = "")
  print(DataInfo(x))
  cat("-- Asymptotic (AsyEst) --\n")
  out <- iNEXT(x, q = c(0, 1, 2), datatype = "abundance",
               size = c(100, 168, 200, 336), se = FALSE)
  print(out$AsyEst)
  cat("-- Size-based iNEXT --\n")
  print(out$iNextEst$size_based)
  cat("-- estimateD size --\n")
  print(estimateD(x, datatype = "abundance", base = "size",
                  level = c(100, 168, 200, 336), conf = NULL))
  cat("-- estimateD coverage C=0.95 --\n")
  print(estimateD(x, datatype = "abundance", base = "coverage",
                  level = 0.95, conf = NULL))
}
