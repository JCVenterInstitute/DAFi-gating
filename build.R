devtools::document()
devtools::check()
devtools::build()


devtools::install_github("JCVenterInstitute/DAFi-gating", build_vignettes = TRUE)

library(DAFi)
SampleData <- system.file("extdata", "sample.fcs",package="DAFi")
inclusionConfig <- system.file("extdata", "inclusion.config",package="DAFi")
exclusionConfig <- system.file("extdata", "exclusion.config",package="DAFi")
