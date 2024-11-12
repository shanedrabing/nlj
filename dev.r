# PACKAGING

# 0) Set working directory
setwd("/Users/Shane/Code/packaging/nlj")

# 1) Create framework
# devtools::create("/Users/Shane/Code/packaging/nlj")

# 2) usethis (run once)
# usethis::use_readme_rmd()
# usethis::use_mit_license()
# usethis::use_vignette("nlj_guide")
# usethis::use_testthat()

# 3) Documentation
devtools::document()
devtools::build_readme()
# devtools::build_vignettes()

# 4) check code
devtools::check()
# devtools::check_man()
# devtools::spell_check()
# devtools::run_examples()

# test code
# devtools::unload()
# devtools::load_all()
# devtools::test()

# build
# devtools::build()


# TESTING


{
    # reset    
    tryCatch(devtools::unload(), error = function(e) {})
    devtools::load_all()
}

{
    detach(mtcars)
    attach(mtcars)
}

formula <- qsec ~ hp * mpg * wt

m <- lm(qsec ~ hp * mpg * wt, mtcars)
gat <- nlj::lm.gat(qsec ~ hp * mpg * wt, mtcars,
                   iterations = 6, penalty = 1e-9, verbose = TRUE)

summary(m)
summary(gat$fit)

df <- data.frame(qsec = mtcars$qsec,
                 simple = m$fitted.values,
                 nonlin = gat$z)

round(df, 2)
sum(abs(log(df$qsec / df$simple)))
sum(abs(log(df$qsec / df$nonlin)))
