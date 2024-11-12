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
