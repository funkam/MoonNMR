## code to prepare `examples` dataset goes here
experiments<-read.csv("C:/Users/funkalexa/Documents/MoonShiny/local/MoonNMR/example_files/Format_Experiments.csv")
usethis::use_data(experiments, overwrite = TRUE)


templateuploader<-read.csv("C:/Users/funkalexa/Documents/MoonShiny/local/MoonNMR/example_files/Format_for_templateuploader.csv")
usethis::use_data(templateuploader, overwrite = TRUE)


