## code to prepare `archive` dataset goes here
archive_backup<-read.csv("C:/Users/funkalexa/Documents/MoonShiny/local/MoonNMR/data-raw/archive/archive.csv")
usethis::use_data(archive, overwrite = TRUE)
usethis::use_data(archive_backup, overwrite = TRUE)

