#' Create an archive template
#' 
#' A simple script that creates an empty archive (+archive_backup) template in the current directory
#'
#' @param archive consists of the columns Date,Name,Projet,Type,Size
#' See `?golem::get_golem_options` for more details.
#' @export

moonnmr_archivecreator<-function(){

  if(file.exists("archive.csv")==FALSE){
empty_archive<-data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Date","Name","Project","Type","Size"))))
write.csv(empty_archive,"archive.csv",row.names=FALSE)
write.csv(empty_archive,"archive_backup.csv",row.names=FALSE)
stop("Archives created!")
} else {
  stop("Archive already exists in current directory!")
}
}

