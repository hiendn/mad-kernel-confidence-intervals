sha256 <- function(paths) {
 vapply(paths,function(p) {
  result<-system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE)
  if(!is.null(attr(result,'status'))||length(result)!=1L)stop('Cannot hash ',p)
  substr(result,1L,64L)
 },character(1))
}
write_json <- function(x,path)jsonlite::write_json(x,path,pretty=TRUE,auto_unbox=TRUE,digits=17,null='null')
