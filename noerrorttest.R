noerrorttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(list(p.value=NA)) else return(obj)
}
