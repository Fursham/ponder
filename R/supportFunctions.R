

## Below are functions to write and possibly print the warnings and information
#   There are packages out there that can do this professionally, but I want the ability
#   for user to decide whether to print these info into console
#   this is where the quiet argument comes in

stopLog <- function(text) {
  
  message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                sprintf("[ERROR] %s", text)))
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

infoLog <- function(text, file, quiet = FALSE) {
  
  message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                sprintf("[INFO] %s", text)))
  
  Sys.sleep(0.3)
}

warnLog <- function(text) {
  
  message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                sprintf("[WARN] %s", text)))
  
  Sys.sleep(0.3)
}

unpack <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}
