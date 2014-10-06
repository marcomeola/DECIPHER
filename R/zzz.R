.onUnload <- function(libpath)
	library.dynam.unload("DECIPHER", libpath)

dbIsValid <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		do.call("dbIsValid", as.list(...), envir=getNamespace("RSQLite"))
	} else {
		do.call("isIdCurrent", as.list(...), envir=getNamespace("RSQLite"))
	}
}

dbBegin <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		do.call("dbBegin", as.list(...), envir=getNamespace("RSQLite"))
	} else {
		do.call("dbBeginTransaction", as.list(...), envir=getNamespace("RSQLite"))
	}
}
