.onUnload <- function(libpath)
	library.dynam.unload("DECIPHER", libpath)

dbIsValid <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		do.call("dbIsValid", list(...), envir=getNamespace("RSQLite"))
	} else {
		do.call("isIdCurrent", list(...), envir=getNamespace("RSQLite"))
	}
}

dbBegin <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		do.call("dbBegin", list(...), envir=getNamespace("RSQLite"))
	} else {
		do.call("dbBeginTransaction", list(...), envir=getNamespace("RSQLite"))
	}
}
