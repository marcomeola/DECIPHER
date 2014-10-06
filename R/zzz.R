.onUnload <- function(libpath)
	library.dynam.unload("DECIPHER", libpath)

dbIsValid <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		DBI::dbIsValid(...)
	} else {
		RSQLite::isIdCurrent(...)
	}
}

dbBegin <- function(...) {
	if (packageVersion("RSQLite") >= package_version("1.0.0")) {
		RSQLite::dbBegin(...)
	} else {
		RSQLite::dbBeginTransaction(...)
	}
}
