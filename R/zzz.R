.onUnload <- function(libpath)
    library.dynam.unload("DECIPHER", libpath)

dbIsValid <- function(...)
{
    if ("dbIsValid" %in% getNamespace)
        RSQLite::dbIsValid(...)
    else
        RSQLite::isIdCurrent(...)
}