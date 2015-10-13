Disambiguate <- function(myXStringSet) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		stop("myXStringSet must be a DNAStringSet or RNAStringSet."))
	
	myXStringSet <- .Call("expandAmbiguities",
		myXStringSet,
		ifelse(type==1,
			"T",
			"U"),
		PACKAGE="DECIPHER")
	
	if (type==1) {
		myXStringSetList <- DNAStringSetList(myXStringSet)
	} else {
		myXStringSetList <- RNAStringSetList(myXStringSet)
	}
	
	return(myXStringSetList)
}
