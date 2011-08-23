##################################################
# Miscellaneous useful functions
# Jonathan Whiteley		R v2.12		2010-01-26
##################################################

strip_empty_dims  <- function( data = NULL, dim = c(1, 2), 
							  rows = NULL, cols = NULL, col.class = NULL ) {
  ## strip whole rows or columns if they are entirely empty (NA or NaN)
  ## `data` is the data frame to strip empty values from.
  ## `rows` is a vector of indices or logical values of 
  ##        which rows to check for empty values (default all)
  ## `cols` see `rows` above (but applied to columns)
  ## `dim`  the dimensions to strip if all empty.  1 = rows, 2 = cols
  ## `col.class` is a vector of classes: columns that match will be included
  ##             This applies to the `cols` argument.

  if (!is.data.frame(data) && !is.matrix(data)) stop(
	  "The first argument must be a \"data.frame\" or \"matrix\"."
     )

  cols <- as.vector(cols)
  if ( is.null(cols) ) {
	cols <-  1:ncol(data)
  }
  if ( !is.null(col.class) ) {
	check.cols <- which( lapply(data[, cols], class) %in% as.vector(col.class) )
	check.cols <- cols[check.cols]
  } else {
	check.cols <- cols
  }
  rows <- as.vector(rows)
  if ( is.null(rows) ) {
	rows <-  1:nrow(data)
  }  
  
  empty.rows <- which( apply( data[, check.cols], 1, function(x) all(is.na(x)) ) )
  empty.cols <- which( apply( data[rows,       ], 2, function(x) all(is.na(x)) ) )

  if (all(dim == 1)) {
	if (length(empty.rows) > 0) data <- data[-empty.rows, ]
  } else if (all(dim == 2))  {
	if (length(empty.cols) > 0) data <- data[, -empty.cols]
  } else if (all(dim %in% c(1, 2))) {
	if (length(empty.rows) > 0 && length(empty.cols) > 0) 
      data <- data[-empty.rows, -empty.cols]
  }
  
  return(data)
  
  if (FALSE) {  # Testing
	test.df <-  data.frame( )
	df.proc <- strip_empty_dims(test.df)
	df.proc
  }
}


