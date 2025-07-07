# Returns TRUE if x is a numeric/integer which also is a whole number
# If x is a vector of numerics/integers, it returns a vector, where element i is
# TRUE if x[i] is a whole number, and FALSE if not.
# (Whole numbers = {0,1,2,3,....}; they do not include negative numbers)
is_whole_number <- function(x){
  # Must be a numeric or an integer:
  stopifnot(
    "is_whole_number() only supports numerics and integers" =
      any ( class(x) %in% c("numeric", "integer") )
  )
  # Must be single number or a vector:
  stopifnot(
    "is_whole_number() does not support vectors of length zero" =
      length(x) >= 1
  )
  res <- suppressWarnings( as.numeric( x[[1]]) %% 1 == 0 &
                             as.numeric( x[[1]]) >= 0 )
  res <- suppressWarnings( as.numeric( x ) %% 1 == 0 &
                             as.numeric( x ) >= 0 )
  if( any( is.na(res) ) ) res[ which( is.na(res) ) ] <- FALSE
  return( res )
} 
