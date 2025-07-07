# Matrix level summarising function - MEDIAN


# INPUT: 
#           1) string with SL ID
#           2) an m x n matrix, where 
#               m = rows; representing probe IDs,
#               n = columns; representing EWAS family member IDs, 
#               and the cells are the individuals' observed measurements for the 
#               specific probes.
#               Must have column names = EWAS family member IDs
#               Must have row names = EWAS probe IDs


# OUTPUT: an 1 x n matrix, where 
#             1 = the number of row; representing the SL,
#             n = columns; representing EWAS family member IDs, 
#             and each cell contains an individual's summarised measurement for 
#             this particular state locus (SL).
#             Row names = SL ID
#             Column names = EWAS family member IDs

# mat <- sl_probe_matrix
# sl_id <- sl

summarising_function <- function( mat, sl_id ){
  
  summarised_mat_vector <- apply(mat, 2, median, na.rm = TRUE) # ignore NAs
  
  # Stop if summarised_mat_vector is not a 1 x m vector:
  stopifnot( is.vector(summarised_mat_vector) &
               length( summarised_mat_vector ) == ncol(mat) )
  
  
  # Add the vector with column means to mat for validation check:
  mat <- rbind( mat, summarised_mat_vector)
  # Set sl_id as row name
  rownames(mat)[nrow(mat)] <- sl_id
  
  # Test if the vector of column means was correctly paired with columns in mat.
  # Pick 5 random columns to check:
  for( k  in sample(1:ncol(mat), 5, replace = FALSE) ){
    stopifnot( 
      median( mat[ -nrow(mat), k] ) == mat[ nrow(mat), k]
    )
    # The sum of the probe rows in mat in column k,
    # divided by the number of probe rows in mat, 
    # should be equal to the column mean found in the newly added bottom row
  }
 
  # Return a matrix consisting of the newly added bottom row in mat
  return( mat[ nrow(mat), ] )
}