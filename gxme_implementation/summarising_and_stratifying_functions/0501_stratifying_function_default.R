# Matrix level stratifying function - DEFAULT (i.e. tertiles)


# INPUT: 
#           1) a 1 x n matrix, where 
#               1 = number of rows; where the row name = SL ID,
#               n = number of columns; where the column names represent the EWAS
#               family member IDs, 
#               and the cells are the individuals' summarised measurements that 
#               were calculated using the summarising function.
#               Must have column names = EWAS family member IDs
#               Must have row names = EWAS probe IDs


# OUTPUT: an 1 x n matrix, where 
#             1 = the number of rows; where the row name = SL ID,
#             n = columns; where the column names represent the EWAS
#             family member IDs,
#             and each cell contains the number of the stratum that the 
#             individual/family has been allocated to based on their summarised 
#             measurement for this particular state locus (SL).
#             Row names = SL ID
#             Column names = EWAS family member IDs

# mat <- sl_summarised_matrix
# sl_id <- sl

stratifying_function <- function( mat ){
  
  # First, find the breakpoints for the tertiles
  
  # Use the summarised values to compute the breaks:
  breakpoints <- quantile( mat, probs = c( 0:3 / 3 ) )
  
  # Then, use the breakpoints to divide the individuals into tertiles
  
  mat <- rbind(
    mat
    , cut(
      # The summarised values to be cut into a factor/strata:
      x = mat[ 1, , drop = TRUE ] 
      # Unique cut points:
      , breaks = breakpoints 
      # Whether to include any values equal to the lowest value in 'breaks':
      , include.lowest = TRUE
      # Make interval closed on the right and open on the left:
      , right = TRUE
      # Intervals will be like this: [min, T1], (T1 ,
      # T2], (T2 , max]
      , labels = FALSE
      # return integer codes (1,2,..) instead of factor
    )
  )
  
  # Stop if there are any NAs in the row with strata integers:
  stopifnot( any( is.na( mat[2, ] )) == FALSE )
  
  # Check if all indivs with strata = 1 have summarised values >= 0 percentile
  # and <= 33.3333 percentile:
  stopifnot( 
    "Individuals placed in top tertile do not have summarised values to match" =
      all( mat[ 1, which( mat[2,] == 1 )] >= min( mat[ 1, ] ) &
             mat[ 1, which( mat[2,] == 1 )] <= breakpoints["33.33333%"] )
  )
  # Check if all indivs with strata = 2 have summarised values > 33.33333
  # percentile and <= 66.66667 percentile:
  stopifnot( 
    "Individuals placed in top tertile do not have summarised values to match" =
      all( mat[ 1, which( mat[2,] == 2 )] > breakpoints["33.33333%"] &
             mat[ 1, which( mat[2,] == 2 )] <= breakpoints["66.66667%"] )
  )
  # Check if all indivs with strata = 3 have summarised values > 66.66667
  # percentile and <= 100 percentile:
  stopifnot( 
    "Individuals placed in top tertile do not have summarised values to match" =
      all( mat[ 1, which( mat[2,] == 3 )] > breakpoints["66.66667%"] &
             mat[ 1, which( mat[2,] == 3 )] <= max( mat[ 1, ] ) )
  )
  
  # We're now satisfied that the newly added row with strata is correct.
  # We will now remove the row with summarised values and name the strata row
  # after the SL ID.
  
  sl_id <- rownames(mat)[1]
  
  # Remove row with summarised values:
  mat <- mat[2, , drop = FALSE ]
  
  # Set SL ID as row name in row with strata numbers:
  rownames(mat)[1] <- sl_id
  
  # Return a matrix consisting of the newly added strata row in mat
  return( mat )
}
