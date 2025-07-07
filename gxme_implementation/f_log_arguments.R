# Input:  - args, a list with arguments
#         - output_path, the path to write the output to
# This function creates a data frame with two columns, where each row contains
# the name of the argument and the argument value. Then it writes the data frame
# to the path provided in `output_path`.

log_arguments <- function( args_list, output_path ){
  
  # args must be list:
  stopifnot( all( class(args_list) %in% "list" ) )
  
  # Make a data frame containing the argument name and value:
  args_df <- data.frame( arg_name = character(0), arg_value = character(0) )
  for( i in 1:length(args_list) ){
    # Add argument name to row i in df:
    args_df[ i, 1] <- names(args_list[i])
    # Add argument value to row i in df:
    args_df[ i, 2] <- as.character(args_list[i])
  }
  
  # Write the data frame to /Reports/Argument_log to log what the specified
  # arguments were at this particular stage
  # (Useful should you find yourself in a situation where you suspect that you
  # might perhaps have altered the arguments during the pipeline)
  
  readr::write_excel_csv2( args_df 
                           , file = output_path
                           , append = FALSE # overwrite
                           , col_names = TRUE)
  # Throw error if writing the file was not successful:
  # (the readr function will probably throw an error if this happens, but just
  # in case)
  status <- file.exists( output_path )
  if( !(status %in% TRUE) ){
    stop( "Writing the argument log to ", output_path, " failed." )
  }
}
