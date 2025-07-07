# Function that sources a file, then assigns the objects in the file to a
# temporary environment, and subsequently returns a list with those objects.
# N.B.! It only works as intended if the source file only contains one list. 
# Returns a list with the elements from the list in the sourced file 

# Based on the function "my_source" provided in
# https://stackoverflow.com/a/58880544

source_arguments <- function( file_path, local = NULL){
  stopifnot( "The given file path does not exist!" = file.exists(file_path) )
  # Create temporary environment:
  tmp <- new.env( parent = parent.frame() )
  
  # Assign objects in sourced file to temp env:
  source( file_path, local = tmp ) 
  
  # Get the names of the objects:
  objs <- names(tmp)
  
  # Stop if there is more than one object in the provided argument file:
  stopifnot( "The provided argument file contains multiple objects" =
               length( objs ) == 1 )
  
  # Stop unless the object in the provided argument file is a list:
  stopifnot( "The provided argument file must contain one list" =
               any( class( get(objs, envir = tmp) ) %in% "list" ) )
  
  # # Assign the objects in the temporary environment to the global environment:
  # for( x in names(tmp) ){
  #   assign( x, tmp[[x]], envir = parent.frame() )
  # }
  
  # Return the objects inside objs (this could be the elements of the params
  # list from the argument file for example):
  get(objs, envir = tmp)
}
