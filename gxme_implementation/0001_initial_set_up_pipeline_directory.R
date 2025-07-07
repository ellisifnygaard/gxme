# setwd("O:/Documents")

message( "Commencing phase 0001...\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1: SET UP PIPELINE DIRECTORY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine path of the pipeline directory:
pipeline_dir <- file.path( pipeline_dir_parent, pipeline_dir_name ) 

# Create pipeline directory if it does not already exist. If it already exists,
# throw error.
if( dir.exists( pipeline_dir ) == FALSE ){
  status <- dir.create( pipeline_dir )
  stopifnot( "Creating pipeline_dir failed" = status == TRUE )
} else{
  stop( "There is already a directory at \""
        , pipeline_dir
        , "\".\nPlease specify a different path/directory name, move the "
        , "existing directory to another location, or rename the existing "
        , "directory." )
}

# Set the working directory to pipeline_dir:
setwd( pipeline_dir )
stopifnot( identical( getwd(), pipeline_dir ) )


# Write a .here file to pipeline_dir so that the here package recognises it as
# the root directory
if( file.exists( file.path( pipeline_dir, ".here")) == FALSE ){
  status <- file.create( file.path( pipeline_dir, ".here") )
  stopifnot( "Creating .here file failed" = status == TRUE )
}

here::i_am( ".here" )
# here::here()
library(here)
# Note: the user may need to tweak this script according to their particular
# local circumstances. 

# Check that here recognises the pipeline directory as the root directory:
here_dir <- here::here()
stopifnot(
  "here_dir is not equal to pipeline_dir!  \n
           *** DO NOT PROCEED BEFORE THIS HAS BEEN FIXED! ***" =
    here_dir == pipeline_dir 
)
message( "`here()` recognises ", here_dir, " as root directory.\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2: SET UP R SUBDIRECTORY 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove any punctuation characters (!"#%&'()*+,-./ etc.) from the end of the
# provided path (in case it is a forward slash that would result in paths
# with "//")
if( grepl( "[[:punct:]]"
           , substring( path_pipeline_scripts, nchar(path_pipeline_scripts) ) )
){
  warning( "`path_pipeline_scripts` ends with a special character.\n"
           , "Removing the last character of the string.")
  # Remove last character from path:
  path_pipeline_scripts <- substr( path_pipeline_scripts
                                   , start = 1
                                   , stop = nchar(path_pipeline_scripts) - 1
  )
}

# Check that path_pipeline_scripts exists:
stopifnot( "`path_pipeline_scripts` does not exist" =
             dir.exists( path_pipeline_scripts )
)

# List files in path_pipeline_scripts that we want to copy to pipeline_dir/R:
# (Only .R and .Rmd files with names starting with 4 digits or "f_" followed by
# an underscore)
path_pipeline_scripts_files <- list.files(
  path_pipeline_scripts 
  , pattern = paste0(
    "^f_.+\\.R|", paste0("^[[:digit:]]{4}_.+\\.", c("R", "Rmd"), collapse = "|")
  )
  , full.names = FALSE
)


# Check that path_pipeline_scripts contains files:
stopifnot( "`path_pipeline_scripts` contains fewer files than required" =
             length( path_pipeline_scripts_files ) > 30
           #xxxxx Update the integer later according to the number of scripts
           # needed to run the pipeline
)

# Create R subdirectory if it does not already exist:
if( dir.exists( here("R") ) == FALSE ){
  status <- dir.create( here("R") )
  stopifnot( "Creating R subdirectory failed" = status == TRUE )
  message( "R subdirectory created at ", here("R"), ".\n")
}

# Copy pipeline scripts to R subdirectory:
copy_status <- file.copy( 
  # Vector with paths to the files we will copy:
  from = file.path( path_pipeline_scripts, path_pipeline_scripts_files)
  # Vector with paths to the copies in the .R subdirectory:
  , to = here( "R", path_pipeline_scripts_files )
  , overwrite = TRUE
  # Preserve file date if possible:
  , copy.date = TRUE 
)

# Stop if copying files to subdirectory was not successful:
stopifnot( all( copy_status == TRUE) )
message("All scripts successfully copied from "
        , path_pipeline_scripts 
        , "\nto "
        , here( "R" )
        , "\n"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3: SET UP REPORT SUBDIRECTORY 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create Reports subdirectory if it does not already exist:
if( dir.exists( here("Reports") ) == FALSE ){
  status <- dir.create( here("Reports") )
  stopifnot( "Creating Reports subdirectory failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Reports"), ".\n")
}

# Create subfolder within /Reports:
if( dir.exists( here("Reports", "Argument_logs") ) == FALSE ){
  status <- dir.create( here("Reports", "Argument_logs") )
  stopifnot( "Creating folder within Reports subdirectory failed" = 
               status == TRUE )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4: SET UP RESULT SUBDIRECTORY 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create Reports subdirectory if it does not already exist:
if( dir.exists( here("Results") ) == FALSE ){
  status <- dir.create( here("Results") )
  stopifnot( "Creating Results subdirectory failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Results"), ".\n")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5: MOVE COPY OF ARGUMENT FILE TO PIPELINE DIRECTORY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check that the file in arguments_path is a R file:
stopifnot( "`arguments_path` does not point to an R file" =
             grepl( "\\.R$", arguments_path )
)

# Check that there exists a file in arguments_path:
stopifnot( "`arguments_path` does not exist" = file.exists( arguments_path ) )

# Copy argument file to pipeline_dir:
copy_status <- file.copy( 
  # Paths to the file we will copy:
  from = arguments_path
  # Path to where the copy will be moved to:
  , to = here( base::basename(arguments_path) )
  , overwrite = TRUE
  # Preserve file date if possible:
  , copy.date = TRUE 
)

# Stop if copying files to subdirectory was not successful:
stopifnot( all( copy_status == TRUE) )
message("R file  successfully copied from "
        , path_pipeline_scripts 
        , "\nto "
        , here( "R" )
        , "\n"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6: LOAD ALL ARGUMENTS SPECIFIED IN R FILE - CHECKS AND LOG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The arguments file should now be located in pipeline_dir

# Source the source_arguments function from R/f_source_arguments.R
if( file.exists( here("R", "f_source_arguments.R") ) ){
  source( here("R", "f_source_arguments.R") )
} else{ stop(  here("R", "f_source_arguments.R"), " does not exist!" )}

# Source the log_arguments function from R/f_log_arguments.R
if( file.exists( here("R", "f_log_arguments.R") ) ){
  source( here("R", "f_log_arguments.R") )
} else{ stop(  here("R", "f_log_arguments.R"), " does not exist!" )}


# Import the arguments:
args <- source_arguments( file_path = here( base::basename(arguments_path) ) )

# Check that .R file with arguments was successfully loaded:
stopifnot( is.list(args) & length(args) > 16 ) 
#xxxxx Update accordingly wrt the number of arguments

# Check args for a couple of arguments:
check_args <- c( "chr_numbers"
                 ,"ewas_fileset_dir"
                 ,"ewas_fileset_name_base_regex"
                 ,"ewas_map_file_dir"
                 ,"ewas_map_file_name_base_regex"
                 ,"ewas_annotation_file_path"
                 ,"gwas_fileset_dir"
                 ,"gwas_fileset_name_base_regex"
                 ,"summarising_function_path"
                 ,"stratifying_function_path"
                 ,"mqtl_file_name_base_regex"
                 ,"mqtl_file_dir"
                 ,"plink_exe_path" )

if( all( check_args %in% names(args) ) == FALSE ){
  stop( here( base::basename(arguments_path) )
              , " is missing one or more arguments." )
}


# Check that chr_numbers are whole numbers:

# Source the is_whole_number function from R/f_is_whole_number.R
if( file.exists( here("R", "f_is_whole_number.R") ) ){
  source( here("R", "f_is_whole_number.R") )
} else{ stop(  here("R", "f_is_whole_number.R"), " does not exist!" )}

stopifnot(
  "`chr_numbers` must consist exclusively of numerics that are whole numbers" =
    all( sapply( args$chr_numbers, is.numeric ) ) &
    all( sapply( args$chr_numbers, is_whole_number ) )
)


# Log arguments:
message( "Logging arguments...\n" )
log_arguments( args_list = args
               , output_path =
                 here( "Reports", "Argument_logs", "0001_args.csv") )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7: COPY SUMMARISING AND STRATIFYING FUNCTIONS TO /R:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# (The paths to these files must be provided as arguments.)

summ_and_strat_function_paths <- c(
  args[["summarising_function_path"]]
  , args[["stratifying_function_path"]]
)

# Stop if any of these files do not exist:
if( any( file.exists( summ_and_strat_function_paths ) == FALSE ) ){
  stop( 
    "Cannot locate these files:  \n\n"
    , paste( summ_and_strat_function_paths[ 
      which( file.exists( summ_and_strat_function_paths ) )
    ], collapse = "\n\n" )
    , "\n\nPlease review the paths in `summarising_function_path` and/or "
    , "`stratifying_function_path`."
  )
}

# Copy files with summarising/stratifying functions to R subdirectory:
copy_status <- file.copy( 
  # Vector with paths to the files we will copy:
  from = summ_and_strat_function_paths
  # Vector with paths to the copies in the .R subdirectory:
  , to = here( "R", basename(summ_and_strat_function_paths) )
  , overwrite = TRUE
  # Preserve file date if possible:
  , copy.date = TRUE 
)

# Stop if copying files to subdirectory was not successful:
stopifnot( all( copy_status == TRUE) )
message( "These files:\n\n  -"
         , paste( summ_and_strat_function_paths, collapse = ",\n  - ")
         , ",\n\nwere succesfully copied to "
         , here( "R" )
         , "\n"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8: IF PERFORMING REPLICATION, COPY REPLICATION MODULE FILES TO R SUBDIRECTORY 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If performing replication analysis, copy R scripts from Replication folder to
# the pipeline_dir/R subdirectory:
if( args[["replication_analysis"]] == TRUE ){
  
  # Stop if /Replication folder does not exist:
  stopifnot( dir.exists( file.path(path_pipeline_scripts, "Replication") ) )
  
  # List the replication module-specific files to be copied to pipeline_dir/R:
  replication_module_files <- list.files(
    file.path(path_pipeline_scripts, "Replication")
    , pattern = paste0(
      "^f_.+\\.R|"
      , paste0("^[[:digit:]]{4}R_.+\\."
               , c("R", "Rmd"), collapse = "|" )
    )
    , full.names = FALSE )                                
  
  # Check that path_pipeline_scripts/Replication contains files:
  stopifnot(
    "`path_pipeline_scripts/Replication` contains fewer files than required" =
      length( replication_module_files ) > 1
  )
  
  # Copy pipeline scripts to R subdirectory:
  copy_status <- file.copy( 
    # Vector with paths to the files we will copy:
    from = file.path( path_pipeline_scripts
                      , "Replication"
                      , replication_module_files )
    # Vector with paths to the copies in the .R subdirectory:
    , to = here( "R", replication_module_files )
    , overwrite = TRUE
    # Preserve file date if possible:
    , copy.date = TRUE 
  )
  
  # Stop if copying files to subdirectory was not successful:
  stopifnot( all( copy_status == TRUE) )
  message("All replication module scripts successfully copied from\n"
          , file.path(path_pipeline_scripts, "Replication") 
          , "\nto "
          , here( "R" )
          , "\n"
  )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE VARIABLES IN GLOBAL ENV THAT WE NO LONGER NEED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(copy_status
   , path_pipeline_scripts
   , path_pipeline_scripts_files
   )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0001 IS NOW COMPLETE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message( "Phase 0001 complete!\n.\n.\n.\n")
