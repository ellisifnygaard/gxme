#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP CHROMOSOME DIRECTORY WITH ONE SUBDIRECTORY PER CHROMOSOME BY RENDERING
# THE RMD REPORT BELONGING TO THIS PHASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message( "Commencing phase 0002...\n")


rmd_file_name <- "0002_set_up_chromosome_directory.Rmd"

# Import the arguments again just in case:
args <- source_arguments( file_path = here( base::basename(arguments_path) ) )


# Create Reports/0002 subdirectory if it does not already exist:
if( dir.exists( here("Reports", "0002") ) == FALSE ){
  status <- dir.create( here("Reports", "0002") )
  stopifnot( "Creating Reports/0002 subdirectory failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Reports", "0002"), ".\n")
} else{
  message( file.path("Reports", "0002") , " already exists.\n")
}


# REMOVE ANY HTML REPORTS ALREADY PRESENT IN REPORTS/0002

# Remove html reports left in /Reports/0002 from previous runs (should not
# occur, but as a safety measure)
if( file.exists( here(
  "Reports", "0002", paste0( tools::file_path_sans_ext(rmd_file_name), ".html" ) 
) ) ){
  
  # Remove the html reports left in /Reports/0002 from previous runs:
  message( "Removing existing versions of "
           , tools::file_path_sans_ext(rmd_file_name)
           , ".html\n" )
  
  file.remove( here( "Reports"
                     , "0002"
                     , paste0( tools::file_path_sans_ext(rmd_file_name)
                               , ".html" ) )
  )
}



# Load rmarkdown library
library(rmarkdown) # Not sure if this will work in parallel render scripts...

# Make temporary directory to serve as intermediate files directory:
tf <- tempfile()
dir.create(tf)
#jjj Are the 2 lines above okay or should tf be created in a specific directory?
# (It is used as intermediates_dir in rmarkdown::render below.)

message("Rendering 0002 report...\n")

render_output <- rmarkdown::render(
  # Specify the path to the .Rmd file we wish to render:
  input = here( "R", rmd_file_name )
  # The working directory inside the Rmd file:
  # , knit_root_dir = ""
  # Intermediate directory:
  , intermediates_dir = tf
  # Parameters:
  , params =
    args[ c( "chr_numbers"
             , "ewas_fileset_dir"
             , "ewas_fileset_name_base_regex"
             , "ewas_map_file_dir"
             , "ewas_map_file_name_base_regex"
             , "gwas_fileset_dir"
             , "gwas_fileset_name_base_regex"
             # , "summarising_function_path"
             # , "stratifying_function_path"
             , "mqtl_file_name_base_regex"
             , "mqtl_file_dir"
             , "ewas_annotation_file_path"
             , "plink_exe_path"
    ) ]
  # File name of html report:
  , output_file = paste0( tools::file_path_sans_ext(rmd_file_name), ".html" )
  # Where to place the html report:
  , output_dir = here( "Reports", "0002" )
  , quiet = TRUE
  # Ensure that script only uses objects defined inside the Rmd file:
  , envir = new.env()
)
# Delete temporary directory
unlink(tf, force = TRUE)



# render_output should now be the path to the report that was produced.
if( file.exists( render_output ) ){
  message(
    "Report successfully rendered and written to:\n", render_output, "\n"
  )
} else{
  stop( "No html report written to ", here( "Reports"), "!\n"
           , "Do not proceed unless you are certain that the code in "
           , rmd_file_name
           , " was run successfully (i.e. that the report was not rendered "
           , "successfully due to something unrelated, such as a pandoc error"
           , " for example).\n"
  )
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOG ARGUMENTS USED IN THIS PHASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message( "Logging arguments...\n" )
log_arguments( args_list = args
               , output_path =
                 here( "Reports", "Argument_logs", "0002_args.csv") )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE VARIABLES IN GLOBAL ENV THAT WE NO LONGER NEED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm( rmd_file_name
   , render_output
)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# AT THIS STAGE, THE USER SHOULD INSPECT THE GENERATED HTML REPORT, AS WELL AS
# THE CONTENTS OF THE CHROMOSOME-SPECIFIC FOLDERS IN pipeline_dir/chromosome
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE SET-UP PHASE IS NOW COMPLETE AND WE CAN START RUNNING THE STAGES OF THE
# PIPELINE, STARTING WITH 0101.

message( "Phase 0002 complete!\n.\n.\n.\n")