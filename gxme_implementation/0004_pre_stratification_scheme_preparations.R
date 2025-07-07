#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# PHASE 0004
#
# PRE-STRATIFICATION SCHEME PREPARATIONS OF GWAS AND EWAS DATA
#
# RUN THE FOLLOWING STAGES: 
#       - 0201, 0202, (0203 if strict mQTL filtering procedure)
#       - 0301
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message( "Commencing phase 0004...\n")
require(here)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE THE FUNCTION THAT WILL RENDER THE REPORTS AND THEIR ARGUMENTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (The code for rendering Rmds in parallel is based on the gist from
# hrbrmstr/do_rpt.r on Github)

# DEFINE FUNCTION TO BE USED WHEN RENDERING THE REPORTS AT EACH SEPARATE STAGE:
render_chr_report <- function(chr){
  # Set wd to pipeline dir, otherwise the here function will use another root
  # directory:
  setwd(pipeline_dir) 
  require(here)
  # Make temporary directory to serve as intermediate files directory:
  tf <- tempfile()
  dir.create(tf)
  rmarkdown::render(
    input = here( "R", rmd_file_name )
    # The working directory inside the Rmd file:
    , knit_root_dir = chr$params$root_dir
    # Intermediate directory:
    , intermediates_dir = tf
    # Parameters:
    , params = chr$params
    # File name of html report:
    , output_file = chr$out
    # Where to place the html report:
    , output_dir = here( "Reports", stage )
    , quiet = TRUE
    # , clean = TRUE
    # Ensure that script only used objects defined inside the Rmd file:
    , envir = new.env()
  )
  # Delete temporary directory
  unlink(tf, force = TRUE)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stage 0201
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# n_cores <- parallel::detectCores() - 2 # (set in 0000_run_entire_pipeline.R)
rmd_file_name <- "0201_preliminary_qc_gwas_data.Rmd"
stage <- "0201"

message( "\n-------------------------\n"
         , "Performing Stage ", stage, "..."
         , "\n-------------------------\n"
)

# Check that we're still in the right root:
stopifnot( here() == pipeline_dir )

# Import the arguments again just in case:
args <- source_arguments( file_path = here( base::basename(arguments_path) ) )


# Create Reports/stage subdirectory if it does not already exist:
if( dir.exists( here("Reports", stage) ) == FALSE ){
  status <- dir.create( here("Reports", stage) )
  stopifnot( "Creating stage subdirectory in /Reports/ failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Reports", stage), ".\n")
}

# Load rmarkdown library
library(rmarkdown) 

# Delete any html reports left in /Reports/stage from previous runs:
files_in_Reports_stage <- list.files( here("Reports", stage)
                                , full.names = TRUE )

if( length(files_in_Reports_stage) > 0 ){
  message( "Deleting files that were already present in "
           , file.path("Reports", stage)
           , "...\n" )
  lapply( files_in_Reports_stage, file.remove )
}

# Stop if there are still html reports left in Reports/stage:
if( length( list.files( here("Reports", stage) ) ) > 0 ){
  stop(
    "There are still html reports left in "
    , here("Reports", stage)
    , "!\nRemove them manually and try again."
  )
}
  


# MAKE LIST WITH ARGUMENTS/PARAMETERS

# Create empty list to store arguments in - one set per chromosome
chr_reports <- list()

# For every chromosome, add chromosome-specific list element with arguments and
# other parameters:
# for( i in as.character( c(1:3, 9, 21) ) ){
for( i in as.character( args$chr_numbers ) ){
  chr_reports[[i]] <- list( 
    # Create file name for the chromosome-specific html report:
    out = gsub("\\.Rmd", sprintf("_chr%02d.html", as.integer(i) ), rmd_file_name)
    # (replace ".Rmd" in file name with _chr$$.html)
    , params = list( 
      chr_number = as.integer(i) 
      # Create root directory for the Rmd script (i.e. working directory for
      # inside the Rmd):
      , root_dir = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
      # PLINK-related:
      , plink_memory_mb = args[["plink_memory_mb"]]
      , plink_timeout = args[["plink_timeout"]]
      # QC-related parameter:
      , maf_minimum = args[["maf_minimum"]]
      # , stage_dir = "0201" # use the default
    ) 
  )
}


# RUN 0201 FOR ALL CHROMOSOMES, IN PARALLEL

library(doParallel)

# doParallel::registerDoParallel( cores = n_cores ) 

# Parallel rendering of chromosome-specific reports:
message( "Performing Stage ", stage, " and rendering reports...\n")
system.time({
  status <- foreach( chr = chr_reports, .combine = c ) %dopar% {
    tryCatch( 
      expr = {
        render_chr_report(chr)
        paste0("Chromosome ", chr$params$chr_number, " successfully rendered.")
      }
      , error = function( e ){
        paste0( "Chromosome ", chr$params$chr_number, " failed to render: "
                , e$message)
      } )
  }
  doParallel::stopImplicitCluster()
})

# STATUS MESSAGES AND CHECKS

# If rendering one or more chromosomes failed for some reason:
if( sum(grepl("successfully", status)) < length(chr_reports) ){
  # if( sum(status) < length(chr_reports) ){
  stop( "Rendering ", stage, " report failed for the following chromosomes:\n\n"
        # , paste( sort( sapply( chr_reports[ which( !( status %in% 1 ) ) ]
        , paste( status[ which( grepl("failed", status) ) ], collapse = "\n" )
  )
} else{ message("Rendering successful for all chromosomes.\n")}


message(
  "Checking if html reports were succesfully created for all chromosomes....\n"
)
report_status <- sapply( chr_reports, function(chr){
  file.exists( here( "Reports", stage, chr$out ) ) 
})

# report_status should now be the path to the report that was produced.
if( all( report_status ) ){
  message( "Stage "
           , stage
           , " reports successfully rendered and written to \n\'"
           , here( "Reports", stage )
           , "\'\nfor all chromosomes: \n\n" 
           , paste( 
             sapply( chr_reports[ which( report_status == TRUE ) ]
                     , function(x) x$params$chr_number )
             , collapse = ", " )
  )
} else{
  stop( "No html report written to "
        , here( "Reports")
        , " for the following chromosomes:\n  
           "
        , paste( 
          sapply( chr_reports[ which( report_status == FALSE ) ]
                  , function(x) x$params$chr_number )
          , collapse = ", ")
        , "\n\n**Do not** proceed unless you are certain that the code in \'"
        , rmd_file_name
        , "\' was run successfully for all chromosomes (i.e. that the report "
        , "was not rendered successfully due to something unrelated, such as"
        , " a pandoc error etc.).\n"
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stage 0202
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# n_cores <- parallel::detectCores() - 2 # (set in 0000_run_entire_pipeline.R)
rmd_file_name <- "0202_make_snp_list.Rmd"
stage <- "0202"

message( "\n-------------------------\n"
         , "Performing Stage ", stage, "..."
         , "\n-------------------------\n"
         )

# Check that we're still in the right root:
stopifnot( here() == pipeline_dir )

# Create Reports/stage subdirectory if it does not already exist:
if( dir.exists( here("Reports", stage) ) == FALSE ){
  status <- dir.create( here("Reports", stage) )
  stopifnot( "Creating stage subdirectory in /Reports/ failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Reports", stage), ".\n")
}

# Delete any html reports left in /Reports/stage from previous runs:
files_in_Reports_stage <- list.files( here("Reports", stage)
                                , full.names = TRUE )

if( length(files_in_Reports_stage) > 0 ){
  message( "Deleting files that were already present in "
           , file.path("Reports", stage)
           , "...\n" )
  lapply( files_in_Reports_stage, file.remove )
}

# Stop if there are still html reports left in Reports/stage:
if( length( list.files( here("Reports", stage) ) ) > 0 ){
  stop(
    "There are still html reports left in "
    , here("Reports", stage)
    , "!\nRemove them manually and try again."
  )
}



# MAKE LIST WITH ARGUMENTS/PARAMETERS

# Create empty list to store arguments in - one set per chromosome
chr_reports <- list()

# For every chromosome, add chromosome-specific list element with arguments and
# other parameters:
# for( i in as.character( c(1:3, 9, 21) ) ){
for( i in as.character( args$chr_numbers ) ){
  chr_reports[[i]] <- list( 
    # Create file name for the chromosome-specific html report:
    out = gsub("\\.Rmd", sprintf("_chr%02d.html", as.integer(i) ), rmd_file_name)
    # (replace ".Rmd" in file name with _chr$$.html)
    , params = list( 
      chr_number = as.integer(i) 
      # Create root directory for the Rmd script (i.e. working directory for
      # inside the Rmd):
      , root_dir = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
      # PLINK-related:
      , plink_memory_mb = args[["plink_memory_mb"]]
      , plink_timeout = args[["plink_timeout"]]
      # MQTL LD:
      , mqtl_ld_r2 = args[["mqtl_ld_r2"]]
      # REPLICATION ANALYSIS STATUS:
      , replication_analysis = args[["replication_analysis"]]
      # , stage_dir = "0202" # use the default
    ) 
  )
}


# RUN 0202 FOR ALL CHROMOSOMES, IN PARALLEL

# doParallel::registerDoParallel( cores = n_cores ) 

# Parallel rendering of chromosome-specific reports:
message( "Performing Stage ", stage, " and rendering reports...\n")
system.time({
  status <- foreach( chr = chr_reports, .combine = c ) %dopar% {
    tryCatch( 
      expr = {
        render_chr_report(chr)
        paste0("Chromosome ", chr$params$chr_number, " successfully rendered.")
      }
      , error = function( e ){
        paste0( "Chromosome ", chr$params$chr_number, " failed to render: "
                , e$message)
      } )
  }
  doParallel::stopImplicitCluster()
})


# STATUS MESSAGES AND CHECKS

# If rendering one or more chromosomes failed for some reason:
if( sum(grepl("successfully", status)) < length(chr_reports) ){
  # if( sum(status) < length(chr_reports) ){
  stop( "Rendering ", stage, " report failed for the following chromosomes:\n\n"
        # , paste( sort( sapply( chr_reports[ which( !( status %in% 1 ) ) ]
        , paste( status[ which( grepl("failed", status) ) ], collapse = "\n" )
  )
} else{ message("Rendering successful for all chromosomes.\n")}


message(
  "Checking if html reports were succesfully created for all chromosomes....\n"
)
report_status <- sapply( chr_reports, function(chr){
  file.exists( here( "Reports", stage, chr$out ) ) 
})

# report_status should now be the path to the report that was produced.
if( all( report_status ) ){
  message( "Stage "
           , stage
           , " reports successfully rendered and written to \n\'"
           , here( "Reports", stage )
           , "\'\nfor all chromosomes: \n\n" 
           , paste( 
             sapply( chr_reports[ which( report_status == TRUE ) ]
                     , function(x) x$params$chr_number )
             , collapse = ", " )
  )
} else{
  stop( "No html report written to "
        , here( "Reports")
        , " for the following chromosomes:\n  
           "
        , paste( 
          sapply( chr_reports[ which( report_status == FALSE ) ]
                  , function(x) x$params$chr_number )
          , collapse = ", ")
        , "\n\n**Do not** proceed unless you are certain that the code in \'"
        , rmd_file_name
        , "\' was run successfully for all chromosomes (i.e. that the report "
        , "was not rendered successfully due to something unrelated, such as"
        , " a pandoc error etc.).\n"
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stage 0301
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rmd_file_name <- "0301_remove_probes_on_chromosomes_without_snps.Rmd"
stage <- "0301"

message( "\n-------------------------\n"
         , "Performing Stage ", stage, "..."
         , "\n-------------------------\n"
)

# Check that we're still in the right root:
stopifnot( here() == pipeline_dir )

# Create Reports/stage subdirectory if it does not already exist:
if( dir.exists( here("Reports", stage) ) == FALSE ){
  status <- dir.create( here("Reports", stage) )
  stopifnot( "Creating stage subdirectory in /Reports/ failed" = status == TRUE )
  message( "Reports subdirectory created at ", here("Reports", stage), ".\n")
}

# Delete any html reports left in /Reports/stage from previous runs:
files_in_Reports_stage <- list.files( here("Reports", stage)
                                      , full.names = TRUE )

if( length(files_in_Reports_stage) > 0 ){
  message( "Deleting files that were already present in "
           , file.path("Reports", stage)
           , "...\n" )
  lapply( files_in_Reports_stage, file.remove )
}

# Stop if there are still html reports left in Reports/stage:
if( length( list.files( here("Reports", stage) ) ) > 0 ){
  stop(
    "There are still html reports left in "
    , here("Reports", stage)
    , "!\nRemove them manually and try again."
  )
}



# MAKE LIST WITH ARGUMENTS/PARAMETERS

# Create empty list to store arguments in - one set per chromosome
chr_reports <- list()

# For every chromosome, add chromosome-specific list element with arguments and
# other parameters:
# for( i in as.character( c(1:3, 9, 21) ) ){
for( i in as.character( args$chr_numbers ) ){
  chr_reports[[i]] <- list( 
    # Create file name for the chromosome-specific html report:
    out = gsub("\\.Rmd", sprintf("_chr%02d.html", as.integer(i) ), rmd_file_name)
    # (replace ".Rmd" in file name with _chr$$.html)
    , params = list( 
      chr_number = as.integer(i) 
      # Create root directory for the Rmd script (i.e. working directory for
      # inside the Rmd):
      , root_dir = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
      # Get name of EWAS map file based on regex:
      , ewas_map_file_name = ifelse( 
        test = length( list.files( 
          path = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
          , pattern = args[["ewas_map_file_name_base_regex"]] ) ) == 1
        # If only one matching file in root_dir:
        , yes = list.files( 
          path = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
          , pattern = args[["ewas_map_file_name_base_regex"]]
          , full.names = FALSE )
        # 1f zero or multiple matching files in root_dir:
        , no = stop(
          " There are zero or mutiple files in "
          , here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
          , " that match \""
          , args[["ewas_map_file_name_base_regex"]]
          , "\"."
        )
      )
      # , stage_dir = "0301" # use the default
    ) 
  )
}


# RUN 0301 FOR ALL CHROMOSOMES, IN PARALLEL

# doParallel::registerDoParallel( cores = n_cores ) 

# Parallel rendering of chromosome-specific reports:
message( "Performing Stage ", stage, " and rendering reports...\n")
system.time({
  status <- foreach( chr = chr_reports, .combine = c ) %dopar% {
    tryCatch( 
      expr = {
        render_chr_report(chr)
        paste0("Chromosome ", chr$params$chr_number, " successfully rendered.")
      }
      , error = function( e ){
        paste0( "Chromosome ", chr$params$chr_number, " failed to render: "
                , e$message)
      } )
  }
  doParallel::stopImplicitCluster()
})


# STATUS MESSAGES AND CHECKS

# If rendering one or more chromosomes failed for some reason:
if( sum(grepl("successfully", status)) < length(chr_reports) ){
  # if( sum(status) < length(chr_reports) ){
  stop( "Rendering ", stage, " report failed for the following chromosomes:\n\n"
        # , paste( sort( sapply( chr_reports[ which( !( status %in% 1 ) ) ]
        , paste( status[ which( grepl("failed", status) ) ], collapse = "\n" )
  )
} else{ message("Rendering successful for all chromosomes.\n")}


message(
  "Checking if html reports were succesfully created for all chromosomes....\n"
)
report_status <- sapply( chr_reports, function(chr){
  file.exists( here( "Reports", stage, chr$out ) ) 
})

# report_status should now be the path to the report that was produced.
if( all( report_status ) ){
  message( "Stage "
           , stage
           , " reports successfully rendered and written to \n\'"
           , here( "Reports", stage )
           , "\'\nfor all chromosomes: \n\n" 
           , paste( 
             sapply( chr_reports[ which( report_status == TRUE ) ]
                     , function(x) x$params$chr_number )
             , collapse = ", " )
  )
} else{
  stop( "No html report written to "
        , here( "Reports")
        , " for the following chromosomes:\n  
           "
        , paste( 
          sapply( chr_reports[ which( report_status == FALSE ) ]
                  , function(x) x$params$chr_number )
          , collapse = ", ")
        , "\n\n**Do not** proceed unless you are certain that the code in \'"
        , rmd_file_name
        , "\' was run successfully for all chromosomes (i.e. that the report "
        , "was not rendered successfully due to something unrelated, such as"
        , " a pandoc error etc.).\n"
  )
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOG ARGUMENTS USED IN THIS PHASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message( "\n-------------------------------\n"
         , "Logging phase 0004 arguments..."
         , "\n-------------------------------\n"
)

log_arguments( args_list = args
               , output_path =
                 here( "Reports", "Argument_logs", "0004_args.csv") )

if( file.exists(  here( "Reports", "Argument_logs", "0004_args.csv") ) ){
  message( "Log exported to "
           , here( "Reports", "Argument_logs", "0004_args.csv")
  )
} else{
  warning( "Log was not succesfully exported to "
           ,  here( "Reports", "Argument_logs", "0004_args.csv") )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE VARIABLES IN GLOBAL ENV THAT WE NO LONGER NEED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rm( rmd_file_name
#     , report_status
# )




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# AT THIS STAGE, THE USER SHOULD INSPECT THE GENERATED HTML REPORT, AS WELL AS
# THE CONTENTS OF THE CHROMOSOME-SPECIFIC FOLDERS IN pipeline_dir/chromosome
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE SET-UP PHASE IS NOW COMPLETE AND WE CAN START RUNNING THE STAGES OF THE
# PIPELINE, STARTING WITH 0101.

message( "Phase 0004 complete!\n.\n.\n.\n")



