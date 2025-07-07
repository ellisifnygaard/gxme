# THE SET-UP PHASE IS NOW COMPLETE AND WE CAN START RUNNING THE STAGES OF THE
# PIPELINE, STARTING WITH 0101.

# THIS STAGE IS IMPORTANT BECAUSE IT PERFORMS A LOT OF TESTS/DIAGNOSTICS TO
# DETERMINE WHETHER THERE ARE, OR COULD POTENTIALLY BE, ANY PROBLEMS WITH THE
# USER PROVIDED ARGUMENTS AND DATA THAT COULD NEGATIVELY IMPACT THE PIPELINE OR
# CAUSE THE PIPELINE TO FAIL.
# ALL ERRORS, WARNINGS AND MESSAGES SHOULD BE EXAMINED, ESPECIALLY THE ERRORS
# AND WARNINGS. 
# ___DO NOT__ PROCEED WITH THE SUBSEQUENT STAGES OF THE PIPELINE UNTIL THIS
# STAGE GIVES THE CONTENTS OF YOUR PIPELINE DIRECTORY "A CLEAN BILL OF HEALTH"
# FOR EACH CHROMOSOME BEING STUDIED. ONE HTML REPORT IS GENERATED PER CHROMOSOME
# - STUDY EACH REPORT BEFORE CONTINUING. IF THERE IS A PROBLEM WITH JUST ONE
# CHROMOSOME, THIS MAY CAUSE TECHNICAL ISSUES DOWNSTREAM IN THE PIPELINE.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PERFORM PRELIMINARY TESTS AND PREPARE DATA FOR PREPROCESSING, SEPARATELY FOR
# EACH CHROMOSOME, BY RENDERING THE RMD REPORT BELONGING TO THIS PHASE:
#     *0101*
#
# THIS STAGE FEATURES PARALLEL PROCESSING
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message( "Commencing phase 0003...\n")

start_time_stage <- Sys.time()

# n_cores <- parallel::detectCores() - 2 # (set in 0000_run_entire_pipeline.R)
rmd_file_name <- "0101_check_files_and_extract_study_population.Rmd"
stage <- "0101"

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
library(rmarkdown) # Not sure if this will work in parallel render scripts...
# library(rmarkdown, lib.loc = "C:/Program Files/R/library") # C-disk library works on SAFE

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE THE FUNCTION THAT WILL RENDER THE REPORTS AND THEIR ARGUMENTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (The code for rendering Rmds in parallel is based on the gist from
# hrbrmstr/do_rpt.r on Github)


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
      
      # , ewas_fileset_name = args[["ewas_fileset_name"]] 
      # "ewas" is the default (i.e. "ewas_env.ffData/RData)
      , ewas_family_member = args[["ewas_family_member"]]
      # Get name of EWAS map file based on regex:
      , ewas_map_file_name = list.files( 
        path = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
        , pattern = args[["ewas_map_file_name_base_regex"]]
        , full.names = FALSE )
      , ewas_annotation_file =
        base::basename( args[["ewas_annotation_file_path"]] )
      
      # , gwas_fileset_name = "gwas" # 0101 uses "gwas" as default.
      
      , plink_memory_mb = args[["plink_memory_mb"]]
      , plink_timeout = args[["plink_timeout"]]
      
      , scheme_states = args[["scheme_states"]]
      , cut_off = args[["cut_off"]]
      
      , sl_summarising_function = args[["sl_summarising_function"]]
      , sl_stratifying_function = args[["sl_stratifying_function"]]
      
      # Get name of mQTL file based on regex:
      , mqtl_file_name = list.files( 
        path = here( "chromosomes", sprintf("Chr%02d", as.integer(i)) )
        , pattern = args[["mqtl_file_name_base_regex"]]
        , full.names = FALSE )
      , mqtl_ld_r2 = args[["mqtl_ld_r2"]]
      # , stage_dir = "0101" # Use default stage directory
    ) 
  )
}


# DEFINE FUNCTION

render_chr_report <- function(chr){
  
  # Set wd to pipeline dir, otherwise the here function will use another root
  # directory:
  setwd(pipeline_dir) 
  
  require(here)
  # require(rmarkdown)
  # library(rmarkdown, lib.loc = "C:/Program Files/R/library") 
  # # (The gist used require, but I need to specify the library on C:/ in order to
  # # ensure that the rendering is successful -> SAFE quirk)
  
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
#
# RUN 0101 FOR ALL CHROMOSOMES, IN PARALLEL
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# RUN ALL CHROMOSOMES IN PARALLEL (I.E. RENDER ALL REPORTS IN PARALLEL)
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


# STATUS CHECK - GET NUMBERS OF FAILED CHROMOSOMES FOR RE-RUN, IF ANY

# Initiate empty vector:
failed_chr_number <- c()
# (Any failed chromosomes will be added to this vector)

if( sum(grepl("successfully", status)) == length(chr_reports) ){
  # IF RENDERING WAS SUCCESSFUL FOR ALL CHROMOSOMES:
  
  # If all chromosomes were rendered according to the status variable, check if
  # html reports were successfully created and stored in the /Reports folder:
  message("Rendering successful for all chromosomes.\n")
  message(
    "Checking if html reports were succesfully created for all chromosomes...\n"
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
             , "\'\nfor all remaining chromosomes: \n\n" 
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
  
  
} else if( sum(grepl("successfully", status)) < length(chr_reports) ){
  # IF RENDERING FAILED FOR ONE OR MORE CHROMOSOMES:
  
  #* WARN USER:
  message( "\nRendering "
           , stage
           , " report failed for the following chromosome(s):\n\n"
           , paste(status[ which( grepl("failed", status) ) ], collapse = "\n")
  )
  
  # Identify the messages in status belonging to the chromosomes that failed:
  status_failed_chr <- status[ which( grepl("failed", status) ) ]
  
  # Get the "Chromosome x failed" part of the elements in the status vector
  # belonging to failed chromosomes:
  failed_chr_status_msg <- regmatches(
    status
    , regexpr("^Chromosome[[:space:]]\\d{1,2}[[:space:]]failed", status )
  )
  # Get the chromosome number from those elements (as character):
  failed_chr_number <- regmatches(
    failed_chr_status_msg
    , regexpr("\\d{1,2}", failed_chr_status_msg)
  )
  
  # Test that the failed chr number matches the corresponding status message in
  # status_failed_chr
  for( i in 1:length(failed_chr_number) ){
    if( !grepl( paste0("Chromosome ", failed_chr_number[i])
                , status_failed_chr[i] ) ){
      stop( "The chromosomes in `failed_chr_number` do not match their "
            , "corresponding element in `status_failed_chr`.\n"
            , "failed_chr_number[i] = ", failed_chr_number[i]
            , "\nstatus_failed_chr[i] = ", substr( status_failed_chr[i], 1, 20)
            , "(...)"
      )
    }
  }
  
  # At this point, `failed_chr_number` contains the chromosome numbers we want
  # to re-run and try to render reports for (if the failures are due to
  # plink.exe failing due to memory issues, otherwise there is no point in
  # retrying)
  
  message( "\n\nAttempting reruns of chromosome(s) "
           , paste( failed_chr_number, collapse = ", ")
           , "...\n"
  )
  
} else{
  # If none of the statements above are TRUE then something has gone very wrong
  stop( "The number of chromsome reports ("
        , length(chr_reports)
        , ") is incongruous with the length of the `status` variable ("
        , length(status)
        , "). This should not occur."
  )
}


# As long as there are failed chromosomes in failed_chr_number, attempt to
# render the 0101 report for the remaining failed chromosomes.
# (But a maximum of x times) #xxxx parameterise x
rerun_limit <- 5
n_reruns <- 0

# AS LONG AS  1) there is one or more failed chr number,
#             2) there is is one or more status containing an error message
#             stating that plink returned a non-zero exit status, and
#             3) we have not surpassed our maximum limit of reruns,
# PERFORM THIS LOOP:
while( length(failed_chr_number) > 0 &
       sum(grepl("plink.*returned non-zero exit status", status)) > 0 &
       n_reruns <= rerun_limit
){
  n_reruns <- n_reruns + 1
  
  message( "\n\nAttempting reruns of chromosome(s) "
           , paste( failed_chr_number, collapse = ", ")
           , "...\n"
  )
  
  # Subset chr_reports containing all params etc.
  chr_reports_failed <- chr_reports[ failed_chr_number ]
  stopifnot( length(chr_reports_failed) == length(failed_chr_number) )
  
  
  # ATTEMPT TO RENDER FOR THE n_rerunsTH TIME:
  message( "\nPerforming a new attempt at Stage ", stage, " for chromosome(s) "
           , paste( failed_chr_number, collapse = ", ")
           , "...\n\n"
           , "(This is attempt number ", n_reruns, ")\n"
  )
  
  system.time({
    status <- foreach( chr = chr_reports_failed, .combine = c ) %dopar% {
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
  
  # STATUS CHECK OF RE-RUN
  
  # IF ALL REMAINING CHROMOSOMES WERE SUCCESSFULLY RENDERED:
  if( sum(grepl("successfully", status)) == length(chr_reports_failed) ){
    # IF RENDERING WAS SUCCESSFUL FOR ALL CHROMOSOMES:
    
    # If all chromosomes were rendered according to the status variable, check
    # if html reports were successfully created and stored in the /Reports
    # folder:
    message("Rendering successful for all remaining chromosomes.\n")
    message(
      "Checking if html reports were succesfully created for all chromosomes...\n"
    )
    report_status <- sapply( chr_reports_failed, function(chr){
      file.exists( here( "Reports", stage, chr$out ) ) 
    })
    
    # report_status should now be the path to the report that was produced.
    if( all( report_status ) ){
      message( "Stage "
               , stage
               , " reports successfully rendered and written to \n\'"
               , here( "Reports", stage )
               , "\'\nfor all remaining chromosomes: \n\n" 
               , paste( 
                 sapply( chr_reports_failed[ which( report_status == TRUE ) ]
                         , function(x) x$params$chr_number )
                 , collapse = ", " )
      )
    } else{
      stop( "No html report written to "
            , here( "Reports")
            , " for the following chromosomes:\n  
           "
            , paste( 
              sapply( chr_reports_failed[ which( report_status == FALSE ) ]
                      , function(x) x$params$chr_number )
              , collapse = ", ")
            , "\n\n**Do not** proceed unless you are certain that the code in \'"
            , rmd_file_name
            , "\' was run successfully for all chromosomes (i.e. that the report "
            , "was not rendered successfully due to something unrelated, such as"
            , " a pandoc error etc.).\n"
      )
    }
    
    # UPDATE failed_chr_number:
    failed_chr_number <- c()
    
    
    # IF SOME CHROMOSOMES FAILED TO RENDER AGAIN DUE TO PLINK PROBLEMS:
  } else if( sum(grepl("successfully", status)) < length(chr_reports_failed) &
             sum(grepl("plink.*returned non-zero exit status", status)) > 0 ){
    # IF RENDERING FAILED FOR ONE OR MORE CHROMOSOMES:
    
    #* ALERT USER:
    message( "\nRendering "
             , stage
             , " report failed for the following chromosome(s):\n\n"
             , paste(status[ which( grepl("failed", status) ) ], collapse = "\n")
    )
    
    
    # Identify the messages in status belonging to the chromosomes that failed:
    status_failed_chr <- status[ which( grepl("failed", status) ) ]
    
    # Get the "Chromosome x failed" part of the elements in the status vector
    # belonging to failed chromosomes:
    failed_chr_status_msg <- regmatches(
      status
      , regexpr("^Chromosome[[:space:]]\\d{1,2}[[:space:]]failed", status )
    )
    # Get the chromosome number from those elements (as character):
    failed_chr_number <- regmatches(
      failed_chr_status_msg
      , regexpr("\\d{1,2}", failed_chr_status_msg)
    )
    
    
    # DEV/BUG-FIXING:
    message( paste( paste0( "`failed_chr_number` = "
                            , failed_chr_number
                            , " <=> `status_failed_chr` =" )
                    , paste( substr( status_failed_chr, 1, 50), "\n" ) ) )
    
    
    # Test that the failed chr number matches the corresponding status message in
    # status_failed_chr
    for( i in 1:length(failed_chr_number) ){
      if( !grepl( paste0("Chromosome ", failed_chr_number[i])
                  , status_failed_chr[i] ) ){
        stop( "The chromosomes in `failed_chr_number` do not match their "
              , "corresponding element in `status_failed_chr`.\n"
              , "failed_chr_number[i] = ", failed_chr_number[i]
              , "\nstatus_failed_chr[i] = ", substr( status_failed_chr[i], 1, 20)
              , "(...)"
        ) } }
    
    
    # IF SOME CHROMOSOMES FAILED TO RENDER AGAIN BUT NOT DUE TO PLINK PROBLEMS:
  } else if( sum(grepl("successfully", status)) < length(chr_reports_failed) &
             sum(grepl("plink.*returned non-zero exit status", status)) == 0 ){
    
    stop( "The remaining chromosomes failed to render, but not due to plink"
          , " returning a non-zero exit status.\n"
          , "Look into/remedy these errors before attempting to run the "
          , "preprocessing phase again.\n"
          , "Stopping preprocessing run..."
    )
    
    # IF NONE OF THE ABOVE:
  } else{
    # If none of the statements above are TRUE then something has gone very wrong
    stop( "The number of chromosome reports ("
          , length(chr_reports_failed)
          , ") is incongruous with the length of the `status` variable ("
          , length(status)
          , "). This should not occur."
    )
  }
  
  # At this point, `failed_chr_number` contains the chromosome numbers we want
  # to re-run and try to render reports for.
}




# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # LOG 0101 EXECUTION TIME
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# message( "Logging execution time...\n" )
# 
# stage_execution_time <- Sys.time() - start_time_stage
# 
# # Create folder for csv file if it does not already exist
# stage_execution_time_dir <- 
#   here( pipeline_dir, "Results", "Preprocessing_stage_execution_times" )
# 
# if( !dir.exists( stage_execution_time_dir ) ){
#   dir.create( stage_execution_time_dir )
# }
# 
# exec_time_df <- data.frame(
#   stage = stage_dir
#   , chr = chr_number
#   , script_execution_time_seconds =
#     as.numeric( script_execution_time, units = "secs" )
#   , script_start_time = start_time
#   , computername = Sys.getenv("COMPUTERNAME")
#   , pipeline_dir = pipeline_dir
#   , cpu_model = benchmarkme::get_cpu()$model_name
#   , no_of_cores = benchmarkme::get_cpu()$no_of_cores
#   , ram_iec_units = print(benchmarkme::get_ram(), unit_system = "iec")
#   , system_memory_total_Mb = ps::ps_system_memory()$total / 1024^2
#   , system_memory_avail_Mb = ps::ps_system_memory()$avail / 1024^2
#   , R_platform = R.version$platform
#   , R_version = R.version$version.string
#   , Platform_GUI = .Platform$GUI
#   , RStudio_version = ifelse( .Platform$GUI == "RStudio"
#                               , yes = as.character(rstudioapi::getVersion())
#                               , no = NA )
# )
# 
# log_arguments( args_list = args
#                , output_path =
#                  here( "Reports", "Argument_logs", "0003_args.csv") )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOG ARGUMENTS USED IN THIS PHASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message( "Logging arguments...\n" )
log_arguments( args_list = args
               , output_path =
                 here( "Reports", "Argument_logs", "0003_args.csv") )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE VARIABLES IN GLOBAL ENV THAT WE NO LONGER NEED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rm( rmd_file_name
#     , render_output
# )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# AT THIS STAGE, THE USER SHOULD INSPECT THE GENERATED HTML REPORT, AS WELL AS
# THE CONTENTS OF THE CHROMOSOME-SPECIFIC FOLDERS IN pipeline_dir/chromosome
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# THE SET-UP PHASE IS NOW COMPLETE AND WE CAN START RUNNING THE STAGES OF THE
# PIPELINE, STARTING WITH 0101.

message( "Phase 0003 complete!\n.\n.\n.\n")
