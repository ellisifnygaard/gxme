# setwd("O:/Documents")
# # Delete entire pipeline directory:
# unlink( pipeline_dir, recursive = TRUE, force = TRUE )
# (Doing this, or using a "fresh" pipeline_dir is sometimes necessary to avoid
# errors )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SELECT ALL AND RUN ENTIRE SCRIPT IN ORDER TO GET ACCURATE TOTAL EXECUTION TIME 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


start_time <- Sys.time() # start timer


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         THE USER MUST EXPLICITLY SPECIFY THESE PARAMETERS:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PIPELINE DIRECTORY

# Fill in the path to the folder where you want the pipeline directory to be
# created and your chosen name for the pipeline directory. 

# The pipeline directory will be the root directory for all processing and
# analysis. It must be created from scratch to ensure that everything runs
# smoothly; the run will be stopped if there already exists a directory with
# the same name in pipeline_dir_parent.


# pipeline_dir_parent <- "S:/Project/SAM/Julia/Ellisif"
pipeline_dir_parent <- "C:/Temp/edj079"
# pipeline_dir_name <- paste0( Sys.Date(), "_pipeline_dir" )
pipeline_dir_name <- paste0( format( Sys.time(), "%Y-%m-%d_%H%M" ), "_pipeline_dir" )
# If you plan on performing multiple separate runs and then plot resulting data,
# it's recommended to choose pipeline directory names that are intuitive and
# easily distinguishable, as the pipeline directory features in several reports

# SCRIPT DIRECTORY 

# Fill in the path to a folder containing all the .Rmd and .R files that
# comprise the pipeline:
path_pipeline_scripts <- "S:/Project/SAM/Julia/Ellisif/Paper 1/paper1"

# ARGUMENT FILE

# Fill in the path to the R file containing all your pipeline arguments:
arguments_path <- "S:/Project/SAM/Julia/Ellisif/Paper 1/paper1/0000_arguments.R"

# IF USING the 'plink_memory_mb' ARGUMENT, REMEMBER TO USE A REASONABLE VALUE
# (see f_get_available_mb_per_process.R if you need help with this)


# NUMBER OF CORES TO USE IN THE PARALLEL PROCESSING

# Number of cores to use in the pre-processing stages with parallel processing
# (the number of chromosomes or less):
n_cores <- 6



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       PHASE 0001 - setting up pipeline directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_0001 <- file.path( path_pipeline_scripts
                         , "0001_initial_set_up_pipeline_directory.R" )

system.time( source( phase_0001 ) )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PHASE 0002 - setting up chromosome directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_0002 <- here( "R", "0002_set_up_chromosome_directory.R" )

system.time( source( phase_0002 ) )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      SET-UP FOR PARALLEL RENDERING OF CHROMOSOME REPORTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Setting up a cluster that will be used throughout the preprocessing stages of
# the pipeline. 

# Close any open connections/active clusters
base::closeAllConnections()

# No point in creating cluster with more cores than chromosomes:
if( n_cores > length( args$chr_numbers ) ){
  warning( "`n_cores` ("
           , n_cores
           , ") is greater than the number of chromosomes ("
           , length( args$chr_numbers )
           , ").\nSetting n_cores = length( args$chr_numbers )..." )
  n_cores <- min( n_cores, length( args$chr_numbers ) )
}

# Create parallel socket cluster
cl <- parallel::makePSOCKcluster( n_cores )

# Register the parallel backend with the foreach package
doParallel::registerDoParallel( cl, cores = n_cores )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PHASE 0003 - perform preliminary tests and prepare data for pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# phase_0003 <- here( "R", "0003_perform_preliminary_tests_and_establish_study_population.R" )
phase_0003 <- here( "R", "0003_perform_preliminary_tests_and_establish_study_population_HANDLE_FAILED_CHROMOSOMES.R" )

system.time( source( phase_0003 ) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PHASE 0004 - pre-stratification scheme preparations of gwas and ewas data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_0004 <- here( "R", "0004_pre_stratification_scheme_preparations.R" )

system.time( source( phase_0004 ) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PHASE 0005 - stratification scheme (pairing of snps and ewas probes)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_0005 <- here( "R", "0005_stratification_scheme.R" )

system.time( source( phase_0005 ) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PHASE 0006 - calculate strata allocations and prepare data for analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_0006 <- here( "R", "0006_calculate_strata_and_make_data_ready_for_haplin.R" )

system.time( source( phase_0006 ) )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      PREPROCESSING STAGES COMPLETE. SHUT DOWN CLUSTER.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Shut down the workers used by foreach:
doParallel::stopImplicitCluster()

# Clean up any remnants from previous clusters
# (Code from answer from Steve Weston on Stack Overflow dated 4/8/2014 )
foreach_env <- foreach:::.foreachGlobals
rm( list = ls( name = foreach_env ), pos = foreach_env )

# Close any open connections (this should not affect the ff data)
base::closeAllConnections()
base::showConnections()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      CALCULATE EXECUTION TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end_time <- Sys.time() # stop timer
execution_time <- end_time -start_time

cat("Total execution time: "
    , as.numeric( execution_time, units = "mins" )
    , "minutes")
