# SNP-WISE ANALYSIS FUNCTION USED FOR TESTING PRIOR TO INITIATING PARALLEL
# EXECUTION

# CONTAINS ADDITIONAL TESTS (SOME OF WHICH WERE CUT FROM f_analyse_snp.R TO
# SPEED UP EXECUTION TIME)

# FUNCTION THAT TAKES SNP ID, SNP x SL DATA FRAME AND FF MATRIX WITH STRATA AS
# INPUT AND OUTPUTS A DATA FRAME/DATA.TABLE WITH THE GXE RESULTS
analyse_snp_for_testing_pre_execution <- function(
    snp_id
    , snp_sl_dt = snp_sl # data.table with all snp x sl pairings
    , n_unique_chr_snps = n_unique_snps # no of unique SNPs on chr being analysed
    , chunk_number = chunk_i # chunk currently being executed
    , max_chunk_snp_number # index of the last SNP in the chunk to be analysed
    , first_last_snp = first_last_snp_in_each_chunk # vector of snp ids whose log file we want to keep
    , strata_ff_matrix # ff_matrix with all SL strata
    , root_dir = root_dir # Path to root directory in rmd
    , stage_dir = stage_dir # "0701"
    , log_dir = log_dir # Path to log folder
    , map_file = map #qqq filter first by snp_id?
    , gwas = gwas_all_snps
    , n_cores = n_cores # Number of workers (used to determine log_limit)
    , ewas_family_member = ewas_family_member
    # Haplin arguments:
    , haplin_args = haplin_args
    # Maximum number of log files to keep in log directory (in addition to one
    # file per core)
    , keep_max_x_log_files = 1000
    # Vector with strings containing the names of p-values from the results of
    # Haplin::gxe() to return:
    # , gxe_pvals = ifelse( haplin_args$haplin_poo == FALSE, yes = "child", no = "poo" ) 
    , gxe_pvals = c( "child", "child.trend" )
){
  # profvis::profvis({
  
  # Should not be necessary due to packages being exported via foreach()
  # # LOAD PACKAGES -----
  # library(data.table)
  # library(dplyr)
  # library(ff)
  # library(Haplin)
  # library(logr)
  # library(readr)
  
  
  # Specify that data.table is only to use one core/thread:
  data.table::setDTthreads( threads = 1 )
  
  
  # LOG FILE -----
  
  # Initialise log file in "stage_dir/log", using Process ID, timestamp, and SNP
  # ID in the file name:
  log_filename <- file.path(
    log_dir
    , paste0(
      # Date + time
      format( Sys.time(), "%Y-%m-%d_%H%M%OS4" ), "_TEST_"
      # , Sys.getpid(), "_"
      # Chunk number
      , "c", chunk_number, "_"
      # SNP number x out of Y
      # (where Y = the maximum snp number in the chunk; indicates how many
      # SNPs are left in the current chunk)
      , "snp_"
      , unique(snp_sl_dt[ snp == snp_id, snp_order_overall])
      # , unique(snp_sl_dt[ snp == snp_id, snp_index])
      , "_out_of_"
      , max_chunk_snp_number, "_"
      # SNP ID:
      , snp_id, "_"
      # Number of SNP x SL pairings
      , snp_sl_dt %>% filter( snp == snp_id ) %>% nrow()
      , "_sls"
    )
  )
  # Open log:
  log_path <- logr::log_open(
    file_name = log_filename
    # Place the log file in a directory named "log" inside the working directory:
    , logdir = TRUE
    # Show traceback if error occurs:
    , traceback = TRUE
    # Disable automatic printing of notes (notes show Log Print Time and Elapsed
    # Time for every entry in the log -> clutter)
    , show_notes = FALSE
  )
  
  # # Log available memory:
  # memuse::Sys.meminfo() %>% logr::put()
  
  # # Start SNP analysis timer:
  # start_analyse_snp_time <- Sys.time()
  
  # LOAD HAPLIN ARGUMENTS:
  # Add all the haplin parameters/arguments to the environment
  for( x in names(haplin_args) ){
    assign( x, haplin_args[[x]], envir = parent.frame() )
  }
  # Stop if assigning the haplin parameters did not succeed:
  stopifnot( all( sapply( names(haplin_args), exists) == TRUE ) )
  
  
  # CHECK THAT THE SPECIFIED P-VALUE COLUMN NAMES IN gxe_pvals ARE IN
  # CONCORDANCE WITH THE GIVEN HAPLIN 'poo' ARGUMENT
  if( haplin_poo == FALSE ){
    # When poo = FALSE, the p-value column name(s) cannot start with "poo":
    if( any( grepl( "^poo", gxe_pvals ) ) ){
      stop( "`gxe_pvals` contains column names starting with 'poo'"
            , ", but `haplin_poo` = ", haplin_poo, " !" )
    }
  } else if( haplin_poo == TRUE ){
    # When poo = TRUE, the p-value column name(s) cannot start with "child":
    if( any( grepl( "^child", gxe_pvals ) ) ){
      stop( "`gxe_pvals` contains column names starting with 'child'"
            , ", but `haplin_poo` = ", haplin_poo, " !" )
    }
  } else{
    stop( "`haplin_poo` = "
          , haplin_poo
          , ", which is an unexpected value."
    )
  }
  
  # CHECK THAT THE SPECIFIED P-VALUE COLUMN NAMES IN gxe_pvals ARE IN
  # CONCORDANCE WITH THE P-VALUES GIVEN IN A GXE RESULTS TABLE
  if( !all( gxe_pvals %in% c( "haplo.freq", "child", "poo"
                              , "haplo.freq.trend", "child.trend", "poo.trend" )
  ) ){
    stop( "`gxe_pvals` contains non-valid names of `gxe.test` p-values!" )
  }
  
  # If the analysis of this SNP is successfully completed, the log file will be
  # deleted permanently.
  # If something goes wrong, however, we want to save the log file so that the
  # user can examine it to determine the cause of the problem.
  # If something went wrong, we want to skip the code deleting the log file,
  # BUT: only if there isn't a lot of log files there already. If all the
  # analyses for all SNPs should fail, that could entail millions of log files,
  # which is not viable.
  
  
  # FILTER LARGE DATA FRAMES/TABLES BY SNP ID -----
  
  # Filter map_file by snp_id:
  map_file <- map_file %>% filter( snp == snp_id )
  
  # Filter snp_sl_dt by snp_id:
  snp_sl_dt <- snp_sl_dt[ snp == snp_id ]
  
  # Stop if snp_sl_dt has zero rows:
  stopifnot( nrow(snp_sl_dt) > 0 )
  
  
  
  # PROCESS DIRECTORY -----
  
  # Path to process directory:
  process_dir <- file.path( root_dir, stage_dir, Sys.getpid() )
  
  # Create process-specific temporary directory unless it already exists:
  if( !dir.exists( process_dir ) ){ dir.create( process_dir ) }
  stopifnot( dir.exists( process_dir ) )
  
  # Delete files left over from previous runs:
  files_in_process_dir <- list.files( process_dir, full.names = TRUE )
  
  if( length(files_in_process_dir) > 0 ){
    # Deleting files that were already present in process_dir:
    unlink( file.path( process_dir, "*" ), recursive = TRUE, force = TRUE )
  }
  
  
  # SUBSET GWAS DATA -----
  
  # Select only GWAS data related to the SNP currently being analysed (snp_id)
  
  # Check that current SNP is in fact in the gwas data:
  if( !( nrow( map_file ) == 1 ) ){
    # If SNP is not in gwas data, add message to log, and stop iteration:
    status_message <- paste0(
      snp_id
      , " could not be located in the .map file "
      , "(i.e. the haplin.data object from ", gwas_fileset_dir , ". "
      , "Please review your genotype data and .map file."
    )
    logr::put( status_message, console = FALSE )
    # If this occurs, something has gone very wrong indeed. Throw error:
    stop( status_message )
  }
  
  # Determine snp_id_haplin_index:
  snp_id_haplin_index <- map_file$haplin_index
  
  # Use genDataGetPart to create a haplin.data object containing only the SNP
  gwas_snp_id <- Haplin::genDataGetPart(
    data.in = gwas
    , design = haplin_design
    , markers = snp_id_haplin_index
    #qqq Do we need the cc and sex arguments to be parametrised as well?
    # Name of file with subset contains SNP ID and Process ID:
    , file.out = paste0( stage_dir
                         , "_gwas_snp_id_"
                         , snp_id
                         , "_"
                         , Sys.getpid() )
    # (File names contain process ID as an extra insurance against different
    # processes attempting to write to the same files.)
    , dir.out = process_dir
    , overwrite = TRUE
  )
  
  
  # If subsetting the haplin.data was unsuccessful, something went very wrong
  # and it it is likely to go wrong in all iterations -> throw error
  if( !( exists("gwas_snp_id") & gwas_snp_id$aux$marker.names == snp_id ) ){
    stop( "exists(\"gwas_snp_id\") = "
          , exists("gwas_snp_id")
          , " and gwas_snp_id$aux$marker.names = "
          , gwas_snp_id$aux$marker.names
          , ".\nHaplin::genDataGetPart failed."
          ,"\nThis should not occur. "
          , "Was your genotype data sucessfully preprocessed?"
    )
  }
  
  
  # EXTRACT SL STRATA -----
  
  logr::put( "Extract SL strata data from ff_matrix..." )
  # cat( "Extract SL strata data from ff_matrix..." )
  
  # Vector with SLs belonging to this SNP (sorted by number in sl_id):
  sls_snp_id <- snp_sl_dt %>%
    mutate( sl_as_number = as.integer( gsub("^SL", "", sl_id) ) ) %>%
    arrange( sl_as_number ) %>%
    pull( sl_id )
  # Stop if sls_snp_id contains duplicates:
  stopifnot( "sls_snp_id contains duplicates" =
               any( duplicated(sls_snp_id) ) == FALSE )
  
  # Extract the strata for the SLs that are paired with the current SNP:
  sls_strata <- all_sl_strata_ff[ sls_snp_id, , drop = FALSE ]
  stopifnot( exists("sls_strata") )
  stopifnot( ncol( sls_strata ) > 0 )
  # Note: subsetting the ff_matrix produces a regular matrix. Does not appear to
  # be a simple way of avoiding this.
  # class(sls_strata)
  
  
  # ADD THE STRATA COLUMN(S) TO HAPLIN DATA -----
  
  # This is a little bit elaborate due to the possibility of family units with
  # multiple offspring and the parameterisation of the family member that the
  # EWAS data is from.
  
  # Use both id.fam and id.ewas_family_member as keys when joining SL strata
  # (split the individual IDs from the SL strata file using "_" as separator and
  # transform "01"-suffix into integer, then character)
  
  # Transpose the matrix and turn the rownames (i.e. individual IDs) into column:
  sls_strata_dt <- data.table::setDT(
    as.data.frame.matrix( t( sls_strata ) )
    , keep.rownames = "id"
  )

  # Check that this did not alter the individuals' assigned strata:
  # Test 20 randomly selected individuals:
  random_indiv_ids <-
    colnames(sls_strata)[ sample( 1:ncol(sls_strata)
                                  , min(20, ncol(sls_strata))
    ) ]
  stopifnot(
    all( sapply( random_indiv_ids, function(i){
      all( sls_strata[, i, drop = TRUE]  ==
             sls_strata_dt[ id == i, -1] %>% unlist() )
    }) )
  )
  # sls_strata_dt is a transposed version of the data from the ff_matrix, and it
  # contains individual strata allocations identical to those in the ff_matrix.
  rm(sls_strata)
  
  
  # ADD STRATA COLUMNS TO $COV.DATA -----
  
  # # Peek at cov.data:
  # gwas_snp_id$cov.data %>% head()
  
  # Join strata columns to $cov.data in gwas_snp_id
  # First, split the id column into two columns using "_" as delimiter:
  sls_strata_dt <- sls_strata_dt %>%
    tidyr::separate_wider_delim( cols = id
                                 , delim = "_"
                                 , names = c("id.fam", "id.c")
    ) %>%
    # Remove leading zero from id.c by converting to integer and then back to
    # character:
    mutate( id.c = as.character( as.integer( id.c ) ) )
  
  
  # Then, create a copy of cov.data with strata columns added using left_join:
  cov.data_with_strata <- gwas_snp_id$cov.data %>%
    as.data.frame.matrix() %>%
    dplyr::left_join( . , sls_strata_dt, by = c("id.fam", "id.c") ) %>%
    as.matrix()
  
  # Check that the copy of the cov.data matrix that was merged with strata data
  # is identical to cov.data when disregarding the strata columns:
  stopifnot( all.equal( cov.data_with_strata[, 1:6], gwas_snp_id$cov.data ) )
  
  # Replace cov.data with the copy that has strata columns:
  gwas_snp_id$cov.data <- cov.data_with_strata
  logr::put( gwas_snp_id, console = FALSE, hide_notes = TRUE )
  
  # # Peek at cov.data:
  # gwas_snp_id$cov.data %>% head()
  
  # PREPARE USING genDataPreprocess -----
  
  # First, create a map file containing only the current SNP:
  # (otherwise the marker in the haplin.ready object will get dummy names)
  map_snp_id_file <-
    #qqq Why is this line a problem all of a sudden?
    readr::read_delim( file.path( root_dir, gwas_snp_id$aux$map.filename )
                       # readr::read_delim( gwas_snp_id$aux$map.filename
                       , show_col_types = FALSE ) %>%
    filter( snp == snp_id )
  
  stopifnot( nrow(map_snp_id_file) == 1 )
  
  # Write map file to process directory:
  map_snp_id_file %>% readr::write_tsv(
    .
    , file = file.path( process_dir, paste0(snp_id, ".map") )
    , col_names = TRUE
    , append = FALSE # overwrite the existing map file w/o header
  )
  
  # Then, run genDataPreprocess:
  
  # Initiate empty object for capturing non-critical warnings from
  # genDataPreprocess():
  genDataPreprocess_warning <- NULL 
  
  suppressWarnings( 
    # suppressWarnings is necessary, otherwise non-critical warnings are not
    # silenced and trigger the creation of an msg file due to the logr package.
    withCallingHandlers({ gwas_snp_id <- Haplin::genDataPreprocess(
      # gwas_snp_id <- Haplin::genDataPreprocess(
      data.in = gwas_snp_id
      , design = haplin_design
      # Overwrite the file containing the haplin.object we got when subsetting:
      , file.out = paste0( stage_dir
                           , "_gwas_snp_id_"
                           , snp_id
                           , "_"
                           , Sys.getpid() )
      # (File names contain process ID as an extra insurance against different
      # processes attempting to write to the same files.)
      , dir.out = process_dir
      , map.file = file.path( process_dir, paste0(snp_id, ".map") )
      , map.header = TRUE # suddenly necessary for some reason
      , overwrite = TRUE
      # Don't specify multiple cores when we're already parallelised; it slows
      # down the code significantly...
      , ncpu = 1
    )
    
    # gwas_snp_id
    }, warning = function(w){
      # Check if warning message is a non-critical warning:
      if(
        grepl("Found family size larger than 3!", w$message) |
        grepl("refer to non-existing individuals and have been set to missing"
              , w$message)
      ){
        # Don't print non-critical warnings to log in order to speed up
        # execution time.
        # # If warning is not critical, write it to the log file as regular text
        # # (thus avoiding triggering the creation of an auxiliary .msg file)
        # paste0( "Haplin::genDataPreprocess issued the following non-critical "
        #         , "warning:\n"
        #         , w$message
        # ) %>% logr::put()
      } else{
        # If another warning than those predefined as non-critical is issued,
        # then reissue warning outside of this withCallingHandlers call. (This
        # will trigger logr to create an .msg file)
        genDataPreprocess_warning <<- w
      }
      # gwas_snp_id
    }
    
    ) # withCallingHandlers end
  ) #suppressWarnings end
  
  # If a possibly critical warning occurred, reissue warning from the
  # withCallingHandlers call above:
  if( length(genDataPreprocess_warning) > 0 ){
    warning( "Haplin::genDataPreprocess issued the following warning(s) "
             , "that might require close attention:\n"
             , paste0( genDataPreprocess_warning, collapse = "\n")
    )
  }
  
  # # Peek at cov.data:
  # gwas_snp_id$cov.data %>% head()
  
  # Stop if genDataPreprocess() did not result in a haplin.ready object:
  stopifnot(
    exists("gwas_snp_id") & class(gwas_snp_id) %in% "haplin.ready"
  )
  
  # Identify the strata column(s) in haplin.ready$cov.data belonging to the
  # ewas_family_member so that we can test them for NA:
  strata_col_ewas_family_member <- which(
    grepl( paste0( sls_snp_id, "\\.", ewas_family_member, collapse = "|" )
           , colnames( gwas_snp_id$cov.data ) ) )
  # Stop if the strata column(s) belonging to the EWAS family members contains
  # any missing values
  if( any( is.na( gwas_snp_id$cov.data[, strata_col_ewas_family_member] ) ) ){
    stop(
      "gwas_snp_id$cov.data contains missing data. This is not allowed as all"
      , " EWAS family members included in the analysis have to be allocated to "
      , "one stratum per state locus." )
    # The error message should be printed to log without using logr::put()
  }
  
  gwas_snp_id %>% logr::put( . , console = FALSE, hide_notes = TRUE )
  
  # Remove superfluous strata columns created by genDataPreprocess
  
  # Identify the strata columns comprised only of NA because they were intended
  # to contain the strata numbers for the "non-EWAS" family members:
  superfluous_cols <- which( grepl(
    paste0("SL[[:digit:]]+\\.", setdiff(c("c", "m", "f"), ewas_family_member)
           , collapse = "|")
    , colnames( gwas_snp_id$cov.data )
  ) )
  
  # colnames( gwas_snp_id$cov.data )[superfluous_cols]
  # Remove the superfluous columns from the haplin.ready data:
  gwas_snp_id$cov.data <- gwas_snp_id$cov.data[, -superfluous_cols]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   WE ARE NOW READY TO RUN HAPLINSTRAT FOR EACH SL
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # RUN HAPLINSTRAT FOR EACH SNP X SL COMBINATION -----
  
  
  # FOR EACH SL PAIRED WITH SNP_ID:
  sl_res <- lapply( seq_along(sls_snp_id), function(i){
    sl <- sls_snp_id[i]
    
    # Additional logging:
    cat( "\n\n\nSNP x state locus pairing", which( sls_snp_id == sl )
         , "out of", length(sls_snp_id)
         , paste0("\n\n", snp_id), "x", sl , "\n")
    
    paste0( "Pairing ", sprintf("%-7.0d", which( sls_snp_id == sl ))
            , ": Analysing ", snp_id, " x ", sl, " with haplinStrat()" ) %>%
      logr::sep( console = FALSE )
    
    paste0( "SL number ", which( sls_snp_id == sl )
            , " out of the ", length(sls_snp_id)
            , " SLs in total that are paired with ", snp_id , ".") %>%
      logr::put( console = FALSE, hide_notes = TRUE )

    
    # Get the column number in cov.data of the column with the SL's strata:
    strata_col_number <- which( grepl(
      paste0( sl, "\\.", ewas_family_member ), colnames( gwas_snp_id$cov.data )
    ) )
    
    # RUN HAPLINSTRAT -----
    
    # List to store warnings generated in by haplinStrat
    hapstrat_warnings <- list()
    # List to store error generated in by haplinStrat
    hapstrat_error <- list()
    
    # Run haplinStrat inside withCallingHandlers (to catch all warnings that are
    # generated), which is run inside tryCatch (for error handling):
    res_sl_hapstrat <- tryCatch(
      expr = withCallingHandlers({
        Haplin::haplinStrat(
          # res_sl_hapstrat <- Haplin::haplinStrat(
          data = gwas_snp_id
          , strata = strata_col_number # column with SL-specific strata
          , markers = "ALL" # hard-code
          , design = haplin_design
          , use.missing = haplin_use.missing
          , x_chrom = haplin_xchrom #qqq hard-code to FALSE for the time being?
          , maternal = haplin_maternal #qqq hard-code to FALSE for the time being?
          , test.maternal = haplin_test.maternal #qqq hard-code to FALSE for the time being?
          , poo = haplin_poo #qqq hard-code to FALSE for the time being? Or allow PoO?
          , scoretest = haplin_scoretest #qqq hard-code to FALSE for the time being?
          , ccvar = haplin_ccvar #qqq hard-code to FALSE for the time being?
          , sex = haplin_sex #qqq hard-code to FALSE for the time being?
          , comb.sex = haplin_comb.sex
          , reference = haplin_reference
          , response = haplin_response
          # N.B. If response = "mult", but reference is not "ref.cat", then
          # reference is changed to "ref.cat" by haplin and a warning is given.
          , threshold = haplin_threshold
          , max.haplos = haplin_max.haplos # not yet implemented in Haplin!
          , haplo.file = NULL # not yet implemented in Haplin!
          , resampling = "no" # hard-coded; see Haplin reference manual.
          , max.EM.iter = haplin_max.EM.iter
          , data.out = haplin_data.out # hard-coded
          , verbose = haplin_verbose # hard-coded
          , printout = haplin_printout # hard-coded
        )
        # # Return result:
        # res_sl_hapstrat
      }
      # withCallingHandlers handling or warnings:
      , warning = function( w ){
        #qqq Skip logging the warning message? If it is registered by its code,
        #then this does not add much informational value once the all the
        #possible error and warning scenarios have been mapped out
        logr::put( paste0(
          "////// Warning from Haplin::haplinStrat when analysing "
          , snp_id, " x ", sl,":" ), console = FALSE )
        logr::put( w, console = FALSE, msg = FALSE )
        # Append warning message to hapstrat_warnings:
        hapstrat_warnings <<- append(hapstrat_warnings, w)
      }
      )
      , error = function( e ){
        # Append error message to hapstrat_error:
        hapstrat_error <<- append(hapstrat_error, e)
        #qqq Skip logging the error message? If it is registered by its code,
        #then this does not add much informational value once the all the
        #possible error and warning scenarios have been mapped out
        # Log the error:
        logr::sep( paste0(
          "\nError produced by Haplin::haplinStrat when analysing "
          , snp_id, " x ", sl,":" ), console = FALSE  )
        logr::put( e, console = FALSE )
        
        # Return a data.table with error code = 9 (unknown error):
        res_sl_gxe_df <- data.table( snp = snp_id
                                     , sl_id = sl
                                     , error_code = "9" # 9 = unknown error
                                     , warning_code = NA
                                     , gxe_execution_allowed = FALSE
        ) 
        res_sl_gxe_df
      }
    ) # tryCatch end
    
    
    # CREATE DATA.TABLE ----
    
    # The data.table where we will place any captured error and warning codes,
    # as well as the p-value from gxe if the results from haplinStrat is cleared
    # for running gxe. Same structure as res_sl_gxe_df.
    
    # If res_sl_hapstrat is a data.table due to an error occurring, set
    # res_sl_gxe_df equal to res_sl_hapstrat.
    # If res_sl_hapstrat is a haplinStrat object (i.e. no error occurred), then
    # create from scratch.
    if( "data.table" %in% class(res_sl_hapstrat) ){
      res_sl_gxe_df <- res_sl_hapstrat
    } else if( "haplinStrat" %in% class(res_sl_hapstrat) ){
      res_sl_gxe_df <- data.table( snp = snp_id
                                   , sl_id = sl
                                   , error_code = NA
                                   , warning_code = NA
                                   , gxe_execution_allowed = NA
      ) 
    } else{
      stop(
        "`res_sl_hapstrat` is neither a data.table nor a haplinStrat object!"
        , " Something has gone wrong.\n"
        , "class(res_sl_hapstrat) = \"", class(res_sl_hapstrat), "\""
      )
    }
    
    
    # FIX NAMES OF ELEMENTS IN HAPLINSTRAT OBJECT ----
    
    # Fix the names of the list elements containing the strata-specific results
    # in res_sl_hapstrat if necessary
    
    # If the "all" 'stratum' has a named list element, but the individual strata
    # do not, fix the names
    if( identical(
      names(res_sl_hapstrat)
      , c("all", rep.int(NA, times = length(res_sl_hapstrat) - 1 ) )
    ) ){

      # Log occurrence:
      logr::put( paste0(
        "names(res_sl_hapstrat) contains NA.\n"
        , "names(res_sl_hapstrat) = "
        , paste0( names(res_sl_hapstrat), collapse = ", " )
        , "\nFixing the list element names..."
      ) )
      
      na_element_name_index <- which( is.na( names(res_sl_hapstrat) ) )
      
      # If there are s strata, then the indexes should be 2, 3, ..., s:
      stopifnot(
        identical( na_element_name_index, 2:( length(res_sl_hapstrat) ) )
      )
      # Get strata 'names' (i.e. 1,2, ... ) by subtracting 1 from each index:
      names(res_sl_hapstrat)[na_element_name_index] <- na_element_name_index - 1
      #xxxxx Perform test executions of entire chromosomes with other
      #stratifying functions/numbers of strata to check that this works in all
      #settings
      
      # Check if renaming was successful:
      stopifnot(
        "Renaming res_sl_hapstrat elements was unsuccessful" =
          identical( names( res_sl_hapstrat )
                     , c( "all", na_element_name_index - 1 ) )
      )
    }
    
    
    # REGISTER ANY WARNINGS OR ERRORS ----
    
    # The errors and warnings registered and added to res_sl_gxe_df will later
    # determine whether res_sl_hapstrat is cleared for analysis with gxe()
    
    # WARNINGS ----
    # REGISTER ANY WARNING SCENARIOS/CODES THAT OCCURRED
    
    # If the hapstrat_warnings list has elements in it, a warning occurred while
    # haplinStrat was running
    if( length( hapstrat_warnings ) > 0 ){
      
      # Start with removing any list elements in hapstrat_warnings that are NULL
      hapstrat_warnings[ which( sapply( hapstrat_warnings, is.null ) ) ] <- NULL
      # (Sometimes $call contains NULL)
      
      # Print all warning messages to the log:
      logr::put( "Contents of `hapstrat_warnings`:  "
                 , console = FALSE )
      sapply( hapstrat_warnings, logr::put )
      
      
      ## WARNING CODE 1 ----
      # The following warning was produced when running haplinStrat():
      # "Something's strange with the frequency count in HWE test!"
      # (The warning stems from Haplin:::f.sel.haplos() )
      if( any( sapply( hapstrat_warnings, function(w){
        grepl( "strange with the frequency count in HWE test", w )
      } ) ) ){
        
        logr::put( "Warning 1 identified.", console = FALSE )
        
        # Add warning code to warning_code column:
        res_sl_gxe_df$warning_code <- paste0( res_sl_gxe_df$warning_code, 1 )
        # Remove "NA" from warning_code column:
        res_sl_gxe_df$warning_code <-
          gsub( "^NA", "", as.character(res_sl_gxe_df$warning_code) )
        
        
        # WARNING CODE 2 ----
        # haplinStrat()'s call of haplin resulted in the following warning
        # "Maximum number of EM iterations reached!
        # Convergence not yet obtained. Setting max.EM.iter higher may help.",
        # in one stratum or in multiple strata.
        # In these cases, we consider the estimates in the haplinStrat object as
        # reliable enough and return a data.table with warning_code = 2.
      } else if( any( sapply( hapstrat_warnings, function(w){
        grepl( "Maximum number of EM iterations reached", w )
      } ) ) ){
        
        
        # Determine in which strata max.EM.iter was reached without convergence:
        if( "haplinStrat" %in% class(res_sl_hapstrat) ){
          # (Only do this if we have a haplinStrat object.)
          failed_strata <- c()
          for( s in seq_along(res_sl_hapstrat) ){
            if( any( res_sl_hapstrat[[s]]$info$estimation$EM.conv %in% FALSE ) ){
              failed_strata <- c(failed_strata, names(res_sl_hapstrat)[s] )
            }
          }
          
          # Log the strata for which convergence was not reached:
          paste0( "Warning scenario 2 identified."
                  , "\nThe maximum number of EM iterations ("
                  , haplin_max.EM.iter
                  , ") was reached in the following strata:\n"
                  , paste0( failed_strata, collapse = ", ")
                  , "\nThe estimates for this pairing will be disregarded since"
                  , " convergence was not reached."
          ) %>% logr::put()
        } else{
          paste0( "Warning 2 identified, but class(res_sl_hapstrat) = "
                  , "\"", class(res_sl_hapstrat) , "\"" ) %>% 
            logr::put()
        }
        
        # Add warning code to warning_code column:
        res_sl_gxe_df$warning_code <- paste0( res_sl_gxe_df$warning_code, 2 )
        # Remove "NA" from warning_code column:
        res_sl_gxe_df$warning_code <-
          gsub( "^NA", "", as.character(res_sl_gxe_df$warning_code) )
        
        # THIS WARNING DISQUALIFIES res_sl_hapstrat FROM GXE ANALYSIS
        # Update gxe_execution_allowed column accordingly:
        res_sl_gxe_df$gxe_execution_allowed <- FALSE
        
        ## WARNING CODE 9 (OTHER/UNSPECIFIED) ----
        # If haplinStrat produced a warning, but not of the types outlined
        # above, let warning_code = 9, meaning "other", unspecified warning
        # scenario
      } else{
        logr::put( "NOTE! Unspecified warning scenario registered.
                   warning_code = 9", console = FALSE )
        
        # Add warning code to warning_code column:
        res_sl_gxe_df$warning_code <- paste0( res_sl_gxe_df$warning_code, 9 )
        # Remove "NA" from warning_code column:
        res_sl_gxe_df$warning_code <-
          gsub( "^NA", "", as.character(res_sl_gxe_df$warning_code) )
      }
    }
    
    
    # ERRORS ----
    
    
    # ERROR CAPTURED BY TRYCATCH WHEN RUNNING HAPLINSTRAT ----
    
    # If haplinStrat threw an error, code the error and update res_sl_gxe_df
    # accordingly.
    
    if(  "data.table" %in% class(res_sl_hapstrat) ){
      # If this happens, it means that res_sl_hapstrat is the empty
      # res_sl_gxe_df data frame that was created in the error part of the
      # tryCatch block.
      
      # CODE THE ERROR
      
      # Start with removing any list elements in hapstrat_warnings that are NULL
      hapstrat_error[ which( sapply( hapstrat_error, is.null ) ) ] <- NULL
      
      
      ## ERROR CODE 1 ----
      # haplinStrat error message = "Less than 2 haplotypes above threshold.
      # Locus may have low information content (or threshold set too high)"
      
      # NB: SNP-WIDE ERROR! WE WANT TO SKIP THE REMAINING PAIRINGS
      # This error is SNP-wide, meaning that if it occurs, it occurs for all SNP
      # x SL pairings of that particular SNP.
      #qqq Is it in fact SNP-wide? Try running it without skipping the remaining pairings and see.
      
      if(
        any( sapply( hapstrat_error, function(e){
          grepl( "Less than 2 haplotypes above threshold", e ) } ) )
      ){
        
        logr::put( "Error 1 identified.", console = FALSE )
        
        # Add error code to error_code column:
        res_sl_gxe_df$error_code <- paste0( res_sl_gxe_df$error_code, 1 )
        # Remove "9" from error_code column:
        res_sl_gxe_df$error_code <-
          gsub( "^9", "", as.character(res_sl_gxe_df$error_code) )
        
        # THIS ERROR DISQUALIFIES res_sl_hapstrat FROM GXE ANALYSIS
        # Update gxe_execution_allowed column accordingly:
        res_sl_gxe_df$gxe_execution_allowed <- FALSE
        
        # ERROR CODE 9 (OTHER/UNSPECIFIED) ----
        # If an error occurred, but not of the type outlined above, let
        # error_code = 9, meaning "other", unspecified error scenario
      } else{
        
        logr::put( "Unspecified error scenario registered. error_code = 9"
                   , console = FALSE )
        
        # Error code 9 was added in tryCatch block, but let's ensure that it's
        # there:
        if( grepl( "9", res_sl_gxe_df$error_code ) ){
          # Add error code to error_code column:
          res_sl_gxe_df$error_code <- paste0( res_sl_gxe_df$error_code, 1 )
          # Remove "NA" from error_code column:
          res_sl_gxe_df$error_code <-
            gsub( "^9$", "", as.character(res_sl_gxe_df$error_code) )
          
          # THIS ERROR DISQUALIFIES res_sl_hapstrat FROM GXE ANALYSIS
          # Update gxe_execution_allowed column accordingly:
          res_sl_gxe_df$gxe_execution_allowed <- FALSE
          
        }
        
      }
    }
    
    
    # ERROR CAPTURED *AFTER* RUNNING HAPLINSTRAT ----
    
    # If haplinStrat produced a haplinStrat object:
    # (At this point, res_sl_hapstrat should be either a data.table or a
    # haplinStrat object. If it's not a haplinStrat object then this is due to
    # errors occurred while running haplinStrat.)
    if( "haplinStrat" %in% class(res_sl_hapstrat) ){
      
      
      ## ERROR CODE 2 ----
      
      # The haplinStrat object contains one (or more) try-error object. I.e., a
      # stratum's element in the haplinStrat object contains a try-error object
      # instead of a haplin object with estimates etc.
      # (This one's a bit weird because haplinStrat does not throw an error when
      # this occurs.) 
      
      if( class( res_sl_hapstrat ) == "haplinStrat" &
          any(sapply( res_sl_hapstrat, class ) == "try-error" )
      ){
        # Get the  index of the list element(s) in res_sl_hapstrat belonging to
        # strata where there were problems:
        failed_strata <- which(sapply( res_sl_hapstrat, class ) == "try-error") 
        
        paste0( "Error 2 identified. Strata: "
                , paste0( names(failed_strata), collapse = ", " )
        ) %>% 
          logr::put()
        
        # Add error code to error_code column:
        res_sl_gxe_df$error_code <- paste0( res_sl_gxe_df$error_code, 2 )
        # Remove "NA" from error_code column:
        res_sl_gxe_df$error_code <-
          gsub( "^NA", "", as.character(res_sl_gxe_df$error_code) )
        
        # THIS ERROR DISQUALIFIES res_sl_hapstrat FROM GXE ANALYSIS
        # Update gxe_execution_allowed column accordingly:
        res_sl_gxe_df$gxe_execution_allowed <- FALSE
        
        
        # ERROR CODE 3 ----
        
        # A special case of Error 2, where it is confirmed that the try-error
        # object (or at least one of the try-error objects) in the haplinStrat
        # object consists of this error message:
        # "Error : Haplotypes * not found in file!"
        # This is a stop message from f.sel.haplos which is given if haplotype *
        # is not present in a stratum at all. (* = an integer such as 1, 2, etc)
        
        # If any failed strata contain the "not found in file" error, code this
        # SNP x SL pairing as error code 3
        if( any(
          sapply( res_sl_hapstrat[ failed_strata ], function(s){
            grepl( "Haplotypes \\d not found in file!", s[[1]]) } )
        ) ){
          
          logr::put( "Error 3 identified.")
          
          # Add error code to error_code column:
          res_sl_gxe_df$error_code <- paste0( res_sl_gxe_df$error_code, 3 )
          # Remove "NA" from error_code column:
          res_sl_gxe_df$error_code <-
            gsub( "^NA", "", as.character(res_sl_gxe_df$error_code) )
          
          # gxe_execution_allowed was updated above when registering error 2.
        }
        
        
        # Print the failed strata name and their corresponding try-error objects
        for( s in failed_strata ){
          paste0( "\ntry-error belonging to stratum "
                  , names(res_sl_hapstrat[s])
                  , ":" ) %>%
            logr::put()
          logr::put( res_sl_hapstrat[s][[1]] )
        }
        
        
      }
    }
    
    
    # INSPECT ERROR AND WARNING CODES ----
    
    # If res_sl_hapstrat is a data.table, then gxe_execution_allowed must be
    # FALSE:
    if(  "data.table" %in% class(res_sl_hapstrat) ){
      stopifnot( res_sl_gxe_df$gxe_execution_allowed == FALSE )
    }
    
    # If warnings and errors that "disqualify" this pairing from gxe() analysis
    # have been registered, then gxe_execution_allowed must be FALSE:
    if( grepl( "1|2|3", res_sl_gxe_df$error_code ) |
        grepl( "2", res_sl_gxe_df$warning_code ) ){
      stopifnot( res_sl_gxe_df$gxe_execution_allowed == FALSE )
    }
    
    # If we have a haplinStrat object and we did not register error 2, error 3
    # or warning 2, then gxe_execution_allowed should be NA, and we will set it
    # to TRUE so that we can run gxe() later.
    if(  "haplinStrat" %in% class(res_sl_hapstrat) &
         !grepl( "1|2|3", res_sl_gxe_df$error_code ) &
         !grepl( "2", res_sl_gxe_df$warning_code ) ){
      
      stopifnot( is.na( res_sl_gxe_df$gxe_execution_allowed ) )
      
      # Set to TRUE:
      res_sl_gxe_df$gxe_execution_allowed <- TRUE
    }
    
    # At this point, gxe_execution_allowed should not be NA:
    stopifnot(  !is.na(res_sl_gxe_df$gxe_execution_allowed) )
    
    
    # IF GXE_EXECUTION_ALLOWED = FALSE, SKIP GXE AND RETURN DATA.TABLE ----
    if(  res_sl_gxe_df$gxe_execution_allowed == FALSE ){
      
      # Add empty p-value column(s) to res_sl_gxe_df:
      res_sl_gxe_df <- cbind( 
        res_sl_gxe_df
        # Create dt with empty pval column(s)
        , data.table::fread(
          paste0( paste0( gxe_pvals, collapse = ";") # headers
                  , "
                                        " # newline
                  # Column values:
                  , paste0( rep("NA", length(gxe_pvals)), collapse = ";") 
          )
        )
      )
      
      # Check that there are columns with the names specified via gxe_pvals:
      stopifnot( all( gxe_pvals %in% colnames(res_sl_gxe_df) ) )
      
      # Return data.table without p-values and skip the rest of the iteration
      return( res_sl_gxe_df )
    }
    
    
    # IF GXE_EXECUTION_ALLOWED = TRUE ----
    
    # At this point, res_sl_hapstrat must be a haplinStrat object and all the
    # elements in res_sl_hapstrat must be a haplin object:
    stopifnot( "haplinStrat" %in% class(res_sl_hapstrat) & 
                 all( sapply(res_sl_hapstrat, class) == "haplin" ) )
    
    
    
    
    
    # RUN GXE()  ----
    
    # Additional logging:
    logr::put("Run Haplin::gxe()...", console = FALSE)
    
    res_sl_gxe <- Haplin::gxe( res_sl_hapstrat )
    res_sl_gxe <- res_sl_gxe[[1]] # The result df is inside a list
    
    # Additional logging:
    res_sl_gxe %>% logr::put( ., console = FALSE)
    
    
    # Check that the gxe results contain all the estimates we expect to find:
    stopifnot(
      "Object returned by Haplin::gxe() does not contain the expected p-value estimate(s)" =
        all( gxe_pvals %in% unique( res_sl_gxe$gxe.test ) )
    )
    
    
    # UPDATE res_sl_gxe_df WITH P-VALUE ----
    
    # Create dt with empty pval column(s)
    pvals_dt <- data.table::fread(
      paste0( paste0( gxe_pvals, collapse = ";") # headers
              , "
                                        " # newline
              # Column values:
              , paste0( rep("NA", length(gxe_pvals)), collapse = ";") 
      )
    )
    
    # Add p-value(s) from the gxe results to the dt with empty pval column(s):
    for( pval_name in colnames(pvals_dt) ){
      pvals_dt[, c(pval_name) :=
                 res_sl_gxe[ res_sl_gxe$gxe.test == pval_name,]$pval ]
      # pvals_dt[, c(pval_name) := res_sl_gxe %>%
      #            filter( gxe.test == pval_name ) %>%
      #            pull( pval )]
    }
    
    # Stop if a p-value column is still empty:
    stopifnot( all(!is.na(pvals_dt) ) )
    
    # Add the p-value column(s) to res_sl_gxe_df:
    # res_sl_gxe_df <- cbind( res_sl_gxe_df, pvals_dt )
    res_sl_gxe_df[,  names(pvals_dt) := pvals_dt ]
    
    # Check that res_sl_gxe_df has columns with the names given in gxe_pvals:
    stopifnot( all( gxe_pvals %in% colnames(res_sl_gxe_df) ) )
    
    
    # Check that each p-value column now contains the correct p-value from the
    # gxe results table:
    for( pval_name in gxe_pvals ){
      stopifnot(
        res_sl_gxe_df %>% pull( !!sym( pval_name ) ) ==
          res_sl_gxe$pval[ res_sl_gxe$gxe.test == pval_name ]
      )
    }
    
    
    # If either the error_code or warning_code column contains NA, replace NA
    # with 0 to indicate that no errors/warnings occurred:
    if( is.na( res_sl_gxe_df$error_code ) ){
      res_sl_gxe_df$error_code <- "0"
    }
    if( is.na( res_sl_gxe_df$warning_code ) ){
      res_sl_gxe_df$warning_code <- "0"
    }
    
    # RETURN DATA.TABLE UPDATED WITH SL-SPECIFC GXE RESULTS
    return( res_sl_gxe_df )
    
  }) # END OF LAPPLY
  
  # Additional logging:
  logr::sep( "All SNP x SL combinations have now been analysed with haplinStrat"
             , console = FALSE )
  
  
  # All SNP x SL combinations have now been analysed with haplinStrat
  
  # COMBINE ALL THE SNP x SL RESULTS INTO ONE DATA FRAME (data.table):
  sl_res <- data.table::rbindlist( sl_res )
  
  # Additional logging:
  sl_res %>% logr::put( ., console = FALSE )
  
  
  # Throw error if sl_res does not contain all SL IDs from sls_snp_id:
  if( !all( sls_snp_id %in% sl_res$sl_id ) ){
    stop( "lapply() did not process all of the SNP x SL pairings!\n\n"
          , "The following SL IDs are missing from the resulting data.table: "
          , paste0( setdiff(sls_snp_id, sl_res$sl_id), collapse = ", " )
    )
  }
  
  # THE OBJECT WE WANT THE analyse_snp FUNCTION TO RETURN, sl_res, IS NOW READY!
  # (A data frame containing one row per SNP x SL pairing for the SNP currently
  # being processed.)
  
  # ALL THAT REMAINS NOW IS SOME FINAL HOUSEKEEPING:
  
  # # LOG SNP EXECUTION TIME -----
  # 
  # # Get SNP analysis end time:
  # end_analyse_snp_time <- Sys.time()
  # 
  # logr::put( paste0(
  #   "SNP EXECUTION TIME: "
  #   , round( as.numeric( end_analyse_snp_time - start_analyse_snp_time
  #                        , units = "secs" ), 8)
  #   , " SECONDS" ) )
  
  # DELETE LOG FILE -----
  
  # Delete the log file if there were no errors registered during the execution
  # of analyse_snp().
  # If an error or warning was registered (using the error_code and
  # warning_code columns in sl_res), then the log file is retained as long
  # as there aren't too many log files in the log folder already. Users can
  # study these log files and use them to figure out why the errors/warnings
  # occurred.
  
  # A warning occurred if sl_res$warning_code contains an integer greater
  # than 0.
  
  log_limit <- n_cores + keep_max_x_log_files
  # log_limit <- n_cores + 500
  
  # (At any point during parallel processing, there will most likely be at least
  # n_cores log files in the log_dir since each process have their own log file.
  # If we choose a log_limit < n_cores then that would result in the log file
  # always being deleted at the end of the processing of the SNP.)
  
  # GET TOTAL NUMBER OF SNPS WHO HAVE LOG FILES IN LOG DIRECTORY 
  
  # Determine how many SNPs there are in the current chunk who have log files in
  # the log directory that are not due to errors/warnings
  # If there are only 1-3, this indicates that these are the last few SNPs
  # being analysed in this particular chunk.
  
  # Get names of all files in log directory:
  all_log_files <- list.files( log_dir )
  # Keep only those marked with the current chunk number;
  chunk_log_files <- 
    all_log_files[ grepl( paste0("_c", chunk_number, "_"), all_log_files ) ]
  
  # We want to discount the log files that are in the directory due to an error
  # occurring. 
  # (Note: There is no 100% fail-safe way of ensuring that we disregard all log
  # files that re in the directory due an error/warning when the number of log
  # files exceeds `log_limit`, but we will probably get a fairly accurate
  # picture most of the time.)
  
  # Create vector of the SNP number that is used in the file names (stems from
  # snp_order_overall)
  # Remove everything up until "snp_1234" part of log/msg file name:
  chunk_log_snp_numbers <- 
    gsub( paste0("^.{10,}_c", chunk_number, "_"), "", chunk_log_files )
  # Remove everything after the "snp_1234" part of log/msg file name:
  chunk_log_snp_numbers <- gsub( "_out_of_\\d+.+$", "", chunk_log_snp_numbers )
  
  
  # Create vector with one logical values per element in chunk_log_snp_numbers
  # indicating whether there is an .msg file whose name contains this SNP
  # number:
  chunk_log_snp_numbers_msg_status <- 
    sapply( chunk_log_snp_numbers, function(s){
      # Get logical value indicating whether there is an .msg file whose name
      # contains SNP number s:
      grepl(  "[[:punct:]]msg$",
              # Create string with regex containing all file names containing
              # SNP number:
              paste0( grep( s, chunk_log_files, value = TRUE  )
                      , collapse = "|" )
      )
    } )
  
  # Create variable with number of SNPs that have "non-error-related" log files
  # in the log directory
  n_chunk_snps_with_non_error_log_files <-
    length(chunk_log_snp_numbers_msg_status[
      which( chunk_log_snp_numbers_msg_status == FALSE )
    ])
  
  
  # Keep the log file in the log directory if 
  #   1) an error or warning was registered, OR
  #   2) at least one of the p-values is NA, AND
  #   3) the number of files already present in the log directory is below the
  #   specified limit
  # 1) and 2) should always coincide, but let's cover all our bases. This
  # ensures that log is retained even though error_code or warning_code does not
  # capture that something went wrong. NA in p-value column => worth
  # investigating.
  if( 
    # 1)
    # An error was registered:
    ( !all( unique( sl_res$error_code ) %in% c(NA, "0") ) |
      # A warning was registered:
      !all( unique( sl_res$warning ) %in% c(NA, "0") ) |
      # 2) 
      # At least one p-value is NA 
      base::anyNA( sl_res[ 
        , .SD
        , .SDcols = names(sl_res) %like% paste0( gxe_pvals, collapse = "|" )
        # sl_res %>% select( matches( paste0( gxe_pvals, collapse = "|" ) ) )
      ] )
    ) &
    # 3)
    # Number of files in log_dir is under specified upper limit:
    length( list.files( log_dir ) ) < log_limit ){
    
    # Print warning so that logr creates a .msg file to accompany the log
    # file, indicating that a warning/error has occurred
    warning("This SNP has one or more pairings that produced missing p-values.")
    
    # logr::put( paste0( "Missing p-values in `sl_res`."
    #                    , "\nSince there are are less than ", log_limit
    #                    , " files in the log folder, refrain from deleting the "
    #                    , "log file so that the user can study it." )
    #            , console = FALSE )
    
    # Close log:
    logr::log_close()
    
    # Keep the log file in the log directory if the SNP is either the first or
    # the last SNP in its respective chunk.
    # Keep the log file so that we can track the time between first and last SNP
    # as well as the time between chunks (and indicator of how long it takes for
    # the cluster to get to business)
  } else if( snp_id %in% first_last_snp ){
    
    logr::put( paste0( "\nKeep this log file as the SNP is either the first"
                       , " or the last SNP in its chunk.\n"
                       , "(Makes it possible to track how much time passes "
                       , "between chunks.)" )
               , console = FALSE )
    
    # NOTE: The last SNP according to snp_order_overall or snp_order_in_chunk
    # will often not be the last SNP being actively processed by the cluster. We
    # will try to catch any "straggler" SNPs in the next else if clause.
  } else if( n_chunk_snps_with_non_error_log_files < 3 ){
    # If there currently is only 1-3 chunk SNPs with log files in the log
    # directory that are not error-related, then it is probable that this SNP is
    # one of the very last SNPs being processed in this chunk.
    # We will therefore keep this SNP's log file so that we can more easily
    # track 1) whether there were "straggler" SNPs with many pairings slowing
    # the chunk execution time down, and 2) the time between when the last chunk
    # SNP was completed and when the next chunk started being processed by the
    # cluster.
    
    logr::put( paste0( "\nKeep this log file as it appeared to be one of the"
                       , " last SNPs to be processed in its chunk.\n"
                       , "(Makes it possible to track how much time passes "
                       , "between chunks, as well keeping an eye out for "
                       , "\"straggler\" SNPs slowing down the chunk execution "
                       , "time." )
               , console = FALSE )
    
  } else{
    # Remove the log file otherwise, i.e. if either 1) no error was registered,
    # or 2) an error was registered, but there are too many existing log files
    # in the log folder.
    
    # Close log:
    logr::log_close()
    
    # Delete .log file
    unlink( log_path, force = TRUE )
    
    # Delete .msg file (if any)
    if( file.exists( gsub( "log$", "msg", log_path) ) ){
      unlink( gsub( "log$", "msg", log_path), force = TRUE )
    }
  }
  # }, interval =  0.005, rerun = TRUE) # ****** PROFVIS end
  
  # Additional printout to console:
  cat( "\nTest of SNP using analyse_snp_for_testing_pre_execution() complete.")
  
  # Return data frame with the results of all the SNP x SL pairings:
  return( sl_res )
} # End of analyse_snp()
