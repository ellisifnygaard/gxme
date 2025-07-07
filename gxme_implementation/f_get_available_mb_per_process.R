# This code shows one way of determining the amount of physical memory per
# process that a Windows computer has available (in MB). The returned value can
# be passed as the `plink_memory_mb` argument.

# DETERMINE AVAILABLE MEMORY FOR PLINK PER PROCESS (plink_memory_mb argument):

# get_available_mb_per_process( n_cores = parallel::detectCores() - 2, buffer_percent = 5 )
# get_available_mb_per_process( n_cores = 22, buffer_percent = 10 )

get_available_mb_per_process <- function( 
    # Number of cores to be used during parallel processing:
  n_cores = parallel::detectCores() - 2
  # How many percent of the total available memory would you like to reserve
  # as a buffer (integer between 0 and 99):
  , buffer_percent = 10
){
  # n_cores must be a single whole number:
  stopifnot( length(n_cores) == 1 )
  source( here( "R", "f_is_whole_number.R") )
  stopifnot( is_whole_number(n_cores) )
  
  # buffer_percent must be whole number in [0,99]:
  stopifnot( is_whole_number( buffer_percent ) & buffer_percent %in% 0:99 )
  
  # Get total physical memory:
  total_physical_memory <-
    system( "wmic ComputerSystem get TotalPhysicalMemory", intern = TRUE)
  total_physical_memory <- as.numeric(
    gsub("[^0-9]", "", grep("[0-9]", total_physical_memory, value = TRUE) )
  )
  # Total physical memory in MB:
  total_physical_memory_mb <- total_physical_memory / (1024^2)
  # Create a `buffer_percent` % buffer to subtract from total in order to get
  # total *available* memory:
  buffer <- ( total_physical_memory_mb / 100 ) * buffer_percent
  
  # Calculate the amount of RAM in MB that PLINK can reserve in each of the
  # processes:
  available_memory_mb_per_process <- 
    floor( (total_physical_memory_mb - buffer) / n_cores )
  message( "Allow PLINK to reserve "
           , available_memory_mb_per_process
           , " MB in each of the processes when using "
           , n_cores
           , " cores."
           )
  
  # Warning if RAM per process is "dangerously" low:
  if( available_memory_mb_per_process < 300 ){
    warning( "'available_memory_mb_per_process' is very low: "
             , available_memory_mb_per_process )
  }
  
  return( available_memory_mb_per_process )
}
