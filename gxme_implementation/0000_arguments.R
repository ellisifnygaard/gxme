args <- list(
  
  #*** DIRECTORIES 
  
  # PIPELINE DIRECTORY
  # Path to pipeline directory (root directory for all processing and analysis)
  # pipeline_dir = "S:/Project/SAM/Julia/Ellisif/pipeline_dir"
  #* Should not be necessary if here is working as intended
  
  # CHROMOSOME DIRECTORY
  # (The root folder which will contain all the chromosome-specific 
  # subdirectories) 
  # chromosome_dir: here::here( "chromosomes")
  #* Just use default = "pipeline_dir/chromosomes"
  

  #*** CHROMOSOMES TO BE STUDIED
  chr_numbers = 1:22
  # chr_numbers = c(1,4,5,7,18,21)

    
  #*** DNAM-RELATED
  
  # DNAM ARRAY DATA
  , ewas_fileset_dir = "path/to/directory/containing/dnam/array/data"
  # (Path to the directory where the DNAm array data filesets in HaplinMethyl's
  # ffarchive file format are located. There must be one fileset per chromosome
  # being studied.)
  , ewas_fileset_name_base_regex = "_dnam_data_chr[[:digit:]]{2}"
  # (Regular expression that matches the basenames of the files with the DNAm
  # array data, *excluding the "_env" suffix that HaplinMethyl automatically
  # adds the filenames. There needs to be one fileset per chromosome.)
  
  , ewas_map_file_dir = "path/to/directory/containing/map/files"
  # (Path to the directory where the map files accompanying the DNAm array data
  # are located. There must be one .feather file per chromosome being studied
  # containing the following columns: 1) `probe_id`, containing character
  # strings with the IDs of the DNAm site probes; 2) `chr` containing the
  # chromosome number in the form of an integer; and 3) `coord`, containing the
  # 1-base coordinate/single nucleotide position of the cytosine base in the
  # DNAm sites, in the form of an integer.)
  , ewas_map_file_name_base_regex = "^ewas_map_hm450k_hg19_chr[[:digit:]]{1,2}"
  # (Regular expression that matches the basenames of the map files with the
  # genomic coordinates of the DNAm sites (basename = the filename excluding the
  # ".feather" file extension.)
  
  # EWAS FAMILY MEMBER
  , ewas_family_member = "c"
  # The implementation currently only supports "c" for child. This was
  # introduced as a parameter to facilitate using parent DNAm measurements as
  # the stratifying layer in future versions (if there is a statistical model
  # forming a basis for this)
  
  # EWAS PROBE ANNOTATION FILE
  , ewas_annotation_file_path =
    "./gxme/annotation_file/hg19_genome_100_segments.bed.gz"
  # (Path to the .bed.gz file containing the annotation file. )

  #*** SNP ARRAY DATA
  
  # SNP ARRAY FILES FILE
  , gwas_fileset_dir = "path/to/directory/containing/snp/array/data"
  # (The path to the directory where the binary PLINK filesets containing the
  # SNP genotype data are located.)
  , gwas_fileset_name_base_regex = "^Chr[[:digit:]]{1,2}_0006"
  # (A regular expression to be passed on as the 'pattern' argument in 
  # base::list.files() when identifying which PLINK filesets to copy
  # to the chromosome-specific subdirectories in the pipeline directory.)
  
  # GWAS QC
  # QC-related parameters:
  , maf_minimum = 0.01 # (for removing monomorphic SNPs)
  
  #*** FUNCTIONS
  
  # PATH(S) TO SUMMARISING FUNCTION(S)
  , summarising_function_path =
    c( "./gxme/summarising_and_stratifying_functions/0501_summarising_function_default.R" 
       , "./gxme/summarising_and_stratifying_functions/0501_summarising_function_median.R")
  # (You can specify multiple paths/R files to copy to the pipeline directory,
  # but you can only use one of them during a pipeline execution.)
  
  # NAME OF R FILE WITH THE SUMMARISING FUNCTION YOU WANT TO USE, OR "default"
  , sl_summarising_function = "default"
  # (If sl_summarising_function = "default" then the pipeline will use
  # 0501_summarising_function_default.R from /summarising_and_stratifying_functions/. Note that you need to
  # include a path to this function in `summarising_function_path`)
  # , sl_summarising_function = "0501_summarising_function_median.R"
  
  # PATH(S) TO STRATIFYING FUNCTION(S)
  , stratifying_function_path =
    "./gxme/summarising_and_stratifying_functions/0501_stratifying_function_default.R"
  # (You can specify multiple paths/R files to copy to the pipeline directory,
  # but you can only use one of them during a pipeline execution.)
  
  # NAME OF R FILE WITH THE STRATIFYING FUNCTION YOU WANT TO USE, OR "default"
  , sl_stratifying_function = "default"
  # (If sl_stratifying_function = "default" then the pipeline will use
  # 0501_stratifying_function_default.R from /summarising_and_stratifying_functions/. Note that you need to
  # include a path to this function in `stratifying_function_path`.)
  
  #*** MQTL
  
  # MQTL FILE(S)
  , mqtl_file_name_base_regex =
    "godmc_mqtl_associations_hg19_chr[[:digit:]]{1,2}"
  , mqtl_file_dir =
    "path/to/directory/containing/mqtl/files"
  # The path to a folder containing all the files with the user's curated list
  # of flagged mQTL association pairs. There must be one feather file per
  # chromosome being studied, i.e., one feather file per chromosome given in
  # chr_numbers.
  
  # MQTL LD LOWER LIMIT
  , mqtl_ld_r2 = 0.9
  # (SNPs with pairwise R2 greater than or equal to mqtl_ld_r2 are defined as 
  # being in high LD with one another)
  
  #*** PLINK
  
  # PATH TO PLINK EXECUTABLE
  # PLINK version: PLINK v1.90b6.26 64-bit (2 Apr 2022)
  , plink_exe_path =
    "./PLINK/plink_win64_20220402/plink.exe"
  
  , plink_memory_mb = 1800 
  # (Can also be NA)
  
  # PLINK timeout 
  # (maximum number of seconds a call to plink.exe can take before it's stopped)
  , plink_timeout = 240
  
  
  #*** STRATIFICATION SCHEME
  
  # All full-stack chromatin states except GapArtf (assembly gaps and alignment
  # artefacts):
  , scheme_states = c( "quies01" = "3_Quies1"
                       , "quies02" = "4_Quies2"
                       , "quies03" = "5_Quies3"
                       , "quies04" = "6_Quies4"
                       , "quies05" = "7_Quies5"
                       , "het01" = "8_HET1"
                       , "het02" = "9_HET2"
                       , "het03" = "10_HET3"
                       , "het04" = "11_HET4"
                       , "het05" = "12_HET5"
                       , "het06" = "13_HET6"
                       , "het07" = "14_HET7"
                       , "het08" = "15_HET8"
                       , "het09" = "16_HET9"
                       , "repr_pol01" = "17_ReprPC1"
                       , "repr_pol02" = "18_ReprPC2"
                       , "repr_pol03" = "19_ReprPC3"
                       , "repr_pol04" = "20_ReprPC4"
                       , "repr_pol05" = "21_ReprPC5"
                       , "repr_pol06" = "22_ReprPC6"
                       , "repr_pol07" = "23_ReprPC7"
                       , "repr_pol08" = "24_ReprPC8"
                       , "repr_pol09" = "25_ReprPC9"
                       , "acet01" = "26_Acet1"
                       , "acet02" = "27_Acet2"
                       , "acet03" = "28_Acet3"
                       , "acet04" = "29_Acet4"
                       , "acet05" = "30_Acet5"
                       , "acet06" = "31_Acet6"
                       , "acet07" = "32_Acet7"
                       , "acet08" = "33_Acet8"
                       , "enh_weak01" = "34_EnhWk1"
                       , "enh_weak02" = "35_EnhWk2"
                       , "enh_weak03" = "36_EnhWk3"
                       , "enh_weak04" = "37_EnhWk4"
                       , "enh_weak05" = "38_EnhWk5"
                       , "enh_weak06" = "39_EnhWk6"
                       , "enh_weak07" = "40_EnhWk7"
                       , "enh_weak08" = "41_EnhWk8"
                       , "enh_act01" = "42_EnhA1"
                       , "enh_act02" = "43_EnhA2"
                       , "enh_act03" = "44_EnhA3"
                       , "enh_act04" = "45_EnhA4"
                       , "enh_act05" = "46_EnhA5"
                       , "enh_act06" = "47_EnhA6"
                       , "enh_act07" = "48_EnhA7"
                       , "enh_act08" = "49_EnhA8"
                       , "enh_act09" = "50_EnhA9"
                       , "enh_act10" = "51_EnhA10"
                       , "enh_act11" = "52_EnhA11"
                       , "enh_act12" = "53_EnhA12"
                       , "enh_act13" = "54_EnhA13"
                       , "enh_act14" = "55_EnhA14"
                       , "enh_act15" = "56_EnhA15"
                       , "enh_act16" = "57_EnhA16"
                       , "enh_act17" = "58_EnhA17"
                       , "enh_act18" = "59_EnhA18"
                       , "enh_act19" = "60_EnhA19"
                       , "enh_act20" = "61_EnhA20"
                       , "tr_enh01" = "62_TxEnh1"
                       , "tr_enh02" = "63_TxEnh2"
                       , "tr_enh03" = "64_TxEnh3"
                       , "tr_enh04" = "65_TxEnh4"
                       , "tr_enh05" = "66_TxEnh5"
                       , "tr_enh06" = "67_TxEnh6"
                       , "tr_enh07" = "68_TxEnh7"
                       , "tr_enh08" = "69_TxEnh8"
                       , "tr_weak01" = "70_TxWk1"
                       , "tr_weak02" = "71_TxWk2"
                       , "tr_str01" = "72_Tx1"
                       , "tr_str02" = "73_Tx2"
                       , "tr_str03" = "74_Tx3"
                       , "tr_str04" = "75_Tx4"
                       , "tr_str05" = "76_Tx5"
                       , "tr_str06" = "77_Tx6"
                       , "tr_str07" = "78_Tx7"
                       , "tr_str08" = "79_Tx8"
                       , "exon01" = "80_TxEx1"
                       , "exon02" = "81_TxEx2"
                       , "exon03" = "82_TxEx3"
                       , "exon04" = "83_TxEx4"
                       , "znf01" = "84_znf1"
                       , "znf02" = "85_znf2"
                       , "dnase" = "86_DNase1"
                       , "prom_biv01" = "87_BivProm1"
                       , "prom_biv02" = "88_BivProm2"
                       , "prom_biv03" = "89_BivProm3"
                       , "prom_biv04" = "90_BivProm4"
                       , "prom_fl01" = "91_PromF1"
                       , "prom_fl02" = "92_PromF2"
                       , "prom_fl03" = "93_PromF3"
                       , "prom_fl04" = "94_PromF4"
                       , "prom_fl05" = "95_PromF5"
                       , "prom_fl06" = "96_PromF6"
                       , "prom_fl07" = "97_PromF7"
                       , "tss01" = "98_TSS1"
                       , "tss02" = "99_TSS2"
  )
  
  # CUT-OFF (PROXIMITY THRESHOLD)
  
  # PROXIMITY THRESHOLD CAN BE TAILORED TO EACH SCHEME
  # BUT THEN EACH SCHEME # (I.E. UNIQUE ANNOTATION/LABEL) MUST BE LISTED
  # SPECIFICALLY WITH A VALID CUT OFF
  # , cut_off = c( "quies01" = 1000000
  #                , "quies02" = 1000000
  #                , "quies03" = 75000
  #                , "quies04" = 500000
  #                , "quies05" = 1000000
  #                , "het01" = 50000
  #                , "het02" = 50000
  #                , ...
  #                , "repr_pol01" = 750000
  #                , "repr_pol02" = 1000000
  #                , ...
  # )
  # OR ONE CAN USE THE SAME PROXIMITY THRESHOLD FOR ALL SCHEMES
  , cut_off = 1000000 # universal cut-off
  
  # HAPLIN ARGUMENTS
  
  # All arguments with the prefix "haplin_" are arguments that will be passed to
  # Haplin::haplin() via Haplin::haplinStrat(). See the Haplin reference manual
  # on CRAN for more information. (The part of the argument name following
  # "haplin" is the name of the haplin argument.)
  
  # Design used in your study (triad, cc, cc.triad; see description of haplin()
  # in the Haplin Reference manual)
  , haplin_design = "triad"
  # Whether to use triads with missing data (uses the EM algorithm):
  , haplin_use.missing = TRUE
  , haplin_xchrom = FALSE # use the default
  , haplin_maternal = FALSE # use the default
  , haplin_test.maternal = FALSE # use the default
  , haplin_poo = FALSE # use the default
  , haplin_scoretest = "no" # use the default
  , haplin_ccvar = NULL # use the default
  , haplin_sex = NULL # use the default
  , haplin_comb.sex = "double" # use the default
  , haplin_reference = "ref.cat" 
  # (when response = "mult", reference must be "ref.cat")
  , haplin_response = "mult" # haplin default is "free"
  , haplin_threshold = 0.01
  , haplin_max.haplos = NULL
  , haplin_haplo.file = NULL
  , haplin_resampling = "no"
  # Maximum number of iterations used by the EM algorithm:
  , haplin_max.EM.iter = 50
  # , haplin_data.out = "no" # Will be hard-coded in the pipeline
  # , haplin_verbose = TRUE # Will be hard-coded
  # , haplin_printout = TRUE # Will be hard-coded in the pipeline
  
  # SELECTED P-VALUES FROM HAPLIN'S GXE RESULTS TABLE:
  , gxe_pvals = c( "child", "child.trend" )
  # Vector with strings containing the names of p-values from the results of
  # Haplin::gxe() to return. 
  # Possible choices when haplin_poo = FALSE: "haplo.freq", "child",
  # "haplo.freq.trend", "child.trend" 
  # Possible choices when haplin_poo = TRUE: "haplo.freq", "poo",
  # "haplo.freq.trend", "poo.trend" 
  # If nothing is specified, the used default will be 
  # gxe_pvals = ifelse( haplin_poo == FALSE, yes = "child", no = "poo" ).
  
  
  # ARGUMENT SIGNALLING WHETHER THE ANALYSES ARE REPLICATION ANALYSES
  , replication_analysis = FALSE
  # (This parameter is only listed here due to plans for future functionality.
  # Do not change from `FALSE`.)
)
