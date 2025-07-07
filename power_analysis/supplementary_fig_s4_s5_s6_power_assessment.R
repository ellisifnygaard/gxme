
library(dplyr)
library(data.table)
library(Haplin)
library(magrittr)
library(ggplot2)
library(scales)
library(grid)

wd <- getwd()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GENERATE POWER CALCULATIONS FOR DIFFERENT SAMPLE SIZES, HAPLO.FREQ AND ALPHA ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SUPPLEMENTARY FIGURES S4, S5 and S6 ALL SHOW RESULTS FROM CALCULATIONS WHERE
# THERE WERE 100 TRIADS PER STRATUM. HOWEVER, THIS CODE PERFORMS CALCULATIONS AT
# SEVERAL DIFFERENT SAMPLE SIZES. 
# REPLACE `n_triads_per_stratum = c(50, 100, 150, 200, 250, 300, 350 )`
# WITH `n_triads_per_stratum = 100 ` IF YOU ONLY WANT TO CALCULATE POWER FOR
# THIS SAMPLE SIZE (THIS GOES MUCH FASTER).



# CREATE DATA TABLE WITH ALL COMBINATIONS OF PARAMETERS

# Create dt with all possible combinations of the number of triads per stratum,
# allele frequency of the most common allele in stratum 1, 2, and 3, and alpha:
n_triads_haplo_freq_alpha_combinations <- base::expand.grid(
  # Constant parameters:
  n_all = 2
  , n_strata = 3
  , RR_stratum_1 = 3
  , RR_stratum_2 = 1
  , RR_stratum_3 = 1/3
  , RR_star_lower = 1
  , RR_star_upper = 1
  , response = "mult"
  # Varying parameters:
  , n_triads_per_stratum = c(50, 100, 150, 200, 250, 300, 350 )
  , haplo_freq_most_common_allele_stratum_1 =
    c(0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975)
  , haplo_freq_most_common_allele_stratum_2 =
    c(0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975)
  , haplo_freq_most_common_allele_stratum_3 =
    c(0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975)
  , alpha = c( 0.05 # "baseline"
               , 0.05*10^(-3) # (1,000 pairings)
               , 0.05*10^(-4) # (10,000 pairings)
               , 0.05*10^(-5) # (100,000 pairings)
               , 0.05*10^(-6) # (1,000,000 pairings)
               , 0.05*10^(-8) # (100,000,000 pairings)
  )
) %>% as.data.table()

dim(n_triads_haplo_freq_alpha_combinations)
summary(n_triads_haplo_freq_alpha_combinations)

# Check that the data set only contains unique combinations:
n_triads_haplo_freq_alpha_combinations %>%
  group_by( n_triads_per_stratum
            , haplo_freq_most_common_allele_stratum_1
            , haplo_freq_most_common_allele_stratum_2
            , haplo_freq_most_common_allele_stratum_3
            , alpha
  ) %>%
  mutate( n_rows_per_unique_combo = n() ) %>%
  ungroup() %>%
  distinct( n_triads_per_stratum
            , haplo_freq_most_common_allele_stratum_1
            , haplo_freq_most_common_allele_stratum_2
            , haplo_freq_most_common_allele_stratum_3
            , alpha
            , n_rows_per_unique_combo
  ) %>%
  count( n_rows_per_unique_combo )




# RUN HAPPOWERASYMP FOR EACH COMBINATION OF PARAMETERS:

power_res <-
  lapply( 1:nrow(n_triads_haplo_freq_alpha_combinations), function(i){

    row_i <- n_triads_haplo_freq_alpha_combinations[i,]

    power_gxe_res <- hapPowerAsymp(
      nall = row_i$n_all
      , n.strata = row_i$n_strata
      , cases = list( c( mfc = row_i$n_triads_per_stratum )
                      , c( mfc = row_i$n_triads_per_stratum )
                      , c( mfc = row_i$n_triads_per_stratum ) )
      , haplo.freq = list(
        c( row_i$haplo_freq_most_common_allele_stratum_1
           , 1 - row_i$haplo_freq_most_common_allele_stratum_1 )
        , c( row_i$haplo_freq_most_common_allele_stratum_2
             , 1 - row_i$haplo_freq_most_common_allele_stratum_2 )
        , c( row_i$haplo_freq_most_common_allele_stratum_3
             , 1 - row_i$haplo_freq_most_common_allele_stratum_3 ) )
      , RR = list( c(1, row_i$RR_stratum_1)
                   , c(1, row_i$RR_stratum_2)
                   , c(1, row_i$RR_stratum_3) )
      , RRstar = c(row_i$RR_star_lower, row_i$RR_star_upper)
      , response = row_i$response
      , alpha = row_i$alpha
    )

    # Add RR.power for Haplotype 2 to row_i:
    row_i$RR.power_Haplotype_2 <- power_gxe_res$haplo.power %>%
      filter(Haplotype == 2) %>%
      pull(RR.power) %>%
      as.numeric()

    return( row_i )

  })

length(power_res)
power_res %<>% dplyr::bind_rows(power_res)
gc()

power_res %>% dim()
power_res %>% summary()

# Export data table:
power_res %>%
  readr::write_csv(
    .
    , file = file.path(
      wd
      , "hapPowerAsymp_results_varying_sample_size_haplo_freq_and_alpha.csv"
    )
  )

power_res



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORT DATA WITH WITHIN-STRATUM MAF THAT VARIES ACROSS STRATA  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat <- fread(
  file.path(
    wd
    , "hapPowerAsymp_results_varying_sample_size_haplo_freq_and_alpha.csv"
  )
)
dat
dat %>% summary()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  PROCESS DATA BEFORE PLOTTING  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dat %<>% 
  mutate(
    maf_stratum_1 = 1 - haplo_freq_most_common_allele_stratum_1
    , maf_stratum_2 = 1 - haplo_freq_most_common_allele_stratum_2
    , maf_stratum_3 = 1 - haplo_freq_most_common_allele_stratum_3
    , min_maf_across_strata = pmin(maf_stratum_1, maf_stratum_2, maf_stratum_3)
    , max_maf_across_strata = pmax(maf_stratum_1, maf_stratum_2, maf_stratum_3)
  )

dat %>% summary()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     SET COMMON GGPLOT ARGUMENTS  (ALPHA AS GROUP) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Palette:
# Okabe–Ito eight-colour palette
okabe_ito <- c(
  "#E69F00"
  , "#56B4E9"
  , "#009E73"
  , "#F0E442"
  , "#CC79A7"
  , "#D55E00"
  , "#0072B2"
  , "#000000"
)

geom_point_size <- 2.5
geom_line_width <- 0.6
caption_width_fig <- 190
caption_size <- 8.5

# Build factor and expression labels for α:
alpha_vals   <- sort(unique(dat$alpha), decreasing = TRUE) # numeric
alpha_breaks <- scientific(alpha_vals, digits = 1) # “5e-02”, …

# Convert to desired math-style expressions:
alpha_lab_chr <- vapply(alpha_vals, function(a) {
  # If alpha = 0.05 :
  if (abs(a - 0.05) < 1e-12) { 
    "0.05"
    # If alpha 0.05·10^k, k = -2, -4, -5, -6, -8 :
  } else { 
    k <- as.integer( round( log10(a / 0.05) ) )
    paste0("0.05%*%10^", k) # e.g. "0.05*10^-3"
  }
}, FUN.VALUE = character(1) )

# Turn into plotmath expressions:
alpha_labels <- parse(text = alpha_lab_chr)   

# Labeller for facet grid:
maf_labeller <- labeller(
  # Stratum 2 label:
  maf_stratum_2 = function(x)
    paste0("MAF in\nstratum 2 = ", x)
  # Stratum 3 label:
  , maf_stratum_3 = function(x)
    paste0("MAF in\nstratum 3 = ", x)
)


# Add column with alpha as factor:
dat %<>% 
  mutate( alpha_f = factor(alpha, levels = alpha_vals, labels = alpha_breaks) )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POWER VS MAF IN STRATUM 1 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GENERATE ONE PLOT FOR EACH SELECTED NUMBER OF TRIADS PER PLOT:

# # If you want one plot per n_triads_per_stratum:
# lapply( unique(dat$n_triads_per_stratum), function(x){
# If just want plots where n_triads_per_stratum = 100
lapply( c(100), function(x){
  
  
  no_of_triads_per_stratum <- x
  
  
  # Labeller for facet grid:
  maf_labeller <- labeller(
    maf_stratum_2 = function(x)
      paste0("MAF in\nstratum 2 = ", x),
    
    maf_stratum_3 = function(x)
      paste0("MAF in\nstratum 3 = ", x)
  )
  
  # Figure title
  fig_title <- paste0(
    "Power at different significance levels (\u03b1) across different combinations of within-stratum minor allele frequency (MAF)"
    , "\nin stratum 1, 2 and 3 ("  
    ,  no_of_triads_per_stratum
    , " triads in each stratum)"
  )
  
  # Figure caption
  fig_caption <- stringr::str_wrap(
    paste0(
      "Asymptotic power for single-SNP GxE analyses as calculated by 
    hapPowerAsymp() when using a diallelic SNP, the case-parent triad design, 
    the multiplicative dose-response model, and three strata."
      , " Sample size held constant at "
      , unique(dat$n_triads_per_stratum) * 3
      , " ("
      , unique(dat$n_triads_per_stratum)
      , " triads per stratum)."
      , " RR equals 3 in stratum 1, 1.0 in stratum 2, and 1/3 in stratum 3. 
    \"Minor allele frequency\" refers to the frequency of the least common 
    allele in all three strata." )
    , width = caption_width_fig
  )
  
  
  # CREATE PLOT:
  p <- dat %>% 
    filter( n_triads_per_stratum == no_of_triads_per_stratum ) %>% 
    filter( !(maf_stratum_1 %in% c(0.25, 0.35) )
            , !(maf_stratum_2 %in% c(0.25, 0.35) )
            , !(maf_stratum_3 %in% c(0.25, 0.35) )
    ) %>% 
    ggplot( .
            , aes( x = maf_stratum_1
                   , y = RR.power_Haplotype_2
                   , colour = alpha_f
                   , linetype = alpha_f
                   , group = alpha_f )
    ) +
    geom_point( 
      alpha = 0.4
      , stroke = 0 
      , size = geom_point_size
    ) +  
    geom_line(linewidth = geom_line_width
              , alpha = 0.8 ) +
    scale_x_continuous(
      name   = "MAF in stratum 1",
      limits = c(0, 0.4),
      breaks = seq(0, 0.4, by = 0.05),
    ) +
    scale_y_continuous(
      name   = "Power",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    
    # Colour scale:
    scale_colour_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = okabe_ito[seq_along(alpha_breaks)],
      name   = expression(alpha)
    ) +
    # Linetype scale:
    scale_linetype_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = c("solid","22","42","46", "13", "1343" ),
      name   = expression(alpha)
    ) +
    
    # Title and caption:
    labs(
      title = fig_title
      , caption = fig_caption
    ) +
    
    # Theme:
    theme_bw(base_size = 11) +
    theme(
      axis.title.x   = element_text(margin = margin(t = 5, b = 7))
      , axis.title.y   = element_text(margin = margin(r = 10))
      
      , legend.box = "vertical"
      , legend.key.width = unit(44, "pt") # wider keys to show linetypes
      , legend.position = "right"
      , legend.title = element_text(face = "bold"
                                    , hjust = 0.5
                                    , size = 12 )
      , legend.margin = margin(0,0,0,0)
      , legend.text = element_text( margin = margin(t = 7, r = 0, b = 7, l = 2 )
                                    , size = 8 )
      
      , axis.text.x = element_text( angle = 90
                                    , vjust = 0.5
                                    , size = 8 )
      , axis.text.y = element_text( size = 8 )
      
      , plot.caption = 
        element_text( hjust = 0
                      , size = caption_size
                      , colour = "gray28"
                      , margin = margin( t = 0, r = 0, b = 0, l = -18, "pt" )
        )
      , plot.title = element_text( hjust = 0
                                   , size = 12
      )
      , strip.text = element_text(size = 9)
    ) +
    facet_grid(
      rows = vars(maf_stratum_2)
      , cols = vars( maf_stratum_3 )
      , labeller = maf_labeller
    )
  p
  
  # EXPORT
  
  p_filename <- paste0(
    "power_vs_maf_stratum_1_in_maf_stratum_2_x_maf_stratum_3_grid_"
    , no_of_triads_per_stratum
    , "_triads_per_stratum"
  )
  
  p %>% 
    ggsave( plot = .
            , filename = paste0( p_filename, ".tiff" )
            , dpi = 300
            , width = 29
            , height = 26
            , units = "cm"
            , device = "tiff"
    )
  
  p %>% 
    ggsave( plot = .
            ,  filename = paste0( p_filename, ".png" )
            , dpi = 600
            , width = 29
            , height = 26
            , units = "cm"
    )
  
} )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POWER VS MAF IN STRATUM 2 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# GENERATE ONE PLOT FOR EACH SELECTED NUMBER OF TRIADS PER PLOT:

# # If you want one plot per n_triads_per_stratum:
# lapply( unique(dat$n_triads_per_stratum), function(x){
# If just want plots where n_triads_per_stratum = 100
lapply( c(100), function(x){
  
  no_of_triads_per_stratum <- x
  
  # Labeller for facet grid:
  maf_labeller <- labeller(
    # Stratum 1 label:
    maf_stratum_1 = function(x)
      paste0("MAF in\nstratum 1 = ", x)
    # Stratum 3 label:
    , maf_stratum_3 = function(x)
      paste0("MAF in\nstratum 3 = ", x)
  )
  
  
  # Figure title
  fig_title <- paste0(
    "Power at different significance levels (\u03b1) across different combinations of within-stratum minor allele frequency (MAF)"
    , "\nin stratum 1, 2 and 3 ("  
    ,  no_of_triads_per_stratum
    , " triads in each stratum)"
  )
  
  # Figure caption
  fig_caption <- stringr::str_wrap(
    paste0(
      "Asymptotic power for single-SNP GxE analyses as calculated by 
    hapPowerAsymp() when using a diallelic SNP, the case-parent triad design, 
    the multiplicative dose-response model, and three strata."
      , " Sample size held constant at "
      , unique(dat$n_triads_per_stratum) * 3
      , " ("
      , unique(dat$n_triads_per_stratum)
      , " triads per stratum)."
      , " RR equals 3 in stratum 1, 1.0 in stratum 2, and 1/3 in stratum 3. 
    \"Minor allele frequency\" refers to the frequency of the least common 
    allele in all three strata." )
    , width = caption_width_fig
  )
  
  
  # CREATE PLOT:
  p <- dat %>% 
    filter( n_triads_per_stratum == no_of_triads_per_stratum ) %>% 
    filter( !(maf_stratum_1 %in% c(0.25, 0.35) )
            , !(maf_stratum_2 %in% c(0.25, 0.35) )
            , !(maf_stratum_3 %in% c(0.25, 0.35) )
    ) %>% 
    ggplot( .
            , aes( x = maf_stratum_2
                   , y = RR.power_Haplotype_2
                   , colour = alpha_f
                   , linetype = alpha_f
                   , group = alpha_f )
    ) +
    geom_point( 
      alpha = 0.4
      , stroke = 0 
      , size = geom_point_size
    ) +  
    geom_line(linewidth = geom_line_width
              , alpha = 0.8 ) +
    
    scale_x_continuous(
      name   = "MAF in stratum 2",
      limits = c(0, 0.4),
      breaks = seq(0, 0.4, by = 0.05),
    ) +
    scale_y_continuous(
      name   = "Power",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    
    # Colour scale:
    scale_colour_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = okabe_ito[seq_along(alpha_breaks)],
      name   = expression(alpha)
    ) +
    # Linetype scale:
    scale_linetype_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = c("solid","22","42","46", "13", "1343" ),
      name   = expression(alpha)
    ) +
    
    # Title and caption:
    labs(
      title = fig_title
      , caption = fig_caption
    ) +
    
    # Theme:
    theme_bw( base_size = 11 ) +
    theme(
      axis.title.x   = element_text(margin = margin(t = 5, b = 7))
      , axis.title.y   = element_text(margin = margin(r = 10))
      
      , legend.box = "vertical"
      , legend.key.width = unit(44, "pt") # wider keys to show linetypes
      , legend.position = "right"
      , legend.title = element_text(face = "bold"
                                    , hjust = 0.5
                                    , size = 12 )
      , legend.margin = margin(0,0,0,0)
      , legend.text = element_text( margin = margin(t = 7, r = 0, b = 7, l = 2 )
                                    , size = 8 )
      
      , axis.text.x = element_text( angle = 90
                                    , vjust = 0.5
                                    , size = 8 )
      , axis.text.y = element_text( size = 8 )
      
      , plot.caption = 
        element_text( hjust = 0
                      , size = caption_size
                      , colour = "gray28"
                      , margin = margin( t = 0, r = 0, b = 0, l = -18, "pt" )
        )
      , plot.title = element_text( hjust = 0
                                   , size = 12
      )
      , strip.text = element_text(size = 9)
    ) +
    facet_grid(
      rows = vars(maf_stratum_1)
      , cols = vars( maf_stratum_3 )
      , labeller = maf_labeller
    )
  p
  
  # EXPORT
  
  p_filename <- paste0(
    "power_vs_maf_stratum_2_in_maf_stratum_1_x_maf_stratum_3_grid_"
    , no_of_triads_per_stratum
    , "_triads_per_stratum"
  )
  
  p %>% 
    ggsave( plot = .
            , filename = paste0( p_filename, ".tiff" )
            , dpi = 300
            , width = 29
            , height = 26
            , units = "cm"
            , device = "tiff"
    )
  
  p %>% 
    ggsave( plot = .
            ,  filename = paste0( p_filename, ".png" )
            , dpi = 600
            , width = 29
            , height = 26
            , units = "cm"
    )
  
} )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POWER VS MAF IN STRATUM 3 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GENERATE ONE PLOT FOR EACH SELECTED NUMBER OF TRIADS PER PLOT:

# # If you want one plot per n_triads_per_stratum:
# lapply( unique(dat$n_triads_per_stratum), function(x){
# If just want plots where n_triads_per_stratum = 100
lapply( c(100), function(x){
  
  no_of_triads_per_stratum <- x
  
  
  # Labeller for facet grid:
  maf_labeller <- labeller(
    # Stratum 2 label:
    maf_stratum_2 = function(x)
      paste0("MAF in\nstratum 2 = ", x)
    # Stratum 1 label:
    , maf_stratum_1 = function(x)
      paste0("MAF in\nstratum 1 = ", x)
  )
  
  # Figure title
  fig_title <- paste0(
    "Power at different significance levels (\u03b1) across different combinations of within-stratum minor allele frequency (MAF)"
    , "\nin stratum 1, 2 and 3 ("  
    ,  no_of_triads_per_stratum
    , " triads in each stratum)"
  )

  # Figure caption
  fig_caption <- stringr::str_wrap(
    paste0(
      "Asymptotic power for single-SNP GxE analyses as calculated by 
    hapPowerAsymp() when using a diallelic SNP, the case-parent triad design, 
    the multiplicative dose-response model, and three strata."
      , " Sample size held constant at "
      , unique(dat$n_triads_per_stratum) * 3
      , " ("
      , unique(dat$n_triads_per_stratum)
      , " triads per stratum)."
      , " RR equals 3 in stratum 1, 1.0 in stratum 2, and 1/3 in stratum 3. 
    \"Minor allele frequency\" refers to the frequency of the least common 
    allele in all three strata." )
    , width = caption_width_fig
  )
  
  
  # CREATE PLOT:
  p <- dat %>% 
    filter( n_triads_per_stratum == no_of_triads_per_stratum ) %>% 
    filter( !(maf_stratum_1 %in% c(0.25, 0.35) )
            , !(maf_stratum_2 %in% c(0.25, 0.35) )
            , !(maf_stratum_3 %in% c(0.25, 0.35) )
    ) %>% 
    ggplot( .
            , aes( x = maf_stratum_3
                   , y = RR.power_Haplotype_2
                   , colour = alpha_f
                   , linetype = alpha_f
                   , group = alpha_f )
    ) +
    geom_point( 
      alpha = 0.4
      , stroke = 0 
      , size = geom_point_size
    ) +  
    geom_line(linewidth = geom_line_width
              , alpha = 0.8 ) +
    
    scale_x_continuous(
      name   = "MAF in stratum 3",
      limits = c(0, 0.4),
      breaks = seq(0, 0.4, by = 0.05),
    ) +
    scale_y_continuous(
      name   = "Power",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    
    # Colour scale:
    scale_colour_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = okabe_ito[seq_along(alpha_breaks)],
      name   = expression(alpha)
    ) +
    # Linetype scale:
    scale_linetype_manual(
      breaks = alpha_breaks,
      labels = alpha_labels,
      values = c("solid","22","42","46", "13", "1343" ),
      name   = expression(alpha)
    ) +
    
    # Title and caption:
    labs(
      title = fig_title
      , caption = fig_caption
    ) +
    
    # Theme:
    theme_bw( base_size = 11 ) +
    theme(
      axis.title.x   = element_text(margin = margin(t = 5, b = 7))
      , axis.title.y   = element_text(margin = margin(r = 10))
      
      , legend.box = "vertical"
      , legend.key.width = unit(44, "pt") # wider keys to show linetypes
      , legend.position = "right"
      , legend.title = element_text(face = "bold"
                                    , hjust = 0.5
                                    , size = 12 )
      , legend.margin = margin(0,0,0,0)
      , legend.text = element_text( margin = margin(t = 7, r = 0, b = 7, l = 2 )
                                    , size = 8 )
      
      , axis.text.x = element_text( angle = 90
                                    , vjust = 0.5
                                    , size = 8 )
      , axis.text.y = element_text( size = 8 )
      
      , plot.caption = 
        element_text( hjust = 0
                      , size = caption_size
                      , colour = "gray28"
                      , margin = margin( t = 0, r = 0, b = 0, l = -18, "pt" )
        )
      , plot.title = element_text( hjust = 0
                                   , size = 12
      )
      , strip.text = element_text(size = 9)
    ) +
    facet_grid(
      rows = vars(maf_stratum_2)
      , cols = vars( maf_stratum_1 )
      , labeller = maf_labeller
    )
  p
  
  # EXPORT
  
  p_filename <- paste0(
    "power_vs_maf_stratum_3_in_maf_stratum_2_x_maf_stratum_1_grid_"
    , no_of_triads_per_stratum
    , "_triads_per_stratum"
  )
  
  p %>% 
    ggsave( plot = .
            , filename = paste0( p_filename, ".tiff" )
            , dpi = 300
            , width = 29
            , height = 26
            , units = "cm"
            , device = "tiff"
    )
  
  p %>% 
    ggsave( plot = .
            ,  filename = paste0( p_filename, ".png" )
            , dpi = 600
            , width = 29
            , height = 26
            , units = "cm"
    )
  
} )

