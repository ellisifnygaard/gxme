
library(dplyr)
library(data.table)
library(Haplin)
library(magrittr)

wd <- getwd()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     GENERATE DATA FOR FIGURE 3 (SUPPLEMENTARY S3) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# VARYING OVER RR, SAMPLE SIZE AND MAF (BUT SAME MAF ACROSS THE 3 STRATA)

# VARY RELATIVE RISK (X AXIS VARIABLE)

# VARY ALLELE FREQUENCY IN ALL THREE STRATA FROM c(0.6, 0.4) TO c(0.95, 0.05) 
# (FACET)

# Create dt with parameters:
fig_3_params <- base::expand.grid(
  # Variables identical for both figures:
  n_all = 2 # diallelic snp
  , n_strata = 3 # 3 strata
  , RR_star_lower = 1 
  , RR_star_upper = 1
  , response = "mult"
  # Facet variable, allele frequency in each stratum:
  , haplo_freq_most_common_allele_stratum_1 = c(0.7, 0.8, 0.85, 0.9, 0.95)
  # Facet variable, number of triads per strata:
  , n_triads_per_stratum = floor( seq(50, 800, by = 50) / 3 )
  # Linetype variable, alpha:
  , alpha = c( 0.05 # "baseline"
               , 0.05*10^(-3) # (1,000 pairings)
               , 0.05*10^(-4) # (10,000 pairings)
               , 0.05*10^(-5) # (100,000 pairings)
               , 0.05*10^(-6) # (1,000,000 pairings)
               , 0.05*10^(-8) # (100,000,000 pairings)
  )
  # X axis variable, "rr":
  , rr = seq( 1, 7, by = 0.3)
) %>% as.data.table()

# Use rr to derive the RR in the three different stratum using the formula Ø
# suggested, and add haplo.freq columns for stratum 2 and 3 (did not add them in
# expand.grid to avoid making one row per unique combination of allele frequency
# in stratum):
fig_3_params %<>% 
  mutate( RR_stratum_1 = rr
          , RR_stratum_2 = 1
          , RR_stratum_3 = 1/rr
          , haplo_freq_most_common_allele_stratum_2 =
            haplo_freq_most_common_allele_stratum_1
          , haplo_freq_most_common_allele_stratum_3 =
            haplo_freq_most_common_allele_stratum_1
  )


power_res_fig_3 <- lapply( 1:nrow(fig_3_params), function(i){
  
  params_i <- fig_3_params[i,]
  
  power_gxe_res <- hapPowerAsymp(
    nall = params_i$n_all 
    , n.strata = params_i$n_strata
    , cases = list( c( mfc = params_i$n_triads_per_stratum )
                    , c( mfc = params_i$n_triads_per_stratum )
                    , c( mfc = params_i$n_triads_per_stratum ) )
    , haplo.freq = list(
      c( params_i$haplo_freq_most_common_allele_stratum_1
         , 1 - params_i$haplo_freq_most_common_allele_stratum_1 )
      , c( params_i$haplo_freq_most_common_allele_stratum_2
           , 1 - params_i$haplo_freq_most_common_allele_stratum_2 )
      , c( params_i$haplo_freq_most_common_allele_stratum_3
           , 1 - params_i$haplo_freq_most_common_allele_stratum_3 ) )
    , RR = list( c(1, params_i$RR_stratum_1)
                 , c(1, params_i$RR_stratum_2)
                 , c(1, params_i$RR_stratum_3) ) 
    , RRstar = c(params_i$RR_star_lower, params_i$RR_star_upper)
    , response = params_i$response
    , alpha = params_i$alpha
  )
  
  # Add RR.power for Haplotype 2 to basis_i:
  params_i$RR.power_Haplotype_2 <- power_gxe_res$haplo.power %>%
    filter(Haplotype == 2) %>%
    pull(RR.power) %>%
    as.numeric()
  
  return( params_i )
  
})

length(power_res_fig_3)
power_res_fig_3 <- dplyr::bind_rows(power_res_fig_3)
power_res_fig_3
gc()

# Export data table:
power_res_fig_3 %>% 
  readr::write_csv( 
    . 
    , file = file.path(
      wd
      , "hapPowerAsymp_results_fig3_varying_alpha_sample_size_and_rr_maf_from_005_to_04.csv"
      
    )
  )



# Any combinations of parameters (i.e., rows) resulting in power = NaN/NA?
power_res_fig_3 %>% count( is.na(RR.power_Haplotype_2) )






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     SET COMMON GGPLOT ARGUMENTS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Libraries:
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

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
caption_size <- 8.5
caption_width_fig3 <- 180

# Build factor and expression labels for α:
alpha_vals   <- sort(unique(power_res_fig_3$alpha), decreasing = TRUE) # numeric
alpha_breaks <- scientific(alpha_vals, digits = 1) # “5e-02”, …

# Convert to desired math-style expressions:
alpha_lab_chr <- vapply(alpha_vals, function(a) {
  # If alpha = 0.05:
  if (abs(a - 0.05) < 1e-12) { 
    "0.05"
    # If alpha 0.05·10^k, k = -2, -4, -5, -6, -8: 
  } else { 
    k <- as.integer( round( log10(a / 0.05) ) )
    paste0("0.05%*%10^", k) # e.g. "0.05*10^-3"
  }
}, FUN.VALUE = character(1) )

# Turn into plotmath expressions:
alpha_labels <- parse(text = alpha_lab_chr)   

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     PLOTTING FIGURE 3 (S3) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Figure 3 caption
fig3_caption <- stringr::str_wrap(
  paste0(
    "Asymptotic power for single-SNP GxE analyses as calculated by
  hapPowerAsymp() when using a diallelic SNP, the
  case-parent triad design, the multiplicative dose-response model, and three
  strata."
    , " Divide the total sample size by three to get the number 
                       of case-parent triads per stratum."
    , " rr is a positive number used to determine the 
    strata-specific RR values: RR equals rr in stratum 1, 1.0 in stratum 2, and 
    1/rr in stratum 3. \"Minor allele frequency\" refers to the frequency of 
    the least common allele in all three strata." )
  , width = caption_width_fig3
)




# Data frame with factors for α and allele freq: 
dat <- power_res_fig_3  %>% 
  mutate(
    alpha_f = factor(alpha, levels = alpha_vals, labels = alpha_breaks)
    
    # Reverse order of faceting variable:
    , haplo_freq_most_common_allele_stratum_1_f = factor(
      haplo_freq_most_common_allele_stratum_1,
      levels = sort(unique(haplo_freq_most_common_allele_stratum_1),
                    decreasing = TRUE)
    )
    # Calculate total sample size (all three strata)
    , sample_size = n_triads_per_stratum*n_strata
  ) %>% 
  # Remove very large effect sizes and tiny sample sizes:
  filter( rr <= 5
          , sample_size%in% c(150, 300, 450, 600, 750) )


# Facet label function:
maf_labeller <- labeller(
  haplo_freq_most_common_allele_stratum_1_f = function(x) {
    f <- as.numeric(as.character(x)) # get the true value
    paste0("Minor allele frequency =\n", formatC(1 - f, digits = 2))
  }
  , sample_size = function(x){
    f <- as.numeric(as.character(x)) # get the true value
    paste0("N triads = ", f)
  }
)


# Plot:
power_vs_risk_and_sample_size <- ggplot(dat
                        , aes(x = rr
                              , y = RR.power_Haplotype_2
                              , colour   = alpha_f
                              , linetype = alpha_f
                              , group    = alpha_f ) ) +
  geom_line(linewidth = geom_line_width
            , alpha = 0.8 ) +
  geom_point(size = geom_point_size
             , shape = 20
             , alpha = 0.4 ) +
  
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
  
  # Axes: 
  scale_x_continuous(
    breaks       = seq(1, max(dat$rr), by = 0.5),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2)
  ) +
  
  # Labels:
  labs(
    x = "rr"
    , y = "Power"
    , title = expression( "Power vs. relative risk (RR) at different significance levels ("
                          * alpha * ") across different sample sizes and minor allele frequencies" )
    , caption = fig3_caption
  ) +
  # Theme: 
  theme_bw(base_size = 11) +
  theme(
    axis.title.x   = element_text(margin = margin(t = 5, b = 7))
    , axis.title.y   = element_text(margin = margin(r = 10))
    , legend.key.width = unit(44, "pt") # wider keys to show linetypes
    , legend.position  = "right"
    , legend.title     = element_text(face = "bold"
                                      , hjust = 0.5
                                      , size = 12 )
    , legend.margin = margin(0,0,0,0)
    , legend.text = element_text( margin = margin(t = 7, r = 0, b = 7, l = 2 )
                                  , size = 8 )
    , plot.subtitle    = element_text(size = 10)
    , plot.title = element_text(size = 12)
    , plot.caption = 
      element_text( hjust = 0
                    , size = caption_size
                    , colour = "gray28"
                    , margin = margin( t = 0, r = 0, b = 0, l = -18, "pt" )
      )
    , strip.text = element_text(size = 9)
  ) +
  facet_grid( cols = vars(haplo_freq_most_common_allele_stratum_1_f)
              , rows = vars(sample_size)
              , labeller = maf_labeller
              )



# EXPORT FIGURE (SUPPLEMENTARY FIGURE S3)

power_vs_risk_and_sample_size %>% 
  ggsave( plot = .
          , filename =
            "power_vs_relative_risk_in_sample_size_x_allele_frequency_grid.tif"
          , dpi = 300
          # , width = 21
          # , height = 29
          , width = 29
          , height = 21
          , units = "cm"
  )

power_vs_risk_and_sample_size %>% 
  ggsave( plot = .
          , filename =
            "power_vs_relative_risk_in_sample_size_x_allele_frequency_grid.png"
          , dpi = 600
          # , width = 21
          # , height = 29
          , width = 29
          , height = 21
          , units = "cm"
  )






