# Description: ----------------------------------------------
#   Make box plots and scatter plots for descriptive purposes. 
#   The concept of the analysis is in ../Doc/mtg_190402.pdf.
#
# 19/ 04/ 03-


# Loading data ------------------------------------------------------------

dir.data           <- "../../Data"
fn.data            <- "AnalysisDataSet_v4.0.RData"

load(sprintf("%s/%s", dir.data, fn.data))

# Settings: -------------------------------------------------

args <- list(TRUE,FALSE)

# args <- commandArgs(trailingOnly = T)



if(length(args)==0){
  extract.primarySS <- FALSE
  extract.secondarySS <- FALSE
}else{
  extract.primarySS   <- args[[1]]
  extract.secondarySS <- args[[2]] 
}


if(is.null(extract.primarySS) & is.null(extract.secondarySS) )
{
  extract.primarySS   <- FALSE
  extract.secondarySS <- FALSE
}

# Setting for fn.suffix of outputs

if(extract.primarySS){
  suffix.df.coef.pooled.res <-
    "primarySS_extracted"
}else{
  if(extract.secondarySS){
    suffix.df.coef.pooled.res <-
      "secondarySS_extracted"
  }else suffix.df.coef.pooled.res <- "total_SS"
}


if(extract.primarySS){
  data <- data %>% filter(
    nijisei == "Primary"
    )
}



# Other settings: -------------------------------------------------

dir.sub               <- "../sub"
fn.require_packages.R <- "require_libraries.R"

dir.data           <- "../../Data"
fn.data            <- "AnalysisDataSet_v4.0.RData"

fn.df.coef.pooled.res <- sprintf(
  "%s/%s_%s.csv",
  dir.output,
  "df.coef.pooled.res.lmrob_MM_init.S_Ro60.mice.data_SS",
  suffix.df.coef.pooled.res
  )



dir.output         <- "../../Output/Final"

var.select <- c(
  'ACA_np', 'Dry_mouth', 'Raynaud_np', 
  "Dry_eye", 'IgG_pn','RF_pn'
  )


# Load subroutines: ----------------------------------------------

Bibtex <- FALSE

source(sprintf("%s/%s", dir.sub, fn.require_packages.R))


# Load coefficients from robust regression --------------------------------


df.coef.pooled.res.lmrob.Ro60.mice.data_SS <- try(
  read.csv(
    sprintf(
      "%s/%s",
      dir.output, 
      "df.coef.pooled.res.lmrob_MM_init.S_Ro60.mice.data_SS_nonNijisei.csv"
      )
    ) %>%
    filter(
      .id %in%
        var.select
      )
  )

df.coef.pooled.res.lmrob.Ro60.mice.data_SS[,'betas'] <- c(
  'b0','b1','b2','b3'
  )

df.spread.coef.pooled.res.lmrob.Ro60.mice.data_SS <- try(
  df.coef.pooled.res.lmrob.Ro60.mice.data_SS %>%
  dplyr::select(
    .id, betas, results.rnd
    ) %>%
  spread(
    betas, results.rnd
    )
  )


ggdata <- ggplot(
  data %>%
    filter(SS=="1")
  )

df.beta <-
  df.spread.coef.pooled.res.lmrob.Ro60.mice.data_SS

# Box Plots: var for y is Ro52 only--------------------------------------------------------------

# ExploratoryDataAnalysis::mf.wrap.boxplot()
#
# "ordering sheet": A dataframe with tidy style <- input.data_boxplot_1
#

input.data_boxplot_1 <- data.frame(
  
  "ind"       = c(1:18),
  "var.x"     = c(
    rep("disease_group", 9),
    rep("Ro60_pn", 9)
    ),
  
  "var.y"     = c(
    rep(
      "Ro52", 
      18
      )
    ),
  
  "scale.var.y"     = c(
    rep(
      "log10", 
      18
      )
    ),
  
  "var.caption"     = c(
    rep(
      "Dr.H.Nakamura, 2019", 
      18
    )
  ),
  
  "size"     = c(
    rep(
      0.2, 
      18
      )
    ),
  
  "var.col"   =
  c(
    rep(
      NA, #"anti_SS_A_ab",
      18
    )
  ),
  "plot.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        18
      )
    ),
  
  "box.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        18
      )
    ),
  
  "str"       = c(
    "Dummy","Ro60_pn","Lip_biopsy_tri","Lip_biopsy_tri","Dry_mouth","Dry_eye","ACA_np","Raynaud_np","HTLV_1"
    ),
  
  "dn.surfix" = sprintf(
    "%s_Fig_%s%s",
    suffix.df.coef.pooled.res,
    formatC(
      rep(1:2, each=9), width=3, flag="0"),  # # of disposition on x-ax = 2
    formatC(
      rep(1:9, each=1), width=3, flag="0")   # # of stratification = 9
    ),
  stringsAsFactors = FALSE
  )



input.data_boxplot_2 <- data.frame(
  
  "ind"       = c(1:18),
  "var.x"     = c(
    rep("disease_group", 9),
    rep("Ro60_pn", 9)
    ),
  
  "var.y"     = c(
    rep(
      "Ro52", 
      18
    )
  ),
  
  "scale.var.y"     = c(
    rep(
      "not_scale", 
      18
      )
    ),
  
  "var.caption"     = c(
    rep(
      "Dr.H.Nakamura, 2019", 
      18
    )
  ),
  
  "size"     = c(
    rep(
      0.2, 
      18
    )
  ),
  
  "var.col"   =
    c(
      rep(
        NA, #"anti_SS_A_ab",
        18
      )
    ),
  "plot.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        18
      )
    ),
  
  "box.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        18
      )
    ),
  
  "str"       = c(
    "Dummy","Ro60_pn","Lip_biopsy_tri","Lip_biopsy_tri","Dry_mouth","Dry_eye","ACA_np","Raynaud_np","HTLV_1"
  ),
  
  "dn.surfix" = sprintf(
    "%s_Fig_%s%s",
    suffix.df.coef.pooled.res, 
    formatC(
      rep(1:2, each=9), width=3, flag="0"),  # # of disposition on x-ax = 2
    formatC(
      rep(1:9, each=1), width=3, flag="0")   # # of stratification = 9
  ),
  stringsAsFactors = FALSE
)

list.plot.boxplot.all_subj <- dlply(
  input.data_boxplot_1 %>%
    mutate(
      dn.surfix=gsub("test", "all_subj",dn.surfix)
      ),
  .(ind),
  #ExploratoryDataAnalysis::
  mf.wrap.boxplot, 
  data, # = data %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data# %>%
      #filter(SS=="1")
    )
  )

list.plot.boxplot.SS_subj <- dlply(
  input.data_boxplot_1 %>%
    mutate(
      dn.surfix=gsub("test", "SS_subj",dn.surfix)
    ),
  .(ind),
  #ExploratoryDataAnalysis::
    mf.wrap.boxplot, 
  data = data %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data %>%
      filter(SS=="1")
    )
  )


list.plot.boxplot.all_subj.y_not_scaled <- dlply(
  input.data_boxplot_2 %>%
    mutate(
      dn.surfix=gsub("test", "all_subj.var.y_notscaled",dn.surfix)
    ),
  .(ind),
  #ExploratoryDataAnalysis::
  mf.wrap.boxplot, 
  data = data, 
  ggdata <- ggplot(
    data
    )
  )


list.plot.boxplot.SS_subj.y_not_scaled <- dlply(
  input.data_boxplot_2 %>%
    mutate(
      dn.surfix=gsub("test", "SS_subj.var.y_notscaled",dn.surfix)
    ),
  .(ind),
  #ExploratoryDataAnalysis::
  mf.wrap.boxplot, 
  data = data %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data %>%
      filter(SS=="1")
  )
)

# Box Plots: stratified by the subtypes anti-Ro52/ anti-60 -------------------------------------------------------------

# ExploratoryDataAnalysis::mf.wrap.boxplot()
#
# "ordering sheet": A dataframe with tidy style <- input.data_boxplot_1
#

data.gathered_by_RoSubtypes <-
  data %>%
  mutate(
    anti_Ro52_60= Ro52,
    Subtypes = "anti-Ro52"
    ) %>%
  bind_rows(
    data %>%
      mutate(
        anti_Ro52_60= anti_SS_A_ab,
        Subtypes = "anti-Ro60"
        )
  )

input.data_boxplot_2 <- data.frame(
  
  "ind"       = c(1:12),
  "var.x"     =
     c(
      "Dummy","Ro60_pn","Lip_biopsy_tri","Dry_mouth","Dry_eye","ACA_np","Raynaud_np","Saxon_test_np", "Schirmer_test_np", "anti_SS_B_ab_pn", "IgG_pn","RF_pn"
      # Schimar test, Saxon test
      ),
  
  "var.y"     = c(
    rep(
      "anti_Ro52_60", 
      12)
    ),
  "var.caption"     = c(
    rep(
      "Dr.H.Nakamura, 2019", 
      12
    )
  ),
  
  "size"     = c(
    rep(
      1.5, 
      12
    )
  ),
  
  "var.col"   =
    c(
      rep(
        NA, #"anti_SS_A_ab",
        12
      )
    ),
  "plot.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        12
      )
    ),
  
  "box.col"   =
    c(
      rep(
        "black", #"anti_SS_A_ab",
        12
      )
    ),
  
  "str"       = 
    c(
      rep(
        "Subtypes", #"anti_SS_A_ab",
        12
      )
    ),
  
  "dn.surfix" = sprintf(
    "Fig_%s%s", 
    formatC(
      rep(1:12, each=1), width=3, flag="0"),  # # of disposition on x-ax = 11
    formatC(
      rep(1, each=1), width=3, flag="0")   # # of stratification = 2
  ),
  stringsAsFactors = FALSE
)


list.plot.str_by_Subtypes_boxplot.SS_subj <- dlply(
  input.data_boxplot_2 %>%
    mutate(
      dn.surfix=gsub("test", "SS_subj",dn.surfix)
    ),
  .(ind),
  #ExploratoryDataAnalysis::
  mf.wrap.boxplot, 
  data.gathered_by_RoSubtypes %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data.gathered_by_RoSubtypes %>%
    filter(SS=="1")
  )
)



# Scatterplot: -------------------------------------------------------------------



input.data_scatter_4 <- data.frame(
  "ind"       = c(1:84),
  "var.x"     = c(
    rep("Age", 14), 
    rep("FS", 14), 
    rep("anti_SS_A_ab", 14), 
    rep("anti_SS_B_ab", 14),
    rep("ESSDAI", 14),
    rep("IgG", 14)
  ),
  
  "var.y"     = c(
    rep(
      "Ro52", 
      84
    )
  ),
  
  "trans.y" = rep("log10",84),
#  "trans.x" = rep("NoScale" ,84),
  "trans.x" = rep("log10" ,84),
  "size"    = rep(2.0 ,84), 
  "var.col"   = c(
    rep(
      NA,
      84
    )
  ),
  "plot.col"   = c(
    rep(
      "black",
      84
    )
  ),
  "line.col"   = c(
    rep(
      "black",
      84
    )
  ),
  "cont.col"   = c(
    rep(
      "black",
      84
    )
  ),
  "str"       = c(
    "Dummy","Ro60_pn","Ro52_ACA","Ro52_Ro60_ACA","Sex","Lip_biopsy",
    "Dry_mouth","Dry_eye","ACA_np","Raynaud_np","HTLV_1",
    "ESSDAI", "IgG_pn", "RF_pn"
  ),
"dn.surfix" = sprintf(
  "%s_Fig_%s%s",
  suffix.df.coef.pooled.res,formatC(
      rep(3:8, each=14), width=3, flag="0"),  # # of disposition on x-ax = 3
    formatC(
      rep(1:14, each=1), width=3, flag="0")   # # of stratification = 14
  ),
  stringsAsFactors = FALSE
)


dlply(
  input.data_scatter_4 %>% filter(
    str %in%
      var.select &
      var.x == 'anti_SS_A_ab'
  # ) %>%
  #   mutate(
  #     dn.surfix=gsub("test", "SS_subj",dn.surfix)
    ),
  .(ind),
  mf.wrap.scatterplot, 
  data = data %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data %>%
      filter(SS=="1")
    ),
  df.beta = df.beta
  )



# Scatter plot with miss_box ----------------------------------------------


input.data_scatter_5 <- data.frame(
  "ind"       = c(1:6),
  "var.x"     = c(
    rep("Age", 1), 
    rep("FS", 1), 
    rep("anti_SS_A_ab", 1), 
    rep("anti_SS_B_ab", 1),
    rep("ESSDAI", 1),
    rep("IgG", 1)
  ),
  
  "var.y"     = c(
    rep(
      "Ro52", 
      6
    )
  ),
  
  "trans.y" = rep("log10",6),
  #  "trans.x" = rep("NoScale" ,84),
  "trans.x" = rep("NoScale" ,6),
  "size"    = rep(1.0 ,6), 
  "var.col"   = c(
    rep(
      NA,
      6
    )
  ),
  "plot.col"   = c(
    rep(
      "black",
      6
    )
  ),
  "line.col"   = c(
    rep(
      "black",
      6
    )
  ),
  "cont.col"   = c(
    rep(
      "black",
      6
    )
  ),
  "str"       = c(
    rep("Dummy", 6)
    ),
  "dn.surfix" = sprintf(
    "test_Fig_x.x_with_miss.box_%s%s", 
    formatC(
      rep(3:8, each=1), width=3, flag="0"),  # # of disposition on x-ax = 3
    formatC(
      rep(1, each=1), width=3, flag="0")   # # of stratification = 14
  ),
  stringsAsFactors = FALSE
)


dlply(
  input.data_scatter_5 %>%
       mutate(
         dn.surfix=gsub("test", "SS_subj",dn.surfix)
         ),
  .(ind),
  mf.wrap.scatterplot.with_missbox, 
  data = data %>% filter(SS=="1"), 
  ggdata <- ggplot(
    data %>%
      filter(SS=="1")
  ),
  df.beta = NULL
)


# File_name_coding --------------------------------------------------------

write.csv(
  input.data_boxplot_1,
  sprintf(
    "%s/%s.csv",
    dir.output,
    "_File_Name_Coding_Boxplot"
    )
  )


write.csv(
  input.data_boxplot_2,
  sprintf(
    "%s/%s.csv",
    dir.output,
    "_File_Name_Coding_Boxplot.Panel_Subtypes"
  )
)

write.csv(
  input.data_scatter_4,
  sprintf(
    "%s/%s.csv",
    dir.output,
    "_File_Name_Coding_Scatterplot"
    )
  )


# For comparison between disease group ------------------------------------

data_tidy.all_pat <- data %>%
  gather(var, val, -SS, -RA, -HC, -SubjID, -varname) %>%
  mutate(
    group = ifelse(
      SS ==1, "SS",
      ifelse(
        RA ==1, "RA",
        ifelse(
          HC ==1, "HC", NA
          )
        )
      )
    ) %>%
  mutate(
    group=factor(group, levels = c("SS", "RA", "HC"))
    )


