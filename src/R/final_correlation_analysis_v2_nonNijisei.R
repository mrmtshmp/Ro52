# Description: ----------------------------------------------
#   . 
#   The concept of the analysis is in ../Doc/mtg_190402.pdf.
#
# 19/ 04/ 21-


# Settings: -------------------------------------------------

set.seed(22)

dir.sub               <- "../sub"
fn.require_packages.R <- "require_libraries.R"

dir.data           <- "../../Data"
fn.data            <- "AnalysisDataSet_v4.0.RData"

load(sprintf("%s/%s", dir.data, fn.data))

data$Ro52_log10 <- log(
  data$Ro52,
  base = 10
  )

data$Ro60_log10 <- log(
  data$anti_SS_A_ab,
  base = 10
  ) 

dir.output         <- "../../Output/Final"




# Load subroutines: ----------------------------------------------

Bibtex <- FALSE

source(
  sprintf(
    "%s/%s", dir.sub, fn.require_packages.R
    )
  )

# Dichotomyze outcomes for Lip_biopsy and ESSDAI  -----------------------------------------------------

data_Dichotomyze <- data %>%
  mutate(
    Ro52_ext.High =   # lately (within this pipe) converted to factor
      ifelse(Ro52 > 500, 1, 0),
    Lip_biopsy_Dimyze_by_0 = 
      ifelse(
        is.na(Lip_biopsy_tri), NA,
             ifelse(
               Lip_biopsy_tri ==
                 0, 0, 1
               )
             ),
    Lip_biopsy_Dimyze_by_3 =
      ifelse(
        is.na(Lip_biopsy_tri), NA,
        ifelse(
          Lip_biopsy_tri %in%
            c(0, 3), 0, 1
          )
        ),
    ESSDAI_Dimyze_by_2 =
      ifelse(
        is.na(ESSDAI), NA,
        ifelse(
          ESSDAI <= 2, 0, 1
          )
        ),
    ESSDAI_Dimyze_by_5 =
      ifelse(
        is.na(ESSDAI), NA,
        ifelse(
          ESSDAI <= 5, 0, 1
          )
        ),
    FS_Dimyze_by_0 =
      ifelse(
        is.na(FS), NA,
        ifelse(
          FS < 1, 0, 1
          )
        ),
    FS_Dimyze_by_2 =
      ifelse(
        is.na(FS), NA,
        ifelse(
          FS < 3, 0, 1
          )
        ),
    FS_Dimyze_by_4 =
      ifelse(
        is.na(FS), NA,
        ifelse(
          FS < 5, 0, 1
          )
        )
    ) %>%
  mutate(Ro52_ext.High=factor(Ro52_ext.High, levels = c(0,1), labels = c("Normal", "High"))) %>%
  filter(SS==1)
  
AECG_component.num <- c(
  "FS",
  "anti_SS_B_ab"
  )

AECG_component <- c(
  "Dry_mouth",       #component
  "Dry_eye",         #component
  "FS",              #component
  "Lip_biopsy_tri",  #component
  "anti_SS_B_ab",    #component
  "anti_SS_B_ab_pn", #component
  "ACA_np",
  "Raynaud_np",
  "RF_pn",
  "IgG_pn",
  "Saxon_test_np",
  "Schirmer_test_np",
  "nijisei"
  )

# Data shaping ----------------------------------------------------------

data_tidy <- data_Dichotomyze %>%
  
  filter(
    disease_group=="SS"
    ) %>%
  
  mutate(
    Ro60 = anti_SS_A_ab
    ) %>%
  
  mutate_if(is.factor, as.numeric) %>%

  gather(
    var, val, 
    -SubjID, -Ro52, -Ro52_log10,-Ro60, -Age,
    -disease_group, -nijisei
    ) %>%
  
  mutate(
    val = as.numeric(val),
    Ro60_log10 = log(Ro60,base = 10)
    ) %>%
  
  filter(
    var %in% 
      c(AECG_component,"Ro52_ext.High")
    ) %>%
  
  gather(
    var_y, val_y, 
    -SubjID, -Age, -var, -val,
    -disease_group,- nijisei
    )



# ANOVA -------------------------------------------------------------------

df.ANOVA <- data_Dichotomyze %>%
  dplyr::select(
    SubjID,
    Ro52_log10,
    Ro60_log10,
    AECG_component
  )

# make list of covariates

list.formula <- dlply(
  data.frame(
    "id" = 1:length(AECG_component[c(1,2,4,6,7,8,9,10,11,12,13)]),
    "var"= AECG_component[c(1,2,4,6,7,8,9,10,11,12,13)] 
    ),
  .(id),
  function(D){
    fmr = sprintf(
      "%s~%s", 
      "Ro52_log10",
      D$var
    )
    return(fmr)
  }
)

res.lmrob.Ro52_log10 <- 
  llply(
    list.formula,
    function(L){
      robustbase::lmrob(
        as.formula(L),
        df.ANOVA,
        method = 'MM',
        #        setting="KS2014",
        control = lmrob.control(maxit.scale = 2000)
      )
    }
  )

# make list of covariates

list.formula_Ro60 <- dlply(
  data.frame(
    "id" = 1:length(AECG_component[c(1,2,4,6,7,8,9,10,11,12)]),
    "var"= AECG_component[c(1,2,4,6,7,8,9,10,11,12)] 
  ),
  .(id),
  function(D){
    fmr = sprintf(
      "%s~%s", 
      "Ro60_log10",
      D$var
    )
    return(fmr)
  }
)

res.lmrob.Ro60_log10 <- 
  llply(
    list.formula_Ro60,
    function(L){
      robustbase::lmrob(
        as.formula(L),
        df.ANOVA,
        method = 'MM',
        #        setting="KS2014",
        control = lmrob.control(maxit.scale = 2000)
      )
    }
  )

# Extract estimated coefficients

coef.lmrob.Ro52_log10 <- llply(
  res.lmrob.Ro52_log10,
  summary
) %>%
  ldply(
    function(L){
      out <- coef(L) %>%
        as.data.frame() %>%
        rownames_to_column("terms")
      return(out)
    }
  )

coef.lmrob.Ro60_log10 <- llply(
  res.lmrob.Ro60_log10,
  summary
) %>%
  ldply(
    function(L){
      out <- coef(L) %>%
        as.data.frame() %>%
        rownames_to_column("terms")
      return(out)
    }
  )


# Boxplot -----------------------------------------------------------------

gg.data_tidy <- data_tidy %>%
  filter(
    var %in% c(
      AECG_component[c(1,2,6,7,8,9,10,11,12)],
      "Ro52_ext.High"
      )
    ) %>%
  mutate(
    thre = ifelse(
      var_y=="Ro52", 10, 
      ifelse(
        var_y=="Ro52_log10", 1, NA
        )
      ),
    Alpha = ifelse(
      var_y %in% c("Ro52", "Ro52_log10"), 0.8, 0
      )
    ) %>%

  ggplot(
    aes(
      x = as.factor(val),
      y = val_y,
      yintercept = thre
      )
    )

plot.boxplot <- plot(
  gg.data_tidy + geom_boxplot(outlier.alpha = 0) + 
    geom_beeswarm(col="black", size=1, alpha=1) + 
#    geom_jitter(col="black", size=0.2, alpha=1, height = 0.01, width = 0.1) +
    geom_hline(aes(yintercept = thre, alpha = Alpha), size=0.5, col="black") +
    facet_grid(
      var + var_y ~ disease_group,# + nijisei,
      scales = "free") +
    theme_bw()
  )


# Scatterplot -----------------------------------------------------------------

gg.data_tidy.scatter <- data_tidy %>%
  filter(
    (var %in% c("FS", "anti_SS_B_ab"))
    ) %>%
  # filter(!is.na(val)) %>%
  ggplot(
    aes(
      x = as.numeric(val),
      y = val_y
    )
  )

plot.scatterplot <- plot(
  gg.data_tidy.scatter + #geom_boxplot(outlier.alpha = 0) + 
    geom_point(alpha=0.6) +
    scale_x_continuous(trans = "log10") +
    facet_grid( var + var_y ~., scales = "free") +
    theme_bw()
  )



gg.data_tidy.scatter <- data_tidy %>%
  filter(
    (var %in% c("FS", "anti_SS_B_ab"))
  ) %>%
  # filter(!is.na(val)) %>%
  ggplot(
    aes(
      x = as.numeric(val),
      y = val_y
    )
  )

plot.scatterplot.Subgroup_pri_sec <- plot(
  gg.data_tidy.scatter + #geom_boxplot(outlier.alpha = 0) + 
    geom_point(alpha=0.6) +
    scale_x_continuous(trans = "log10") +
    facet_wrap( var + var_y ~ disease_group + nijisei, scales = "free") +
    theme_bw()
  )

plot.scatterplot <- plot(
  gg.data_tidy.scatter + #geom_boxplot(outlier.alpha = 0) + 
    geom_point(alpha=0.6) +
    scale_x_continuous(trans = "log10") +
    facet_wrap( ~ var + var_y, scales = "free", ncol=1) +
    theme_bw()
)

# Missingness

gg.data_tidy.scatter_missing <- data_tidy %>%
  filter(
    (var %in% c("FS", "anti_SS_B_ab"))
  ) %>%
  mutate(
    flg.na = factor(
      ifelse(is.na(val), "missing", "observed")
      )
    ) %>%
  ggplot(
    aes(
      x = as.numeric(flg.na),
      y = val_y,
      group = flg.na
      )
    )

plot.boxplot.scatter_missing.Subgroup_pri_sec <- plot(
  gg.data_tidy.scatter_missing + 
    geom_boxplot(
      color ="black",
      outlier.alpha = 0
      ) + 
    geom_jitter(
      color  ="black",
      width = 0.2, alpha=0.8, size=0.75
      ) +
    facet_grid(
      var + var_y ~ disease_group + nijisei, scales = "free"
      ) +
    theme_bw() +
    scale_x_discrete()
  )


plot.boxplot.scatter_missing <- plot(
  gg.data_tidy.scatter_missing + 
    geom_boxplot(color ="black", outlier.alpha = 0) + 
    geom_jitter(color  ="black", width = 0.2, alpha=0.8, size=0.75) +
    facet_wrap( ~ var + var_y, ncol = 1, scales = "free") +
    theme_bw() +
    scale_x_discrete()
  )


# Mutual Information ------------------------------------------------------

MIPermute_Ro52_self <- ExploratoryDataAnalysis::MIPermute(
  #mutinformation(
  X=data_Dichotomyze$Ro52_log10, 
  Y=data_Dichotomyze$Ro52_log10,
  method = "shrink",
  n.sim = 500
  )[1,2]

MIPermute_Ro60_self <- ExploratoryDataAnalysis::MIPermute(
  #mutinformation(
  X=data_Dichotomyze$Ro60_log10, 
  Y=data_Dichotomyze$Ro60_log10,
  method = "shrink",
  n.sim = 500
  )[1,2]


AECG_component

pdf(
  sprintf(
    "%s/%s.pdf",
    dir.output,
    "hist.mutinfo_Ro52_SS.Primary"
    )
  )
for(i in 1:length(AECG_component.num)){
  
  res.MIPermute <- ExploratoryDataAnalysis::MIPermute(
    #mutinformation(
    X=data_Dichotomyze$Ro52_log10[data_Dichotomyze$nijisei=="Primary"], 
    Y=data_Dichotomyze[data_Dichotomyze$nijisei=="Primary",AECG_component.num[i]],
    n.sim = 10000,
    method="MIC",
    disc.X = "none",
    disc.Y = "none",
    use = 'pairwise.complete.obs'
    )
  
  q.95 <-  quantile(res.MIPermute$V1, 0.95)

    assign(
      sprintf(
        "MIPermute_Ro52_%s", AECG_component.num[i]
        ),
      res.MIPermute
      )
  hist(
    res.MIPermute$V1,
    breaks='FD',
    main = AECG_component.num[i]
    )
  abline(
    v=res.MIPermute[res.MIPermute$i==1, 'V1'], 
    col='red'
    )
  abline(
    v=q.95, 
    col='red',
    lty=2
    )
  print(AECG_component.num[i])
  }
dev.off()

pdf(
  sprintf(
    "%s/%s.pdf",
    dir.output,
    "hist.mutinfo_Ro60_SS.Primary"
  )
)
for(i in 1:length(AECG_component.num)){
  
  res.MIPermute <- ExploratoryDataAnalysis::MIPermute(
    #mutinformation(
    X=data_Dichotomyze$Ro60_log10, 
    Y=data_Dichotomyze[,AECG_component.num[i]],
    n.sim = 10000,
    method="MIC",
    use = 'pairwise.complete.obs'
    )
  
  q.95 <-  quantile(
    res.MIPermute$V1, 0.95
    )
  
  assign(
    sprintf(
      "MIPermute_Ro60_%s", AECG_component.num[i]
      ),
    res.MIPermute
    )
  hist(
    res.MIPermute$V1,
    breaks='FD',
    main = AECG_component.num[i]
    )
  
  abline(
    v=res.MIPermute[res.MIPermute$i==1, 'V1'], 
    col='red'
    )
  abline(
    v=q.95, 
    col='red',
    lty=2
    )
    print(AECG_component.num[i])
}
dev.off()


# Tabulate the results from permutation test of MIC analysis --------------

MIPermute_Ro52_anti_SS_B_ab$pval <-
  1 - rank(MIPermute_Ro52_anti_SS_B_ab$V1)/
  nrow(MIPermute_Ro52_anti_SS_B_ab) 

MIPermute_Ro52_anti_SS_B_ab$dataname <-
  "MIPermute_Ro52_anti_SS_B_ab"

MIPermute_Ro60_anti_SS_B_ab$pval <-
  1- rank(MIPermute_Ro60_anti_SS_B_ab$V1)/
  nrow(MIPermute_Ro60_anti_SS_B_ab) 

MIPermute_Ro60_anti_SS_B_ab$dataname <-
  "MIPermute_Ro60_anti_SS_B_ab"

  
MIPermute_Ro52_FS$pval <-
  1 - rank(MIPermute_Ro52_FS$V1)/
  nrow(MIPermute_Ro52_FS) 

MIPermute_Ro52_FS$dataname <-
  "MIPermute_Ro52_FS"

MIPermute_Ro60_FS$pval <-
  1 - rank(MIPermute_Ro60_FS$V1)/
  nrow(MIPermute_Ro60_FS) 

MIPermute_Ro60_FS$dataname <-
  "MIPermute_Ro60_FS"

MIPermute <- MIPermute_Ro52_anti_SS_B_ab %>%
  rbind(MIPermute_Ro52_FS) %>%
  rbind(MIPermute_Ro60_anti_SS_B_ab) %>%
  rbind(MIPermute_Ro60_FS) %>%
  filter(i==1) %>%
  dplyr::select(pval, dataname)

# Output ------------------------------------------------------------------


g1 <- ggplotGrob(plot.scatterplot)
g2 <- ggplotGrob(plot.boxplot.scatter_missing)

pdf(
  sprintf(
    "%s/%s.pdf",
    dir.output,
    "scatterplot_with_miss.box"
  ),
  height = 56,
  width = 10
) 
plot_grid(
  g1,g2,
  align = "h",axis = "l", ncol = 2, rel_widths = c(5/7, 2/7)#, 1/8)
)
dev.off()

pdf(
  file=sprintf(
    "%s/%s.pdf",
    dir.output,
    "boxplot"
    ),
#  type = "pdf",
#  device = dev.cur(),
#  dpi = 300,
  width = 10,
  height= 180
  )
plot.boxplot
dev.off()

pdf(
  file=sprintf(
    "%s/%s.pdf",
    dir.output,
    "scatterplot"
  ),
  #  type = "pdf",
  #  device = dev.cur(),
  #  dpi = 300,
  width = 70/8,
  height= 70
)
plot.scatterplot
dev.off()



pdf(
  file=sprintf(
    "%s/%s.pdf",
    dir.output,
    "scatterplot.Subgroup_pri_sec"
  ),
  #  type = "pdf",
  #  device = dev.cur(),
  #  dpi = 300,
  width = 70/4,
  height= 70
)
plot.scatterplot.Subgroup_pri_sec
dev.off()


pdf(
  file=sprintf(
    "%s/%s.pdf",
    dir.output,
    "boxplot.scatter_missing"
  ),
  #  type = "pdf",
  #  device = dev.cur(),
  #  dpi = 300,
  width = 7,
  height= 70
)
plot.boxplot.scatter_missing
dev.off()


pdf(
  file=sprintf(
    "%s/%s.pdf",
    dir.output,
    "boxplot.scatter_missing.Subgroup_pri_sec"
  ),
  #  type = "pdf",
  #  device = dev.cur(),
  #  dpi = 300,
  width = 5,
  height= 70
)
plot.boxplot.scatter_missing.Subgroup_pri_sec
dev.off()

write.csv(
  file = sprintf(
    "%s/%s.csv",
    dir.output, "p_value.MIPermute_SS.Primary"
    ),
  MIPermute
  )

write.csv(
  file = sprintf(
    "%s/%s.csv",
    dir.output, "coef_lmrob_Ro52.AECGcomponent.total_SS"
    ),
  coef.lmrob.Ro52_log10 %>%
    mutate_if(
      is.numeric, function(x)round(x, 3)
    )
  )
write.csv(
  file = sprintf(
    "%s/%s.csv",
    dir.output, "coef_lmrob_Ro60.AECGcomponent.total_SS"
  ),
  coef.lmrob.Ro60_log10 %>%
    mutate_if(
      is.numeric, function(x)round(x, 3)
    )
  )


# Endrant -----------------------------------------------------------------


