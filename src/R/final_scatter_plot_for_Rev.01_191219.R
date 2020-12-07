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
fn.data            <- "AnalysisDataSet_v5.0.RData"

load(sprintf("%s/%s", dir.data, fn.data))

data <- 
  data %>%
  mutate(
    disease_group_2 = ifelse(
      disease_group == 'SS',
      sprintf('%s SS', nijisei),
      as.character(disease_group)
    )
  ) %>%
  mutate(
    disease_group_2 =
      factor(
        disease_group_2,
        c('Primary SS', 'Secondary SS', 'non_SS', 'RA', 'HC'),
        c('Primary SS', 'Secondary SS', 'non SS', 'RA', 'HC'))
    )

data$Ro52_log10 <- log(
  data$Ro52,
  base = 10
  )

data$Ro60_log10 <- log(
  data$anti_SS_A_ab,
  base = 10
  ) 

dir.output         <- "../../Output/Final_Rev01"

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
      ifelse(Ro52 < 500, 0, 1),
    Lip_biopsy_Dimyze_by_0 = 
      ifelse(
        is.na(Lip_biopsy_tri), NA,
             ifelse(
               Lip_biopsy_tri ==
                 0, 0, 1
               )
             )
    ) %>%
  mutate(Ro52_ext.High=factor(Ro52_ext.High, levels = c(0,1), labels = c("Normal", "High")))
  

ESSDAI_component <- c(
  "Articular", "Biological", "Cutaneous", "Glandular", "Hematological",   
  "Muscular", "Pulmonary", "Renal",
  "nijisei"
  )

# Data shaping ----------------------------------------------------------

data_tidy <- data_Dichotomyze %>%
  
  
  mutate(
    Ro60 = anti_SS_A_ab,
    ACA_np = as.numeric(ACA_np),
    Ro60_pn_num = as.numeric(Ro60_pn_num), 
    Ro52_ext.High= as.numeric(Ro52_ext.High)
    ) %>%
  

  gather(
    var, val, 
    -SubjID, -Ro52, -Ro52_log10,-Ro60, -Age,
    -disease_group_2, -nijisei
    ) %>%
  
  mutate(
    val = as.numeric(val),
    Ro60_log10 = log(Ro60,base = 10)
    ) %>%
  
  filter(
    var %in% 
      c(ESSDAI_component, "ACA_np", "Ro60_pn_num", 'Ro52_ext.High')
    ) %>%
  
  gather(
    var_y, val_y, 
    -SubjID, -Age, -var, -val,
    -disease_group_2,- nijisei
    ) %>%
  mutate(
    val = ifelse(var%in% ESSDAI_component & is.na(val), 0, 1)
    ) # For ESSDAI, 


# Boxplots for Fig.2 and Fig.Sx -------------------------------------------

#' Fig.2 is for 
#' Fig.Sx for display distribution of levels in antiRo52



# Boxplot for Fig.2 -------------------------------------------------------

#' Fig.2.1: plot for all subjects, and var.x is 'disease_group_2'
#' Fig.2.2: plot for SS subjects, and var.x are 'ACA_np', 'high and normal in levels of antiRo60' and 'ExtremeHigh in antiRo52'.


# ver.1: Primary and Secondary SS are separated.

data_for_fig2.ver.1 <- data_tidy %>% 
  mutate(
    dummy=1
  ) %>%
  spread(var, val) %>%
  spread(var_y, val_y) %>%
  mutate(
    ACA_np = factor(ACA_np, c(1,2), c('(-)', '(+)')),
    Ro52_ext.High = factor(Ro52_ext.High , c(1,2), c('<=500', '>500')),
    Ro60_pn_num = factor(Ro60_pn_num, c(0,1), c('(-)', '(+)'))
  )

data_for_fig2.ver.1[ ,'disease_group'] <- gsub('^(Primary |Secondary )(SS)', '\\2',data_for_fig2.ver.1$disease_group_2)


# ver.2: Secondary SS plot is added in the side of plot for all SS.

data_for_fig2.ver.2 <- data_for_fig2.ver.1 %>%
  mutate(
    disease_group_2 = ifelse(
      as.character(
        disease_group_2
      ) %in%
        c('Primary SS', 'Secondary SS'),
      'SS',
      disease_group_2
    )
  ) %>%
  rbind(
    data_for_fig2.ver.1 %>%
      dplyr::filter(
        disease_group_2 == 'Secondary SS'
      )
  ) %>%
  mutate(
    disease_group_2 = factor(
      disease_group_2, 
      c('SS', 'Secondary SS', 3, 4, 5),
      c('SS', 'Secondary SS', 'non SS', 'RA', 'HC')
    )
  )

# Boxplot for Fig.2.1 
#
# var.x is 'disease_group_2'


df.input_mf.boxplot_fig2.1 <-
  data.frame(
    expand.grid(
      var.x=c('disease_group_2'),
      var.y=c("Ro52", "Ro60"),
      scale.var.y = c('not_scale', 'log10'),
      stringsAsFactors = FALSE
    )
  ) %>%
  mutate(
    str = ifelse(var.x=='disease_group_2','dummy', 'disease_group'),
    thre = ifelse(
      var.y=="Ro52", 10, NA
      ),
    Alpha = ifelse(
      var.y %in% c("Ro52"), 0.8, 0
    ),
    x.text.angle = ifelse(
      var.x=='disease_group_2', 45, 0
      ),
    dn.surfix = 'Distribution_of_antiRo52'
  )



# Primary and Secondary SS are separated.

dlply(
  df.input_mf.boxplot_fig2.1,
  .(var.x, var.y, scale.var.y),
  function(df){
    mf.boxplot(
      data = data_for_fig2.ver.1,
#      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      scale.var.y = df$scale.var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5,
            face="bold",
            colour="black"
            ),
        axis.text.x = 
          element_text(
            angle =
              df$x.text.angle,
            vjust = 1, hjust = 1
            ),
        axis.title.x =
          element_text(
            size=0,
            lineheight=1.5, 
            face="bold",
            colour="black")
        ),

      dn.surfix = sprintf('%s_%s_%s',df$dn.surfix,df$var.x,df$scale.var.y),
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y
      )
    }
  )


# Secondary SS plot is added in the side of plot for all SS.

dlply(
  df.input_mf.boxplot_fig2.1,
  .(var.y, scale.var.y),
  function(df){
    mf.boxplot(
      data = data_for_fig2.ver.2,
      #      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      scale.var.y = df$scale.var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5,
            face="bold",
            colour="black"
          ),
        axis.text.x = 
          element_text(
            angle = 45,vjust = 1, hjust = 1),
        axis.title.x =
          element_text(
            size=0,
            lineheight=1.5, 
            face="bold",
            colour="black")
      ),
      
      dn.surfix = sprintf('%s_extSecoSS_%s',df$dn.surfix,df$scale.var.y),
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y,
      beeswarm = FALSE
    )
  }
)


# Boxplot for Fig.2.2 
#
# var.x is  "ACA_np", "Ro60_pn_num", 'Ro52_ext.High'


df.input_mf.boxplot_fig2.2 <-
  data.frame(
    expand.grid(
      var.x=c("ACA_np", "Ro60_pn_num", 'ACA_Ro60','Ro52_ext.High'),
      var.y=c("Ro52", "Ro60"),
      scale.var.y = c('not_scale', 'log10'),
      stringsAsFactors = FALSE
    )
  ) %>%
  mutate(
    str = 'dummy',
    thre = ifelse(
      var.y=="Ro52", 10, NA
    ),
    Alpha = ifelse(
      var.y %in% c("Ro52"), 0.8, 0
    ),
    x.text.angle =  0,
    dn.surfix = 'Distribution_of_antiRo52_SSonly_'
  )



# Primary and Secondary SS are separated.

dlply(
  df.input_mf.boxplot_fig2.2,
  .(var.x, var.y, scale.var.y),
  function(df){
    mf.boxplot(
      data = data_for_fig2.ver.1 %>%
        filter(as.character(disease_group)=='SS') %>%
        mutate(
          ACA_Ro60 =
            ifelse(
              ACA_np=='(-)' & 
                Ro60_pn_num=='(-)',
              '(-/-)',
            ifelse(
              ACA_np=='(-)' & 
                Ro60_pn_num=='(+)',
              '(-/+)',
              ifelse(
                ACA_np=='(+)' & 
                  Ro60_pn_num=='(-)',
                '(+/-)',
                ifelse(
                  ACA_np=='(+)' & 
                    Ro60_pn_num=='(+)',
                  '(+/+)', NA
                  )
                )
              )
            )
            ),
      #      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      scale.var.y = df$scale.var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5,
            face="bold",
            colour="black"
          ),
        axis.text.x = 
          element_text(
            angle =
              df$x.text.angle,
            vjust = 1
          ),
        axis.title.x =
          element_text(
            size=12,
            lineheight=1.5, 
            face="bold",
            colour="black")
      ),
      
      dn.surfix = sprintf('%s_%s_%s',df$dn.surfix,df$var.x,df$scale.var.y),
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y
    )
  }
)


# Secondary SS plot is added in the side of plot for all SS.

dlply(
  df.input_mf.boxplot_fig2.2,
  .(var.x, var.y, scale.var.y),
  function(df){
    mf.boxplot(
      data = data_for_fig2.ver.2 %>%
        filter(as.character(disease_group)=='SS') %>%
        mutate(
          ACA_Ro60 =
            ifelse(
              ACA_np=='(-)' & 
                Ro60_pn_num=='(-)',
              '(-/-)',
              ifelse(
                ACA_np=='(-)' & 
                  Ro60_pn_num=='(+)',
                '(-/+)',
                ifelse(
                  ACA_np=='(+)' & 
                    Ro60_pn_num=='(-)',
                  '(+/-)',
                  ifelse(
                    ACA_np=='(+)' & 
                      Ro60_pn_num=='(+)',
                    '(+/+)', NA
                  )
                )
              )
            )
        ),
      #      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      scale.var.y = df$scale.var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5,
            face="bold",
            colour="black"
          ),
        axis.text.x = 
          element_text(
            angle = df$x.text.angle,
            vjust = 1, hjust = 1
            ),
        axis.title.x =
          element_text(
            size=12,
            lineheight=1.5, 
            face="bold",
            colour="black")
      ),
      
      dn.surfix = sprintf('%s_extSecoSS_%s',df$dn.surfix,df$scale.var.y),
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y,
      beeswarm = FALSE
    )
  }
)



# Boxplot (separated files) -----------------------------------------------------------------

df.input_mf.boxplot <-
  data.frame(
    expand.grid(
      var.x=ESSDAI_component,
      var.y=c("Ro52", "Ro52_log10","Ro60", "Ro60_log10"),
      stringsAsFactors = FALSE
    )
  ) %>%
  mutate(
    str = 'dummy',
    thre = ifelse(
      var.y=="Ro52", 10, 
      ifelse(
        var.y=="Ro52_log10", 1, NA
      )
    ),
    Alpha = ifelse(
      var.y %in% c("Ro52", "Ro52_log10"), 0.8, 0
    ),
    dn.surfix = 'EESDAI_subcmp'
  )


df.input_mf.boxplot_figS4and5 <-
  data.frame(
    expand.grid(
      var.x=ESSDAI_component,
      var.y=c("Ro52", "Ro60"),
      scale.var.y = c('not_scale', 'log10'),
      stringsAsFactors = FALSE
    )
  ) %>%
  mutate(
    str = 'dummy',
    thre = ifelse(
      var.y=="Ro52", 10, NA
    ),
    Alpha = ifelse(
      var.y %in% c("Ro52"), 0.8, 0
    ),
    x.text.angle =  0,
    dn.surfix = 'Distribution_of_antiRo52_SSonly_ESSDAIdomains_'
  )


dlply(
  df.input_mf.boxplot_figS4and5,
  .(var.x, var.y, scale.var.y),
  function(df){
    mf.boxplot(
      data = data_for_fig2.ver.2 %>%
        filter(as.character(disease_group)=='SS'),
      #      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      scale.var.y = df$scale.var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5,
            face="bold",
            colour="black"
          ),
        axis.text.x = 
          element_text(
            angle = df$x.text.angle,
            vjust = 1, hjust = 1
          ),
        axis.title.x =
          element_text(
            size=12,
            lineheight=1.5, 
            face="bold",
            colour="black")
      ),
      
      dn.surfix = sprintf('%s_extSecoSS_%s',df$dn.surfix,df$scale.var.y),
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y,
      beeswarm = FALSE
    )
  }
)



gg.data_tidy <- data_tidy %>%
  spread(var, val) %>%
  ggplot()

dlply(
  df.input_mf.boxplot,
  .(var.x,var.y),
  function(df){
    mf.boxplot(
      data = data_tidy %>% 
        mutate(
          dummy=1,
          val = factor(val, c(0, 1), c('0', '>0'))
        ) %>%
        spread(var, val) %>%
        spread(var_y, val_y),
      #      ggdata = gg.data_tidy,
      var.x = df$var.x, 
      var.y = df$var.y,
      str = df$str, 
      theme.input = theme(
        strip.background = element_blank(), 
        strip.text.x = element_text(size=0, angle=-90, face="bold"),
        axis.text = 
          element_text(
            size=12, 
            lineheight=1.5, 
            face="bold",
            colour="black"
          ),
        axis.title.x =
          element_text(
            size=16, 
            lineheight=1.5, 
            face="bold",
            colour="black")
      ),
      
      dn.surfix = df$dn.surfix,
      plot.y_intcpt.alpha = 0.8,
      plot.y_intcpt = df$thre,
      scale.var.y = 'not_scale',
      var.caption = '',
      ax.lab.x = df$var.x,
      ax.lab.y = df$var.y
    )
  }
)


# Endrant -----------------------------------------------------------------


