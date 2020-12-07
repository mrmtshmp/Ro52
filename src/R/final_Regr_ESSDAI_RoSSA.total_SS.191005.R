# Description: ----------------------------------------------
#   . 
#   The concept of the analysis is in ../Doc/mtg_190402.pdf.
#
# 19/ 04/ 21-
dir.data           <- "../../Data"
.fn.data            <- "AnalysisDataSet_v4.0.RData"

load(sprintf("%s/%s", dir.data, .fn.data))

subset_SS <- "total_SS"

# Settings: -------------------------------------------------

set.seed(333)

require('tibble')
require(tidyverse)



dir.sub               <- "../sub"
fn.require_packages.R <- "require_libraries.R"


data$Ro52_log10 <- log(
  data$Ro52,
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

source(
  sprintf(
    "%s/%s", dir.sub, 'setting_plot.R'
    )
  )


# Settings ----------------------------------------------------------------

theme_set(theme_classic(base_size = 36, base_family = "serif"))

ExploratoryDataAnalysis::setting_plot()

# ordered logistic reg -----------------------------------------------------

data <- data %>%
  mutate(
    ESSDAI_ana = ifelse(
      ESSDAI >=4, 4, ESSDAI
      )
    )

table(data$ESSDAI_ana, data$ESSDAI)
table(data$ESSDAI_ana, data$ESSDAI, data$nijisei)


res.GLM <- lmrob(
  Ro52_log10 ~ ESSDAI_ana,
  data=data %>% 
    filter(
      disease_group=='SS'
      ),
  method = 'MM'
  ) 

confint(res.GLM)

data$predict <- predict(res.GLM, data)

plot(data$predict, data$Ro52_log10)



data_tidy <- data %>%
  filter(disease_group=='SS') %>%
  mutate(
    Ro60_raw = anti_SS_A_ab,
    Ro52_raw = Ro52
    ) %>%
  mutate(Ro60_log10 = log(anti_SS_A_ab, 10)) %>%
  dplyr::select(
    Ro52_log10, Ro60_log10, Ro52_raw, Ro60_raw, ESSDAI_ana,  ESSDAI #, nijisei
    ) %>%
  gather(
    var, val, -ESSDAI_ana, -ESSDAI #, -nijisei
    )


mf.robustbase_smoothing <- function(D){
  
  gg_data <- ggplot(
    D,
    aes(x = factor(ESSDAI_ana, levels = seq(0,4,by = 1),labels = c('0','1','2','3','≥4')), y=val)
    )
  
  str.ylab <- unique(D$var)
  str.xlab <- 'ESSDAI'
  
  res.GLM <- lmrob(
    val ~ ESSDAI_ana,
    data=D,
    method = 'MM'
    ) 
  
  gg_geom_1   <- geom_point(
    col="black", #as.factor(ESSDAI_ana)), 
    size=1.5
    )
  
  gg_geom_box   <- geom_boxplot(
    col="black"
    )
  
  geom_regress <- geom_abline(
    intercept = coef(res.GLM)[1], 
    slope = coef(res.GLM)[2],
    size=2
    )
  
  geom_regress_confint_UL <- geom_abline(
    intercept = confint(res.GLM)[1,2], 
    slope     = confint(res.GLM)[2,2],
    color = 'gray',
    size=2
    )
  
  geom_regress_confint_LL <- geom_abline(
    intercept = confint(res.GLM)[1,1], 
    slope     = confint(res.GLM)[2,1],
    color = 'gray',
    size=2
    )
  
  gg_facet <- facet_grid(
    var~ .,#nijisei, 
    scale='free_y'
    )
  
  output <- list(
    plot(
      gg_data + gg_geom_box + gg_geom_1 +  
        geom_regress + geom_regress_confint_UL + geom_regress_confint_LL + 
        gg_facet + summ_plot_theme + xlab(str.xlab) + ylab(str.ylab)
    ),
    res.GLM
    )
  names(output) <- sprintf('%s.SS', unique(D$var))#, unique(D$nijisei))
  return(
   output
    )
  }



mf.robustbase_smoothing.loess_fit <- function(D){
  
  gg_data <- ggplot(
    D,
    aes(x = as.numeric(ESSDAI), y=val)
  )
  
  str.ylab <- unique(D$var)
  str.xlab <- 'ESSDAI'
  

  
  gg_geom_1   <- geom_point(
    aes(x = factor(ESSDAI, levels = seq(0,17)), y=as.numeric(val)),
    col="black", #as.factor(ESSDAI_ana)), 
    size=1.5
  )
  
  gg_geom_box   <- geom_boxplot(
    aes(x = factor(ESSDAI, levels = seq(0,17)), y=as.numeric(val)),
    col="black"
  )
  
  geom_regress <- geom_smooth(
    method = 'loess',
    col='black',
    size=2
  )
  

  gg_facet <- facet_grid(
    var~ .,#nijisei, 
    scale='free_y'
  )
  
  output <- list(
    plot(
      gg_data + gg_geom_box + gg_geom_1 +  
        geom_regress +
        gg_facet + summ_plot_theme + xlab(str.xlab) + ylab(str.ylab)
    ),
    res.GLM
  )
  names(output) <- sprintf('%s.SS', unique(D$var))#, unique(D$nijisei))
  return(
    output
  )
}


setting_plot()

print(unique(data_tidy$var))

main.output <- dlply(
  data_tidy,
  .(var),#, nijisei),
  mf.robustbase_smoothing
  )

main.output.wo_over_4 <- dlply(
  data_tidy %>% filter(ESSDAI<=4),
  .(var),#, nijisei),
  mf.robustbase_smoothing
)

loess.output <- dlply(
  data_tidy,
  .(var),#, nijisei),
  mf.robustbase_smoothing.loess_fit
)

# Output — plots ----------------------------------------------------------

quartz(
  type = 'pdf',
  file = sprintf(
    '%s/%s_%s.%s.pdf',
    dir.output, 
    'regressline.lmrob_ESSDAI_by_log10Ro',
    .fn.data,
    subset_SS
    ),
  pointsize = 10,
  height = 7,
  width = 7
  )
par(family='serif')
lapply(main.output, FUN = function(L){L[[1]]})
dev.off()



quartz(
  type = 'pdf',
  file = sprintf(
    '%s/%s_%s.%s.pdf',
    dir.output, 
    'regressline.loess_ESSDAI_by_log10Ro',
    .fn.data,
    subset_SS
  ),
  pointsize = 10,
  height = 7,
  width = 7
)
par(family='serif')
lapply(loess.output, FUN = function(L){L[[1]]})
dev.off()


# Output — Results of the robust regressions --------------------------

for(
  i in 1: length(
    names(main.output)
    )
  ){
  xi = data.frame(
    round(
      summary(main.output[[i]][[2]])$coefficients[,c(1,2)], 2
      ),
    'p.val' = format.pval(
      summary(main.output[[i]][[2]])$coefficients[,c(4)], eps = .001, digits = 2)
    ) %>% 
    rownames_to_column('terms')
  
  xi$Title <- names(main.output)[i]
  print(xi)
  
  if(i==1){output=xi} else{ output <- rbind(output, xi)}

  }

write_csv(
  output, 
  path =
    sprintf(
      '%s/%s_%s_%s.%s.csv', 
      dir.output, 
      'result.regressline.lmrob_ESSDAI_by_',
      'all',
      .fn.data,
      subset_SS
      #        names(main.output)[i],
    )
  )

# Endrant -----------------------------------------------------------------


