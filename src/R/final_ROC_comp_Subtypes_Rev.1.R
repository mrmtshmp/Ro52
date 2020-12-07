# Description: ----------------------------------------------
#   . 
#   The concept of the analysis is in ../Doc/mtg_190402.pdf.
#
# 19/ 04/ 21-


# Settings: -------------------------------------------------

set.seed(333)


require('tibble')
require(tidyverse)


dir.sub               <- "../sub"
fn.require_packages.R <- "require_libraries.R"

dir.data           <- "../../Data"
fn.data            <- "AnalysisDataSet_v5.0.RData"

load(sprintf("%s/%s", dir.data, fn.data))

data$Ro52_log10 <- log(
  data$Ro52,
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


# Dichotomyze outcomes for Lip_biopsy and ESSDA  -----------------------------------------------------

data_Dichotomyze <- data %>%
  mutate(
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
    ESSDAI_Dimyze_by_0 =
      ifelse(
        is.na(ESSDAI), NA,
        ifelse(
          ESSDAI < 1, 0, 1
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
  mutate(
    Lip_biopsy_Dimyze_by_0 = factor(Lip_biopsy_Dimyze_by_0, levels = c(0,1), labels = c("Negative", "Positive")),
    Lip_biopsy_Dimyze_by_3 = factor(Lip_biopsy_Dimyze_by_3, levels = c(0,1), labels = c("Negative", "Positive")),
    ESSDAI_Dimyze_by_0 = factor(ESSDAI_Dimyze_by_0, levels = c(0,1), labels = c("Negative", "Positive")),
    ESSDAI_Dimyze_by_2 = factor(ESSDAI_Dimyze_by_2, levels = c(0,1), labels = c("Negative", "Positive")),
    ESSDAI_Dimyze_by_5 = factor(ESSDAI_Dimyze_by_5, levels = c(0,1), labels = c("Negative", "Positive")),
    FS_Dimyze_by_0 = factor(FS_Dimyze_by_0, levels = c(0,1), labels = c("Negative", "Positive")),
    FS_Dimyze_by_2 = factor(FS_Dimyze_by_2, levels = c(0,1), labels = c("Negative", "Positive")),
    FS_Dimyze_by_4 = factor(FS_Dimyze_by_4, levels = c(0,1), labels = c("Negative", "Positive"))
    )
  


# Data shaping ----------------------------------------------------------

data_tidy <- data_Dichotomyze %>%
  
  filter(
    disease_group=="SS"
    ) %>%
  
  
  gather(
    var, val, 
    -SubjID, -anti_SS_A_ab, -Ro52_log10, -Ro52, -Age
    ) %>%
  
  
  mutate(
    Ro60_log10 = log(anti_SS_A_ab,base = 10)
    ) %>%
  
  filter(
    var %in% 
      c(
        #"Sex", 
        "ESSDAI_Dimyze_by_0", 'ESSDAI_Dimyze_by_2', 'ESSDAI_Dimyze_by_5', 
        "FS_Dimyze_by_0","FS_Dimyze_by_2","FS_Dimyze_by_4", 
        "HTLV_1",
        "Lip_biopsy_Dimyze_by_0", "Lip_biopsy_Dimyze_by_3",
        "Dry_mouth", "Dry_eye", "ACA_np", "Raynaud_np", "anti_SS_B_ab_pn",
        "RF_pn", "IgG_pn", "Saxon_test_np", "Schirmer_test_np"
        )
    ) %>%
  
  mutate(
    raw = ifelse(
      var %in% c(
        'ESSDAI_Dimyze_by_0', 'ESSDAI_Dimyze_by_2','ESSDAI_Dimyze_by_5'
        ),
      "ESSDAI",
      ifelse(
        var %in% c(
          'Lip_biopsy_Dimyze_by_0','Lip_biopsy_Dimyze_by_3'
          ),
        'Lip_biopsy',
        ifelse(
          var %in% c(
            'FS_Dimyze_by_0','FS_Dimyze_by_2','FS_Dimyze_by_4' 
          ),
          'FS',
          ifelse(
            var %in% c('RF_pn'),
            'RF',
            var
            )
          )
        )
      )
    ) %>%
  
  left_join(
    data %>%
      mutate(
        RF = ifelse(
          is.na(RF), NA,
          ifelse(
            grep(pattern='<', RF),
            0, as.numeric(RF)
            )
        )
      ) %>%
      gather(raw, raw_val, -SubjID, -nijisei),
    by=c('SubjID', 'raw')
  ) %>%
  mutate(raw_val = as.numeric(raw_val))


# Distribution of raw values within dichotomizes levels --------------------------------------------------


ls.ggdata <- dlply(
  data_tidy%>%
    filter(
      var %in% 
        c(
          #"Sex", 
          "ESSDAI_Dimyze_by_0", 'ESSDAI_Dimyze_by_2', 'ESSDAI_Dimyze_by_5', 
          "FS_Dimyze_by_0","FS_Dimyze_by_2","FS_Dimyze_by_4", 
          "Lip_biopsy_Dimyze_by_0", "Lip_biopsy_Dimyze_by_3",
          "RF_pn", "IgG_pn"
        )
    ),
  .(var),
  function(D){
    var_name = unique(D$var)
    ggdata <- ggplot(
      D ,
      aes(x = raw_val)
      )
    plot(
      ggdata + 
        geom_histogram() + 
        facet_grid(nijisei~val, scales = 'fixed') +
        ggtitle(label = var_name)
      )
    }
  )



pdf(
  file = sprintf(
    '%s/%s.pdf',
    dir.output, 
    'Distribution_of_raw_values_within_DicLevels'
    ),
  pointsize = 16
  )
ls.ggdata
dev.off()


# GLM -----------------------------------------------------

sub.pri_sec <- FALSE

print(sub.pri_sec)

if(sub.pri_sec==TRUE){
  data_tidy$Subgroup <- data_tidy$nijisei
  }else{
  data_tidy$Subgroup <- 1
  }


data_tidy.roc <- 
  data_tidy %>%
  mutate(
    val=as.factor(val)
    )


res.GLM <- data_tidy.roc %>%
    
  dlply(
    .(var, Subgroup),
    function(D){
      .var = D$var[1]
      print(.var)
      .formula.Ro60 = formula(
        "val ~ anti_SS_A_ab"
        )
      .formula.Ro60_log10 = formula(
        "val ~ Ro60_log10"
      )
      .formula.Ro52_log10 = formula(
        "val ~ Ro52_log10"
        )
      .formula.BothRo_log10 = formula(
        "val ~ Ro60_log10 + Ro52_log10"
        )
      
      res.Ro60 <- glm(
        .formula.Ro60,
        data = D,
        family = binomial(link="logit")
        ) 
      res.Ro60_log10 <- glm(
        .formula.Ro60_log10,
        data = D,
        family = binomial(link="logit")
        ) 
      res.Ro52_log10 <- glm(
        .formula.Ro52_log10,
        data = D,
        family = binomial(link="logit")
        )
      res.BothRo_log10 <- glm(
        .formula.BothRo_log10,
        data = D,
        family = binomial(link="logit")
        )
      
      D$res.predict_Ro60       <- predict(res.Ro60, D)
      D$res.predict_Ro60_log10 <- predict(res.Ro60_log10, D)
      D$res.predict_Ro52_log10 <- predict(res.Ro52_log10, D)
      D$res.predict_BothRo_log10 <- predict(res.BothRo_log10, D)
      
      init_thre <- D[
        which( abs(D$Ro52-500) == min(abs(D$Ro52-500)) ),
        ]
      
      D$res.predict_Ro60.init_thre       <- predict(res.Ro60, init_thre)
      D$res.predict_Ro60_log10.init_thre <- predict(res.Ro60_log10, init_thre)
      D$res.predict_Ro52_log10.init_thre <- predict(res.Ro52_log10, init_thre)
      D$res.predict_BothRo_log10.init_thre <- predict(res.BothRo_log10, init_thre)
      
      
      res.coef_Ro60    <- res.Ro60 %>%
        summary() %>%
        coef() %>%
        data.frame() %>%
        rownames_to_column("term")
      
      res.coef_Ro60_log10    <- res.Ro60_log10 %>%
        summary() %>%
        coef() %>%
        data.frame() %>%
        rownames_to_column("term")
      
      res.coef_Ro52_log10   <- res.Ro52_log10 %>%
        summary() %>%
        coef() %>%
        data.frame() %>%
        rownames_to_column("term")
      
      res.coef_BothRo_log10   <- res.BothRo_log10 %>%
        summary() %>%
        coef() %>%
        data.frame() %>%
        rownames_to_column("term")
      
      print(head(D))
      
      return(
      list(
        D, 
        res.coef_Ro60, 
        res.coef_Ro60_log10, 
        res.coef_Ro52_log10,
        res.coef_Ro52_log10
        )
      )
      }
    )


# ROC ---------------------------------------------------------------------

list.res.GLM.predict  <- lapply(
  res.GLM, 
  function(L){
    L[[1]]      # list.res.GLM.predict is a list of the lists with [[1]]s are data.frame with column of predicted values
    }
  ) #%>% llply()


list.ROC_Ro52_log10 <- llply(
  list.res.GLM.predict,
  function(D){
    res.ROC <-
      pROC::roc(
        D$val~D$res.predict_Ro52_log10,
        ci=TRUE, of='thresholds'
        )
    
    res.ROC_Ro60 <-
      pROC::roc(
        D$val~D$res.predict_Ro60_log10,
        ci=TRUE, of='thresholds'
        )
    
    res.roc.test <- roc.test(
      roc1 = res.ROC, 
      roc2 = res.ROC_Ro60,
      method = "bootstrap",
      boot.n = 2000
      )
    
    # ci.thresholds_Ro52 <- 
    #   ci.thresholds(
    #     res.ROC,
    #     conf.level=0.95, 
    #     boot.n=2000,
    #     boot.stratified=TRUE, 
    #     thresholds = "local maximas",
    #     progress=getOption("pROCProgress")$name, 
    #     parallel=FALSE) 
    # 
    # ci.thresholds_Ro60 <- 
    #   ci.thresholds(
    #     res.ROC_Ro60,
    #     conf.level=0.95, 
    #     boot.n=2000,
    #     boot.stratified=TRUE, 
    #     thresholds = "local maximas",
    #     progress=getOption("pROCProgress")$name, 
    #     parallel=FALSE)
    
    init_thre.pred.Ro52 <- unique(D$res.predict_Ro52_log10.init_thre)    
    init_thre.pred.Ro60 <- unique(D$res.predict_Ro60_log10.init_thre)
    
    var <- unique(D$var)
    return(
      list(
        res.ROC, var, 
        res.ROC_Ro60, res.roc.test,
        init_thre.pred.Ro52, init_thre.pred.Ro60#,
#        ci.thresholds_Ro52, ci.thresholds_Ro60
        )
      )
    }
  )

if(sub.pri_sec==TRUE){
  save(
    list.ROC_Ro52_log10,
    file = sprintf(
      '%s/%s_%s.pdf', dir.output, 
      'list.ROC_nonNijisei',
      fn.data
      ),
    pointsize = 16
  )
}else{
  save(
    list.ROC_Ro52_log10,
    file = sprintf(
      '%s/%s_%s.pdf', dir.output, 
      'list.ROC',
      fn.data
    )
  ) 
  }



# Output ------------------------------------------------------------------

quartz(
  type = 'pdf',
  file = sprintf(
    '%s/%s_%s.pdf', dir.output, 
    'list.ROC_Ro52_log10',
    fn.data
    ),
  pointsize = 25
  )
llply(
  list.ROC_Ro52_log10,
  function(L){
    plot(
      L[[1]], main=L[[2]],  #family = 'arial',
      # print.thres=L[[5]], print.thres.cex=2,
      # print.thres.col='royalblue',
      # print.thres.pattern=sprintf("%s%s%s",'','',''),
      ci=TRUE,
      ci.type='bars',
      ci.col= "black",
      segments.lwd= 0.5,
      lwd=4,
      identity.col="black",
      cex.main = 0.5,
      cex.axix = 20,
      cex.lab = 0.5
      )
    # plot(
    #   L[[7]], ci.type='shape'
    #   )
    
    plot(L[[3]], lty = 2, add=TRUE,  #family = 'serif',
         # print.thres=L[[6]], print.thres.cex=0,
         # print.thres.col='royalblue',
         # print.thres.pattern=sprintf("%s%s%s",'','',''),
         ci=TRUE,
         ci.type='bars',
         ci.col= "red",
         lwd=4
         )
    # plot(
    #   L[[8]], ci.type='shape'
    #   )
    
    legend(
      x = 0.65, y=0.4, cex = 0.9,
      lwd = c(2,2,0), lty = 1:2,
      legend = c(
        sprintf("AUC = %s", round(auc(L[[4]]$roc1),2)),
        sprintf("AUC = %s", round(auc(L[[4]]$roc2),2)),
        sprintf("p = %s", round(L[[4]]$p.value,3))
        ),
        bty = "n"
      )
    }
  )
dev.off()


# Endrant -----------------------------------------------------------------
