# Description: ----------------------------------------------
#   Summary statistics. 
#   The concept of the analysis is in ../Doc/mtg_190402.pdf.
#
# 19/ 04/ 21-

dir.data           <- "../Data"
fn.data            <- "AnalysisDataSet_v4.0.R"

load(sprintf("%s/%s", dir.data, fn.data))


# Settings: -------------------------------------------------

set.seed(332)

require(VIM)
require("mice")
require('tibble')

dir.sub               <- "./sub"
fn.require_packages.R <- "require_libraries.R"


data$Ro52_log10 <- log(data$Ro52, base = 10) 

dir.output         <- "../Output/Final/test"


# Load subroutines: ----------------------------------------------

Bibtex <- FALSE

source(sprintf("%s/%s", dir.sub, fn.require_packages.R))



# Multiple Imputation by Chain Equation -----------------------------------

data_SS <- data %>% 
  filter(SS=="1") %>%
  dplyr::select(
    Age,
    # disease_group,
    Ro52,
    Ro52_pn,
    # SS,
    anti_SS_A_ab, # confirmed in mail: 2019/4/7 00:45 that this item is "Ro60"
    Ro60_pn,
    Sex,
    FS,　　　
    Lip_biopsy,
    Dry_mouth,
    Dry_eye,
    ACA_np,
    Raynaud_np,
    RF,
    RF_pn,
    IgG,
    IgG_pn,
    Saxon_test_np,
    Schirmer_test_np,
    ESSDAI,
    HTLV_1,
    nijisei
    ) %>%
  mutate(
    ESSDAI_bin = ifelse(
      as.numeric(as.character(ESSDAI)) < 1,
      0,
      1
    )
  )
  
test.mice_SS <- mice(
  data_SS %>%
    mutate(
      RF_pn = as.factor(RF_pn),
      IgG_pn= as.factor(IgG_pn),
      Saxon_test_np=as.factor(Saxon_test_np),
      Schirmer_test_np=as.factor(Schirmer_test_np),
      ESSDAI= as.numeric(ESSDAI),
      ESSDAI_bin=as.factor(ESSDAI_bin),
      Sex.int = (as.numeric(Sex)-1) * as.numeric(scale(log(Ro52, 10))),
      FS.int = as.numeric(scale(as.numeric(FS))) * as.numeric(scale(log(Ro52, 10))),
      Lip_biopsy.int = (as.numeric(Lip_biopsy)-1) * as.numeric(scale(log(Ro52, 10))),
      Dry_mouth.int = (as.numeric(Dry_mouth)-1) * as.numeric(scale(log(Ro52, 10))),
      Dry_eye.int = (as.numeric(Dry_eye)-1) * as.numeric(scale(log(Ro52, 10))),
      ACA_np.int = (as.numeric(ACA_np)-1) * as.numeric(scale(log(Ro52, 10))),
      Raynaud_np.int = (as.numeric(Raynaud_np)-1) * as.numeric(scale(log(Ro52, 10))),
      RF.int = as.numeric(scale(as.numeric(RF))) * as.numeric(scale(log(Ro52, 10))),
      RF_pn.int = (as.numeric(RF_pn)-1) * as.numeric(scale(log(Ro52, 10))),
      IgG.int = as.numeric(scale(as.numeric(IgG))) * as.numeric(scale(log(Ro52, 10))),
      IgG_pn.int = (as.numeric(IgG_pn)-1) * as.numeric(scale(log(Ro52, 10))),
      Saxon_test_np.int = (as.numeric(Saxon_test_np)-1) * as.numeric(scale(log(Ro52, 10))),
      Schirmer_test_np.int = (as.numeric(Schirmer_test_np)-1) * as.numeric(scale(log(Ro52, 10))),
      ESSDAI.int = (as.numeric(ESSDAI)-1) * as.numeric(scale(log(Ro52, 10))),
      HTLV_1.int = (as.numeric(HTLV_1)-1) * as.numeric(scale(log(Ro52, 10)))
      ),
  method = c(
    'Age'          = "pmm",
    'Ro52'         = "pmm",
    'Ro52_pn'      = "logreg",
    'anti_SS_A_ab' = "pmm",
    'Ro60_pn'="",
    'Sex'="logreg",
    'FS' ="pmm",
    'Lip_biopsy'="polyreg",
    'Dry_mouth'="logreg",
    'Dry_eye'  ="logreg",
    'ACA_np'="",  
    'Raynaud_np'="logreg",
    'RF'  = 'pmm',
    'RF_pn'="logreg",
    'IgG'  = 'pmm',
    'IgG_pn'          ="logreg",
    'Saxon_test_np'   ="logreg",
    'Schirmer_test_np'="logreg",
    'ESSDAI'  = '',
    'ESSDAI_bin'  = '',
    'HTLV_1'="logreg",
    'Sex.int'  = 'pmm',
    'FS.int' ="pmm",
    'Lip_biopsy.int'  = 'pmm',
    'Dry_mouth.int'  = 'pmm',
    'Dry_eye.int'   = 'pmm',
    'ACA_np.int'  = 'pmm',  
    'Raynaud_np.int'  = 'pmm',
    'RF.int'  = 'pmm',
    'RF_pn.int'  = 'pmm',
    'IgG.int'  = 'pmm',
    'IgG_pn.int'  = 'pmm',
    'Saxon_test_np.int' = 'pmm',
    'Schirmer_test_np.int'  = 'pmm',
    'ESSDAI.int'  = 'pmm',
    'HTLV_1.int'  = 'pmm',
    'nijisei'     = 'logreg'
    ),
  m     = 500, # Carpenter and Kenward (2013)， Royston and White (2011) 
  maxit = 10
  )

plot(test.mice_SS)

test.mice_SS$meth


# Imputation of “lip_biopsy_tri” by imputed “lip_biopsy" ------------------

# Lip biopsy

test.mice_SS$imp$Lip_biopsy_tri <- ddply(
  test.mice_SS$imp$Lip_biopsy %>% 
    rownames_to_column("id"),
  .(id),
  function(df){
    res <- ifelse(
      as.numeric(as.character(df)) -1 < 3,
      0,
      as.numeric(as.character(df)) -1
      )
    print(data.frame("res"=t(res)[1:10],"df"=t(df)[1:10]))
    return( t(res))
    }
  ) %>% 
  dplyr::select(-`1`) %>%
  mutate(id=as.numeric(id)) %>%
  arrange(id)%>%
  column_to_rownames("id")


# Object for strip plot ---------------------------------------------------

imp.mice.data_SS.pdf <- dlply(
  data.frame(
    "meth.imp" = test.mice_SS$meth
    ) %>%
    rownames_to_column("varname") %>%
    filter(meth.imp!=""),
  .(varname),
  function(D){
    .varname = as.character(D$varname)
    print(.varname)
    res <- stripplot(
      test.mice_SS,
      data = formula(sprintf("%s~.imp",.varname)),
      pch  = 20, cex=2
      )
    return(res)
    }
  )

# Robust linear regression (robustbase package) ---------------------------
sink("test_sink.list.res.lmrob.Ro60.mice.data_SS_nonNijisei.")

long1 <- mice::complete(
  test.mice_SS, action='long', include=TRUE
  )
long1 <- subset(long1, nijisei = 0)
test.mice_SS_2 <- as.mids(long1)


list.res.lmrob.Ro60.mice.data_SS <- 
  dlply(
    data.frame(
      "meth.imp" = test.mice_SS$meth
    ) %>%
      rownames_to_column("varname") %>%
      filter(
        varname %in% 
          c('nijisei',
            'Age',
            'Ro52',
            'Ro52_pn',
            'anti_SS_A_ab',
            'Ro60_pn',
            'Sex',
            'FS',
            'Lip_biopsy',
            'Dry_mouth',
            'Dry_eye',
            'ACA_np',  
            'Raynaud_np',
            'RF',
            'RF_pn',
            'IgG',
            'IgG_pn',
            'Saxon_test_np',
            'Schirmer_test_np',
            'ESSDAI',
            'ESSDAI_bin',
            'HTLV_1'
            )
          ),
    .(varname),
    function(D){
      .varname = as.character(D$varname)
      print(.varname)
      #      if(is.numeric(test.mice_SS$data[,.varname]))
      # ind.family  =  "gaussian(link='identity')"
      #      else ind.family = "Gamma(link='log')"
      res <- try(
        with(
          test.mice_SS_2,
          lmrob(
            formula(
              sprintf("log(Ro52,10) ~ log(anti_SS_A_ab, 10) * %s",.varname)
              ),
             control = lmrob.control(
               max.it = 500, maxit.scale = 500,
               refine.tol = 1e-3, rel.tol = 1e-7, scale.tol = 1e-10, solve.tol = 1e-7,
               trace.lev = 4,
               # method = "KS2014",
               seed = 1
               ),
             method = "MM"#, init="S"
            
          )
        )
      )
      if(class(res)=="try-error"){
        print(class(res))
        print(res)
        res <- NULL
      }
      print(class(res))
      return(res)
    }
  )

sink()


# Extract regression coefficients (betas) and vcov (vars) matrix -----------------------------------------------

betas<-llply(
  list.res.lmrob.Ro60.mice.data_SS[
    !sapply(                       # remove [[]] of NULL from the list
      list.res.lmrob.Ro60.mice.data_SS, 
      is.null
      )
    ],
  function(L){
    print(names(L))
    res <- MIextract(
      L[[4]],
      fun=coef
    )
    return(res)
  }
)

vars <-llply(
  list.res.lmrob.Ro60.mice.data_SS[
    !sapply(                       # remove [[]] of NULL from the list
      list.res.lmrob.Ro60.mice.data_SS, 
      is.null
      )
    ],
  function(L){
    print(names(L))
    res <- MIextract(
      L[[4]],
      fun = vcov
      )
    return(res)
  }
)  


# Combine the result from multiple imputation -----------------------------

list.pooled.res.lmrob.Ro60.mice.data_SS <- list()
Publist.pooled.res.lmrob.Ro60.mice.data_SS <- list()

for( 
  i in names( betas)
  ){
  res <-MIcombine(
    betas[[i]],
    vars[[i]]
    )
  print(class(res))
  list.pooled.res.lmrob.Ro60.mice.data_SS[[i]] <- res
  }

df.coef.pooled.res.lmrob.Ro60.mice.data_SS <- 
  ldply(
    list.pooled.res.lmrob.Ro60.mice.data_SS,
    function(D){
      res <- summary(D) %>%
        rownames_to_column("terms")

      res$pval <- format(x = pchisq(
        (res$results)**2/
          (res$se)**2, 
        nrow(data_SS)
        ),
        digits = 4
        )
      return(res)
      }
    )


df.coef.pooled.res.lmrob.Ro60.mice.data_SS <-
  df.coef.pooled.res.lmrob.Ro60.mice.data_SS %>%
  mutate(
    results.rnd = round(results, 2),
    se.rnd = round(se, 2),
    lower.rnd = round(`(lower`, 2),
    upper.rnd = round(`upper)`, 2)
    )

 save.image("lmrob.191005_nonNijisei.RData")
 # load("lmrob.190824_nonNijisei.RData")



# Other plots -------------------------------------------------------------

# 001002



# Outputs -----------------------------------------------------------------

# Main results

write.csv(
  df.coef.pooled.res.lmrob.Ro60.mice.data_SS, 
  sprintf(
    "%s/%s",
    dir.output, "df.coef.pooled.res.lmrob_MM_init.S_Ro60.mice.data_SS_nonNijisei.csv"
    )
  )



# Multiple imputations:

pdf(
  sprintf(
    "%s/%s",
    dir.output,"imp.mice_convergence.data_SS.pdf"
    )
  )
plot(test.mice_SS)
dev.off()


pdf(
  sprintf(
    "%s/%s",
    dir.output,
    "imp.mice.data_SS.pdf"
    )
  )
imp.mice.data_SS.pdf
dev.off()



plot(
  log(
    data_SS[data_SS$Dry_mouth==1,'anti_SS_A_ab'],10
    ) %>% 
    apply(1,as.numeric),
  log(
    data_SS[data_SS$Dry_mouth==1,'Ro52'],10
    ) %>%
    apply(1,as.numeric)
  )

curve(
  a2 * x + b2,
  add=T,
  col=2,lty=2,lwd=3
  )

plot(log(data_SS[data_SS$Dry_mouth==0,'anti_SS_A_ab'],10) %>% apply(1,as.numeric),log(data_SS[data_SS$Dry_mouth==0,'Ro52'],10) %>% apply(1,as.numeric))
curve(
  a2*x + b2,
  add=T,
  col=2,lty=2,lwd=3)

plot(log(data_SS[data_SS$Raynaud_np==0,'anti_SS_A_ab'],10) %>% apply(1,as.numeric),log(data_SS[data_SS$Raynaud_np==0,'Ro52'],10) %>% apply(1,as.numeric))
curve(
  a2*x + b2,
  add=T,
  col=2,lty=2,lwd=3)

plot(log(data_SS[data_SS$Raynaud_np==1,'anti_SS_A_ab'],10) %>% apply(1,as.numeric),log(data_SS[data_SS$Raynaud_np==1,'Ro52'],10) %>% apply(1,as.numeric))
curve(
  a2*x + b2,
  add=T,
  col=2,lty=2,lwd=3)

# Endrant -----------------------------------------------------------------


