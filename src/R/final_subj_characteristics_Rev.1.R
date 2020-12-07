#' Table for subjects' characteristics.
#' 
#' For response to the revierw's comment below:
#' 
#'   5)	Please indicate if the healthy and RA disease control subjects were age, race,
#'    and sex matched to the SS group?
#'    If not, please confirm that there are no significant differences in the distributions
#'    of these parameters between the groups.
#'    

# Settings: -------------------------------------------------


require('tibble')
require(tidyverse)
require(coin)

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

fn.table1_csv <- "tab.dist.age_within_diseaseGroup.csv"
fn.boxplot    <- "plot.dist.age_within_diseaseGroup.pdf"



# Load subroutines: ----------------------------------------------

Bibtex <- FALSE

source(
  sprintf(
    "%s/%s", dir.sub, fn.require_packages.R
    )
  )



# Create analysis dataset -------------------------------------------------

ADS <- data[,c('SS', 'RA', 'HC', 'Sex', 'Age', 'Ro52')] %>%
  gather(group, flag, -Sex, -Age, -Ro52) %>%
  dplyr::filter(flag==1) %>%
  dplyr::select(-flag) %>%
  mutate(
    Group = factor(group, c("SS", "RA", "HC")),
    Sex   = factor(Sex, c("F", "M"))
    ) %>%
  dplyr::select(-group)


# Create boxplot of the distribution of year-old age ----------------------

ggADS <- ggplot(ADS, aes(x=Group, y=Age))

plot.boxplot <- plot(
  ggADS + 
    geom_boxplot() + 
    geom_point(aes(col=log(Ro52,10)),position = position_quasirandom(width = 0.3)) + 
    facet_grid(Sex~.) +
    theme_bw() +
    scale_color_gradient(low = "blue", high = "orange") +
    theme(
      axis.title.x = element_text(size=20), 
      axis.title.y = element_text(size=20),
      axis.text.x =  element_text(size=20, colour = "black"),
      axis.text.y =  element_text(size=20, colour = "black"),
      strip.text = element_text(size = 20),strip.background = element_rect(fill = "white")
      )
  )

# Create table of subjects' characteristics -------------------------------

if(!require(tableone)){
  install.packages("tableone")
  require(tableone)
}

#Create a variable list which we want in Table 1
listVars <- colnames(ADS)

#Define categorical variables
catVars <- listVars[-which(listVars %in% c("Age", "Ro52"))]

#Total Population
table1 <- 
  CreateTableOne(
    vars = listVars, data = ADS, factorVars = catVars, strata = "Group",includeNA = TRUE,)

print.table1.contn <- print(summary(table1$ContTable,nonNormal=c("Age", "Ro52")), quote = TRUE)
print.table1.categ <- print(summary(table1$CatTable,nonNormal=c("Sex")), quote = TRUE)

#' Test for differences between disease groups.
#'  1. RA versus SS
#'  2. HC versus SS
#'  
ADS.RA_SS <- ADS %>% 
  dplyr::filter(
    as.character(Group) %in% 
      c("RA","SS")
    )
ADS.HC_SS <- ADS %>% 
  dplyr::filter(
    as.character(Group) %in% 
      c("HC","SS")
  )

res.wilcox_test.Age_Disease <- sprintf(
  "The p-values for group differences of age \n between \n 1. SS vs HC and  \n 2. SS vs RA \n\n 1. %s \n2. %s\n<Wilcoxon rank-sum test>\n\n\n",
        format.pval(pvalue(coin::wilcox_test(Age~Group,ADS.HC_SS,distribution="exact")), eps = .001, digits = 2),
        format.pval(pvalue(coin::wilcox_test(Age~Group,ADS.RA_SS,distribution="exact")), eps = .001, digits = 2)
        )
res.fisher_test.Sex_Disease <- sprintf(
  "The p-values for group differences of sex \n between  \n 1. SS vs HC and  \n 2. SS vs RA \n\n 1. %s \n2. %s\n<Fisher exact test>\n",
  format.pval(
      fisher.test(y = as.character(ADS.HC_SS$Sex), x = as.character(ADS.HC_SS$Group))$p.value,
    eps = .001, digits = 2
    ),
  format.pval(
      fisher.test(y = as.character(ADS.RA_SS$Sex), x = as.character(ADS.RA_SS$Group))$p.value,
    eps = .001, digits = 2
    )
  )

cat(res.wilcox_test.Age_Disease, file = "result.txt",append = FALSE)
cat(res.fisher_test.Sex_Disease, file = "result.txt",append = TRUE)

print(fisher.test(y = as.character(ADS.HC_SS$Sex), x = as.character(ADS.HC_SS$Group),ADS.HC_SS))
print(fisher.test(y = as.character(ADS.RA_SS$Sex), x = as.character(ADS.RA_SS$Group),ADS.HC_SS))


# Output ------------------------------------------------------------------

print.table1 <- print(table1)

write.csv(print.table1.contn, sprintf("%s/%s_contn.csv", dir.output, fn.table1_csv))
write.csv(print.table1.categ, sprintf("%s/%s_categ.csv", dir.output, fn.table1_csv))

pdf(file = sprintf("%s/%s", dir.output, fn.boxplot))
plot.boxplot
dev.off()

# Endrant -----------------------------------------------------------------
