library('wordbankr')
library('reshape2')
library('dotwhisker')
library('MASS')
library('tidyverse')
library('rcompanion') # for cramerV
library('fastDummies') # for correlation table
library('corrplot')
library('gplots') # for balloon plots
source('SM_functions.R')

#read in data
elssp <- read.csv("data/ELSSP_SubjectInfo_07242020.csv", stringsAsFactors=F) %>% 
  filter(VisitNumber==1) %>%
  mutate(SubjectNumber = substr(VIHI_ID, 1, 6), 
         anycomorbid = ifelse(VisionLoss == 1 |
                                DevelopmentalConcerns == 1 |
                                HealthIssues == 1 |
                                IsPremature == 1, "1",
                              "0"),
         HLworse_cat = as.factor(case_when(is.na(HLworse) ~ NA_character_,
                                           HLworse < 40 ~ "mild", 
                                           HLworse>70 ~ "severe_profound",
                                           TRUE ~  "moderate")), 
         SPM_cat = as.factor(case_when(is.na(ServicesReceivedPerMonth) ~ NA_character_, 
                                       ServicesReceivedPerMonth < 4 ~ "0-3", 
                                       ServicesReceivedPerMonth>10 ~ ">10", 
                                       TRUE ~ "4-10"))) %>% 
  mutate_if(is.character, as.factor) 

comorbid <- read.csv("data/elssp_comorbidities.csv") %>% 
  mutate(anycomorbid = ifelse(VisionLoss==1|DevelopmentalConcerns==1|HealthIssues==1|IsPremature==1, 1,
                              0), 
         extremelypremature = ifelse(WeeksGestation > 33, 0, 1))

comorbid_long <- comorbid %>% pivot_longer(cols = c(ANSD, IsPremature:extremelypremature, -WeeksGestation), 
                                           names_to = "condition", 
                                           values_to = "n") %>% 
  mutate(condition = as.factor(condition))

condition_tallies <- aggregate(n ~ condition, data=comorbid_long, FUN = sum) %>% 
  remove_rownames() %>% 
  column_to_rownames(var="condition")
#make.names(condition_tallies$condition)


#add months-delay to dataframes
WG_elssp <- prepare_elssp_df('WG', constants, verbose = F)
WS_elssp <- prepare_elssp_df('WS', constants, verbose = F)

#make dataframes easier to call
full_elssp <- bind_rows(WG_elssp$elssp_df, WS_elssp$elssp_df) %>% 
  mutate_if(is.character, as.factor) 
#(full_elssp is different from elssp in that it has the growth curve values)

elssp_cat <- elssp %>% 
  select(Gender, HealthIssues, DevelopmentalConcerns, 
              IsPremature, PrimaryLanguage, HLworse_cat, SPM_cat, 
              Communication, Meets136, Laterality, Amplification, Etiology)

combos <- combn(ncol(elssp_cat),2)

chi_sq_all <- plyr::adply(combos, 2, function(x) {
  subset_elssp <- elssp_cat %>% 
    filter(elssp_cat[x[1]] !='' & elssp_cat[x[2]] !='')
  
  test <- chisq.test(subset_elssp[, x[1]], subset_elssp[, x[2]])

  out <- data.frame("Var1" = colnames(elssp_cat)[x[1]], 
                    "Var2" = colnames(elssp_cat[x[2]]), 
                    "Chi.Square" = round(test$statistic,3), 
                    "df"= test$parameter, 
                    "p.value" = round(test$p.value, 6),
                    "eff.size" = cramerV(subset_elssp[, x[1]], subset_elssp[, x[2]])) %>% 
    mutate(sig = case_when(p.value>.05 ~ "ns", 
                            p.value>0.0007575758 ~ "sig", 
                           TRUE ~ "survivesbc"))
  return(out)

})

nz_balloons <- function(Var1, Var2){ #balloon plots without NA cells
  nz_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  bp <- balloonplot(table(nz_elssp[[Var1]], nz_elssp[[Var2]]),
                    xlab = Var1,
                    ylab = Var2,
                    main = glue("{Var1} by {Var2}"))
}

chisq_output <- function(Var1, Var2) {
  subset_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  test <- chisq.test(table(subset_elssp[[Var1]], subset_elssp[[Var2]]))
  chisq_output = paste("($X^2$ (", round(test$parameter,2),
                    ", N = ", sum((test$observed)),
                    ") = ", round(test$statistic, 2),
                    ", p = ", format.pval(test$p.value, digits = 2),
                    ")", sep='')
  chisq_output
}

beta_output <- function(model, predictornum) {
  beta_output = paste("(ß = ", round(summary(model)$coefficients[predictornum,1], 2),
                       ", p = ",format.pval(summary(model)$coefficients[predictornum,4], digits=2),
                       ")", sep='')
  beta_output
}

corr_prep <- dummy_cols(elssp, select_columns = c("Gender", "Meets136",
                                                  "Laterality", "Communication", 
                                                  "PrimaryLanguage", "DevelopmentalConcerns")) %>%
  dplyr::select(ServicesReceivedPerMonth, HLworse, Gender_male,
                PrimaryLanguage_English, Communication_spoken, DevelopmentalConcerns, 
                HealthIssues, IsPremature, Meets136_yes) 

elssp_corr <- cor(corr_prep, use="pairwise.complete.obs")

p_corr <- cor.mtest(corr_prep) 

dimnames(p_corr$p) <- dimnames(elssp_corr)

n_condition <- function(condition){
  as.numeric(condition_tallies[paste(enexpr(condition)),"n"])
}

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
