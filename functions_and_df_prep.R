library('wordbankr')
library('plyr')
library('reshape2')
library('dotwhisker')
library('MASS')
library('tidyverse')
library('rcompanion')
library('fastDummies')
library('corrplot')
library('gplots')
source('SM_functions.R')

#read in data
elssp <- read.csv("data/ELSSP_SubjectInfo_07242020.csv", stringsAsFactors=F) %>% 
  filter(VisitNumber==1) %>%
  mutate(subject_id = unlist(lapply(strsplit(SubjectNumber,'_'), function(x){as.numeric(x[2])}))) %>%
  mutate(anycomorbid = ifelse(VisionLoss == 1 |
                                DevelopmentalConcerns == 1 |
                                HealthIssues == 1 |
                                IsPremature == 1, "1",
                              "0")) %>%
  mutate(HLworse_cat = as.factor(ifelse(HLworse=='NA', 'NA', 
                                        ifelse(HLworse < 40, "mild", 
                                               ifelse(HLworse>70, "severe_profound", 
                                                      "moderate"))))) %>%
  mutate(SPM_cat = as.factor(ifelse(ServicesReceivedPerMonth=='NA', 'NA', 
                                    ifelse(ServicesReceivedPerMonth < 4, "0-3", 
                                           ifelse(ServicesReceivedPerMonth>10, ">10", 
                                                  "4-10")))))
#add months-delay to dataframes
elssp_datasets <- lapply(c('WG','WS'), function(x){
  prepare_elssp_df(x, constants, verbose=T)
})

#make dataframes easier to call
full_elssp <- bind_rows(elssp_datasets[[1]]$elssp_df, elssp_datasets[[2]]$elssp_df) %>% 
  mutate_if(is.character, as.factor) 
#(full_elssp is different from elssp in that it has the growth curve values)

elssp_cat <- dplyr::select(elssp, 
              Gender, HealthIssues, DevelopmentalConcerns, 
              IsPremature, PrimaryLanguage, HLworse_cat, SPM_cat, 
              Communication, Meets136, Laterality, Amplification, Etiology)
combos <- combn(ncol(elssp_cat),2)

chi_sq_all <- adply(combos, 2, function(x) {
  subset_elssp <- elssp_cat %>% filter((elssp_cat)[x[1]] !='' & (elssp_cat[x[2]] !=''))
  test <- chisq.test(subset_elssp[, x[1]], subset_elssp[, x[2]])
  eff.size <- cramerV(subset_elssp[, x[1]], subset_elssp[, x[2]])

  out <- data.frame("Var1" = colnames(elssp_cat)[x[1]], 
                    "Var2" = colnames(elssp_cat[x[2]]), 
                    "Chi.Square" = round(test$statistic,3), 
                    "df"= test$parameter, 
                    "p.value" = round(test$p.value, 6),
                    "eff.size" = eff.size) %>% 
    mutate(sig = ifelse(p.value>.05, "ns", 
                            ifelse(p.value>0.0007575758, "sig", "survivesbc")))
  return(out)

})

nz_balloons <- function(Var1, Var2){ #balloon plots without NA cells
  nz_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  bp <- balloonplot(table(nz_elssp[[Var1]], nz_elssp[[Var2]]),
                    xlab = Var1,
                    ylab = Var2,
                    main = glue("{Var1} by {Var2}"))
}

# lm_output <- function(data) {
#   print_res = paste("ß=", round(data$Estimate,2),
#                       ", SE=", round(data$Est.Error,2),
#                       ", 95% CI [", printnum(round(data$l.95..CI,2), format = "f"),
#                                   ",", round(data$u.95..CI,2),"]", sep='')
#   print_res
# }

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


comorbid <- read.csv("data/elssp_comorbidities.csv") %>% 
  mutate(anycomorbid = ifelse(VisionLoss==1|DevelopmentalConcerns==1|HealthIssues==1|IsPremature==1, 1,
                                                                                     0), 
         extremelypremature = ifelse(WeeksGestation > 33, 0, 1))
comorbid_long <- comorbid %>% pivot_longer(cols = c(ANSD, IsPremature:extremelypremature, -WeeksGestation), 
                                           names_to = "condition", 
                                           values_to = "n") %>% 
  mutate(condition = as.factor(condition))
condition_tallies <- aggregate(n ~ condition, data=comorbid_long, FUN = sum) %>% 
  remove_rownames %>% 
  column_to_rownames(var="condition")
#make.names(condition_tallies$condition)

n_condition <- function(condition){
  as.numeric(condition_tallies[paste(enexpr(condition)),"n"])
}

hearing_means <- function(col, uni_bi, amp){ifelse(is.null(amp),
  mean((elssp %>% filter(Laterality==enexpr(uni_bi)))$col, na.rm = TRUE),
  mean((elssp %>% filter(Laterality==enexpr(uni_bi) & Amplification==enexpr(amp)))$col)
)
}
