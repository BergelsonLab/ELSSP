library('wordbankr')
library('reshape2')
library('dotwhisker')
library('MASS')
library('tidyverse')
library('rcompanion') # for cramerV
library('fastDummies') # for correlation table
library('corrplot')
library('gplots') # for balloon plots
library(reshape2)

source('SM_functions.R')

#read in data
elssp <- read.csv("data/ELSSP_SubjectInfo_07242020.csv", stringsAsFactors=F, na.strings=c(""," ","NA")) %>% 
  filter(VisitNumber==1) %>%
  mutate(SubjectNumber = substr(VIHI_ID, 1, 6), 
         anycomorbid = ifelse(VisionLoss == "yes" |
                                DevelopmentalConcerns == "yes" |
                                HealthIssues == "yes" |
                                IsPremature == "yes", "yes",
                              "no"),
         HLworse_cat = as.factor(case_when(is.na(HLworse) ~ NA_character_,
                                           HLworse < 40 ~ "mild", 
                                           HLworse>70 ~ "severe_profound",
                                           TRUE ~  "moderate")), 
         SPM_cat = as.factor(case_when(is.na(ServicesReceivedPerMonth) ~ NA_character_, 
                                       ServicesReceivedPerMonth < 4 ~ "0-3", 
                                       ServicesReceivedPerMonth>10 ~ ">10", 
                                       TRUE ~ "4-10"))) %>% 
  mutate_if(is.character, as.factor) 
#add months-delay to dataframes
elssp_datasets <- lapply(c('WG','WS'), function(x){
  prepare_elssp_df(x, constants, verbose=T)
})
#make dataframes easier to call
full_elssp <- bind_rows(elssp_datasets[[1]]$elssp_df, elssp_datasets[[2]]$elssp_df) %>% 
  mutate_if(is.character, as.factor) 
#(full_elssp is different from elssp in that it has the growth curve values)

comorbid <- read.csv("data/elssp_comorbidities.csv") %>% 
  mutate(anycomorbid = ifelse(VisionLoss==1|DevelopmentalConcerns==1|HealthIssues==1|IsPremature==1, 1,
                              0), 
         extremelypremature = ifelse(WeeksGestation > 33, 0, 1),
         otherchromissues = ifelse(MissingChrom7==1|chromosomalproblems14and15==1, 1, 0),
         ventricular_issues = ifelse(intraventricularhemorrhage==1|periventricularleukomalacia==1, 1, 0),
         heart_issues = ifelse(HeartValveProblems==1|heartabnormalities==1|supraventriculartachycardia==1|
                                 bradycardia==1|coarctationofaorta==1|patentductusarteriosus|perferatedheart==1|
                                 pulmonaryarterystenosis, 1, 0),
         lung_issues = ifelse(asthma==1|bronchopulmonarydysplasia==1|chroniclungdisease|tracheotomy==1|
                                trachoesophagealfistula==1|esophagealatresia==1, 1, 0),
  feeding_issues = ifelse(Gtube==1|reflux==1|dysphagia==1|feedingproblems_nos|esophagealatresia==1, 1, 0),
  illness = ifelse(chroniclungdisease==1|MRSA==1|meningitis==1|frequentfevers==1|Sepsis==1|jaundice==1|
                     pneumonia==1|cytomegalovirus==1, 1, 0),
  preg_birth_comp = ifelse(preeclampsia==1|gestationaldiabetes==1|otherpregnancycomplications==1|jaundice==1|
                             respiratorydistress==1|fetaldistress==1|meconiumaspiration==1|intrauterinegrowthretardation==1, 1, 0),
  mus_skel = ifelse(hypertonia==1|torticollis==1|oralmotorweakness==1|missingribs==1|vertebraldefects==1|
                      bonefracture==1|limbabnormalities==1|PDAscoliosis==1|CHARGEorVATER==1|VACTERL==1|
                      treachercollins==1|sagittalsyntosis==1, 1, 0),
  cleft = ifelse(cleftlip==1|cleftpalate==1,1,0),
  health_issues_other = ifelse(esotropia==1|iriscoloboma==1|dysplastickidney==1|CHARGEorVATER==1|apnea==1|
                                 hydronephrosis==1|seizures==1|laryngomalacia==1|tonguetie==1|midlinedefect==1|hypoglycemia==1|
                                 anemia==1|mitochondrialdisorder==1|skinbreakdown==1|VACTERL==1|analatresia==1|
                                 renalabnormalities==1|enlargedadenoids==1|liptie==1|sacraldimple==1|enlargedliver==1|
                                 petachialrash==1|thrombocytopenia==1, 1, 0),
  pos_acqHL = ifelse(anemia==1|MRSA==1|meningitis==1|frequentfevers==1|Sepsis==1|jaundice==1|pneumonia==1|seizures==1|
                       hydrocephalus==1|cytomegalovirus==1, 1, 0))

comorbid_long <- comorbid %>% pivot_longer(cols = c(ANSD, IsPremature:pos_acqHL, -WeeksGestation), 
                                           names_to = "condition", 
                                           values_to = "n") %>% 
  mutate(condition = as.factor(condition))

condition_tallies <- aggregate(n ~ condition, data=comorbid_long, FUN = sum) %>% 
  remove_rownames() %>% 
  column_to_rownames(var="condition")
#make.names(condition_tallies$condition)



elssp_cat <- elssp %>% 
  dplyr::select(Amplification, Communication, DevelopmentalConcerns, 
                Etiology, Gender, HealthIssues, HLworse_cat, 
                IsPremature, Laterality, Meets136, PrimaryLanguage, SPM_cat)

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
    mutate(sig = case_when(p.value>.05 ~ "not significant", 
                            p.value>0.0007575758 ~ "p<.05", 
                           TRUE ~ "survives bonferroni correction"))
  return(out)

})

nz_balloons <- function(Var1, Var2){ #balloon plots without NA cells
  nz_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  bp <- balloonplot(table(nz_elssp[[Var1]], nz_elssp[[Var2]]),
                    xlab = Var1,
                    ylab = Var2,
                    main = glue("{Var1} by {Var2}"))
}

gg_balloons <- function(Var1, Var2){ #balloon plots without NA cells
  nz_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  nz_table <- melt(table(nz_elssp[[Var1]], nz_elssp[[Var2]])) 
  balloons <- ggplot(nz_table, aes(x = Var1, y = Var2)) +
    geom_point(aes(size=value))+
    theme(panel.background=element_blank(), 
          panel.border = element_rect(fill=NA, size=1)) +
    labs(x={Var1}, y={Var2})
  print(balloons)
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
                                                  "PrimaryLanguage", "DevelopmentalConcerns", 
                                                  "HealthIssues", "IsPremature")) %>%
  dplyr::select(ServicesReceivedPerMonth, HLworse, Gender_male,
                PrimaryLanguage_English, Communication_spoken, DevelopmentalConcerns_yes, 
                IsPremature_yes, HealthIssues_yes, Meets136_yes) 

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

#this is the function that i use in the paper
bp_simple <- function(x_col, plottitle) {
  full_elssp %>% 
    drop_na(x_col) %>%
  ggplot(aes_string(x = x_col, y="diff_age_from_expected")) +
    geom_boxplot(color = "mediumpurple1", fill = "mediumpurple1", alpha = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.2, color = "mediumpurple1", fill = "mediumpurple1", alpha = .8,
                aes(shape = CDIversion)) +
    xlab("") +
    ylab("") +
    ggtitle(plottitle) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none") +
    theme(plot.title=element_text(size=16))
}

#these are two plotting functions that i don't use in the paper. 
#the first one separates WG/WS by creating separate boxes on the same graph
#the second one (bp_facet) has separate facets for WG/WS
bp_double <- function(x_col, plottitle) {
  full_elssp_dropped <- full_elssp %>% 
    drop_na(x_col) 
  
    ggplot(data = full_elssp_dropped, aes_string(x = x_col, y="diff_age_from_expected")) +
    geom_boxplot(data = (full_elssp_dropped %>% filter(CDIversion=='WG')),color = "mediumpurple1", fill = "mediumpurple1", alpha = 0.2, outlier.shape = NA) +
    geom_jitter(data = (full_elssp_dropped %>% filter(CDIversion=='WG')),width = 0.2, color = "mediumpurple1", fill = "mediumpurple1", alpha = .8,
                aes(shape = CDIversion)) +
    geom_boxplot(data = (full_elssp_dropped %>% filter(CDIversion=='WS')),color = "lightskyblue", fill = "lightskyblue", alpha = 0.2, outlier.shape = NA) +
    geom_jitter(data = (full_elssp_dropped %>% filter(CDIversion=='WS')),width = 0.2, color = "lightskyblue", fill = "lightskyblue", alpha = .8,
                aes(shape = CDIversion)) +
    xlab("") +
    ylab("") +
    ggtitle(plottitle) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none") +
    theme(plot.title=element_text(size=24))
}

bp_facet <- function(x_col, plottitle) {
  full_elssp %>% 
    drop_na(x_col) %>%
    ggplot(aes_string(x = x_col, y="diff_age_from_expected")) +
    geom_boxplot(color = "mediumpurple1", fill = "mediumpurple1", alpha = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.2, color = "mediumpurple1", fill = "mediumpurple1", alpha = .8,
                aes(shape = CDIversion)) +
    xlab("") +
    ylab("") +
    ggtitle(plottitle) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none") +
    theme(plot.title=element_text(size=24)) + 
    facet_wrap(CDIversion ~ .)
}