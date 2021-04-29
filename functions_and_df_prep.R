library('wordbankr')
library('dotwhisker')
library('MASS')
library('rcompanion') # for cramerV
library('corrplot')
library('gplots') # for balloon plots
library('extraoperators')
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
                                       ServicesReceivedPerMonth < 3 ~ "0-2", 
                                       ServicesReceivedPerMonth>7 ~ ">7", 
                                       TRUE ~ "3-6"))) %>% 
  mutate(HLworse_cat = fct_relevel(HLworse_cat, "mild", "moderate", "severe_profound"),
  SPM_cat = (fct_relevel(SPM_cat, "0-2", "3-6", ">7"))) %>%
  mutate_if(is.character, as.factor) 
#add months-delay to dataframes
elssp_eng <- filter(elssp, PrimaryLanguage=="English")
elssp_span <- filter(elssp, PrimaryLanguage=="Spanish")

prepare_elssp_df_eng("WG", constants_eng, verbose=T)
prepare_elssp_df_eng("WS", constants_eng, verbose=T)
prepare_elssp_df_span("WG", constants_span, verbose=T)
prepare_elssp_df_span("WS", constants_span, verbose=T)

elssp_curves <- rbind(WG_elssp_eng, WS_elssp_eng, WG_elssp_span, WS_elssp_span) %>% 
  filter((Age>8 & CDIversion=='WG')|(Age>16 & CDIversion=='WS')) %>%
  mutate(diff_age_from_expected = case_when( PrimaryLanguage == 'English' & ProductionCDI==0 ~ (Age - 9),
                                             PrimaryLanguage == 'Spanish' & ProductionCDI==0 ~ (Age - 8),
                                             TRUE ~ diff_age_from_expected)) %>%
  mutate(Amplification=fct_relevel(Amplification, "none"))

#make dataframes easier to call
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
  rename('Developmental Delay' = DevelopmentalConcerns,
         'Health Issues' = HealthIssues,
         'Language Background' = Monolingual_English,
         Degree = HLworse_cat,
         '1-3-6' = Meets136,
         Prematurity = IsPremature, 
         'Services per Month' = SPM_cat) %>%
  dplyr::select('1-3-6', Amplification, Communication, Degree, 'Developmental Delay', 
                Etiology, Gender, 'Health Issues', 
               'Language Background', Laterality, Prematurity, 'Services per Month') %>%
  mutate(Etiology = na_if(Etiology, "Mixed"))


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
                           TRUE ~ "survives Bonferroni correction"))
  return(out)

})

chisq_output <- function(Var1, Var2) {
  subset_elssp <- elssp %>% filter(elssp[[Var1]]!='' & elssp[[Var2]]!='')
  test <- chisq.test(table(subset_elssp[[Var1]], subset_elssp[[Var2]]))
  chisq_output = paste("($X^2$ (", round(test$parameter,2),
                    ", N = ", sum((test$observed)),
                    ") = ", round(test$statistic, 2),
                    ", p ", printp(test$p.value, digits = 4, add_equals = TRUE),
                    ")", sep='')
  chisq_output
}

beta_output <- function(model, predictornum) {
  beta_output = paste("(ß = ", round(summary(model)$coefficients[predictornum,1], 2),
                       ", p ",printp(summary(model)$coefficients[predictornum,4], digits=2, add_equals=TRUE),
                       ")", sep='')
  beta_output
}

n_condition <- function(condition){
  as.numeric(condition_tallies[paste(enexpr(condition)),"n"])
}

lm_pvalue <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(printp(p, add_equals = TRUE))
}