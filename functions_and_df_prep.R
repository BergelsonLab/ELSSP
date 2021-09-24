library('tidyverse')
library('wordbankr')
library('reshape2')
library('dotwhisker')
library('MASS')
library('rcompanion') # for cramerV

#inverse logit function
inv_logit <- function(x){
  return(1 / (exp(-x) + 1))
} 

#get score for a given age using logit function
getScoreForAge = function(lm, age, lang, num_items){ #function that takes lm, child's age, and the number of possible words (from WG/WS)
  # just predict
  prop = inv_logit(predict.glm(lm, newdata = data.frame(age=age)))
  return(prop*num_items)
}

#get age for a given score using logit function
getAgeForScore = function(lm, score, num_items){
  proportion = (score + .000001) / num_items #added point .000001 to avoid getting inf delay when score is 0
  # http://www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
  b0 = lm$coefficients[1]
  b1 = lm$coefficients[2]
  predicted_age = (log(proportion / (1-proportion)) - b0)/ b1
  return(predicted_age)
}

# English constants for growth curve
constants_eng = list() #create a list called constants
constants_eng[['WG']] = list() #add a WG section to the list
constants_eng[['WG']]$lowest_num_id = 33 #sets the lowest_num_id to 33 (the first question on WG asking 'does your child know X?')
constants_eng[['WG']]$highest_num_id = 430 #sets the highest_num_id to 430 (the last question on WG asking 'does your child know X?')
constants_eng[['WG']]$num_items = constants_eng[['WG']]$highest_num_id - constants_eng[['WG']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WG
constants_eng[['WS']] = list() #add a WS section to the list
constants_eng[['WS']]$lowest_num_id = 1 #sets the lowest_num_id to 1 (the first question on WS asking 'does your child know X?')
constants_eng[['WS']]$highest_num_id = 680 #sets the highest_num_id to 680 (the last question on WS asking 'does your child know X?')
constants_eng[['WS']]$num_items = constants_eng[['WS']]$highest_num_id - constants_eng[['WS']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WS

# Spanish constants for growth curve
constants_span = list() #create a list called constants
constants_span[['WG']] = list() #add a WG section to the list
constants_span[['WG']]$lowest_num_id = 1 #sets the lowest_num_id to 33 (the first question on WG asking 'does your child know X?')
constants_span[['WG']]$highest_num_id = 428 #sets the highest_num_id to 430 (the last question on WG asking 'does your child know X?')
constants_span[['WG']]$num_items = 428
constants_span[['WS']] = list() #add a WS section to the list
constants_span[['WS']]$lowest_num_id = 1 #sets the lowest_num_id to 1 (the first question on WS asking 'does your child know X?')
constants_span[['WS']]$highest_num_id = 680 #sets the highest_num_id to 680 (the last question on WS asking 'does your child know X?')
constants_span[['WS']]$num_items = 680

#English growth curves
prepare_elssp_df_eng = function(cdi_form, constants, verbose=F){
  print(paste('Processing ',cdi_form,'...', sep=''))	#Prints "Processing WG..." or WS message
  
  num_items = constants_eng[[cdi_form]][['num_items']]	
  print('Number of items:')
  print(num_items) #Prints number of items possible given CDI version
  
  eng_data <- get_instrument_data(language = "English (American)", 
                                  form = cdi_form, administrations = TRUE) 	%>% filter(norming==TRUE)
  eng_words = subset(eng_data, num_item_id < constants_eng[[cdi_form]][['highest_num_id']] & num_item_id > constants_eng[[cdi_form]][['lowest_num_id']]) #takes the subset of columns related to 'does your child know X word?'
  
  print('Computing counts...') #Prints "Computing counts" message
  counts <- eng_words %>% #creates counts df
    dplyr::filter(!is.na(.data$age)) %>% #filters entries without an age 
    dplyr::mutate(produces = !is.na(.data$value) & .data$value == "produces",
                  understands = !is.na(.data$value) &
                    (.data$value == "understands" | .data$value == "produces")) %>%
    dplyr::select(-.data$value) %>%
    tidyr::gather("measure_name", "value", .data$produces, .data$understands) %>%
    dplyr::filter(.data$measure_name == "produces") %>%
    dplyr::group_by(.data$age, .data$data_id) %>%
    dplyr::summarise(num_true = sum(.data$value),
                     num_false = n() - .data$num_true)
  print(counts)
  
  print('Fitting model...') #prints "Fitting model" message
  model <- stats::glm(cbind(num_true, num_false) ~ age, counts,
                      family = "binomial")
  
  if (verbose){ #if verbose argument is TRUE, prints a summary of the model
    print(summary(model))
  }
  
  
  new_data = data.frame(age = (seq(2*30.5,48*30.5,by=1) / 30.5), cdi_form)
  
  print('Getting scores...')
  
  new_scores = cbind(new_data, data.frame(predict(model, new_data, type='response', se.fit=T)))
  print(names(new_scores))
  new_scores$scores = new_scores$fit * num_items
  new_scores$se_high = new_scores$fit + new_scores$se.fit
  new_scores$se_low = new_scores$fit - new_scores$se.fit
  new_scores$se_high = new_scores$se_high * num_items
  new_scores$se_low = new_scores$se_low * num_items
  new_scores$predict_ages = new_scores$age
  
  print('Getting Wordbank norms...')
  print(num_items)
  wordbank_norms = read.csv(paste("data/vocabulary_norms_table_",cdi_form,"_Prod",".csv", sep=""),
                            stringsAsFactors=F)	
  
  wordbank_norms_melted = melt(wordbank_norms, id.vars = c("language", "form", "measure", "age", "identity"))
  
  elssp_for_form <- subset(elssp_eng, CDIversion == cdi_form)
  elssp_for_form$ProductionCDI_no = num_items - elssp_for_form$ProductionCDI 
  elssp_for_form$expected_score_at_chron_age = sapply(elssp_for_form$Age_in_months,
                                                      function(age){getScoreForAge(model, age, num_items=num_items)})
  elssp_for_form$expected_age_for_score = sapply(elssp_for_form$ProductionCDI,
                                                 function(score){getAgeForScore(model, score, num_items)})
  
  print('Computing differences...')	
  elssp_for_form$diff_score_from_expected = -1 * (elssp_for_form$ProductionCDI - elssp_for_form$expected_score_at_chron_age)
  # more negative, more baf
  elssp_for_form$diff_age_from_expected = elssp_for_form$Age_in_months - elssp_for_form$expected_age_for_score
  
  if (verbose){
    print(elssp_for_form[,c('SubjectNumber','Age_in_months', 
                            'ProductionCDI', 'expected_score_at_chron_age', 
                            'diff_age_from_expected','diff_score_from_expected')])
  }
  
  print(head(elssp_for_form))
  rlist = list()		
  rlist[['normative_growth_curve_model']] = model
  rlist[['samples_from_growth_curve_model']] = new_scores
  rlist[['elssp_df']] = elssp_for_form
  rlist[['wordbank_norms_melted']] = wordbank_norms_melted
  rlist[['cdi_form']] = cdi_form
  assign(paste(cdi_form, "elssp_eng", sep = "_"), elssp_for_form, envir = .GlobalEnv)
  assign(paste(cdi_form, "elssp_eng_curves", sep = "_"), wordbank_norms_melted, envir = .GlobalEnv)
  assign(paste(cdi_form, "elssp_eng_gcurve", sep = "_"), new_scores, envir = .GlobalEnv)
}

#Spanish growth curves
prepare_elssp_df_span = function(cdi_form, constants, verbose=F){
  print(paste('Processing ',cdi_form,'...', sep=''))	#Prints "Processing WG..." or WS message
  
  num_items = constants_span[[cdi_form]][['num_items']] 
  print('Number of items:')
  print(num_items) #Prints number of items possible given CDI version
  
  span_data <- wordbankr::get_instrument_data(language = "Spanish (Mexican)", 
                                              form = cdi_form, administrations = TRUE)	%>% filter(norming==TRUE)
  span_words = subset(span_data, num_item_id < constants_span[[cdi_form]][['highest_num_id']] & num_item_id > constants_eng[[cdi_form]][['lowest_num_id']]) #takes the subset of columns related to 'does your child know X word?'
  
  print('Computing counts...') #Prints "Computing counts" message
  counts <- span_words %>% #creates counts df
    dplyr::filter(!is.na(.data$age)) %>% #filters entries without an age 
    dplyr::mutate(produces = !is.na(.data$value) & .data$value == "produces",
                  understands = !is.na(.data$value) &
                    (.data$value == "understands" | .data$value == "produces")) %>%
    dplyr::select(-.data$value) %>%
    tidyr::gather("measure_name", "value", .data$produces, .data$understands) %>%
    dplyr::filter(.data$measure_name == "produces") %>%
    dplyr::group_by(.data$age, .data$data_id) %>%
    dplyr::summarise(num_true = sum(.data$value),
                     num_false = n() - .data$num_true)
  print(counts)
  
  print('Fitting model...') #prints "Fitting model" message
  model <- stats::glm(cbind(num_true, num_false) ~ age, counts,
                      family = "binomial")
  
  if (verbose){ #if verbose argument is TRUE, prints a summary of the model
    print(summary(model))
  }
  
  
  new_data = data.frame(age = (seq(2*30.5,48*30.5,by=1) / 30.5), cdi_form)
  
  print('Getting scores...')
  
  new_scores = cbind(new_data, data.frame(predict(model, new_data, type='response', se.fit=T)))
  print(names(new_scores))
  new_scores$scores = new_scores$fit * num_items
  new_scores$se_high = new_scores$fit + new_scores$se.fit
  new_scores$se_low = new_scores$fit - new_scores$se.fit
  new_scores$se_high = new_scores$se_high * num_items
  new_scores$se_low = new_scores$se_low * num_items
  new_scores$predict_ages = new_scores$age
  
  print('Getting Wordbank norms...')
  print(num_items)
  wordbank_norms = read.csv(paste("data/vocabulary_norms_table_",cdi_form,"_Prod", "_mexspan",".csv", sep=""),
                            stringsAsFactors=F)	
  
  wordbank_norms_melted = melt(wordbank_norms, id.vars = c("language", "form", "measure", "age", "identity"))
  
  elssp_for_form <- subset(elssp_span, CDIversion == cdi_form)
  elssp_for_form$ProductionCDI_no = num_items - elssp_for_form$ProductionCDI 
  elssp_for_form$expected_score_at_chron_age = sapply(elssp_for_form$Age_in_months,
                                                      function(age){getScoreForAge(model, age, num_items=num_items)})
  elssp_for_form$expected_age_for_score = sapply(elssp_for_form$ProductionCDI,
                                                 function(score){getAgeForScore(model, score, num_items=num_items)})
  
  print('Computing differences...')	
  elssp_for_form$diff_score_from_expected = -1 * (elssp_for_form$ProductionCDI - elssp_for_form$expected_score_at_chron_age)
  # more negative, more baf
  elssp_for_form$diff_age_from_expected = elssp_for_form$Age_in_months - elssp_for_form$expected_age_for_score
  
  if (verbose){
    print(elssp_for_form[,c('SubjectNumber','Age_in_months', 
                            'ProductionCDI', 'expected_score_at_chron_age', 
                            'diff_age_from_expected','diff_score_from_expected')])
  }
  
  head(elssp_for_form, 3)
  
  rlist = list()		
  rlist[['normative_growth_curve_model']] = model
  rlist[['samples_from_growth_curve_model']] = new_scores
  rlist[['elssp_df']] = elssp_for_form
  rlist[['wordbank_norms_melted']] = wordbank_norms_melted
  rlist[['cdi_form']] = cdi_form
  assign(paste(cdi_form, "elssp_span", sep = "_"), elssp_for_form, envir = .GlobalEnv)
  assign(paste(cdi_form, "elssp_span_curves", sep = "_"), wordbank_norms_melted, envir = .GlobalEnv)
  assign(paste(cdi_form, "elssp_span_gcurve", sep = "_"), new_scores, envir = .GlobalEnv)
}



#read in data
elssp <- read.csv("data/ELSSP_SubjectInformation.csv", stringsAsFactors=F, na.strings=c(""," ","NA")) %>% 
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
  filter((Age_in_months>=8 & CDIversion=='WG')|(Age_in_months>=16 & CDIversion=='WS')) %>%
  mutate(diff_age_from_expected = case_when( PrimaryLanguage == 'English' & ProductionCDI==0 ~ (Age_in_months - 9),
                                             PrimaryLanguage == 'Spanish' & ProductionCDI==0 ~ (Age_in_months - 8),
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