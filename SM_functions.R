library('tidyverse')
library('wordbankr')
library('reshape2')

inv_logit <- function(x){
  return(1 / (exp(-x) + 1))
} 

getScoreForAge = function(lm, age, lang, num_items){ #function that takes lm, child's age, and the number of possible words (from WG/WS)
  # just predict
  prop = inv_logit(predict.glm(lm, newdata = data.frame(age=age)))
  return(prop*num_items)
}

getAgeForScore = function(lm, score, num_items){
  proportion = (score + .000001) / num_items #added point .000001 to avoid getting inf delay when score is 0
  # http://www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
  b0 = lm$coefficients[1]
  b1 = lm$coefficients[2]
  predicted_age = (log(proportion / (1-proportion)) - b0)/ b1
  return(predicted_age)
}

constants_eng = list() #create a list called constants
constants_eng[['WG']] = list() #add a WG section to the list
constants_eng[['WG']]$lowest_num_id = 33 #sets the lowest_num_id to 33 (the first question on WG asking 'does your child know X?')
constants_eng[['WG']]$highest_num_id = 430 #sets the highest_num_id to 430 (the last question on WG asking 'does your child know X?')
constants_eng[['WG']]$num_items = constants_eng[['WG']]$highest_num_id - constants_eng[['WG']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WG
constants_eng[['WS']] = list() #add a WS section to the list
constants_eng[['WS']]$lowest_num_id = 1 #sets the lowest_num_id to 1 (the first question on WS asking 'does your child know X?')
constants_eng[['WS']]$highest_num_id = 680 #sets the highest_num_id to 680 (the last question on WS asking 'does your child know X?')
constants_eng[['WS']]$num_items = constants_eng[['WS']]$highest_num_id - constants_eng[['WS']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WS

constants_span = list() #create a list called constants
constants_span[['WG']] = list() #add a WG section to the list
constants_span[['WG']]$lowest_num_id = 1 #sets the lowest_num_id to 33 (the first question on WG asking 'does your child know X?')
constants_span[['WG']]$highest_num_id = 428 #sets the highest_num_id to 430 (the last question on WG asking 'does your child know X?')
constants_span[['WG']]$num_items = 428
constants_span[['WS']] = list() #add a WS section to the list
constants_span[['WS']]$lowest_num_id = 1 #sets the lowest_num_id to 1 (the first question on WS asking 'does your child know X?')
constants_span[['WS']]$highest_num_id = 680 #sets the highest_num_id to 680 (the last question on WS asking 'does your child know X?')
constants_span[['WS']]$num_items = 680


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
  elssp_for_form$expected_score_at_chron_age = sapply(elssp_for_form$Age,
                                                      function(age){getScoreForAge(model, age, num_items=num_items)})
  elssp_for_form$expected_age_for_score = sapply(elssp_for_form$ProductionCDI,
                                                 function(score){getAgeForScore(model, score, num_items)})
  
  print('Computing differences...')	
  elssp_for_form$diff_score_from_expected = -1 * (elssp_for_form$ProductionCDI - elssp_for_form$expected_score_at_chron_age)
  # more negative, more baf
  elssp_for_form$diff_age_from_expected = elssp_for_form$Age - elssp_for_form$expected_age_for_score
  
  if (verbose){
    print(elssp_for_form[,c('SubjectNumber','Age', 
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
  elssp_for_form$expected_score_at_chron_age = sapply(elssp_for_form$Age,
                                                      function(age){getScoreForAge(model, age, num_items=num_items)})
  elssp_for_form$expected_age_for_score = sapply(elssp_for_form$ProductionCDI,
                                                 function(score){getAgeForScore(model, score, num_items=num_items)})
  
  print('Computing differences...')	
  elssp_for_form$diff_score_from_expected = -1 * (elssp_for_form$ProductionCDI - elssp_for_form$expected_score_at_chron_age)
  # more negative, more baf
  elssp_for_form$diff_age_from_expected = elssp_for_form$Age - elssp_for_form$expected_age_for_score
  
  if (verbose){
    print(elssp_for_form[,c('SubjectNumber','Age', 
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