library(tidyverse)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(wordbankr)
library(reshape2)
library(quantregGrowth)

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


plot_elssp_df = function(elssp_dataset, split=NULL, save=T) {
  
  wordbank_norms_melted = elssp_dataset$wordbank_norms_melted
  samples_from_growth_curve_model = elssp_dataset$samples_from_growth_curve_model
  elssp_df = elssp_dataset$elssp_df
  if (!is.null(split)){
    elssp_df = elssp_df[!is.na(elssp_df[[split]]) & elssp_df[[split]] != '',]
  }
  cdi_form = elssp_dataset$cdi_form
  
  if (!is.null(split)){ #if there's a variable listed in split position
    elssp_df$split  = as.factor(elssp_df[[split]]) # make the split var in the df on the fly
    print(unique(elssp_df$split)) #print all of the factor levels for the grouping variable
    p1 = ggplot()  + 
      geom_line(data=wordbank_norms_melted, aes(x=age, y=value,  group=variable), linetype= 'dashed', alpha= .5) + 
      geom_line(data = samples_from_growth_curve_model, aes(x=predict_ages, y=scores),colour='black', alpha= .5) + 
      theme_classic() + ylab(paste0("CDI Score (",cdi_form,")")) + 
      xlab("Age in Months") + 
      coord_cartesian(xlim=c(0,40)) + 
      geom_point(data = subset(elssp_df, !is.na(ProductionCDI)), aes(x= Age, y= ProductionCDI, shape=split, color=split)) +
      labs(colour=split, shape=split)
    
    p2 = make_rainplot(elssp_df, split, 'diff_score_from_expected', 'CDI Score Deficit') + 
      ggtitle(paste0(cdi_form,': ', split)) + 
      coord_cartesian(ylim = c(-2, 25)) + 
      coord_flip()
    
    p3 = make_rainplot(elssp_df, split, 'diff_age_from_expected', 'Delay in Months') + 
      coord_cartesian(ylim = c(-2, 25)) + 
      coord_flip() + 
      ggtitle(paste0(cdi_form,': ', split))
    
    # Save to disk
    p1 = p1 + theme_classic(base_size=16)
    ggsave(plot=p1, paste('figures/', elssp_dataset$cdi_form,'_', split, '_trajectory.pdf', sep=''), width=6, height=4)
    p2 = p2 + theme_classic(base_size=16)
    ggsave(plot=p2, paste('figures/', elssp_dataset$cdi_form,'_', split, '_deficit.pdf', sep=''), width=6, height=4)
    p3 = p3 + theme_classic(base_size=16)
    ggsave(plot=p3, paste('figures/', elssp_dataset$cdi_form,'_', split, '_delay.pdf', sep=''), width=6, height=4)
    
    graph <- p3
    return(graph)
    
  } else { #if there's NO variable in the split position
    
    p1 = ggplot() + geom_point (data=wordbank_norms_melted, aes(x=age, y=value, colour=variable)) +
      theme_classic() + geom_line(data=wordbank_norms_melted, aes(x=age, y=value,
                                                                  colour=variable)) + 
      geom_line(data = samples_from_growth_curve_model, aes(x=predict_ages, y=scores), colour='black') + 
      geom_point(data = subset(elssp_df, !is.na(ProductionCDI)),aes(x= Age, y= ProductionCDI), shape=17) + 
      ggtitle(paste("Normative", cdi_form, "vs. ELSSP"))
    
    seq_and_labels = seq(from=0,to=48,by=6)
    p1 = p1 + scale_x_continuous(breaks = seq_and_labels, labels = seq_and_labels)
    p1 = p1 + ylab(paste0("CDI Score (",cdi_form,")"))
    p1 = p1 + xlab("Age in Months")
    p1 = p1 + coord_cartesian(xlim=c(0,40))
    
    # Jupyter preview
    options(repr.plot.width=4, repr.plot.height=3)
    print(p1)
    # Save to disk
    p1 = p1 + theme_classic(base_size=16) + theme(legend.position = "none")
    ggsave(paste('figures/', elssp_dataset$cdi_form,'_trajectories.pdf', sep=''), width=6, height=4)
    graph <- p1
  }
  return(graph)
}

filter_to_longitudinal_admins = function(elssp_df){
  
  admins_per_subject_number = aggregate(admin_id ~ subject_id, elssp_df, length)
  subjects_with_multiple_admins = subset(admins_per_subject_number,  admin_id > 1)$subject_id
  return(subset(elssp_df, subject_id %in% subjects_with_multiple_admins))
}


lms_elssp_df = function(elssp_dataset) {
  
  wordbank_norms_melted = elssp_dataset$wordbank_norms_melted
  samples_from_growth_curve_model = elssp_dataset$samples_from_growth_curve_model
  elssp_df = elssp_dataset$elssp_df
  cdi_form = elssp_dataset$cdi_form
  
  # models go here
}

raincloud_theme <- theme( #create theme for raincloud plots
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))

string_to_var = function(str){
  return(eval(parse(text = str)))
}

make_rainplot = function(data_df, xvar, yvar,ylabel){
  
  lb <- function(x) mean(x) - sd(x)
  ub <- function(x) mean(x) + sd(x)
  
  data_df = data_df[!is.na(data_df[[xvar]]),]
  
  sum_df = do.call('data.frame', aggregate(as.formula(paste(yvar,'~', xvar)), data_df, FUN = function(x){
    c(mean = mean(x), median=median(x), lb = lb(x), ub = ub(x), n =length(x))}))
  
  names(sum_df) = gsub(paste(yvar,'.',sep=''),'',names(sum_df))
  print(sum_df)
  sum_df[[xvar]] = as.factor(sum_df[[xvar]])
  data_df[[xvar]] = as.factor(data_df[[xvar]])
  
  #note the reversal
  p1 <- ggplot(data = data_df,
               aes_string(x = xvar, y = yvar), colour = 'green3', fill='palegreen1') +
    geom_point(aes_string(y = yvar),
               position = position_jitter(width = .15), size = .5, alpha = 0.8, color = 'green3') +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5, colour='green3', fill='palegreen1') +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    ylab(ylabel) +
    theme_bw() + geom_hline(yintercept=0, colour='black',linetype='dashed') + theme(aspect.ratio=2/4)
  
  if (xvar == 'diff_age_from_expected'){
    p1 = p1 + ylab(ylabel)
  } else if (xvar == 'diff_score_from_expected'){
    p1 = p1 + ylab(ylabel)
  }
  return(p1)
}

get_empirical_percentile_rank = function(score, age, monthly_scores_for_instrument){
  stop('not implemented')
}


get_monthly_scores_for_instrument = function(cdi_form){
  
  eng_data <- get_instrument_data(language = "English (American)",
                                  form = cdi_form, administrations = TRUE)
  eng_words = subset(eng_data, num_item_id < constants_eng[[cdi_form]][['highest_num_id']] & num_item_id > constants_eng[[cdi_form]][['lowest_num_id']])
  
  print('Computing counts...')
  counts = eng_words %>%
    dplyr::filter(!is.na(.data$age)) %>%
    dplyr::mutate(produces = !is.na(.data$value) & .data$value == "produces",
                  understands = !is.na(.data$value) &
                    (.data$value == "understands" | .data$value == "produces")) %>%
    dplyr::select(-.data$value) %>%
    tidyr::gather("measure_name", "value", .data$produces, .data$understands) %>%
    dplyr::filter(.data$measure_name == "produces") %>%
    dplyr::group_by(.data$age, .data$data_id) %>%
    dplyr::summarise(num_true = sum(.data$value),
                     num_false = n() - .data$num_true)
  return(counts)
}


dropIfNA = function(df, vars){
  rdf = df
  for (colname in vars){
    rdf = rdf[(!is.na(rdf[[colname]]) & !is.nan(rdf[[colname]])),]
  }
  return(rdf[,vars])
}


all.subsets <- function(set) {
  n <- length(set)
  bin <- expand.grid(rlply(n, c(F, T)))
  sets = mlply(bin, function(...) { set[c(...)] })
  unname(sets[sapply(sets, function(set){length(set) > 0})])
}


## Mika's quantile regression
#' Fit quantiles to vocabulary sizes using quantile regression
#'
#' @param vocab_data A data frame returned by \code{get_administration_data}.
#' @param measure A column of \code{vocab_data} with vocabulary values
#'   (\code{production} or \code{comprehension}).
#' @param group (Optional) A column of \code{vocab_data} to group by.
#' @param quantiles Either one of "standard" (default), "deciles", "quintiles",
#'   "quartiles", "median", or a numeric vector of quantile values.
#'
#' @importFrom quantregGrowth ps
#' @importFrom rlang ":="
#'
#' @return A data frame with the columns "language", "form", "age", \code{group}
#'   (if specified), "quantile", and \code{measure}, where \code{measure} is the
#'   fit vocabulary value for that quantile at that age.
#' @export
#'
#' @examples
#' \dontrun{
#' eng_ws <- get_administration_data("English (American)", "WS")
#' fit_vocab_quantiles(eng_ws, production)
#' fit_vocab_quantiles(eng_ws, production, sex)
#' fit_vocab_quantiles(eng_ws, production, quantiles = "quartiles")
#' }
fit_vocab_quantiles <- function(vocab_data, measure, group = NULL,
                                quantiles = "standard") {

  quo_measure <- rlang::enquo(measure)
  quo_group <- rlang::enquo(group)

  quantile_opts <- list(
    standard = c(0.10, 0.25, 0.50, 0.75, 0.90),
    deciles = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90),
    quintiles = c(0.20, 0.40, 0.60, 0.80),
    quartiles = c(0.25, 0.50, 0.75),
    median = c(0.5)
  )
  if (is.numeric(quantiles)) {
    if (any(quantiles >= 1, quantiles <= 0))
      stop("Numeric quantiles must be between 0 and 1")
    num_quantiles <- quantiles
  } else if (is.character(quantiles) & length(quantiles) == 1) {
    if (!(quantiles %in% names(quantile_opts)))
      stop("Character quantiles must be one of ",
           paste(names(quantile_opts), collapse = ", "))
    num_quantiles <- quantile_opts[[quantiles]]
  } else {
    stop("Quantiles must be numeric vector or character vector of length 1")
  }

  vocab_data <- vocab_data %>% dplyr::group_by(.data$language, .data$form)

  if (!rlang::quo_is_null(quo_group)) {
    vocab_data <- vocab_data %>%
      dplyr::filter(!is.na(!!quo_group)) %>%
      dplyr::group_by(!!quo_group, add = TRUE)
  }

  vocab_models <- vocab_data %>%
    dplyr::rename(vocab = !!quo_measure) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = purrr::pmap(
      list(.data$language, .data$form, .data$data),
      function(lang, frm, df) {
        tryCatch(
          suppressWarnings(
            quantregGrowth::gcrq(vocab ~ ps(age, monotone = 1, lambda = 1000),
                                 tau = num_quantiles, data = df)
          ),
          error = function(e) {
            message(sprintf("Unable to fit model for %s %s", lang, frm))
            return(NULL)
          })
      })) %>%
    dplyr::filter(purrr::map_lgl(.data$model, ~!is.null(.)))

  ages <- data.frame(age = (seq((min(vocab_data$age)*30.5),(max(vocab_data$age))*30.5,by=1) / 30.5))
  get_predicted <- function(vocab_model) {
    vocab_fits <- stats::predict(vocab_model, newdata = ages)
    if (length(vocab_model$taus) == 1)
      vocab_fits <- rlang::set_names(list(vocab_fits), vocab_model$taus)
    vocab_fits %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(age = ages$age) %>%
      tidyr::gather("quantile", "predicted", -.data$age)
  }
  
  vocab_fits <- vocab_models %>%
    dplyr::mutate(predicted = purrr::map(.data$model, get_predicted)) %>%
    dplyr::select(-.data$data, -.data$model) %>%
    tidyr::unnest(cols = .data$predicted) %>%
    dplyr::rename(!!quo_measure := .data$predicted) %>%
    dplyr::mutate(quantile = factor(.data$quantile))
  return(vocab_fits)

}
