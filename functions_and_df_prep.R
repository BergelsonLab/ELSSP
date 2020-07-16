inv_logit <- function(x){
  return(1 / (exp(-x) + 1))
} 

getScoreForAge = function(lm, age, num_items){ #function that takes lm, child's age, and the number of possible words (from WG/WS)
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

constants = list() #create a list called constants
constants[['WG']] = list() #add a WG section to the list
constants[['WG']]$lowest_num_id = 33 #sets the lowest_num_id to 33 (the first question on WG asking 'does your child know X?')
constants[['WG']]$highest_num_id = 430 #sets the highest_num_id to 430 (the last question on WG asking 'does your child know X?')
constants[['WG']]$num_items = constants[['WG']]$highest_num_id - constants[['WG']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WG

constants[['WS']] = list() #add a WS section to the list
constants[['WS']]$lowest_num_id = 1 #sets the lowest_num_id to 1 (the first question on WS asking 'does your child know X?')
constants[['WS']]$highest_num_id = 680 #sets the highest_num_id to 680 (the last question on WS asking 'does your child know X?')
constants[['WS']]$num_items = constants[['WS']]$highest_num_id - constants[['WS']]$lowest_num_id + 1 #substracts lowest_  and highest_num_ids and adds 1, to get highest possible score on WS



prepare_elssp_df = function(cdi_form, constants, verbose=F){
  print(paste('Processing ',cdi_form,'...', sep=''))	#Prints "Processing WG..." or WS message
  
  num_items = (constants[[cdi_form]][['highest_num_id']] - constants[[cdi_form]][['lowest_num_id']]) + 1	
  print('Number of items:')
  print(num_items) #Prints number of items possible given CDI version
  
  eng_data <- get_instrument_data(language = "English (American)", 
                                  form = cdi_form, administrations = TRUE)	
  eng_words = subset(eng_data, num_item_id < constants[[cdi_form]][['highest_num_id']] & num_item_id > constants[[cdi_form]][['lowest_num_id']]) #takes the subset of columns related to 'does your child know X word?'
  
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
  wordbank_norms = read.csv(paste("data/vocabulary_norms_table_",cdi_form,"_Prod.csv", sep=""),
                            stringsAsFactors=F)	
  
  wordbank_norms_melted = melt(wordbank_norms, id.vars = c("language", "form", "measure", "age", "identity"))
  
  elssp_for_form = subset(elssp, CDIversion == cdi_form)
  # elssp_for_form$SubjectNumber = as.factor(elssp_for_form$SubjectNumber)
  # elssp_for_form$DevelopmentalConcerns[elssp_for_form$DevelopmentalConcerns == 1] = 'yes'
  # elssp_for_form$DevelopmentalConcerns[elssp_for_form$DevelopmentalConcerns == 0] = 'no'
  # elssp_for_form$DevelopmentalConcerns = as.factor(elssp_for_form$DevelopmentalConcerns)
  elssp_for_form$ProductionCDI_no = num_items - elssp_for_form$ProductionCDI 
  elssp_for_form$expected_score_at_chron_age = sapply(elssp_for_form$AgeAtEvaluationMonths,
                                                      function(age){getScoreForAge(model, age, num_items)})
  elssp_for_form$expected_age_for_score = sapply(elssp_for_form$ProductionCDI,
                                                 function(score){getAgeForScore(model, score, num_items)})
  
  print('Computing differences...')	
  elssp_for_form$diff_score_from_expected = -1 * (elssp_for_form$ProductionCDI - elssp_for_form$expected_score_at_chron_age)
  # more negative, more baf
  elssp_for_form$diff_age_from_expected = elssp_for_form$AgeAtEvaluationMonths - elssp_for_form$expected_age_for_score
  
  if (verbose){
    print(elssp_for_form[,c('SubjectNumber','AgeAtEvaluationMonths', 
                            'ProductionCDI', 'expected_score_at_chron_age', 
                            'diff_age_from_expected','diff_score_from_expected')])
  }
  
  
  rlist = list()		
  rlist[['normative_growth_curve_model']] = model
  rlist[['samples_from_growth_curve_model']] = new_scores
  rlist[['elssp_df']] = elssp_for_form
  rlist[['wordbank_norms_melted']] = wordbank_norms_melted
  rlist[['cdi_form']] = cdi_form
  return(rlist)	
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
      geom_point(data = subset(elssp_df, !is.na(ProductionCDI)), aes(x= AgeAtEvaluationMonths, y= ProductionCDI, shape=split, color=split)) +
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
      geom_point(data = subset(elssp_df, !is.na(ProductionCDI)),aes(x= AgeAtEvaluationMonths, y= ProductionCDI), shape=17) + 
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

longitudinal_plot_elssp_df = function(elssp_dataset, colorize_by=NULL) {

  wordbank_norms_melted = elssp_dataset$wordbank_norms_melted
  samples_from_growth_curve_model = elssp_dataset$samples_from_growth_curve_model
  elssp_df = elssp_dataset$elssp_df
  if (!is.null(colorize_by)){
    elssp_df = elssp_df[!is.na(elssp_df[[colorize_by]]),]
  }
  cdi_form = elssp_dataset$cdi_form

  options(repr.plot.width=8, repr.plot.height=5)

  longitudinal_children = filter_to_longitudinal_admins(elssp_df)
  print(paste('Found', length(longitudinal_children), 'longitudinal_children.'))


  p1 = ggplot() + geom_point (data=wordbank_norms_melted, aes(x=age, y=value, group=variable), color ='black', alpha=.25) +
    theme_classic() + geom_line(data=wordbank_norms_melted, aes(x=age, y=value, group=variable), color='black', alpha=.25)

  if (!is.null(colorize_by)){
    p1 = p1 + geom_point(data = subset(longitudinal_children, !is.na(ProductionCDI)),
                         aes_string(x= 'AgeAtEvaluationMonths', y= 'ProductionCDI', colour = colorize_by), shape=17
    ) + ggtitle(paste("Normative", cdi_form, "vs. ELSSP")) + geom_line(data = subset(longitudinal_children, !is.na(ProductionCDI)),
                                                                       aes_string(x= 'AgeAtEvaluationMonths', y= 'ProductionCDI', colour = colorize_by, group = 'subject_id'))

  } else {
    p1 = p1 + geom_point(data = subset(longitudinal_children, !is.na(ProductionCDI)),
                         aes_string(x= 'AgeAtEvaluationMonths', y= 'ProductionCDI'), shape=17) + ggtitle(paste("Normative", cdi_form, "vs. ELSSP")) + geom_line(data = subset(longitudinal_children, !is.na(ProductionCDI)),
                                                                                                                                                                aes_string(x= 'AgeAtEvaluationMonths', y= 'ProductionCDI',  group = 'subject_id'))
  }
  p1 = p1 + ylab(paste0("CDI Score (",cdi_form,")"))
  p1 = p1 + xlab("Age in Months")
  seq_and_labels = seq(from=0,to=48,by=6)
  p1 = p1 + scale_x_continuous(breaks = seq_and_labels, labels = seq_and_labels)

  print(p1)
  p1 = p1 + theme_classic(base_size=16)
  if (!is.null(colorize_by)){
    ggsave(plot=p1, paste('figures/', elssp_dataset$cdi_form,'_', colorize_by, '_longitudinal.pdf', sep=''), width=6, height=4)
  } else {
    ggsave(plot=p1, paste('figures/', elssp_dataset$cdi_form,'_longitudinal.pdf', sep=''), width=6, height=4)
  }
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
  eng_words = subset(eng_data, num_item_id < constants[[cdi_form]][['highest_num_id']] & num_item_id > constants[[cdi_form]][['lowest_num_id']])

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

library('wordbankr')
library('plyr')
library('reshape2')
library('dotwhisker')
library('MASS')
library('tidyverse')
library('rcompanion')
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library('fastDummies')
library('corrplot')

#read in data
elssp = read.csv("data/ELSSP_SubjectInfo_02042020.csv", stringsAsFactors=F) %>% mutate(InSample=as.factor(InSample)) %>% filter(VisitNumber==1)
elssp$AgeAtEvaluationMonths <- elssp$Age #for clarity NEED TO DO THIS TIDYVERSE WAY
elssp$subject_id = unlist(lapply(strsplit(elssp$SubjectNumber,'_'), function(x){as.numeric(x[2])}))
elssp$admin_id = elssp$SubjectNumber

elssp$anycomorbid <- ifelse(elssp$VisionLoss == 1 |
                             elssp$DevelopmentalConcerns == 1 |
                             elssp$HealthIssues == 1 |
                             elssp$IsPremature == 1, "1",
                           "0") #create a column called AnyComorbid, which is a "1" if child has Vision Loss, Developmental Concerns, Health Issues, or Prematurity

#add months-delay to dataframes
elssp_datasets = lapply(c('WG','WS'), function(x){
  prepare_elssp_df(x, constants, verbose=T)
})

#make dataframes easier to call
wg_elssp <- elssp_datasets[[1]]$elssp_df %>% mutate(Gender=as.factor(Gender), Laterality=as.factor(Laterality),
                                                    Meets136=as.factor(Meets136), meets13=as.factor(meets13),
                                                    meets6=as.factor(meets6), Etiology=as.factor(Etiology),
                                                    Side=as.factor(Side), ANSD=as.factor(ANSD),
                                                    Amplification=as.factor(Amplification),
                                                    Communication=as.factor(Communication),
                                                    IsPremature=as.factor(IsPremature),
                                                    HealthIssues=as.factor(HealthIssues),
                                                    CDIversion=as.factor(CDIversion),
                                                    Monolingual_English=as.factor(Monolingual_English),
                                                    InSample=as.factor(InSample)) %>% filter(VisitNumber==1) #change this to strings as factors??


ws_elssp <- elssp_datasets[[2]]$elssp_df %>% mutate(Gender=as.factor(Gender), Laterality=as.factor(Laterality),
                                                    Meets136=as.factor(Meets136), meets13=as.factor(meets13),
                                                    meets6=as.factor(meets6), Etiology=as.factor(Etiology),
                                                    Side=as.factor(Side), ANSD=as.factor(ANSD),
                                                    Amplification=as.factor(Amplification),
                                                    Communication=as.factor(Communication),
                                                    IsPremature=as.factor(IsPremature),
                                                    HealthIssues=as.factor(HealthIssues),
                                                    CDIversion=as.factor(CDIversion),
                                                    Monolingual_English=as.factor(Monolingual_English),
                                                    InSample=as.factor(InSample)) %>% filter(VisitNumber==1) #change to strings as factors??

full_elssp <- rbind(wg_elssp, ws_elssp) #bind together wg_elssp and ws_elssp 
#(full_elssp is different from elssp in that it has the growth curve values)


#rainbow cdi plots
wg_rainbow <- plot_elssp_df(elssp_datasets[[1]])
ws_rainbow <- plot_elssp_df(elssp_datasets[[2]])

#categorical contasts

##gender
wg_gender <- plot_elssp_df(elssp_datasets[[1]], 'Gender')
ws_gender <- plot_elssp_df(elssp_datasets[[2]], 'Gender')

##laterality
wg_laterality <- plot_elssp_df(elssp_datasets[[1]], 'Laterality')
ws_laterality <- plot_elssp_df(elssp_datasets[[2]], 'Laterality')

#amplification
wg_amplification <- plot_elssp_df(elssp_datasets[[1]], 'Amplification')
ws_amplification <- plot_elssp_df(elssp_datasets[[2]], 'Amplification')

#healthissues
wg_healthissues <- plot_elssp_df(elssp_datasets[[1]], 'HealthIssues')
ws_healthissues <- plot_elssp_df(elssp_datasets[[2]], 'HealthIssues')

#devconcerns
wg_devconcerns <- plot_elssp_df(elssp_datasets[[1]], 'DevelopmentalConcerns')
ws_devconcerns <- plot_elssp_df(elssp_datasets[[2]], 'DevelopmentalConcerns')

#prematurity
wg_prematurity <- plot_elssp_df(elssp_datasets[[1]], 'IsPremature')
ws_prematurity <- plot_elssp_df(elssp_datasets[[2]], 'IsPremature')

#meets136
wg_meets136 <- plot_elssp_df(elssp_datasets[[1]], 'Meets136')
ws_meets136 <- plot_elssp_df(elssp_datasets[[2]], 'Meets136')

#communication
wg_communication <- plot_elssp_df(elssp_datasets[[1]], 'Communication')
ws_communication <- plot_elssp_df(elssp_datasets[[2]], 'Communication')

#etiology
wg_etiology <- plot_elssp_df(elssp_datasets[[1]], 'Etiology')
ws_etiology <- plot_elssp_df(elssp_datasets[[2]], 'Etiology')

#degree
wg_degree <- ggplot(data=elssp_datasets[[1]]$elssp_df, aes(x=HLworse, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')
ws_degree <- ggplot(data=elssp_datasets[[2]]$elssp_df, aes(x=HLworse, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')

#services/month
wg_services <- ggplot(data=elssp_datasets[[1]]$elssp_df, aes(x=ServicesReceivedPerMonth, y=diff_age_from_expected))+
  geom_point() + 
   geom_smooth(method='lm')
ws_services <- ggplot(data=elssp_datasets[[2]]$elssp_df, aes(x=ServicesReceivedPerMonth, y=diff_age_from_expected))+
  geom_point() + 
  geom_smooth(method='lm')

#interactions among variables
hm_prep <- data.frame(var1=character(),
                      var2=character(),
                      eff_size=double(),
                      pval=double(),
                      sig=factor())

prep_hm_df <- function(df) {
  hm_list <- list(
    "var1" = as.character(word(deparse(substitute(df)), 1, sep = fixed("_"))),
    "var2" = as.character(word(deparse(substitute(df)), 2, sep = fixed("_"))),
    "eff_size" = round(unname(cramerV(df)), digits=3),
    "pval" = (chisq.test(df))$p.value)
}

# prep_hm_catcont <- function(cat, cont) {
#   hm_list <- list(
#     "var1" = cat,
#     "var2" = substitute(cont),
#     "eff_size" <- unname(ifelse(nlevels(cat)==2, (wilcox.test(cont~cat))$statistic, 
#                          (kruskal.test(cont~cat, data=df))$statistic)),
#     "pval" = unname(ifelse(nlevels(cat)==2, (wilcox.test(cont~cat, data=df))$p.value, 
#                     (kruskal.test(cont~cat, data=df))$p.value)))
#   hm_list
# }
# 
# hm_list <- list(
#   "var1" = cat,
#   "var2" = substitute(cont),
#   "eff_size" <- unname(ifelse(nlevels(cat)==2, (wilcox.test(cont~cat))$statistic, 
#                               (kruskal.test(cont~cat, data=df))$statistic)),
#   "pval" = unname(ifelse(nlevels(cat)==2, (wilcox.test(cont~cat, data=df))$p.value, 
#                          (kruskal.test(cont~cat, data=df))$p.value)))
# 

##should functionalize this
Gender_DevConcerns <- table(full_elssp$Gender, full_elssp$DevelopmentalConcerns)
Gender_DevConcerns_list <- prep_hm_df(Gender_DevConcerns)
hm_prep <- rbind(hm_prep, Gender_DevConcerns_list)

i <- sapply(hm_prep, is.factor) #prevent var2 from turning into factor??
hm_prep[i] <- lapply(hm_prep[i], as.character)

HealthIssues_Gender <- table(full_elssp$Gender, full_elssp$HealthIssues)
HealthIssues_Gender_list <- prep_hm_df(HealthIssues_Gender)
hm_prep <- rbind(hm_prep, HealthIssues_Gender_list)

Prematurity_Gender <- table(full_elssp$Gender, full_elssp$IsPremature)
Prematurity_Gender_list <- prep_hm_df(Prematurity_Gender)
hm_prep <- rbind(hm_prep, Prematurity_Gender_list)

Laterality_Gender <- table(full_elssp$Gender, full_elssp$Laterality)
Laterality_Gender <- Laterality_Gender[, c("Bilateral", "Unilateral")]
Laterality_Gender_list <- prep_hm_df(Laterality_Gender)
hm_prep <- rbind(hm_prep, Laterality_Gender_list)


Meets136_Gender <- table(full_elssp$Gender, full_elssp$Meets136)
Meets136_Gender <- Meets136_Gender[, -1]
Meets136_Gender_list <- prep_hm_df(Meets136_Gender)
hm_prep <- rbind(hm_prep, Meets136_Gender_list)


Gender_Amplification <- table(full_elssp$Gender, full_elssp$Amplification)
Gender_Amplification <- Gender_Amplification[, -1]
Gender_Amplification_list <- prep_hm_df(Gender_Amplification)
hm_prep <- rbind(hm_prep, Gender_Amplification_list)


Gender_Etiology <- table(full_elssp$Gender, full_elssp$Etiology)
Gender_Etiology <- Gender_Etiology[, -1]
Gender_Etiology_list <- prep_hm_df(Gender_Etiology)
hm_prep <- rbind(hm_prep, Gender_Etiology_list)


Gender_Communication <- table(full_elssp$Gender, full_elssp$Communication)
Gender_Communication <- Gender_Communication[, -1]
Gender_Communication_list <- prep_hm_df(Gender_Communication)
hm_prep <- rbind(hm_prep, Gender_Communication_list)


Prematurity_HealthIssues <- table(full_elssp$HealthIssues, full_elssp$IsPremature)
Prematurity_HealthIssues_list <- prep_hm_df(Prematurity_HealthIssues)
hm_prep <- rbind(hm_prep, Prematurity_HealthIssues_list)


Laterality_HealthIssues <- table(full_elssp$HealthIssues, full_elssp$Laterality)
Laterality_HealthIssues <- Laterality_HealthIssues[, c("Bilateral", "Unilateral")]
Laterality_HealthIssues_list <- prep_hm_df(Laterality_HealthIssues)
hm_prep <- rbind(hm_prep, Laterality_HealthIssues_list)

Meets136_HealthIssues <- table(full_elssp$HealthIssues, full_elssp$Meets136)
Meets136_HealthIssues <- Meets136_HealthIssues[, -1]
Meets136_HealthIssues_list <- prep_hm_df(Meets136_HealthIssues)
hm_prep <- rbind(hm_prep, Meets136_HealthIssues_list)

HealthIssues_DevConcerns <- table(full_elssp$DevelopmentalConcerns, full_elssp$HealthIssues)
HealthIssues_DevConcerns_list <- prep_hm_df(HealthIssues_DevConcerns)
hm_prep <- rbind(hm_prep, HealthIssues_DevConcerns_list)

HealthIssues_Amplification <- table(full_elssp$HealthIssues, full_elssp$Amplification)
HealthIssues_Amplification <- HealthIssues_Amplification[, -1]
HealthIssues_Amplification_list <- prep_hm_df(HealthIssues_Amplification)
hm_prep <- rbind(hm_prep, HealthIssues_Amplification_list)

HealthIssues_Etiology <- table(full_elssp$HealthIssues, full_elssp$Etiology)
HealthIssues_Etiology <- HealthIssues_Etiology[, -1]
HealthIssues_Etiology_list <- prep_hm_df(HealthIssues_Etiology)
hm_prep <- rbind(hm_prep, HealthIssues_Etiology_list)

HealthIssues_Communication <- table(full_elssp$HealthIssues, full_elssp$Communication)
HealthIssues_Communication <- HealthIssues_Communication[, -1]
HealthIssues_Communication_list <- prep_hm_df(HealthIssues_Communication)
hm_prep <- rbind(hm_prep, HealthIssues_Communication_list)

Prematurity_Laterality <- table(full_elssp$IsPremature, full_elssp$Laterality)
Prematurity_Laterality <- Prematurity_Laterality[, c("Bilateral", "Unilateral")]
Prematurity_Laterality_list <- prep_hm_df(Prematurity_Laterality)
hm_prep <- rbind(hm_prep, Prematurity_Laterality_list)

Prematurity_Meets136 <- table(full_elssp$IsPremature, full_elssp$Meets136)
Prematurity_Meets136 <- Prematurity_Meets136[, -1]
Prematurity_Meets136_list <- prep_hm_df(Prematurity_Meets136)
hm_prep <- rbind(hm_prep, Prematurity_Meets136_list)

Prematurity_DevConcerns <- table(full_elssp$IsPremature, full_elssp$DevelopmentalConcerns)
Prematurity_DevConcerns
Prematurity_DevConcerns_list <- prep_hm_df(Prematurity_DevConcerns)
hm_prep <- rbind(hm_prep, Prematurity_DevConcerns_list)

Prematurity_Amplification <- table(full_elssp$IsPremature, full_elssp$Amplification)
Prematurity_Amplification <- Prematurity_Amplification[, -1]
Prematurity_Amplification_list <- prep_hm_df(Prematurity_Amplification)
hm_prep <- rbind(hm_prep, Prematurity_Amplification_list)

Prematurity_Etiology <- table(full_elssp$IsPremature, full_elssp$Etiology)
Prematurity_Etiology <- Prematurity_Etiology[, -1]
Prematurity_Etiology_list <- prep_hm_df(Prematurity_Etiology)
hm_prep <- rbind(hm_prep, Prematurity_Etiology_list)

Prematurity_Communication <- table(full_elssp$IsPremature, full_elssp$Communication)
Prematurity_Communication <- Prematurity_Communication[, -1]
Prematurity_Communication_list <- prep_hm_df(Prematurity_Communication)
hm_prep <- rbind(hm_prep, Prematurity_Communication_list)

Meets136_Laterality <- table(full_elssp$Meets136, full_elssp$Laterality)
Meets136_Laterality <- Meets136_Laterality[, c("Bilateral", "Unilateral")]
Meets136_Laterality_list <- prep_hm_df(Meets136_Laterality)
hm_prep <- rbind(hm_prep, Meets136_Laterality_list)

Laterality_Amplification <- table(full_elssp$Amplification, full_elssp$Laterality)
Laterality_Amplification <- Laterality_Amplification[-1, c("Bilateral", "Unilateral")]
Laterality_Amplification_list <- prep_hm_df(Laterality_Amplification)
hm_prep <- rbind(hm_prep, Laterality_Amplification_list)

Laterality_DevConcerns <- table(full_elssp$DevelopmentalConcerns, full_elssp$Laterality)
Laterality_DevConcerns <- Laterality_DevConcerns[, c("Bilateral", "Unilateral")]
Laterality_DevConcerns_list <- prep_hm_df(Laterality_DevConcerns)
hm_prep <- rbind(hm_prep, Laterality_DevConcerns_list)

Laterality_Etiology <- table(full_elssp$Etiology, full_elssp$Laterality)
Laterality_Etiology <- Laterality_Etiology[-1, c("Bilateral", "Unilateral")]
Laterality_Etiology_list <- prep_hm_df(Laterality_Etiology)
hm_prep <- rbind(hm_prep, Laterality_Etiology_list)

Laterality_Communication <- table(full_elssp$Communication, full_elssp$Laterality)
Laterality_Communication <- Laterality_Communication[-1, c("Bilateral", "Unilateral")]
Laterality_Communication_list <- prep_hm_df(Laterality_Communication)
hm_prep <- rbind(hm_prep, Laterality_Communication_list)

Meets136_DevConcerns <- table(full_elssp$Meets136, full_elssp$DevelopmentalConcerns)
Meets136_DevConcerns <- Meets136_DevConcerns[-1, ]
Meets136_DevConcerns_list <- prep_hm_df(Meets136_DevConcerns)
hm_prep <- rbind(hm_prep, Meets136_DevConcerns_list)

Meets136_Amplification <- table(full_elssp$Meets136, full_elssp$Amplification)
Meets136_Amplification <- Meets136_Amplification[-1, -1]
Meets136_Amplification_list <- prep_hm_df(Meets136_Amplification)
hm_prep <- rbind(hm_prep, Meets136_Amplification_list)

Meets136_Etiology <- table(full_elssp$Meets136, full_elssp$Etiology)
Meets136_Etiology <- Meets136_Etiology[-1, -1]
Meets136_Etiology_list <- prep_hm_df(Meets136_Etiology)
hm_prep <- rbind(hm_prep, Meets136_Etiology_list)

Meets136_Communication <- table(full_elssp$Meets136, full_elssp$Communication)
Meets136_Communication <- Meets136_Communication[-1, -1]
Meets136_Communication_list <- prep_hm_df(Meets136_Communication)
hm_prep <- rbind(hm_prep, Meets136_Communication_list)

DevConcerns_Amplification <- table(full_elssp$Amplification, full_elssp$DevelopmentalConcerns)
DevConcerns_Amplification <- DevConcerns_Amplification[-1, ]
DevConcerns_Amplification_list <- prep_hm_df(DevConcerns_Amplification)
hm_prep <- rbind(hm_prep, DevConcerns_Amplification_list)

Etiology_Amplification <- table(full_elssp$Amplification, full_elssp$Etiology)
Etiology_Amplification <- Etiology_Amplification[-1, -1]
Etiology_Amplification_list <- prep_hm_df(Etiology_Amplification)
hm_prep <- rbind(hm_prep, Etiology_Amplification_list)

Communication_Amplification <- table(full_elssp$Amplification, full_elssp$Communication)
Communication_Amplification <- Communication_Amplification[-1, -1]
Communication_Amplification_list <- prep_hm_df(Communication_Amplification)
hm_prep <- rbind(hm_prep, Communication_Amplification_list)

Etiology_DevConcerns <- table(full_elssp$DevelopmentalConcerns, full_elssp$Etiology)
Etiology_DevConcerns <- Etiology_DevConcerns[, -1]
Etiology_DevConcerns_list <- prep_hm_df(Etiology_DevConcerns)
hm_prep <- rbind(hm_prep, Etiology_DevConcerns_list)

DevConcerns_Communication <- table(full_elssp$DevelopmentalConcerns, full_elssp$Communication)
DevConcerns_Communication <- DevConcerns_Communication[, -1]
DevConcerns_Communication_list <- prep_hm_df(DevConcerns_Communication)
hm_prep <- rbind(hm_prep, DevConcerns_Communication_list)

Etiology_Communication <- table(full_elssp$Communication, full_elssp$Etiology)
Etiology_Communication <- Etiology_Communication[-1, -1]
Etiology_Communication_list <- prep_hm_df(Etiology_Communication)
hm_prep <- rbind(hm_prep, Etiology_Communication_list)

Language_DevConcerns <- table(full_elssp$PrimaryLanguage, full_elssp$DevelopmentalConcerns)
Language_DevConcerns_list <- prep_hm_df(Language_DevConcerns)
hm_prep <- rbind(hm_prep, Language_DevConcerns_list)

Language_HealthIssues <- table(full_elssp$PrimaryLanguage, full_elssp$HealthIssues)
Language_HealthIssues_list <- prep_hm_df(Language_HealthIssues)
hm_prep <- rbind(hm_prep, Language_HealthIssues_list)

Prematurity_Language <- table(full_elssp$PrimaryLanguage, full_elssp$IsPremature)
Prematurity_Language_list <- prep_hm_df(Prematurity_Language)
hm_prep <- rbind(hm_prep, Prematurity_Language_list)

Laterality_Language <- table(full_elssp$PrimaryLanguage, full_elssp$Laterality)
Laterality_Language <- Laterality_Language[, c("Bilateral", "Unilateral")]
Laterality_Language_list <- prep_hm_df(Laterality_Language)
hm_prep <- rbind(hm_prep, Laterality_Language_list)


Meets136_Language <- table(full_elssp$PrimaryLanguage, full_elssp$Meets136)
Meets136_Language <- Meets136_Language[, -1]
Meets136_Language_list <- prep_hm_df(Meets136_Language)
hm_prep <- rbind(hm_prep, Meets136_Language_list)


Language_Amplification <- table(full_elssp$PrimaryLanguage, full_elssp$Amplification)
Language_Amplification <- Language_Amplification[, -1]
Language_Amplification_list <- prep_hm_df(Language_Amplification)
hm_prep <- rbind(hm_prep, Language_Amplification_list)


Language_Etiology <- table(full_elssp$PrimaryLanguage, full_elssp$Etiology)
Language_Etiology <- Language_Etiology[, -1]
Language_Etiology_list <- prep_hm_df(Language_Etiology)
hm_prep <- rbind(hm_prep, Language_Etiology_list)


Language_Communication <- table(full_elssp$PrimaryLanguage, full_elssp$Communication)
Language_Communication <- Language_Communication[, -1]
Language_Communication_list <- prep_hm_df(Language_Communication)
hm_prep <- rbind(hm_prep, Language_Communication_list)

Language_Gender <- table(full_elssp$PrimaryLanguage, full_elssp$Gender)
Language_Gender_list <- prep_hm_df(Language_Gender)
hm_prep <- rbind(hm_prep, Language_Gender_list)





Degree_Amplification <- kruskal.test(HLworse ~ Amplification,
                                     data=(full_elssp %>% filter(full_elssp$Amplification=="CI"|full_elssp$Amplification=="HA"|full_elssp$Amplification=="none")))
Degree_Amplification_list <- list(
  "var1" = "Degree",
  "var2" = "Amplification",
  "eff_size" = unname(Degree_Amplification$statistic),
  "pval" = Degree_Amplification$p.value
)
hm_prep <- rbind(hm_prep, Degree_Amplification_list)

Services_DevConcerns <- wilcox.test(ServicesReceivedPerMonth ~ DevelopmentalConcerns,
                                 data=(full_elssp %>% filter(full_elssp$DevelopmentalConcerns==1|full_elssp$DevelopmentalConcerns==0)))
# Services_DevConcerns_list <- list(
#   "var1" = "Services",
#   "var2" = "DevConcerns",
#   "eff_size" = unname(Services_DevConcerns$statistic),
#   "pval" = Services_DevConcerns$p.value
# )

hm_prep <- hm_prep %>% mutate(sig = ifelse(pval>.05, "ns", 
                                           ifelse(pval>0.0007575758, "sig", "
                                                  survivesbc")))


Services_Communication <- kruskal.test(ServicesReceivedPerMonth ~ Communication,
                                       data=(full_elssp %>% filter(full_elssp$Communication=="total communication"|full_elssp$Communication=="spoken")))

hm_plot <- ggplot(hm_prep, aes(y=var1, x=var2)) +
  geom_tile(aes(fill = sig)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = "RdPu", direction=-1) +
  xlab(NULL) + ylab(NULL) +
  geom_text(aes(label=ifelse(sig=='survives_bc',eff_size,'')))
hm_plot

# lm_output <- function(data) {
#   print_res = paste("ß=", round(data$Estimate,2),
#                       ", SE=", round(data$Est.Error,2),
#                       ", 95% CI [", printnum(round(data$l.95..CI,2), format = "f"),
#                                   ",", round(data$u.95..CI,2),"]", sep='')
#   print_res
# }

chisq_output <- function(data) {
  chisq_output = paste("($X^2$ (", round((chisq.test(data))$parameter,2),
                    ", N = ", sum(((chisq.test(data))$observed)),
                    ") = ", round(chisq.test(data)$statistic, 2),
                    ", p = ", format.pval(chisq.test(data)$p.value, digits = 2),
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
                                                       "Etiology", "Laterality", "Amplification",
                                                       "Communication", "PrimaryLanguage", "DevelopmentalConcerns")) %>%
  dplyr::select(ServicesReceivedPerMonth, HLworse, Gender_male,
                Amplification_CI, Amplification_HA, Amplification_none,
                PrimaryLanguage_Spanish, PrimaryLanguage_English, PrimaryLanguage_Hindi,
                Etiology_SNHL, Etiology_Mixed, Etiology_Conductive,
                Communication_spoken, DevelopmentalConcerns, HealthIssues, IsPremature,
                VisionLoss) %>%
  mutate()
elssp_corr <- cor(corr_prep, use="pairwise.complete.obs")

comorbid <- read.csv("data/elssp_comorbidities.csv") %>% mutate(anycomorbid = ifelse(VisionLoss==1|
                                                                                       DevelopmentalConcerns==1|
                                                                                       HealthIssues==1|
                                                                                       IsPremature==1, 1,
                                                                                     0), 
                                                                extremelypremature = ifelse(WeeksGestation > 33, 0,
                                                                                            1))
comorbid_long <- comorbid %>% pivot_longer(cols = c(ANSD, IsPremature:extremelypremature, -WeeksGestation), 
                                           names_to = "condition", 
                                           values_to = "n") %>% 
  mutate(condition = as.factor(condition))
condition_tallies <- aggregate(n ~ condition, data=comorbid_long, FUN = sum) %>% remove_rownames %>% column_to_rownames(var="condition")


#make.names(condition_tallies$condition)

n_condition <- function(condition){
  as.numeric(condition_tallies[paste(enexpr(condition)),"n"])
}

hearing_means <- function(col, uni_bi, amp){ifelse(is.null(amp),
  mean((elssp %>% filter(Laterality==enexpr(uni_bi)))$col, na.rm = TRUE),
  mean((elssp %>% filter(Laterality==enexpr(uni_bi) & Amplification==enexpr(amp)))$col)
)
}
