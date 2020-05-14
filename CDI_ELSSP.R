inv_logit <- function(x){
	return(1 / (exp(-x) + 1))
} 

getScoreForAge = function(lm, age, num_items){
# just predict
  prop = inv_logit(predict.glm(lm, newdata = data.frame(age=age)))
  return(prop*num_items)
}

getAgeForScore = function(lm, score, num_items){
  proportion = score / num_items
# http://www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
  b0 = lm$coefficients[1]
  b1 = lm$coefficients[2]
  predicted_age = (log(proportion / (1-proportion)) - b0)/ b1
  return(predicted_age)
}

constants = list()
constants[['WG']] = list()
constants[['WG']]$lowest_num_id = 33
constants[['WG']]$highest_num_id = 430
constants[['WG']]$num_items = constants[['WG']]$highest_num_id - constants[['WG']]$lowest_num_id + 1

constants[['WS']] = list()
constants[['WS']]$lowest_num_id = 1
constants[['WS']]$highest_num_id = 680
constants[['WS']]$num_items = constants[['WS']]$highest_num_id - constants[['WS']]$lowest_num_id + 1



prepare_elssp_df = function(cdi_form, constants, verbose=F){
	print(paste('Processing ',cdi_form,'...', sep=''))	

	num_items = (constants[[cdi_form]][['highest_num_id']] - constants[[cdi_form]][['lowest_num_id']]) + 1	
	print('Number of items:')
	print(num_items)
	
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
	print(counts)
	
	print('Fitting model...')
	model <- stats::glm(cbind(num_true, num_false) ~ age, counts,
		family = "binomial")
		
	if (verbose){
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
	elssp_for_form$SubjectNumber = as.factor(elssp_for_form$SubjectNumber)
	elssp_for_form$DevelopmentalConcerns[elssp_for_form$DevelopmentalConcerns == 1] = 'yes'
	elssp_for_form$DevelopmentalConcerns[elssp_for_form$DevelopmentalConcerns == 0] = 'no'
	elssp_for_form$DevelopmentalConcerns = as.factor(elssp_for_form$DevelopmentalConcerns)
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
		print(elssp_for_form[,c('SubjectNumber','AgeAtEvaluationMonths', 'ProductionCDI', 'expected_score_at_chron_age', 'diff_age_from_expected','diff_score_from_expected')])
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

	if (!is.null(split)){
		elssp_df$split  = as.factor(elssp_df[[split]]) # make the split var in the df on the fly		
		print(unique(elssp_df$split))
		p1 = ggplot()  + geom_line(data=wordbank_norms_melted, aes(x=age, y=value,  group=variable), linetype= 'dashed', alpha= .5) + geom_line(data = samples_from_growth_curve_model, aes(x=predict_ages, y=scores), 
		colour='black', alpha= .5) + theme_classic() 
		p1 = p1 + ylab(paste0("CDI Score (",cdi_form,")"))
		p1 = p1 + xlab("Age in Months") 
		p1 = p1 + coord_cartesian(xlim=c(0,40)) + geom_point(data = subset(elssp_df, !is.na(ProductionCDI)), aes(x= AgeAtEvaluationMonths, y= ProductionCDI, shape=split, color=split)) +labs(colour=split, shape=split)

		p2 = make_rainplot(elssp_df, split, 'diff_score_from_expected', 'CDI Score Deficit')
		p2 = p2 + ggtitle(paste0(cdi_form,': ', split)) + coord_cartesian(ylim = c(-2, 25)) + coord_flip()
		p3 = make_rainplot(elssp_df, split, 'diff_age_from_expected', 'Delay in Months') + coord_cartesian(ylim = c(-2, 25)) + coord_flip()
		p3 = p3 + ggtitle(paste0(cdi_form,': ', split))

		# Jupyter preview
		#options(repr.plot.width=4, repr.plot.height=3)
		#print(p1)
		#print(p2)
		#print(p3)
        
        
       
		# Save to disk
		p1 = p1 + theme_classic(base_size=16)
		ggsave(plot=p1, paste('figures/', elssp_dataset$cdi_form,'_', split, '_trajectory.pdf', sep=''), width=6, height=4)
		p2 = p2 + theme_classic(base_size=16)
		ggsave(plot=p2, paste('figures/', elssp_dataset$cdi_form,'_', split, '_deficit.pdf', sep=''), width=6, height=4)
		p3 = p3 + theme_classic(base_size=16)
		ggsave(plot=p3, paste('figures/', elssp_dataset$cdi_form,'_', split, '_delay.pdf', sep=''), width=6, height=4)

        graph <- p3
        return(graph)

	} else {

		p1 = ggplot() + geom_point (data=wordbank_norms_melted, aes(x=age, y=value, colour=variable)) + 
		theme_classic() + geom_line(data=wordbank_norms_melted, aes(x=age, y=value, 
		colour=variable)) + geom_line(data = samples_from_growth_curve_model, aes(x=predict_ages, y=scores), 
		colour='black') + geom_point(data = subset(elssp_df, !is.na(ProductionCDI)), 
		                             aes(x= AgeAtEvaluationMonths, y= ProductionCDI), shape=17
		) + ggtitle(paste("Normative", cdi_form, "vs. ELSSP"))	
		
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

raincloud_theme <- theme(
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


cv_model_search = function(response_vars, predictor_vars, covar, df, parallel=T, mixed_effects=F){
	# enumerate the set of all predictor_vars
	
	# generate the powerset of model predictors
	predictor_sets = all.subsets(predictor_vars)
	if (length(predictor_sets) > 1024){
		print('Requested evaluation of too many models... reduce the set of predictor vars')
	} 

	# parallel call to cv_model_eval
	if (parallel){
		print(paste('Requested parallel evaluation of', length(predictor_sets),'models...'))
		model_results = mclapply(predictor_sets, function(predictor_set){
			return(cv_model_eval(response_vars, predictor_set, covar, df, mixed_effects))
		}, mc.cores= detectCores())
	} else {
		print(paste('Requested serial evaluation of', length(predictor_sets),'models...'))
		model_results = lapply(predictor_sets, function(predictor_set){
			return(cv_model_eval(response_vars, predictor_set, covar, df, mixed_effects))
		})
	}

	# make a df of the scores
	scores = do.call('rbind', lapply(model_results, function(x){
		x$score
	}))	
	scores$index = 1:nrow(scores)
	scores = scores[order(scores$aic),]	
	
	#return both the scores and the models
	rlist = list()
	rlist[['model_scores']] = scores
	rlist[['models']] = lapply(model_results, function(x){
		x$model
	})

	if (nrow(rlist[['model_scores']]) != length(rlist[['models']])){
		stop('Mismatch in the number of model scores and models')
	}
	
	return(rlist)
}



cv_model_eval = function(response_vars, predictor_vars, covar, df, mixed_effects=F){
	
	# generate a model name
	model_name = paste(predictor_vars, collapse=' + ') 

	# build LHS, which can be reused for the null model
	if (length(response_vars) == 2){
		LHS = paste("cbind(",response_vars[1],",",response_vars[2],") ~ ", sep="")
	} else if (length(response_vars) == 1){
		LHS = paste(response_vars," ~ ")
	} else {
		stop('cv_model_search only handles 1 or 2 reponse_vars (t/f obs)')
	}

	# build RHS for the real model
	if (mixed_effects){
		RHS = paste("(", covar, " | SubjectNumber) + ", covar, "* (", paste(predictor_vars, collapse=" + "), ")")
	} else {
		RHS = paste(covar, "* (", paste(predictor_vars, collapse=" + "), ")")
	}
	
	# fit it
	model_eq = paste(LHS, RHS, sep='')
	print(paste('Model eq: ', model_eq))
	
	if (mixed_effects){
		fitted_model = lme4::glmer(as.formula(model_eq), df, family = "binomial")
	} else {
		fitted_model = stats::glm(as.formula(model_eq), df, family = "binomial")
	}
	#print('Finished fitting!')
	#print(fitted_model)

	# get McFadden's pseudo R-squared
	# fit the null model (intercept_only)
	if (mixed_effects){
		null_model_eq = paste(LHS, "(1 | SubjectNumber)")	
		null_model = lme4::glmer(as.formula(null_model_eq), df, family = "binomial")
	} else {
		null_model_eq = paste(LHS, "1")	
		null_model = stats::glm(as.formula(null_model_eq), df, family = "binomial")
	}
	print(paste('null model eq: ', null_model_eq))
	

	pseudoR2 = 1 - (logLik(fitted_model) / logLik(null_model))

	# do H1O CV and get 1 - delta
	df_clean = dropIfNA(subset(df, type=='real'), c(predictor_vars, response_vars, covar))
	if (mixed_effects){
		stop("Not implemented. Need to compute by hand with Mika's approach")

	} else {
		cv_scores = cv.glm(df_clean, fitted_model)
		mspe = cv_scores$delta[1]
		mspe_adjusted = cv_scores$delta[2]	
	}	
	
	aic = AIC(fitted_model)
	logLik = logLik(fitted_model)

	# return parameter settings plus the above two measures
	rlist = list()
	rlist[['score']] = data.frame(model_name=RHS, mspe, mspe_adjusted, pseudoR2, aic)
	rlist[['model']] = fitted_model
	return(rlist)
}


predictFromModel = function(model, participant_properties, num_items, range='full'){
    # build a dataframe with predictions over the range of ages
    if (range == 'full'){
    	AgeAtEvaluationMonths = seq(from=8, to=36, by=.1)
    } else if (range == 'last'){
    	AgeAtEvaluationMonths = c(36)
    }	
    rdf = data.frame(AgeAtEvaluationMonths)
    for (property in names(participant_properties)){
        rdf[[property]] = participant_properties[[property]]
    }
    
    predicted_probs = data.frame(predict(model, rdf, type='response', se.fit=T))
    if (range == 'full'){    	
	    predicted_probs$se_high = predicted_probs$fit + predicted_probs$se.fit
	    predicted_probs$se_low = predicted_probs$fit - predicted_probs$se.fit
	    rdf$score = predicted_probs$fit * num_items
	    rdf$se_high = predicted_probs$se_high * num_items
	    rdf$se_low = predicted_probs$se_low * num_items
	    return(rdf)
	} else if (range == 'last'){
		score = predicted_probs$fit * num_items		
    	return(score)	
    }

}


get_all_obs_levels = function(model, df, covar){
	# enumerate all combinations of factor levels, 
	# and max and min for each continuous variable
	var_inventory = list()	
	for (column in names(df)){
		if (column %in% colnames(attr(model$terms, 'factors')) & 
			column != covar
			){				
			if (is.factor(df[[column]]) || is.character(df[[column]])){
				var_inventory[[column]] = levels(as.factor(df[[column]]))
			} else if (is.numeric(df[[column]])){
				var_inventory[[column]] = c(min(df[[column]], na.rm=T), max(df[[column]], na.rm=T))
			}	
		}
	} 

	print(var_inventory)
	all_obs = do.call(expand.grid, var_inventory)
	return(all_obs)
}

# auto-discover the best and worst trajectories
get_trajectory_plot = function(model, df, instrument, covar=c()){
	
	all_obs = get_all_obs_levels(model, df, covar)
	
	all_obs$score = sapply(1:nrow(all_obs), function(i){
		return(predictFromModel(model, all_obs[i,], constants[[instrument]]$num_items,
		 range='last'))		
	})
	all_obs = all_obs[order(all_obs$score, decreasing=T),]
	
	print(all_obs)

	best_performing_properties =  all_obs[1,]
	worst_performing_properties =  all_obs[nrow(all_obs),]


	best_curve = predictFromModel(model, best_performing_properties, 
		constants[[instrument]]$num_items, range='full')
	worst_curve = predictFromModel(model, worst_performing_properties, 
		constants[[instrument]]$num_items, range='full')

	seq_and_labels = seq(from=0,to=48,by=6)
	options(repr.plot.idth=6, repr.plot.height=4)
	
	p1 = ggplot(worst_curve) + geom_line(aes(x=AgeAtEvaluationMonths, y=score), color='red') + geom_errorbar(aes(x=AgeAtEvaluationMonths, ymin=se_low,
ymax = se_high), alpha =.25, color='red') + geom_line(data= best_curve, 
aes(x=AgeAtEvaluationMonths, y=score), color='green') + geom_errorbar(data= best_curve, aes(x=AgeAtEvaluationMonths, ymin=se_low,
ymax = se_high), alpha =.25, color='green') + theme_classic() + ylab(paste('Score on', instrument)) + xlab('Age at Evaluation') + geom_hline(yintercept=constants[[instrument]]$num_items, color='black',
    alpha=.2)


	if (instrument  == 'WG'){
		p1 =  p1 + geom_errorbar(data= subset(
elssp_datasets[[1]]$samples_from_growth_curve_model, predict_ages < 36 & predict_ages > 8),
aes(x=predict_ages, ymin=se_low,
ymax = se_high), alpha =.75, color='black') + theme_classic() + ylab('Score on WG') + xlab('Age at Evaluation') + scale_x_continuous(
breaks = seq_and_labels, labels = seq_and_labels)
	} else if (instrument == 'WS'){
		p1 =  p1 + geom_errorbar(data= subset(
elssp_datasets[[2]]$samples_from_growth_curve_model, predict_ages < 36 & predict_ages > 8),
aes(x=predict_ages, ymin=se_low,
ymax = se_high), alpha =.75, color='black') + theme_classic() + ylab('Score on WS') + xlab('Age at Evaluation') + scale_x_continuous(breaks = seq_and_labels, labels = seq_and_labels) 
	}

	rlist = list()
	rlist[['plot']] = p1
	rlist[['best_performing_properties']] = best_performing_properties
	rlist[['worst_performing_properties']] = worst_performing_properties
	return(rlist)
}


augment_data = function(df, instrument){

	# get the first instance of each participant
	print(paste('Starting with ', nrow(df),' observations...'))
	df = df[order(df$VisitNumber),]
	df$type = 'real'
	df$ProductionCDI_no = constants[[instrument]]$num_items - df$ProductionCDI

	df_synthetic = df[!duplicated(df$SubjectNumber),] 	
	# set the months to 6 
	df_synthetic$AgeAtEvaluationMonths = 8
	df_synthetic$ProductionCDI = 0 
	df_synthetic$ProductionCDI_no = constants[[instrument]]$num_items
	df_synthetic$type = 'synthetic'
	# set the score to 0

	# bind the observations	
	augmented_df = rbind.fill(df, df_synthetic)
	print(paste('Now ', nrow(augmented_df),' observations!'))
	return(augmented_df)
}	