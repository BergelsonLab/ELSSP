library('wordbankr')
library('plyr')
library('dplyr')
library('ggplot2')
library('reshape2')
library('dotwhisker')
library('gplots')
library('MASS')
source('CDI_ELSSP.R')
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


#read in data
elssp = read.csv("data/ELSSP_SubjectInfo_02042020.csv", stringsAsFactors=F) %>% mutate(InSample=as.factor(InSample)) %>% filter(VisitNumber==1)	
elssp$AgeAtEvaluationMonths = elssp$Age #for clarity
elssp$subject_id = unlist(lapply(strsplit(elssp$SubjectNumber,'_'), function(x){as.numeric(x[2])}))
elssp$admin_id = elssp$SubjectNumber
summary(elssp)

#add months-delay to dataframes
source('CDI_ELSSP.R')
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
                                                    InSample=as.factor(InSample)) %>% filter(VisitNumber==1)


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
                                                    InSample=as.factor(InSample)) %>% filter(VisitNumber==1)

full_elssp <- rbind(wg_elssp, ws_elssp)

#rainbow cdi plots
plot_elssp_df(elssp_datasets[[1]])
plot_elssp_df(elssp_datasets[[2]])

#categorical contasts

##gender
wilcox.test(diff_age_from_expected ~ Gender, data=wg_elssp)
wg_gender <- plot_elssp_df(elssp_datasets[[1]], 'Gender')

wilcox.test(diff_age_from_expected ~ Gender, data=ws_elssp)
ws_gender <- plot_elssp_df(elssp_datasets[[2]], 'Gender')

##laterality
wilcox.test(diff_age_from_expected ~ Laterality, data=wg_elssp)
wg_laterality <- plot_elssp_df(elssp_datasets[[1]], 'Laterality')

wilcox.test(diff_age_from_expected ~ Laterality, data=ws_elssp)
ws_laterality <- plot_elssp_df(elssp_datasets[[2]], 'Laterality')

#amplification
kruskal.test(diff_age_from_expected ~ Amplification, data=wg_elssp)
wg_amplification <- plot_elssp_df(elssp_datasets[[1]], 'Amplification')

kruskal.test(diff_age_from_expected ~ Laterality, data=ws_elssp)
ws_amplification <- plot_elssp_df(elssp_datasets[[2]], 'Amplification')

#healthissues
wilcox.test(diff_age_from_expected ~ HealthIssues, data=wg_elssp)
wg_healthissues <- plot_elssp_df(elssp_datasets[[1]], 'HealthIssues')

wilcox.test(diff_age_from_expected ~ HealthIssues, data=ws_elssp)
ws_healthissues <- plot_elssp_df(elssp_datasets[[2]], 'HealthIssues')

#devconcerns
wilcox.test(diff_age_from_expected ~ DevelopmentalConcerns, data=wg_elssp)
wg_devconcerns <- plot_elssp_df(elssp_datasets[[1]], 'DevelopmentalConcerns')

wilcox.test(diff_age_from_expected ~ DevelopmentalConcerns, data=ws_elssp)
ws_devconcerns <- plot_elssp_df(elssp_datasets[[2]], 'DevelopmentalConcerns')

#prematurity
wilcox.test(diff_age_from_expected ~ IsPremature, data=wg_elssp)
wg_prematurity <- plot_elssp_df(elssp_datasets[[1]], 'IsPremature')

wilcox.test(diff_age_from_expected ~ IsPremature, data=ws_elssp)
ws_prematurity <- plot_elssp_df(elssp_datasets[[2]], 'IsPremature')

#meets136
wilcox.test(diff_age_from_expected ~ Meets136, data=wg_elssp)
wg_meets136 <- plot_elssp_df(elssp_datasets[[1]], 'Meets136')

wilcox.test(diff_age_from_expected ~ Meets136, data=ws_elssp)
ws_meets136 <- plot_elssp_df(elssp_datasets[[2]], 'Meets136')

#communication
wilcox.test(diff_age_from_expected ~ Communication, data=wg_elssp%>%
              filter(Communication=="spoken"|Communication=="total communication"))
wg_communication <- plot_elssp_df(elssp_datasets[[1]], 'Communication')

wilcox.test(diff_age_from_expected ~ Communication, data=ws_elssp%>%
              filter(Communication=="spoken"|Communication=="total communication"))
ws_communication <- plot_elssp_df(elssp_datasets[[2]], 'Communication')

#etiology
kruskal.test(diff_age_from_expected ~ Etiology, data=wg_elssp)
wg_etiology <- plot_elssp_df(elssp_datasets[[1]], 'Etiology')

kruskal.test(diff_age_from_expected ~ Etiology, data=ws_elssp)
ws_etiology <- plot_elssp_df(elssp_datasets[[2]], 'Etiology')

#degree
cor.test(wg_elssp$HLworse, y=wg_elssp$diff_age_from_expected, method="kendall")
cor.test(ws_elssp$HLworse, y=ws_elssp$diff_age_from_expected, method="kendall")

wg_degree <- ggplot(data=elssp_datasets[[1]]$elssp_df, aes(x=HLworse, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')
ws_degree <- ggplot(data=elssp_datasets[[2]]$elssp_df, aes(x=HLworse, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')

#services/month
cor.test(wg_elssp$ServicesReceivedPerMonth, y=wg_elssp$diff_age_from_expected, method="kendall")
cor.test(ws_elssp$ServicesReceivedPerMonth, y=ws_elssp$diff_age_from_expected, method="kendall")

wg_services <- ggplot(data=elssp_datasets[[1]]$elssp_df, aes(x=ServicesReceivedPerMonth, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')
ws_services <- ggplot(data=elssp_datasets[[2]]$elssp_df, aes(x=ServicesReceivedPerMonth, y=diff_age_from_expected))+
  geom_point() + geom_smooth(method='lm')

#meets136 predictors
##diagnosis
diagnosis_full <- lm(IdentificationOfHLMonths ~ Gender + Etiology +HLbetter + 
                       HLworse + Laterality + IsPremature + HealthIssues + 
                       DevelopmentalConcerns + Monolingual_English + VisionLoss + ANSD, (data=full_elssp %>% 
                                                                                           filter(DevelopmentalConcerns=="yes"|DevelopmentalConcerns=="no") 
                                                                                         %>% filter(HLworse!="NA")))
summary(diagnosis_full)
stepAIC(diagnosis_full)
diagnosis_aic <- lm(formula = IdentificationOfHLMonths ~ HealthIssues + Monolingual_English, 
                    data = full_elssp)
summary(diagnosis_aic)
dw_diagnosis = dwplot(diagnosis_aic)
dw_diagnosis = dw_diagnosis + geom_vline(xintercept=0, linetype = 'dashed') + theme_classic(
)+ theme(legend.position="none")
options(repr.plot.width=6, repr.plot.height=4)
print(dw_diagnosis)

##services
services_full <- lm(AgeStartedServices ~ Gender + Etiology +HLbetter + 
                      HLworse + Laterality + IsPremature + HealthIssues + 
                      DevelopmentalConcerns + Monolingual_English + Communication+
                      VisionLoss + ANSD + IdentificationOfHLMonths, (data=full_elssp %>% 
                                                                       filter(HLworse!="NA")))
summary(services_full)
stepAIC(services_full)
services_aic <- lm(formula = AgeStartedServices ~ HLworse + Laterality + IsPremature + 
                       IdentificationOfHLMonths, data = (data = full_elssp %>% filter(HLworse != 
                                                                                        "NA")))
summary(services_aic)
dw_services = dwplot(services_aic)
dw_services = dw_services + geom_vline(xintercept=0, linetype = 'dashed') + theme_classic(
)+ theme(legend.position="none")
options(repr.plot.width=6, repr.plot.height=4)
print(dw_services)

##event plot
event <- c('Diagnosis','Intervention','Amplification', 'Implantation')
eventmin <- c(min(full_elssp$IdentificationOfHLMonths, na.rm = TRUE), min(full_elssp$AgeStartedServices, na.rm = TRUE), 
              min(full_elssp$AgeAmplification, na.rm = TRUE), min(full_elssp$AgeImplantation, na.rm = TRUE))
eventmax <- c(max(full_elssp$IdentificationOfHLMonths, na.rm = TRUE), max(full_elssp$AgeStartedServices, na.rm = TRUE), 
              max(full_elssp$AgeAmplification, na.rm = TRUE), max(full_elssp$AgeImplantation, na.rm = TRUE))
eventmean <- c(mean(full_elssp$IdentificationOfHLMonths, na.rm = TRUE), mean(full_elssp$AgeStartedServices, na.rm = TRUE), 
               mean(full_elssp$AgeAmplification, na.rm = TRUE), mean(full_elssp$AgeImplantation, na.rm = TRUE))
eventsd <- c(sd(full_elssp$IdentificationOfHLMonths, na.rm = TRUE), sd(full_elssp$AgeStartedServices, na.rm = TRUE), 
             sd(full_elssp$AgeAmplification, na.rm = TRUE), sd(full_elssp$AgeImplantation, na.rm = TRUE))
event_data <- data.frame(event, eventmin, eventmax, eventmean, eventsd)

summary(event_data)

eventplot <- ggplot(data = event_data, 
                    mapping = aes(x = reorder(event, eventmean), 
                                  y = eventmean, color=event)) + 
  geom_pointrange(mapping = aes(ymin = eventmean - eventsd,
                                ymax = eventmean + eventsd)) +
  geom_pointrange(mapping = aes(ymin = eventmin,
                                ymax = eventmax), width=2, alpha=0.5) +
  labs(x= "", y= "Age Months") + 
  coord_flip() +
  theme(legend.position = "none") +
  theme_classic()

print(eventplot)

#interactions among variables
##should functionalize this


##test plots
hm_test <- read.csv("data/heat_map_test.csv")
hm_test_plot <- ggplot(hm_test, aes(x=var1, y=var2)) +
  geom_tile(aes(fill = p)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_distiller(palette = "RdPu") +
  xlab(NULL) + ylab(NULL) 

hm_test_plot2 <- ggplot(hm_test, aes(x=var1, y=var2)) +
  geom_tile(aes(fill = sig)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab(NULL) + ylab(NULL) 