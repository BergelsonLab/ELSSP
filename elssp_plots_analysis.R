library('wordbankr')
library('plyr')
library('dplyr')
library('ggplot2')
library('reshape2')
library('dotwhisker')
library('gplots')
library('MASS')
source('../CDI_ELSSP.R')
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


#read in data
elssp = read.csv("../data/ELSSP_SubjectInfo_02042020.csv", stringsAsFactors=F) %>% mutate(InSample=as.factor(InSample)) %>% filter(VisitNumber==1)	
elssp$AgeAtEvaluationMonths = elssp$Age #for clarity
elssp$subject_id = unlist(lapply(strsplit(elssp$SubjectNumber,'_'), function(x){as.numeric(x[2])}))
elssp$admin_id = elssp$SubjectNumber
summary(elssp)

#add months-delay to dataframes
source('../CDI_ELSSP.R')
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
wg_rainbow <- plot_elssp_df(elssp_datasets[[1]])
ws_rainbow <- plot_elssp_df(elssp_datasets[[2]])

#categorical contasts

##gender
wilcox.test(diff_age_from_expected ~ Gender, data=wg_elssp)
wg_gender <- plot_elssp_df(elssp_datasets[[1]], 'Gender')

wilcox.test(diff_age_from_expected ~ Gender, data=ws_elssp)
ws_gender <- plot_elssp_df(elssp_datasets[[2]], 'Gender')

##laterality
wilcox.test(diff_age_from_expected ~ Laterality, data=(wg_elssp %>% 
                                                         dplyr::filter(Laterality=='Unilateral'|Laterality=='Bilateral')))
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
wilcox.test(diff_age_from_expected ~ Meets136, data=(wg_elssp %>% filter(Meets136=='yes'|Meets136=='no')))
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
Gender_DevConcerns <- table(full_elssp$Gender, full_elssp$DevelopmentalConcerns)
sum((chisq.test(Gender_DevConcerns))$observed)

Gender_HealthIssues <- table(full_elssp$Gender, full_elssp$HealthIssues)
chisq.test(Gender_HealthIssues)

Gender_Premature <- table(full_elssp$Gender, full_elssp$IsPremature)
chisq.test(Gender_Premature)

Gender_Laterality <- table(full_elssp$Gender, full_elssp$Laterality)
Gender_Laterality <- Gender_Laterality[, c("Bilateral", "Unilateral")]
chisq.test(Gender_Laterality)

Gender_Meets136 <- table(full_elssp$Gender, full_elssp$Meets136)
Gender_Meets136 <- Gender_Meets136[, -1]
Gender_Meets136
chisq.test(Gender_Meets136)

Gender_Amp <- table(full_elssp$Gender, full_elssp$Amplification)
Gender_Amp <- Gender_Amp[, -1]
chisq.test(Gender_Amp)

Gender_Etiology <- table(full_elssp$Gender, full_elssp$Etiology)
Gender_Etiology <- Gender_Etiology[, -1]
Gender_Etiology
chisq.test(Gender_Etiology)

Gender_Communication <- table(full_elssp$Gender, full_elssp$Communication)
Gender_Communication <- Gender_Communication[, -1]
Gender_Communication
chisq.test(Gender_Communication)

HealthIssues_Premature <- table(full_elssp$HealthIssues, full_elssp$IsPremature)
chisq.test(HealthIssues_Premature)

HealthIssues_Laterality <- table(full_elssp$HealthIssues, full_elssp$Laterality)
HealthIssues_Laterality <- HealthIssues_Laterality[, c("Bilateral", "Unilateral")]
chisq.test(HealthIssues_Laterality)

HealthIssues_Meets136 <- table(full_elssp$HealthIssues, full_elssp$Meets136)
HealthIssues_Meets136 <- HealthIssues_Meets136[, -1]
chisq.test(HealthIssues_Meets136)

HealthIssues_DevDelay <- table(full_elssp$DevelopmentalConcerns, full_elssp$IsPremature)
chisq.test(HealthIssues_DevDelay)

HealthIssues_Amplification <- table(full_elssp$HealthIssues, full_elssp$Amplification)
HealthIssues_Amplification <- HealthIssues_Amplification[, -1]
chisq.test(HealthIssues_Amplification)

HealthIssues_Etiology <- table(full_elssp$HealthIssues, full_elssp$Etiology)
HealthIssues_Etiology <- HealthIssues_Etiology[, -1]
chisq.test(HealthIssues_Etiology)
#balloonplot(HealthIssues_Etiology)

HealthIssues_Communication <- table(full_elssp$HealthIssues, full_elssp$Communication)
HealthIssues_Communication <- HealthIssues_Communication[, -1]
chisq.test(HealthIssues_Communication)
#balloonplot(HealthIssues_Communication)

Prematurity_Laterality <- table(full_elssp$IsPremature, full_elssp$Laterality)
Prematurity_Laterality <- Prematurity_Laterality[, c("Bilateral", "Unilateral")]
chisq.test(Prematurity_Laterality)
#balloonplot(Prematurity_Laterality)

Prematurity_Meets136 <- table(full_elssp$IsPremature, full_elssp$Meets136)
Prematurity_Meets136 <- Prematurity_Meets136[, -1]
chisq.test(Prematurity_Meets136)
#balloonplot(Prematurity_Meets136)

Prematurity_DevConcerns <- table(full_elssp$IsPremature, full_elssp$DevelopmentalConcerns)
Prematurity_DevConcerns
chisq.test(Prematurity_DevConcerns)
#balloonplot(Prematurity_DevConcerns)

Prematurity_Amplification <- table(full_elssp$IsPremature, full_elssp$Amplification)
Prematurity_Amplification <- Prematurity_Amplification[, -1]
chisq.test(Prematurity_Amplification)
#balloonplot(Prematurity_Amplification)

Prematurity_Etiology <- table(full_elssp$IsPremature, full_elssp$Etiology)
Prematurity_Etiology <- Prematurity_Etiology[, -1]
chisq.test(Prematurity_Etiology)
#balloonplot(Prematurity_Etiology)

Prematurity_Communication <- table(full_elssp$IsPremature, full_elssp$Communication)
Prematurity_Communication <- Prematurity_Communication[, -1]
chisq.test(Prematurity_Communication)
#balloonplot(Prematurity_Communication)

Meets136_Laterality <- table(full_elssp$Meets136, full_elssp$Laterality)
Meets136_Laterality <- Meets136_Laterality[, c("Bilateral", "Unilateral")]
chisq.test(Meets136_Laterality)
#balloonplot(Meets136_Laterality)

Amplification_Laterality <- table(full_elssp$Amplification, full_elssp$Laterality)
Amplification_Laterality <- Amplification_Laterality[-1, c("Bilateral", "Unilateral")]
chisq.test(Amplification_Laterality)
#balloonplot(Amplification_Laterality)

DevConcerns_Laterality <- table(full_elssp$DevelopmentalConcerns, full_elssp$Laterality)
DevConcerns_Laterality <- DevConcerns_Laterality[, c("Bilateral", "Unilateral")]
chisq.test(DevConcerns_Laterality)
#balloonplot(DevConcerns_Laterality)

Etiology_Laterality <- table(full_elssp$Etiology, full_elssp$Laterality)
Etiology_Laterality <- Etiology_Laterality[-1, c("Bilateral", "Unilateral")]
chisq.test(Etiology_Laterality)
#balloonplot(Etiology_Laterality)

Communication_Laterality <- table(full_elssp$Communication, full_elssp$Laterality)
Communication_Laterality <- Communication_Laterality[-1, c("Bilateral", "Unilateral")]
chisq.test(Communication_Laterality)
#balloonplot(Communication_Laterality)

Meets136_DevDelay <- table(full_elssp$Meets136, full_elssp$DevelopmentalConcerns)
Meets136_DevDelay <- Meets136_DevDelay[-1, ]
chisq.test(Meets136_DevDelay)
#balloonplot(Meets136_DevDelay)

Meets136_Amplification <- table(full_elssp$Meets136, full_elssp$Amplification)
Meets136_Amplification <- Meets136_Amplification[-1, -1]
chisq.test(Meets136_Amplification)
#balloonplot(Meets136_Amplification)

Meets136_Etiology <- table(full_elssp$Meets136, full_elssp$Etiology)
Meets136_Etiology <- Meets136_Etiology[-1, -1]
chisq.test(Meets136_Etiology)
#balloonplot(Meets136_Etiology)

Meets136_Communication <- table(full_elssp$Meets136, full_elssp$Communication)
Meets136_Communication <- Meets136_Communication[-1, -1]
chisq.test(Meets136_Communication)
#balloonplot(Meets136_Communication)

Amplification_DevConcerns <- table(full_elssp$Amplification, full_elssp$DevelopmentalConcerns)
Amplification_DevConcerns <- Amplification_DevConcerns[-1, ]
chisq.test(Amplification_DevConcerns)
#balloonplot(Amplification_DevConcerns)

Amplification_Etiology <- table(full_elssp$Amplification, full_elssp$Etiology)
Amplification_Etiology <- Amplification_Etiology[-1, -1]
chisq.test(Amplification_Etiology)
#balloonplot(Amplification_Etiology)

Amplification_Communication <- table(full_elssp$Amplification, full_elssp$Communication)
Amplification_Communication <- Amplification_Communication[-1, -1]
chisq.test(Amplification_Communication)
#balloonplot(Amplification_Communication)

DevConcerns_Etiology <- table(full_elssp$DevelopmentalConcerns, full_elssp$Etiology)
DevConcerns_Etiology <- DevConcerns_Etiology[, -1]
chisq.test(DevConcerns_Etiology)
#balloonplot(DevConcerns_Etiology)

DevConcerns_Communication <- table(full_elssp$DevelopmentalConcerns, full_elssp$Communication)
DevConcerns_Communication <- DevConcerns_Communication[, -1]
chisq.test(DevConcerns_Communication)
#balloonplot(DevConcerns_Communication)

Communication_Etiology <- table(full_elssp$Communication, full_elssp$Etiology)
Communication_Etiology <- Communication_Etiology[-1, -1]
chisq.test(Communication_Etiology)
#balloonplot(Communication_Etiology)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ Gender, data=full_elssp)
#ggplot(data=full_elssp, aes(Gender, HLworse, fill=Gender))+
  #geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
wilcox.test(ServicesReceivedPerMonth ~ Gender, data=full_elssp)
#ggplot(data=full_elssp, aes(Gender, ServicesReceivedPerMonth, fill=Gender))+
 # geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ HealthIssues, data=full_elssp)
# ggplot(data=full_elssp, aes(HealthIssues, HLworse, fill=HealthIssues))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
wilcox.test(ServicesReceivedPerMonth ~ HealthIssues, data=full_elssp)
# ggplot(data=full_elssp, aes(HealthIssues, ServicesReceivedPerMonth, fill=HealthIssues))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ IsPremature, data=full_elssp)
# ggplot(data=full_elssp, aes(IsPremature, HLworse, fill=IsPremature))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
wilcox.test(ServicesReceivedPerMonth ~ IsPremature, data=full_elssp)
# ggplot(data=full_elssp, aes(IsPremature, ServicesReceivedPerMonth, fill=IsPremature))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ Laterality, data=full_elssp)
# ggplot(data=full_elssp, aes(Laterality, HLworse, fill=Laterality))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
wilcox.test(ServicesReceivedPerMonth ~ Laterality, data=full_elssp)
# ggplot(data=full_elssp, aes(Laterality, ServicesReceivedPerMonth, fill=Laterality))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ Meets136, 
            data=(full_elssp %>% filter(full_elssp$Meets136=="yes"|full_elssp$Meets136=="no")))
# ggplot(data=full_elssp, aes(Meets136, HLworse, fill=Meets136))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
wilcox.test(ServicesReceivedPerMonth ~ Meets136, 
            data=(full_elssp %>% filter(full_elssp$Meets136=="yes"|full_elssp$Meets136=="no")))
# ggplot(data=full_elssp, aes(Meets136, ServicesReceivedPerMonth, fill=Meets136))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
wilcox.test(HLworse ~ DevelopmentalConcerns, 
            data=(full_elssp %>% filter(full_elssp$DevelopmentalConcerns=="yes"|full_elssp$DevelopmentalConcerns=="no")))
# ggplot(data=full_elssp, aes(DevelopmentalConcerns, HLworse, fill=DevelopmentalConcerns))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
Services_DevDelay <- wilcox.test(ServicesReceivedPerMonth ~ DevelopmentalConcerns, 
            data=(full_elssp %>% filter(full_elssp$DevelopmentalConcerns=="yes"|full_elssp$DevelopmentalConcerns=="no")))
# ggplot(data=full_elssp, aes(DevelopmentalConcerns, ServicesReceivedPerMonth, fill=DevelopmentalConcerns))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
Degree_Amplification <- kruskal.test(HLworse ~ Amplification, 
             data=(full_elssp %>% filter(full_elssp$Amplification=="CI"|full_elssp$Amplification=="HA"|full_elssp$Amplification=="none")))
# ggplot(data=full_elssp, aes(Amplification, HLworse, fill=Amplification))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
kruskal.test(ServicesReceivedPerMonth ~ Amplification, 
             data=(full_elssp %>% filter(full_elssp$Amplification=="CI"|full_elssp$Amplification=="HA"|full_elssp$Amplification=="none")))
# ggplot(data=full_elssp, aes(Amplification, ServicesReceivedPerMonth, fill=Amplification))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
kruskal.test(HLworse ~ Etiology, 
             data=(full_elssp %>% filter(full_elssp$Etiology=="Conductive"|full_elssp$Etiology=="Mixed"|full_elssp$Etiology=="SNHL")))
# ggplot(data=full_elssp, aes(Etiology, HLworse, fill=Etiology))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
kruskal.test(ServicesReceivedPerMonth ~ Etiology, 
             data=(full_elssp %>% filter(full_elssp$Etiology=="Conductive"|full_elssp$Etiology=="Mixed"|full_elssp$Etiology=="SNHL")))
# ggplot(data=full_elssp, aes(Etiology, ServicesReceivedPerMonth, fill=Etiology))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$HLworse)
kruskal.test(HLworse ~ Communication, 
             data=(full_elssp %>% filter(full_elssp$Communication=="total communication"|full_elssp$Communication=="spoken")))
# ggplot(data=full_elssp, aes(Communication, HLworse, fill=Communication))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

shapiro.test(full_elssp$ServicesReceivedPerMonth)
Services_Communication <- kruskal.test(ServicesReceivedPerMonth ~ Communication, 
             data=(full_elssp %>% filter(full_elssp$Communication=="total communication"|full_elssp$Communication=="spoken")))
# ggplot(data=full_elssp, aes(Communication, ServicesReceivedPerMonth, fill=Communication))+
#   geom_violin() + geom_jitter(height = 0, width = 0.1)

cor.test(full_elssp$HLworse, full_elssp$ServicesReceivedPerMonth, method="kendall")

##test plots
hm_test <- read.csv("../data/heat_map_test.csv")
hm_test_plot <- ggplot(hm_test, aes(x=var1, y=var2)) +
  geom_tile(aes(fill = p)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_distiller(palette = "RdPu") +
  xlab(NULL) + ylab(NULL) 


hm_test_plot2 <- ggplot(hm_test, aes(x=var1, y=var2)) +
  geom_tile(aes(fill = sig)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab(NULL) + ylab(NULL) 