This repo contains the data for Characterizing North Carolina's Deaf/Hard-of-Hearing Infants and Toddlers: Predictors of Vocabulary, Diagnosis, and Intervention.

**elssp_comorbidities.csv**: contains information about the diagnoses experienced by different children in the sample. Columns:
	- VIHI_ID: subject id
	- Age: age in months
	- WeeksGestation: weeks gestation at which child was born. if child was reported full-term, 40 is listed in column
	- other columns: different diagnoses with 1/0 to indicate whether child has/doesn't have that condition
**ELSSP_SubjectInformation**: contains information about the children in the sample. Columns:
	- VIHI_ID: subject id with age identifier
	- SubjectNumber: subject id without age identifier
	- VisitNumber: whether this was the first evaluation for the child
	- Age: in months
	- Gender: male/female
	- Meets136: whether child received diagnosis by 3 months and intervention by 6 months
	- Meets13: whether child received diagnosis by 3 months
	- Meets6: whether child received diagnosis by 6 months
	- NBHS: whether child received newborn hearing screening (yes/no/no_info)
	- IdentificationOfHLMonths: age at which hearing loss was identified
	- HLbetter: degree of hearing loss in better ear
	- HLworse: degree of hearing loss in worse ear
	- Hlaverage: average of HLbetter, HLworse
	- Etiology: conductive/SNHL/mixed
	- Laterality: of hearing loss, unilateral/bilateral
	- Side: for unilateral hearing loss, which ear has hearing loss
	- ANSD: 0/1 for auditory neuropathy spectrum disorder
	- AgeAmplification: age at amplification in months
	- AgeImplantation: age at implantation in months
	- Amplification: hearing aid, cochlear implant, or none (HA/CI/none)
	- AgeStartedServices: age at which intervention began
	- ServicesReceivedPerMonth: # of clinical services received per month
	- Communication: spoken/total communication/cued speech
	- IsPremature: whether child was born >3 months premature
	- HealthIssues: whether the child has other medical diagnoses
	- DevelopmentalConcerns: whether child has cognitive/developmental diagnoses or documented concerns about development
	- VisionLoss: whether child has a visual impairment
	- CDIversion: version of cdi received (words & gestures [WG] or words & sentences [WS])
	- ProductionCDI: words produced on CDI
	- ComprehensionCDI: words understood on CDI (NA for WS)
	- Monolingual_English: whether household is monolingual & english-speaking
	- LanguageBackground: language spoken most in home
**vocabulary_norms_table_XXX.csv**: for Words & Gestures (WG) or Words & Sentences (WS) either for American English or Mexican Spanish. Downloaded from wordbank.stanford.edu

Last updated: 4/30/2021
Maintained by: Erin Campbell (erin.e.campbell@duke.edu)
