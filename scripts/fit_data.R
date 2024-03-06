# 12/22/2023: Modified from fit_data.ipynb by Malin Pinsky to produce site-specific log likelihoods.
# assumes that the working directory has been set to the CamotesROMS directory
# the log-likelihoods are from the kernel fit against the biophysical data
# doesn't really mean much. we instead want the LL froms the biophysical prediction of the parentage data

Packages <- c("dplyr",  "nleqslv", "broom","cubature", "geosphere", "data.table",  "ggplot2", "bbmle", "stringr",  "lubridate", "RColorBrewer")
invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
source("scripts/ll_kt_both_bbmle.R")
source("scripts/ll_kt_both_bbmle_per_site.R") # new as of 12/2023
source("scripts/ll_kt_both_grid_search.R")
source("scripts/LL_biophys.R") # not used? loads the LL_biophys() function
source("scripts/GenGausKernInt_sum0.5.R") #integrate_kernel_sum1
source("scripts/GenGausKernInt_sum1.R")
source("scripts/cdf_solve.R") #median
source("scripts/cdf_solve90.R") #dist 90% retained

set.seed(7)
#read in the kernel fitting summary
kernels <- fread(file="script_output/KernelFits/summary_tables/kernel_fitting_summary.csv")
kernel2012_14 <- fread(file="empirical_data/genetics/GenKernelsForROMSComp2012-14.csv")

#read in the centroids adjusted for the simulation, so the Magbangons combined 
#centroids <- fread(file="script_output/SurveyData/SimulationCentroids.csv")
Centroids <- fread(file="empirical_data/site_centroids_SimTest.csv") 
Centroids$site <- gsub("_full", "", Centroids$site, fixed=TRUE)
Centroids$site <- gsub("_ten_per_cover", "", Centroids$site, fixed=TRUE)

Centroids <- Centroids %>%
    arrange(site)

#setorder(Centroids, site)#warning! This sets order based on site, and then lat/lon. So the table is not alphabetical by site, but that's fine as long as all of the "sampled_reef" vectors reflect this, so that reef_sizes, distance, and sampled_reefs match up by row/col index 
#read in the table with number of recruits sampled at each site for each year
AnnualRecsSamp <- fread(file="script_output/SurveyData/AnnualRecruitsSampled.csv")
#read in the table of the proportion of anemones sampled at each site for each year
PropSamp <- unique(fread(file="script_output/SurveyData/ProportionHabitatSampled.csv")[
    , .(site, year=end_year, prop_anem_samp=total_prop_hab_sampled_anems_tidied)][ #select and rename columns with the tideied data to use
    site %like% "Magbangon", site := "Magbangon"][ #collapse Magbangon values
    , prop_anem_samp := sum(prop_anem_samp), by=c("site", "year")], by=c("site", "year"))[ #collapse magbangons to match ROMS data
    site=="Sitio Lonas" & year %in% c(2012, 2013, 2014), prop_anem_samp :=1][site=="Caridad Proper" & year %in% c(2013, 2014), prop_anem_samp :=1]
                   
                   
### List of source locations
SitesSource <- Centroids

### List of destination locations
SitesDest <- Centroids

DistMatm <- distm(SitesSource[,c('lon','lat')], SitesSource[,c('lon','lat')], fun=distVincentyEllipsoid)
Distances <- DistMatm*10^-3
#read in the reef areas for the kernel fitting
Area <- fread("empirical_data/site_area_header_nonsurveyed_simulation_kernels_test.csv") %>%
    arrange(site) %>%
    filter(site %!in% c("near_north_full1", "near_north_full2", "near_north_full3", "near_south_full1", "near_south_full2", "near_south_full3")) %>%
    mutate(kmsq=msq*10^-6)# %>%
    #select(kmsq) #need to uncomment for functions to work
Area$site <- gsub("_ten_per_cover", "", Area$site, fixed=TRUE)

reef_sizes <- as.matrix(Area$kmsq)

#make a site index table, use this for Sampled_reefs input in kernel fitting
SiteIndex <- unique(Centroids %>% arrange(site), by="site")[, index := .I] #add the row number as the unique site index, leave CAI in if fitting a kernel 
SiteIndexBioPhys <- unique(Centroids %>% arrange(site), by="site")[site != "CAI" ][, index := .I] #add the row number as the unique site index, take CAI out for biophysical likelihood function

#make a table with the survey information for each site (how many fish sampled, prop anems sampled, total number of anems at site)
SurveyData <- AnnualRecsSamp[PropSamp, on=.(year, site)][#join the sampling tables together
    is.na(n_offs_gen), n_offs_gen := 0]#change NA's to 0
#setnames(SurveyData, c("PropAnemSamp", "TotalAnems"), c("prop_anem_samp", "total_anems"))
#setkey(SurveyData, site)
#check all sites are represented in centroids and area (and indirectly distances, which comes from centroids)
#Area[site %!in% centroids$site] #should be nothing

#Allison's abundance time series data 
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Script_outputs/females_df_F.RData?raw=true", destfile = "empirical_data/genetics/females_df_F.RData")
load("empirical_data/genetics/females_df_F.RData")
Abundance <- as.data.table(females_df_F)
setnames(Abundance, "nF", "num_females")
Abundance <- unique(Abundance[site %like% "Magbangon", site := "Magbangon"][ #collapse Magbangon values
            , num_females := sum(num_females), by=c("site", "year")], by=c("site", "year"))
#join the survey sampling tables together
SurveyData <- AnnualRecsSamp[PropSamp, on=.(year, site)][
    is.na(n_offs_gen), n_offs_gen := 0]#change NA's to 0


SurveyData <- Abundance[, c("year", "site", "num_females")][SurveyData, on=.(year, site)]#join in Allison's estimate of female abundance. There are NA values, but that's okay we can figure those out when we start thinking about incorporating uncertainty in this
#quick check that all components are in the same, alphabetical order
sum(which(SiteIndex$site==Area$site)==FALSE) #needs to be 0!! sites have to be in the same order
sum(which(Area$site==Centroids$site)==FALSE) #needs to be 0!! sites have to be in the same order

##read in genetic parentage matrices for each time frame MATCHING DIMENSIONS WITH BIOPHYSICAL FOR BIOPHYSICAL LIKELIHOOD- NOT KERNEL FITTING
#GenMat2012_4 <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrix2012-14ForROMSComp.csv"))
#GenMat2012 <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrix2012ForROMSComp.csv"))
#GenMat2013 <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrix2013ForROMSComp.csv"))
#GenMat2014 <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrix2014ForROMSComp.csv"))
#GenMatNEM <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrixNEM2012-14ForROMSComp.csv"))
#GenMatSWM <- as.matrix(fread(file="script_output/SurveyData/20210625_ParentageMatrixSWM2012-14ForROMSComp.csv"))

#read in genetic parentage matrices for each time frame FOR KERNEL FITTING
KernelGenMat2012_4 <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrix2012-14ForROMSComp.csv"))
KernelGenMat2012 <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrix2012ForROMSComp.csv"))
KernelGenMat2013 <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrix2013ForROMSComp.csv"))
KernelGenMat2014 <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrix2014ForROMSComp.csv"))
KernelGenMatNEM <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrixNEM2012-14ForROMSComp.csv"))
KernelGenMatSWM <- as.matrix(fread(file="script_output/SurveyData/20210701_KernelParentageMatrixSWM2012-14ForROMSComp.csv"))

##read in biophysical data for KERNEL FITTING
FullBiophysMat <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrix2012-14ForROMSComp08DayPLD.csv"))
#fwrite(biophys_parentage_matrices$file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrix2012-14ForROMSComp08DayPLD.csv")

AnnualBiophysMat2012 <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrix2012ForROMSComp08DayPLD.csv"))
AnnualBiophysMat2013 <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrix2013ForROMSComp08DayPLD.csv"))
AnnualBiophysMat2014 <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrix2014ForROMSComp08DayPLD.csv"))
MonsoonBiophysMatNEM <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrixNEMForROMSComp08DayPLD.csv"))
MonsoonBiophysMatSWM <- as.matrix(fread(file="script_output/ROMSDataTables/2022-01-20_BioPhysParentageMatrixSWMForROMSComp08DayPLD.csv"))

sum(AnnualBiophysMat2012)+sum(AnnualBiophysMat2013)+sum(AnnualBiophysMat2014)==sum(FullBiophysMat)#should be true

sum(MonsoonBiophysMatNEM)
sum(MonsoonBiophysMatSWM)


#fit the kernels, get the biophysical data together
biophys_par_data2012 <- list(Distances=Distances, Assignments=AnnualBiophysMat2012, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(AnnualBiophysMat2012), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(AnnualBiophysMat2012), ncol=1, 1))
Sim2012Fit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-3, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_data2012, control=list(maxit=500)))

biophys_par_data2013 <- list(Distances=Distances, Assignments=AnnualBiophysMat2013, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(AnnualBiophysMat2013), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(AnnualBiophysMat2013), ncol=1, 1))
Sim2013Fit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-3, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_data2013, control=list(maxit=500)))

biophys_par_data2014 <- list(Distances=Distances, Assignments=AnnualBiophysMat2014, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(AnnualBiophysMat2014), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(AnnualBiophysMat2014), ncol=1, 1))
Sim2014Fit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-3, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_data2014, control=list(maxit=500)))

biophys_par_data2012_4 <- list(Distances=Distances, Assignments=FullBiophysMat, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(FullBiophysMat), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(FullBiophysMat), ncol=1, 1))
Sim2012_4Fit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-3, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_data2012_4, control=list(maxit=500)))

biophys_par_dataNEM <- list(Distances=Distances, Assignments=MonsoonBiophysMatNEM, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(FullBiophysMat), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(MonsoonBiophysMatNEM), ncol=1, 1))
SimNEMFit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-.3, theta=.6), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_dataNEM, control=list(maxit=500)))

biophys_par_dataSWM <- list(Distances=Distances, Assignments=MonsoonBiophysMatSWM, Sampled_reefs=t(as.matrix(SiteIndex[site %in% colnames(FullBiophysMat), .(index)])), 
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(MonsoonBiophysMatSWM), ncol=1, 1))
SimSWMFit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-.3, theta=.6), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=biophys_par_dataSWM, control=list(maxit=500)))


# extract the site-by-site log-likelihoods.
LL_by_site_year <- LL_kt_bbmle_per_site(k=Sim2012Fit@coef[1], theta=Sim2012Fit@coef[2], 
	Assignments= biophys_par_data2012$Assignments, 
	Sampled_reefs= biophys_par_data2012$Sampled_reefs, 
	Distances= biophys_par_data2012$Distances, 
	Reef_sizes= biophys_par_data2012$Reef_sizes, 
	Adult_sample_proportions= biophys_par_data2012$Adult_sample_proportions,
	LLname = '2012')
LL_by_site_year <- merge(LL_by_site_year, 
	LL_kt_bbmle_per_site(k=Sim2013Fit@coef[1], theta=Sim2013Fit@coef[2], 
		Assignments= biophys_par_data2013$Assignments, 
		Sampled_reefs= biophys_par_data2013$Sampled_reefs, 
		Distances= biophys_par_data2013$Distances, 
		Reef_sizes= biophys_par_data2013$Reef_sizes, 
		Adult_sample_proportions= biophys_par_data2013$Adult_sample_proportions,
		LLname = '2013'))
LL_by_site_year <- merge(LL_by_site_year, 
	LL_kt_bbmle_per_site(k=Sim2014Fit@coef[1], theta=Sim2014Fit@coef[2], 
		Assignments= biophys_par_data2014$Assignments, 
		Sampled_reefs= biophys_par_data2014$Sampled_reefs, 
		Distances= biophys_par_data2014$Distances, 
		Reef_sizes= biophys_par_data2014$Reef_sizes, 
		Adult_sample_proportions= biophys_par_data2014$Adult_sample_proportions,
		LLname = '2014'))
LL_by_site_year <- merge(LL_by_site_year, 
	LL_kt_bbmle_per_site(k=SimNEMFit@coef[1], theta=SimNEMFit@coef[2], 
		Assignments= biophys_par_dataNEM$Assignments, 
		Sampled_reefs= biophys_par_dataNEM$Sampled_reefs, 
		Distances= biophys_par_dataNEM$Distances, 
		Reef_sizes= biophys_par_dataNEM$Reef_sizes, 
		Adult_sample_proportions= biophys_par_dataNEM$Adult_sample_proportions,
		LLname = 'NEM'))
LL_by_site_year <- merge(LL_by_site_year, 
	LL_kt_bbmle_per_site(k=SimSWMFit@coef[1], theta=SimSWMFit@coef[2], 
		Assignments= biophys_par_dataSWM$Assignments, 
		Sampled_reefs= biophys_par_dataSWM$Sampled_reefs, 
		Distances= biophys_par_dataSWM$Distances, 
		Reef_sizes= biophys_par_dataSWM$Reef_sizes, 
		Adult_sample_proportions= biophys_par_dataSWM$Adult_sample_proportions,
		LLname = 'SWM'))

LL_by_site_year # examine the results