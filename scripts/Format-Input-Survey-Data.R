Packages <- c("dplyr",  "nleqslv", "broom","cubature", "geosphere", "data.table",  "ggplot2", "bbmle", "stringr",  "lubridate", "RColorBrewer", "viridis")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

Centroids <- fread(file="empirical_data/site_centroids_SimTest.csv") %>%
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
#add in the numbers of particles seeded at each site
SeededParticles <- fread("ROMS/data/Particles_Per_Release_Site_Renamed.csv")
setnames(SeededParticles,c("source", "daily_particles_released")) 
#DateJoin <- SeededParticles[DateJoin, on="source"][, particles_released_daily := as.numeric(particles_released_daily)] 

#make vectors defining sites we didn't sample, but that are in the model, and the sandflats specifically 
unsampled_sites <- c("SF1", "SF2", "SF3", "SF4", "SF5", "SF6", "Pangasugan", "Other", "CAI") 
sand_flats <- c("SF1", "SF2", "SF3", "SF4", "SF5", "SF6") 
unrealistic_sources <- c("SF1", "SF2", "SF3", "SF4", "SF5", "SF6", "Pangasugan") 
#make the constant inputs for the kernel fitting function
#distance matrix using the centroids with combined Magbangon
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
#setorder(Area, site)
reef_sizes <- as.matrix(Area$kmsq)

#make a site index table, use this for Sampled_reefs input in kernel fitting
SiteIndex <- unique(Centroids %>% arrange(site), by="site")[, index := .I] #add the row number as the unique site index, leave CAI in if fitting a kernel 
SiteIndexBioPhys <- unique(Centroids %>% arrange(site), by="site")[site != "CAI"][, index := .I] #add the row number as the unique site index, take CAI out for biophysical likelihood function

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

save(Centroids,  AnnualRecsSamp,  PropSamp, SeededParticles,  unsampled_sites,  sand_flats, unrealistic_sources,  SitesSource,  SitesDest,  Distances, Area, reef_sizes, SiteIndex, SiteIndexBioPhys, SurveyData, Abundance, SurveyData, file= "script_output/SurveyData/for_likelihood_functions/2021-11-04_InputSurveyData.Rdata")