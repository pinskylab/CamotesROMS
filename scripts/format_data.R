Packages <- c("dplyr",  "nleqslv", "broom","cubature", "geosphere", "data.table",  "ggplot2", "bbmle", "stringr",  "lubridate", "RColorBrewer")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/oceanography')
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_bbmle.R")
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_grid_search.R")
source("~/oceanography/scripts/neg_LL_biophys.R")
#source("~/oceanography/scripts/PredictedProportions.R")

#read in the kernel fitting summary
kernels <- fread(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv")
kernel2012_14 <- fread(file="~/oceanography/empirical_data/genetics/GenKernelsForROMSComp2012-14.csv")

#read in the centroids adjusted for the simulation, so the Magbangons combined 
#centroids <- fread(file="~/oceanography/script_output/SurveyData/SimulationCentroids.csv")
Centroids <- fread(file="~/oceanography/empirical_data/site_centroids_SimTest.csv")
setorder(Centroids, site)
#read in the table with number of recruits sampled at each site for each year
AnnualRecsSamp <- fread(file="~/oceanography/script_output/SurveyData/AnnualRecruitsSampled.csv")
#read in the table of the proportion of anemones sampled at each site for each year
PropSamp <- fread(file="~/oceanography/script_output/SurveyData/ProportionHabitatSampled.csv")
setnames(PropSamp, c("PropAnemSamp", "TotalAnems"), c("prop_anem_samp", "total_anems"))
#read in the ROMS simulation connectivity table with metadata, not yet subsetted (*but check this)
#SimConn <- fread(file="~/oceanography/script_output/ROMSDataTables/SimConnectivityTableWithMetaLongForm.csv")

#add in the numbers of particles seeded at each site
SeededParticles <- fread("~/oceanography/ROMS/data/Particles_Per_Release_Site_Renamed.csv")
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
Area <- fread("~/oceanography/empirical_data/site_area_header_nonsurveyed_simulation_kernels_test.csv") %>%
    arrange(site) %>%
    filter(site %!in% c("near_north_full1", "near_north_full2", "near_north_full3", "near_south_full1", "near_south_full2", "near_south_full3")) %>%
    mutate(kmsq=msq*10^-6)# %>%
    #select(kmsq) #need to uncomment for functions to work
setorder(Area, site)
reef_sizes <- as.matrix(Area$kmsq)

#make a site index table, use this for Sampled_reefs input in kernel fitting
SiteIndex <- unique(Centroids, by="site")[, index := .I] #add the row number as the unique site index, leave CAI in if fitting a kernel 
SiteIndexBioPhys <- unique(Centroids, by="site")[site != "CAI"][, index := .I] #add the row number as the unique site index, take CAI out for biophysical likelihood function

#make a table with the survey information for each site (how many fish sampled, prop anems sampled, total number of anems at site)
SurveyData <- AnnualRecsSamp[PropSamp, on=.(year=end_year, site)][#join the sampling tables together
    is.na(n_offs_gen), n_offs_gen := 0][,#change NA's to 0
    -"time_frame"]#drop the time_frame column, we can key with end_year
#setnames(SurveyData, c("PropAnemSamp", "TotalAnems"), c("prop_anem_samp", "total_anems"))
#setkey(SurveyData, site)
#check all sites are represented in centroids and area (and indirectly distances, which comes from centroids)
#Area[site %!in% centroids$site] #should be nothing

#Allison's abundance time series data 
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Script_outputs/females_df_F.RData?raw=true", destfile = "~/oceanography/empirical_data/genetics/females_df_F.RData")
load("~/oceanography/empirical_data/genetics/females_df_F.RData")
Abundance <- as.data.table(females_df_F)
setnames(Abundance, "nF", "num_females")
Abundance <- unique(Abundance[site %like% "Magbangon", site := "Magbangon"][ #collapse Magbangon values
            , num_females := sum(num_females), by=c("site", "year")], by=c("site", "year"))
#join the survey sampling tables together
SurveyData <- AnnualRecsSamp[PropSamp, on=.(year=end_year, site)][
    is.na(n_offs_gen), n_offs_gen := 0][,#change NA's to 0
    -"time_frame"]#drop the time_frame

SurveyData <- Abundance[, c("year", "site", "num_females")][SurveyData, on=.(year, site)]#join in Allison's estimate of female abundance. There are NA values, but that's okay we can figure those out when we start thinking about incorporating uncertainty in this



#make a table of the dates of release for simulations, for calculating the number of particles released in each time frame
season1 <- data.table(date=seq(as.Date("2010-10-01"), as.Date("2011-05-31"), by="days"))

season2 <- data.table(date=seq(as.Date("2011-10-01"), as.Date("2012-05-31"), by="days"))

season3 <- data.table(date=seq(as.Date("2012-10-01"), as.Date("2013-05-31"), by="days"))

season4 <- data.table(date=seq(as.Date("2013-10-01"), as.Date("2014-04-18"), by="days"))

AllDates <- rbind(season1, season2, season3, season4)

#mark the monsoon seasons, based on the same criteria I used for the parentage indirectly through the growth estimates
NEM <- c(11, 12, 1, 2, 3, 4, 5, 6)
SWM <- c(7, 8, 9, 10)

AllDates[,date := ymd(date)][, #format as ymd
             sim_monsoon := ifelse(month(date) %in% NEM, "NEM", "SWM")][,#mark monsoon season based on month
             sim_year:=year(date)][,#add year column
            year_sampled:= ifelse(date %in% season1$date, 2011, ifelse(date %in% season2$date, 2012, ifelse(date %in% season3$date, 2013, 2014)))]#and then add a year_sampled for the empircal sampling season of that particle

ReleaseDays <- AllDates[, .(num_release_days_seasonal=.N), by=c("year_sampled", "sim_monsoon")][, num_release_days_annual:= sum(num_release_days_seasonal), by=year_sampled]

total_release_days <- AllDates[year_sampled %in% c(2012, 2013, 2014), .N]#for the all year kernel- how many days of the simulation conincide with our particle sampling?
total_release_days #should be 687


##prep biophysical connectivity matrix
##outside of the loop, trim this to only be the destinations we sampled
#SourceJoin <- SurveyData[SimConn, on = .(site = source, year=year_sampled)]
#setnames(SourceJoin, skip_absent=TRUE, c("site", "n_offs_gen", "prop_anem_samp", "total_anems", "num_females"), c("source", "source_num_rec_sampled_annual",  "source_prop_anem_samp", "source_total_anems", "source_num_females"))
#DestJoin <- SurveyData[SourceJoin, on = .(site = dest, year)]
#setnames(DestJoin, skip_absent=TRUE, c("site", "n_offs_gen", "prop_anem_samp", "total_anems", "num_females"), c("dest", "dest_num_rec_sampled_annual",  "dest_prop_samp", "dest_total_anems", "dest_num_females"))
#
#SimConn <- DestJoin[source %!in% unrealistic_sources & dest %!in% unrealistic_sources & year %in% c(2012, 2013, 2014)][#sand flats and Pangasugan are not realistic source or destination sites because there's almost no habitat. Safe to drop, but keep the rest of the possibilities so we can subsample iteratively all possibilities.
#    , daily_particles_released := as.numeric(daily_particles_released)] #change from integer to numeric
#SimConn <- ReleaseDays[SimConn, on=.(year_sampled=year, sim_monsoon)]#join in the info for number of release days in the time frame
#SimConn <- kernels[Year %in% c("2012", "2013", "2014")][, year:=as.integer(Year)][,c("year", "NumParentageMatches")][SimConn, on=.(year=year_sampled)]#add in a column for the observed number of parentage matches
##rename the monsoon column in the full table for consistency
#setnames(SimConn, c("sim_monsoon", "NumParentageMatches"), c("monsoon", "num_route_parentage_matches")) #get rid of upper case and inconsistent naming
#setcolorder(SimConn, c("particle_id", "source", "dest", "year", "monsoon", "date"))
#
##at this point, we can make the raw number assignment matrix, but we want to make a normalized version that is num assigned from a source to a destination/ num released from that source
##fwrite(SimConn, file="~/oceanography/script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm.csv")

SimConn <- fread(file="~/oceanography/script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm.csv")[dest != "CAI"] #filter out CAI as a destination for now, not very well spatially defined


#each year will require a different set of survey data, so make a list of each and index by site for fast look up
SampledTable <- SurveyData[prop_anem_samp >0, c("year", "site")]#previously named PropSampTable

#make sure all sampled sites are represented when joining the survey data to the sampled simulation- this chunk has the tables to add to a subsampled particle table. no need for the full
SampTable <- rbind(SurveyData[prop_anem_samp >0 & year %in% c(2012, 2013, 2014), c("year", "site")][, .(source=site, dest=site, year=year)][, #will join to the simulated sampling table by source and dest, so make those each a column from site and preserve the year variable as a key
     c("year", "source", "dest")][, monsoon := "NEM"], SurveyData[prop_anem_samp >0 & year %in% c(2012, 2013, 2014), c("year", "site")][, .(source=site, dest=site, year=year)][, #will join to the simulated sampling table by source and dest, so make those each a column from site and preserve the year variable as a key
     c("year", "source", "dest")][, monsoon := "SWM"])

UnqSurvey <- unique(SampTable, by=c("source", "dest", "year", "monsoon"))#add in the diff Monsoon seasons so there are complete parentage matrices later

AddDestSim <- rbind(rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "NEM"],
    rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "SWM"])



#make a parentage matrix for the whole biophysical results- CAI included with Other as "unknown", use this to fit a kernel
FullBiophysMat <- as.matrix(rbind(dcast(SimConn[source != "Other" & source != "CAI" , .(source, dest)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    , parentage :=1][ #mark each row as a parentage match, because at this point I'm using all particles as matches for the simulations
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="parentage", fun.aggregate = sum)[#use sum to count the matches for each id variable combo, that populated the cells of the matrix
    ,-"source"], #remove the source column after casting
      dcast(SimConn[source == "Other" | source == "CAI", .(source, dest)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          , parentage :=1][order(source, dest)][, source := "unknown"], source ~ dest, value.var="parentage", fun.aggregate = sum)[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(FullBiophysMat)

x <- list(Distances=Distances, Assignments=FullBiophysMat, Sampled_reefs=t(SiteIndex[site %in% SurveyData[, site], index]), #if CAI is it's own site- site %in% AllYearsRec[, dest]
                  Reef_sizes=reef_sizes, Adult_sample_proportions=matrix(nrow=ncol(FullBiophysMat), ncol=1, 1)) #put inputs into a list because that's the bbmle format
Sim2012_4Fit <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=-3, theta=1), lower=c(-10, 0.15), upper=c(10, 8), method="L-BFGS-B", data=x, control=list(maxit=500)))
Sim2012_4Fit
#Next, do grid search to get the likelihood profile to compare to genetics

# use a grid search to find k that minimizes the log likelihood
k_eval <- c(seq(from=-10, to=10, by=0.01))
theta_eval <- c(seq(from=0.1, to=5, by=.01))
nll_matrix <- matrix(data=NA, nrow=length(k_eval), ncol=length(theta_eval), 
                     dimnames=list(k_eval, theta_eval))
pb <- txtProgressBar(min = 0, max = length(k_eval), style = 3)#7 years in each interation

#Begin grid search: i <- 1; j <- 1
for(i in 1:length(k_eval)){
  k <- k_eval[i]
  for(j in 1:length(theta_eval)){
    theta <- theta_eval[j]
    nll_matrix[i,j] <- LL_kt_grid(k=k, theta=theta, Distances=x$Distances, Assignments=x$Assignments, Sampled_reefs=x$Sampled_reefs, Reef_sizes=x$Reef_sizes, Adult_sample_proportions=x$Adult_sample_proportions)
  }
    setTxtProgressBar(pb, i)

}


close(pb)

max(nll_matrix, na.rm = T)
min(nll_matrix, na.rm = T)

sum(is.na(nll_matrix)==T)
best_params_index <- which(nll_matrix == min(nll_matrix, na.rm = T), arr.ind=TRUE)
best_params <- c(k_eval[best_params_index[1]], theta_eval[best_params_index[2]])
best_params


#write profile results
#write.csv(nll_matrix, file="~/oceanography/script_output/KernelFits/LikelihoodProfileBiophysical2012-4.csv", row.names=T, quote=FALSE)



#at this point, we can make the raw number assignment matrix, but we want to make a normalized version that is num assigned from a source to a destination/ num released from that source
total_release_days <- 687

AllYearsRec <- SimConn[ , .(total_particles_rec = .N), by= c("source","dest")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released)*as.numeric(total_release_days)), by= c("source")], by="source") #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

#make sure all possible routes are represented
AddDestAllYearsSim <- unique(AddDestSim, by=c("source","dest"))[, -c("year","monsoon")] 

AllYearsRecInt <- rbind(AddDestAllYearsSim[!AllYearsRec, on =.(source, dest)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , total_particles_rec:=0 ], AllYearsRec) 

AllYearsReleaseInt <- rbind(AddDestAllYearsSim[!AllYearsRelease, on =.(source)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , total_particles_released:=0 ][, -"dest"], AllYearsRelease) 

#collapse CAI and "Other" together into "unknown" before making the normalized matrices
AllYearsRec <- unique(AllYearsRecInt[source== "Other" | source == "CAI", source := "unknown"][
    , total_particles_rec := sum(total_particles_rec),  by=c("source", "dest")], by=c("source", "dest"))
AllYearsRelease <- unique(AllYearsReleaseInt[source== "Other" | source == "CAI", source := "unknown"][
    , total_particles_released := sum(total_particles_released), by="source"], by="source")

#join recruited and released tables together and make a column for the normalized values
AllYearsNormConn <-  AllYearsRelease[AllYearsRec, on="source"][
    , source_norm_rec := total_particles_rec/total_particles_released]

#check that they sum to =< 1
#AllYearsNormConn[,sum(source_norm_rec), by="source"]#nothing should be greater than 1. It isn't- great

#make sure all possible routes are represented!!*
#cast into wide format
FullBiophysMatNorm <- as.matrix(rbind(dcast(AllYearsNormConn[source !="unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AllYearsNormConn[source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(FullBiophysMatNorm )
FullBiophysMatNorm[is.na(FullBiophysMatNorm)] <- 0 #change NAs to zeros


#make annual matrices with all of the particle data
AnnualRec <- SimConn[ , .(annual_particles_rec = .N), by= c("source","dest", "year")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
AnnualRelease <- unique(SimConn[, .(annual_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_annual)), by= c("source", "year")], by= c("source", "year")) #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

#make sure all possible routes are represented
AddDestAnnualSim <- unique(AddDestSim, by=c("source","dest", "year"))[, -"monsoon"] 

AnnualRecInt <- rbind(AddDestAnnualSim[!AnnualRec, on =.(source, dest, year)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , annual_particles_rec:=0 ], AnnualRec) 

AnnualReleaseInt <- rbind(AddDestAnnualSim[!AnnualRelease, on =.(source, year)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , annual_particles_released:=0 ][, -"dest"], AnnualRelease) 

#collapse CAI and "Other" together into "unknown" before making the normalized matrices
AnnualRec <- unique(AnnualRecInt[source== "Other" | source == "CAI", source := "unknown"][
    , annual_particles_rec := sum(annual_particles_rec), by=c("source", "dest", "year")], by=c("source", "dest", "year"))
AnnualRelease <- unique(AnnualReleaseInt[source== "Other" | source == "CAI", source := "unknown"][
    , annual_particles_released := sum(annual_particles_released), by=c("source", "year")], by=c("source", "year"))

#join recruited and released tables together and make a column for the normalized values
AnnualNormConn <-  AnnualRelease[AnnualRec, on=c("source", "year")][
    , source_norm_rec := annual_particles_rec/annual_particles_released]
#check that they sum to =< 1
#summary(AnnualNormConn[, .(sum=sum(source_norm_rec)), by=c("year", "source")][, sum])#nothing should exceed 1, it doesn't- great

#cast into wide format for each year

AnnualBiophysMatNorm2012 <- as.matrix(rbind(dcast(AnnualNormConn[year==2012 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2012][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(AnnualBiophysMatNorm2012)
AnnualBiophysMatNorm2012[is.na(AnnualBiophysMatNorm2012)] <- 0 #change NAs to zeros

AnnualBiophysMatNorm2013 <- as.matrix(rbind(dcast(AnnualNormConn[year==2013 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2013][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(AnnualBiophysMatNorm2013)
AnnualBiophysMatNorm2013[is.na(AnnualBiophysMatNorm2013)] <- 0 #change NAs to zeros

AnnualBiophysMatNorm2014 <- as.matrix(rbind(dcast(AnnualNormConn[year==2014 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2014][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(AnnualBiophysMatNorm2014)
AnnualBiophysMatNorm2014[is.na(AnnualBiophysMatNorm2014)] <- 0 #change NAs to zeros

#make monsoon matrices with all of the particle data
MonsoonRec <- SimConn[ , .(monsoon_particles_rec = .N), by= c("source","dest", "monsoon")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_seasonal)), by= c("source", "monsoon")], by= c("source", "monsoon")) #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

#make sure all possible routes are represented

AddDestMonsoonSim <- unique(AddDestSim[, -"year"], by=c("source", "dest", "monsoon"))
MonsoonRec <- SimConn[ , .(monsoon_particles_rec = .N), by= c("source","dest", "monsoon")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_seasonal)), by= c("source", "monsoon")], by= c("source", "monsoon")) #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

MonsoonRecInt <- rbind(AddDestMonsoonSim[!MonsoonRec, on =.(source, dest, monsoon)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , monsoon_particles_rec:=0 ], MonsoonRec) 

MonsoonReleaseInt <- rbind(AddDestMonsoonSim[!MonsoonRelease, on =.(source, monsoon)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , monsoon_particles_released:=0 ][, -"dest"], MonsoonRelease) 

#collapse CAI and "Other" together into "unknown" before making the normalized matrices
MonsoonRec <- unique(MonsoonRecInt[source== "Other" | source == "CAI", source := "unknown"][
    , monsoon_particles_rec := sum(monsoon_particles_rec), by=c("source", "dest", "monsoon")], by=c("source", "dest", "monsoon"))
MonsoonRelease <- unique(MonsoonReleaseInt[source== "Other" | source == "CAI", source := "unknown"][
    , monsoon_particles_released := sum(monsoon_particles_released), by=c("source", "monsoon")], by=c("source", "monsoon"))

#join recruited and released tables together and make a column for the normalized values
MonsoonNormConn <-  MonsoonRelease[MonsoonRec, on=c("source", "monsoon")][
    , source_norm_rec := monsoon_particles_rec/monsoon_particles_released]
#check that they sum to =< 1
#summary(MonsoonNormConn[, .(sum=sum(source_norm_rec)), by=c("monsoon", "source")][, sum])#nothing should exceed 1, it doesn't- great

MonsoonBiophysMatNormNEM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon=="NEM" & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(MonsoonNormConn[monsoon=="NEM"][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(MonsoonBiophysMatNormNEM)
MonsoonBiophysMatNormNEM[is.na(MonsoonBiophysMatNormNEM)] <- 0 #change NAs to zeros

MonsoonBiophysMatNormSWM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon=="SWM" & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(MonsoonNormConn[monsoon=="SWM"][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
dim(MonsoonBiophysMatNormSWM)
MonsoonBiophysMatNormSWM[is.na(MonsoonBiophysMatNormSWM)] <- 0 #change NAs to zeros

#read in the genetic parentage data and format for comparison
TotalPar2012_4 <- fread(file="~/oceanography/empirical_data/genetics/parentage_table_2012-14.csv")
TotalPar2012_4 <- unique(TotalPar2012_4[offs_site %like% "Magbangon", offs_site := "Magbangon"][par_site %like% "Magbangon", par_site := "Magbangon"][ #collapse Magbangon values
            , n_matches := sum(n_matches), by=c("offs_site", "par_site", "year")], by=c("offs_site", "par_site", "year"))

#add in all the sampled sites and numbers of recruits sampled at each site
TotalParInt <- unique(TotalPar2012_4[, num_matches := sum(n_matches), by=c("year", "offs_site")][, -"par_site"], by=c("year", "offs_site"))[, -"n_matches"]
#sum(TotalParInt$num_matches) #should be 37

TotalUnassigned2012_4 <- TotalParInt[SurveyData[year %in% c(2012, 2013, 2014)], on=.(offs_site=site, year)]
TotalUnassigned2012_4 [is.na(TotalUnassigned2012_4 )] <- 0

#sum(TotalUnassigned2012_4 $n_offs_gen) #should be 394
#sum(TotalUnassigned2012_4 $num_matches) #should be 37
#nrow(unique(TotalUnassigned2012_4 , by="offs_site"))#should be 18, so that every site is represented so that the years are all 18*18 sites
TotalUnassigned2012_4 <- TotalUnassigned2012_4[, num_unassigned := n_offs_gen-num_matches, by=c("year", "offs_site")][
    , .(offs_site, year, num_unassigned)]
#sum(TotalUnassigned2012_4$num_unassigned) #should be 357

#add in sites that aren't represented in every year
AddDestGen <- rbindlist(list(unique(cbind(SurveyData[year %in% c(2012, 2013, 2014) ][, .(offs_site=site)], 
                  SurveyData[year %in% c(2012, 2013, 2014)][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
unique(cbind(SurveyData[year %in% c(2012, 2013, 2014) ][, .(offs_site=site)], 
                  SurveyData[year %in% c(2012, 2013, 2014)][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
unique(cbind(SurveyData[year %in% c(2012, 2013, 2014) ][, .(offs_site=site)], 
                  SurveyData[year %in% c(2012, 2013, 2014)][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))

TotalPar2012_4 <- rbind(AddDestGen[!TotalPar2012_4, on =.(par_site, offs_site, year)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , n_matches:=0 ], TotalPar2012_4[,-"num_matches"])  #add the parentage column, add back into the parentage table but drop the num_matches column that's a summary column I used to make the unassigned table 
#sum(TotalPar2012_4$n_matches)#should be 37 still

#format genetic parentage matrices for each year and all years combined
GenMat2012 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2012, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2012][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
GenMat2012[is.na(GenMat2012)] <- 0

GenMat2013 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2013, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2013][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
GenMat2013[is.na(GenMat2013)] <- 0

GenMat2014 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2014, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2014][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
GenMat2014[is.na(GenMat2014)] <- 0

GenMat2012_4 <- as.matrix(rbind(dcast(TotalPar2012_4[, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches", fun.aggregate = sum)[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[, .(num_unassigned=sum(num_unassigned)), by="offs_site"][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
GenMat2012_4[is.na(GenMat2012_4)] <- 0

#format seasonal genetic parentage
NEMParentageMat <- fread(file="~/parentage/kernel_fitting/1340_loci/input/20200624_parentage_matrix_NEM2012-14ForROMSComp.csv")
NEMParentageMat[, par_site := c(as.character(colnames(NEMParentageMat)), "unassigned")]

#make into a table format to make sure Magbangon names get corrected and all routes are represented
NEMParentage <- melt(NEMParentageMat, id.vars="par_site",  variable.name="offs_site", value.name="n_matches")[, offs_site := as.character(offs_site)]
NEMParentage <- unique(NEMParentage[offs_site %like% "Magbangon", offs_site := "Magbangon"][par_site %like% "Magbangon", par_site := "Magbangon"][ #collapse Magbangon values
            , n_matches := sum(n_matches), by=c("offs_site", "par_site")], by=c("offs_site", "par_site"))

#turn back into a matrix format
GenMatNEM <- rbind(dcast(NEMParentage[par_site != "unassigned"][order(offs_site, par_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches", fun.aggregate= sum)[
    ,-"par_site"],
t(as.matrix(NEMParentage[par_site =="unassigned"][order(offs_site)][, .(n_matches)])), use.names=FALSE)

SWMParentageMat <- fread(file="~/parentage/kernel_fitting/1340_loci/input/20200624_parentage_matrix_SWM2012-14ForROMSComp.csv")
SWMParentageMat[, par_site := c(as.character(colnames(SWMParentageMat)), "unassigned")]

#make into a table format to make sure Magbangon names get corrected and all routes are represented
SWMParentage <- melt(SWMParentageMat, id.vars="par_site",  variable.name="offs_site", value.name="n_matches")[, offs_site := as.character(offs_site)]
SWMParentage <- unique(SWMParentage[offs_site %like% "Magbangon", offs_site := "Magbangon"][par_site %like% "Magbangon", par_site := "Magbangon"][ #collapse Magbangon values
            , n_matches := sum(n_matches), by=c("offs_site", "par_site")], by=c("offs_site", "par_site"))

#turn back into a matrix format
GenMatSWM <- rbind(dcast(SWMParentage[par_site != "unassigned"][order(offs_site, par_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches", fun.aggregate= sum)[
    ,-"par_site"],
t(as.matrix(SWMParentage[par_site =="unassigned"][order(offs_site)][, .(n_matches)])), use.names=FALSE)

fwrite(GenMat2012_4, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrix2012-14ForROMSComp.csv")
fwrite(GenMat2012, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrix2012ForROMSComp.csv")
fwrite(GenMat2013, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrix2013ForROMSComp.csv")
fwrite(GenMat2014, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrix2014ForROMSComp.csv")
fwrite(GenMatNEM, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrixNEM2012-14ForROMSComp.csv")
fwrite(GenMatSWM, file="~/oceanography/script_output/SurveyData/20210625_ParentageMatrixSWM2012-14ForROMSComp.csv")



#to compare, keeping N/S Magbangons as separate sites







