# script to calculate McFadden 1977's pseudo-R2. A subset of the code in bootstrap_data.ipynb

## Load functions
require(data.table)
require(bit64)
require(dplyr)
source("scripts/format_genetic_parentage_matrix.R")
source("scripts/format_genetic_kernel_parentage_matrix.R")
source("scripts/format_biophys_normalized_matrix.R")
source("scripts/format_biophys_parentage_matrix.R")
source("scripts/LL_biophys_bysite.R") # load the LL_biophys_bysite() function

## Load data
load(file= "script_output/SurveyData/2022-01-20_AddDestTables.Rdata")
load(file= "script_output/ROMSDataTables/2021-11-04_SeededParticleInfo.Rdata")
load(file= "script_output/SurveyData/2021-11-04_InputSurveyData.Rdata")
AllRecruitsGen <- fread("empirical_data/genetics/AllFishObsWithPar.csv")
MonsoonRecSampPar <- fread(file="empirical_data/genetics/RecruitsByMonsoon2012-14ForROMSComp.csv")
SimConn <- fread(cmd="unzip -cq script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm08DayPLDCorrectSeasonalReleaseDays.csv.zip", fill=TRUE) # adds two extraneous rows of NAs
SimConn <- SimConn[1:(nrow(SimConn)-2),]
nrow(SimConn) # should be 955263

## Format data
biophys_norm_matrices <- format_biophys_normalized_matrix(SimConn=SimConn, total_release_days=687, AddDestSim=AddDestSim)
biophys_parentage_matrices <- format_biophys_parentage_matrix(SimConn=SimConn, total_release_days=687, AddDestSim=AddDestSim)
sum(biophys_parentage_matrices$MonsoonBiophysMatNEM)+sum(biophys_parentage_matrices$MonsoonBiophysMatSWM) # should be 955263
sum(biophys_parentage_matrices$AnnualBiophysMat2012)+sum(biophys_parentage_matrices$AnnualBiophysMat2013)+sum(biophys_parentage_matrices$AnnualBiophysMat2014)==sum(biophys_parentage_matrices$FullBiophysMat) # should be TRUE
annual_rec_per_dest <- SimConn[,.(annual_rec_per_dest=.N), by=c("dest", "year")]
monsoonal_rec_per_dest <- SimConn[,.(monsoonal_rec_per_dest= .N), by=c("dest", "monsoon")]
SimConn <- left_join(SimConn, annual_rec_per_dest, by=c("dest", "year"))
SimConn <- left_join(SimConn, monsoonal_rec_per_dest, by=c("dest", "monsoon"))

sampled_reefs_vec <- as.matrix(SiteIndexBioPhys[site %in% SurveyData[, site], .(index)])
AllRecruitsGenWithMonsoon <- left_join(AllRecruitsGen, MonsoonRecSampPar[, .(fish_indiv, year, monsoon=season)], by=c("fish_indiv", "year"))[offs_site %like% "Magbangon", offs_site := "Magbangon"][par_site %like% "Magbangon", par_site := "Magbangon"][site %like% "Magbangon", site := "Magbangon"][year %in% c(2012, 2013, 2014)]
genetic_parentage_matrices <- suppressMessages(format_genetic_parentage_matrix(AllRecruitsGenWithMonsoon, SurveyData))

## Set up data vs. model comparisons for each year, monsoon, and multi-period averages (full data)
Data_par2012_bio2012 <- list(BioPhysMat=as.matrix(biophys_norm_matrices$AnnualBiophysMatNorm2012[1:nrow(sampled_reefs_vec),]), 
                             Assignments=genetic_parentage_matrices$GenMat2012[1:nrow(genetic_parentage_matrices$GenMat2012)-1,], 
                             pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                             sampled_reefs_vec=sampled_reefs_vec, 
                             prop_samp_vec=as.matrix(SurveyData[year == 2012, .(prop_anem_samp)]), 
                             unassigned_vec=as.matrix(genetic_parentage_matrices$GenMat2012[nrow(genetic_parentage_matrices$GenMat2012),]))
Data_par2013_bio2013 <- list(BioPhysMat=as.matrix(biophys_norm_matrices$AnnualBiophysMatNorm2013[1:nrow(sampled_reefs_vec),]), 
                             Assignments=genetic_parentage_matrices$GenMat2013[1:nrow(genetic_parentage_matrices$GenMat2013)-1,], 
                             pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                             sampled_reefs_vec=sampled_reefs_vec, 
                             prop_samp_vec=as.matrix(SurveyData[year == 2013,  .(prop_anem_samp)]), 
                             unassigned_vec=as.matrix(genetic_parentage_matrices$GenMat2013[nrow(genetic_parentage_matrices$GenMat2013),]))
Data_par2014_bio2014 <- list(BioPhysMat=as.matrix(biophys_norm_matrices$AnnualBiophysMatNorm2014[1:nrow(sampled_reefs_vec),]), 
                             Assignments=genetic_parentage_matrices$GenMat2014[1:nrow(genetic_parentage_matrices$GenMat2014)-1,], 
                             pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                             sampled_reefs_vec=sampled_reefs_vec, 
                             prop_samp_vec=as.matrix(SurveyData[year == 2014,  .(prop_anem_samp)]), 
                             unassigned_vec=as.matrix(genetic_parentage_matrices$GenMat2014[nrow(genetic_parentage_matrices$GenMat2014),]))
Data_par2012_4_bio2012_4 <- list(BioPhysMat=as.matrix(biophys_norm_matrices$FullBiophysMatNorm[1:nrow(sampled_reefs_vec),]), 
                                 Assignments=genetic_parentage_matrices$GenMat2012_4[1:nrow(genetic_parentage_matrices$GenMat2012_4)-1,], 
                                 pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                                 sampled_reefs_vec=sampled_reefs_vec, 
                                 prop_samp_vec=as.matrix(SurveyData[year == 2014,  .(prop_anem_samp)]), 
                                 unassigned_vec=as.matrix(genetic_parentage_matrices$GenMat2012_4[nrow(genetic_parentage_matrices$GenMat2012_4),]))
Data_parNEM_bioNEM <- list(BioPhysMat=as.matrix(biophys_norm_matrices$MonsoonBiophysMatNormNEM[1:nrow(sampled_reefs_vec),]), 
                           Assignments=as.matrix(genetic_parentage_matrices$GenMatNEM[1:nrow(genetic_parentage_matrices$GenMatNEM)-1,]), 
                           pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                           sampled_reefs_vec=sampled_reefs_vec, 
                           prop_samp_vec=as.matrix(SurveyData[year == 2014,  .(prop_anem_samp)]), 
                           unassigned_vec=as.matrix(genetic_parentage_matrices$GenMatNEM[nrow(genetic_parentage_matrices$GenMatNEM),]))
Data_parSWM_bioSWM <- list(BioPhysMat=as.matrix(biophys_norm_matrices$MonsoonBiophysMatNormSWM[1:nrow(sampled_reefs_vec),]), 
                           Assignments=as.matrix(genetic_parentage_matrices$GenMatSWM[1:nrow(genetic_parentage_matrices$GenMatSWM)-1,]), 
                           pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                           sampled_reefs_vec=sampled_reefs_vec, 
                           prop_samp_vec=as.matrix(SurveyData[year == 2014,  .(prop_anem_samp)]), 
                           unassigned_vec=as.matrix(genetic_parentage_matrices$GenMatSWM[nrow(genetic_parentage_matrices$GenMatSWM),]))
Data_parNEMSWM_bio2012_4 <- list(BioPhysMat=as.matrix(biophys_norm_matrices$FullBiophysMatNorm[1:nrow(sampled_reefs_vec),]), 
                                 Assignments=(as.matrix(genetic_parentage_matrices$GenMatNEM[1:nrow(genetic_parentage_matrices$GenMatNEM)-1,])+ as.matrix(genetic_parentage_matrices$GenMatSWM[1:nrow(genetic_parentage_matrices$GenMatSWM)-1,])), 
                                 pop_size_vec=as.matrix(SurveyData[,.(avg_num_females=mean(num_females, na.rm = TRUE)), by=site][order(site)][, .(avg_num_females)]), 
                                 sampled_reefs_vec=sampled_reefs_vec, 
                                 prop_samp_vec=as.matrix(SurveyData[year == 2014,  .(prop_anem_samp)]), 
                                 unassigned_vec=as.matrix(genetic_parentage_matrices$GenMatNEM[nrow(genetic_parentage_matrices$GenMatNEM),])+as.matrix(genetic_parentage_matrices$GenMatSWM[nrow(genetic_parentage_matrices$GenMatSWM),]))


## Log likelihoods for the full data, separately by site
LL_annual_annual_2012 <- LL_biophys_bysite(Data_par2012_bio2012) # year-specific
LL_annual_annual_2013 <- LL_biophys_bysite(Data_par2013_bio2013) # year-specific
LL_annual_annual_2014 <- LL_biophys_bysite(Data_par2014_bio2014) # year-specific
LL_total_total <- LL_biophys_bysite(Data_par2012_4_bio2012_4) # multi-annual average
LL_monsoonal_monsoonal_NEM <- LL_biophys_bysite(Data_parNEM_bioNEM)
LL_monsoonal_monsoonal_SWM <- LL_biophys_bysite(Data_parSWM_bioSWM) # monsoon-specific
LL_monsoonal_total <- LL_biophys_bysite(Data_parNEMSWM_bio2012_4) # cross-monsoon average

## LL differences by site
LL_byyear = LL_annual_annual_2012 + LL_annual_annual_2013 + LL_annual_annual_2014
LLdiff_ave_vs_annual = LL_byyear - LL_total_total # positive values indidate sites where the year-specific sims are better than the average
LLdiff_ave_vs_monsoon = LL_monsoonal_monsoonal_NEM + LL_monsoonal_monsoonal_SWM- LL_monsoonal_total # positive values indidate sites where the monsoon-specific sims are better than the average


## barplot of negative log-likelihoods. lower is better fit. normalized to min 0
toplot <- matrix(data = c(LLdiff_ave_vs_annual,
                          LLdiff_ave_vs_monsoon), 
                 ncol=18, byrow = TRUE, dimnames = list(c('Annual', 'Monsoonal'), SiteIndexBioPhys$site[sampled_reefs_vec]))
toplot = toplot[,c('Palanas', 'Wangag', 'Magbangon', 'Cabatoan', 'Caridad Cemetery', 'Caridad Proper', 'Hicgop South', 'Sitio Tugas', 'Elementary School', 'Sitio Lonas', 'San Agustin', 'Poroc San Flower', 'Poroc Rose', 'Visca', 'Gabas', 'Tamakin Dacot', 'Haina', 'Sitio Baybayon')]

# plot 
png(filename = 'script_output/plots/loglikelihood_differences_bysite.png', width=6, height = 6, units = 'in', res = 300)
par(mar=c(10,4,2,2)); barplot(toplot, beside = TRUE, las = 2, ylab = 'Difference in log likelihood', col = c('gray90', 'gray30'))
legend(1,1200, legend = c('Annual', 'Monsoonal'), fill = c('gray90', 'gray30'))
dev.off()