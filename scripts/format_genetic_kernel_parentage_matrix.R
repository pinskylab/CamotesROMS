format_genetic_kernel_parentage_matrix <- function(TotalPar2012_4, MonsoonRecSampPar, SurveyData){
	
#add in all the sampled sites and numbers of recruits sampled at each site
TotalParInt <- unique(TotalPar2012_4[, num_matches := sum(n_matches), by=c("year", "offs_site")][, -"par_site"], by=c("year", "offs_site"))[, -"n_matches"]
#sum(TotalParInt$num_matches) #should be 37

TotalUnassigned2012_4 <- TotalParInt[SurveyData[year %in% c(2012, 2013, 2014) & prop_anem_samp > 0 ], on=.(offs_site=site, year)]
TotalUnassigned2012_4 [is.na(TotalUnassigned2012_4 )] <- 0

#sum(TotalUnassigned2012_4 $n_offs_gen) #should be 394
#sum(TotalUnassigned2012_4 $num_matches) #should be 37
#nrow(unique(TotalUnassigned2012_4 , by="offs_site"))#should be 18, so that every site is represented so that the years are all 18*18 sites
TotalUnassigned2012_4 <- TotalUnassigned2012_4[, num_unassigned := n_offs_gen-num_matches, by=c("year", "offs_site")][
    , .(offs_site, year, num_unassigned)]
#sum(TotalUnassigned2012_4$num_unassigned) #should be 357

#add in sites that aren't represented in every year
AddDestGen <- rbindlist(list(unique(cbind(SurveyData[year==2012 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2012 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
unique(cbind(SurveyData[year==2013 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2013 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
unique(cbind(SurveyData[year==2014 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2014 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))

TotalPar2012_4 <- rbind(AddDestGen[!TotalPar2012_4, on =.(par_site, offs_site, year)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , n_matches:=0 ], TotalPar2012_4[,-"num_matches"])  #add the parentage column, add back into the parentage table but drop the num_matches column that's a summary column I used to make the unassigned table 
#sum(TotalPar2012_4$n_matches)#should be 37 still

#format genetic parentage matrices for each year and all years combined
KernelGenMat2012 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2012, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2012][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
KernelGenMat2012[is.na(KernelGenMat2012)] <- 0

KernelGenMat2013 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2013, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2013][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
KernelGenMat2013[is.na(KernelGenMat2013)] <- 0
#for monsoon seasons

#sum(unique(MonsoonRecSamp, by=c("offs_site","monsoon"))[, num_gen]) #should be 256MonsoonRecSamp <- unique(MonsoonRecSamp, by=c("site", "par_site", "monsoon", "year"))[
MonsoonRecSampPar <- unique(MonsoonRecSampPar, by=c("offs_site", "par_site", "monsoon"))[
    is.na(par_site), n_matches := 0, by=c("offs_site", "par_site", "monsoon")][
    , .(offs_site, par_site, monsoon, num_gen, n_matches)]
#sum(unique(MonsoonRecSampPar, by=c("offs_site","monsoon"))[, num_gen]) #should be 256
#sum(unique(MonsoonRecSampPar, by=c("offs_site", "par_site","monsoon"))[, n_matches]) #should be 17

#total recruits genotyped in each year and season at each site
MonsoonRecSamp <- unique(MonsoonRecSampPar, by=c("offs_site","monsoon"))[,-c("par_site", "n_matches")] #should be 256

#total parentage matches along each route in each year and season
ParMonsoon <- MonsoonRecSampPar[n_matches >0, -"num_gen"]
#print a summary table to check 
#ParMonsoon[, .(n_matches=sum(n_matches)) , by=c("monsoon")]
#sum(AnnualParMonsoon$n_matches) #should be 17

MonsoonParInt <- ParMonsoon[, .(num_matches = sum(n_matches)), by=c("monsoon", "offs_site")]#sum(MonsoonParInt$num_matches) #should be 17

#join all genotyped and assigned to get unassigned

MonsoonUnassigned <- MonsoonParInt[MonsoonRecSamp, on=.(offs_site, monsoon)]
MonsoonUnassigned[is.na(MonsoonUnassigned)] <- 0
MonsoonUnassigned[, num_unassigned := num_gen-num_matches]
#sum(MonsoonUnassigned$num_gen) #should be 256
#sum(MonsoonUnassigned$num_matches) #should be 17
MonsoonUnassigned <- MonsoonUnassigned[, .(offs_site,monsoon, num_unassigned)]

#AnnualParMonsoon #for assigned in matrix
#MonsoonUnassigned #for unassigned in matrix

#add in sites that aren't represented in every year
AddDestGen <- cbind(PropSamp[prop_anem_samp > 0 & year==2014, .(site)][, .(offs_site=site)], PropSamp[prop_anem_samp > 0 & year==2014, .(site)][, .(par_site=site)])

AddDestGenMonsoon <- rbind(cbind(AddDestGen, data.table(monsoon = rep("NEM",nrow(AddDestGen)))), cbind(AddDestGen, data.table(monsoon = rep("SWM",nrow(AddDestGen)))))

ParMonsoon <- rbind(AddDestGenMonsoon[!ParMonsoon, on =.(par_site, offs_site, monsoon)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , n_matches:=0 ], ParMonsoon)  #add the parentage column, add back into the parentage table but drop the num_matches column that's a summary column I used to make the unassigned table 
#ParMonsoon <- unique(ParMonsoon[, n_matches := sum(n_matches), by=c("par_site", "offs_site", "monsoon")], by=c("par_site", "offs_site", "monsoon"))
#sum(ParMonsoon$n_matches)#should be 17 still

MonsoonUnassigned <- rbind(AddDestGenMonsoon[!MonsoonUnassigned, on =.( offs_site, monsoon)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
    , num_unassigned:=0 ][, -"par_site"], MonsoonUnassigned)  #add the parentage column, add back into the parentage table but drop the num_matches column that's a summary column I used to make the unassigned table 
#MonsoonUnassigned <- unique(MonsoonUnassigned[, num_unassigned := sum(num_unassigned), by=c("offs_site", "monsoon")], by=c("offs_site", "monsoon"))

#sum(MonsoonUnassigned$num_unassigned)#should be 256-17=239

#make into parentage matrices
KernelGenMatNEM <- rbind(dcast(ParMonsoon[monsoon=="NEM"][, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(MonsoonUnassigned[monsoon=="NEM"][order(offs_site)][, .(num_unassigned)])), use.names=FALSE)
KernelGenMatNEM[is.na(KernelGenMatNEM)] <- 0
#dim(KernelGenMatNEM)
#sum(KernelGenMatNEM)

KernelGenMatSWM <- rbind(dcast(ParMonsoon[monsoon=="SWM"][, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(MonsoonUnassigned[monsoon=="SWM"][order(offs_site)][, .(num_unassigned)])), use.names=FALSE)
KernelGenMatSWM[is.na(KernelGenMatSWM)] <- 0
KernelGenMat2014 <- as.matrix(rbind(dcast(TotalPar2012_4[year==2014, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches")[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[year==2014][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
KernelGenMat2014[is.na(KernelGenMat2014)] <- 0

KernelGenMat2012_4 <- as.matrix(rbind(dcast(TotalPar2012_4[, .(par_site, offs_site, n_matches)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(par_site, offs_site)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , par_site ~ offs_site, value.var="n_matches", fun.aggregate = sum)[
    ,-"par_site"], #remove the source column after casting
      t(as.matrix(TotalUnassigned2012_4[, .(num_unassigned=sum(num_unassigned)), by="offs_site"][order(offs_site)][, .(num_unassigned)])), use.names=FALSE))
KernelGenMat2012_4[is.na(KernelGenMat2012_4)] <- 0

kernel_matrices <- list(KernelGenMat2012=KernelGenMat2012, KernelGenMat2013=KernelGenMat2013, KernelGenMat2014=KernelGenMat2014, KernelGenMatNEM=KernelGenMatNEM, KernelGenMatSWM=KernelGenMatSWM, KernelGenMat2012_4=KernelGenMat2012_4)
return(kernel_matrices)
}
#fwrite(KernelGenMat2012_4, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2012-14ForROMSComp.csv")
#fwrite(KernelGenMat2012, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2012ForROMSComp.csv")
#fwrite(KernelGenMat2013, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2013ForROMSComp.csv")
#fwrite(KernelGenMat2014, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2014ForROMSComp.csv")
#fwrite(KernelGenMatNEM, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrixNEM2012-14ForROMSComp.csv")
#fwrite(KernelGenMatSWM, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrixSWM2012-14ForROMSComp.csv")
