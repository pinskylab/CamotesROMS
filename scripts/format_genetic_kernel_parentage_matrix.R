format_genetic_kernel_parentage_matrix <- function(AllRecruitsGenWithMonsoon, SurveyData){
	

#add in sites that aren't represented in every year
AddDestGen <- rbindlist(list(unique(cbind(SurveyData[year==2012 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2012 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
unique(cbind(SurveyData[year==2013 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2013 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
unique(cbind(SurveyData[year==2014 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2014 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))


#format genetic parentage matrices for each year and all years combined
KernelGenMat2012 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2012 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2012, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2012 & matched_offs=="N", .N, by=site], AddDestGen[year==2012, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMat2012[is.na(KernelGenMat2012)] <- 0

KernelGenMat2013 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2013 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2013, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2013 & matched_offs=="N", .N, by=site], AddDestGen[year==2013, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMat2013[is.na(KernelGenMat2013)] <- 0

KernelGenMat2014 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2014 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2014 & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMat2014[is.na(KernelGenMat2014)] <- 0

#all year
KernelGenMat2012_4 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[ matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMat2012_4[is.na(KernelGenMat2012_4)] <- 0



#for monsoon seasons
KernelGenMatSWM <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[monsoon=="SWM" & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[monsoon=="SWM" & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMatSWM[is.na(KernelGenMatSWM)] <- 0

KernelGenMatNEM <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[monsoon=="NEM" & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
t(as.matrix(full_join(AllRecruitsGenWithMonsoon[monsoon=="NEM" & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

KernelGenMatNEM[is.na(KernelGenMatNEM)] <- 0

kernel_matrices <- list(KernelGenMat2012=KernelGenMat2012, KernelGenMat2013=KernelGenMat2013, KernelGenMat2014=KernelGenMat2014, KernelGenMatNEM=KernelGenMatNEM, KernelGenMatSWM=KernelGenMatSWM, KernelGenMat2012_4=KernelGenMat2012_4)
return(kernel_matrices)
}
#fwrite(KernelGenMat2012_4, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2012-14ForROMSComp.csv")
#fwrite(KernelGenMat2012, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2012ForROMSComp.csv")
#fwrite(KernelGenMat2013, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2013ForROMSComp.csv")
#fwrite(KernelGenMat2014, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrix2014ForROMSComp.csv")
#fwrite(KernelGenMatNEM, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrixNEM2012-14ForROMSComp.csv")
#fwrite(KernelGenMatSWM, file="~/oceanography/script_output/SurveyData/20210701_KernelParentageMatrixSWM2012-14ForROMSComp.csv")
