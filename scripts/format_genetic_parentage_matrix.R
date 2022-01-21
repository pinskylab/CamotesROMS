format_genetic_parentage_matrix <- function(AllRecruitsGenWithMonsoon, SurveyData){
	
	#add in sites that aren't represented in every year
	AddDestGen <- rbindlist(list(unique(cbind(SurveyData[year==2012 & prop_anem_samp >0 ][, .(offs_site=site)], 
	                  SurveyData[year==2012 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
	unique(cbind(SurveyData[year==2013 & prop_anem_samp >0 ][, .(offs_site=site)], 
	                  SurveyData[year==2013 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
	unique(cbind(SurveyData[year==2014 ][, .(offs_site=site)], 
	                  SurveyData[year==2014 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))


	#format genetic parentage matrices for each year and all years combined
	GenMat2012 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2012 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2012 & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMat2012[is.na(GenMat2012)] <- 0

	GenMat2013 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2013 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2013 & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMat2013[is.na(GenMat2013)] <- 0

	GenMat2014 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[year==2014 & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[year==2014 & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMat2014[is.na(GenMat2014)] <- 0

	#all year
	GenMat2012_4 <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[ matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMat2012_4[is.na(GenMat2012_4)] <- 0



	#for monsoon seasons
	GenMatSWM <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[monsoon=="SWM" & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014][, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[monsoon=="SWM" & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMatSWM[is.na(GenMatSWM)] <- 0

	GenMatNEM <- rbind(as.matrix(dcast(full_join(AllRecruitsGenWithMonsoon[monsoon=="NEM" & matched_offs=="Y", .N, by=c("par_site", "offs_site")], AddDestGen[year==2014][, -"year"])[order(par_site, offs_site)], par_site ~ offs_site, value.var="N")[, -"par_site"]),  #[is.na(N), N :=0]
	t(as.matrix(full_join(AllRecruitsGenWithMonsoon[monsoon=="NEM" & matched_offs=="N", .N, by=site], AddDestGen[year==2014, .(site=offs_site)])[order(site)][, -"site"])))

	GenMatNEM[is.na(GenMatNEM)] <- 0
	
	
matrices <- list(GenMat2012=GenMat2012, GenMat2013=GenMat2013, GenMat2014=GenMat2014, GenMatNEM=GenMatNEM, GenMatSWM=GenMatSWM, GenMat2012_4=GenMat2012_4)
return(matrices)
}
#fwrite(GenMat2012_4, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrix2012-14ForROMSComp.csv")
#fwrite(GenMat2012, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrix2012ForROMSComp.csv")
#fwrite(GenMat2013, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrix2013ForROMSComp.csv")
#fwrite(GenMat2014, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrix2014ForROMSComp.csv")
#fwrite(GenMatNEM, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrixNEM2012-14ForROMSComp.csv")
#fwrite(GenMatSWM, file="~/oceanography/script_output/SurveyData/20210701_ParentageMatrixSWM2012-14ForROMSComp.csv")
