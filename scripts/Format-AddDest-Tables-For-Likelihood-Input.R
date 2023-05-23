SimConn <- fread(file="script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm08DayPLD.csv")[year %in% c(2012, 2013, 2014) & dest != "CAI"][#filter out CAI as a destination for now, not very well spatially defined
   ,`:=`(daily_particles_rec=as.numeric(daily_particles_rec),daily_particles_released=as.numeric(daily_particles_released))][
    source == "Other" | source == "CAI", source:="unknown"][ #aggregate the CAI and other sites as "unknown for likelihood functions"
    , daily_particles_rec := sum(daily_particles_rec), by=c("date", "source", "dest")][
   , daily_particles_released :=sum(as.numeric(daily_particles_released)), by=c("date", "source")]
   
   
AddDestSim <- rbind(rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "NEM"],
    rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "SWM"])

##add in sites that aren't represented in every year by the genetics for the biophysical PARENTAGE matrix used to fit the biophysical kernel
AddDestGenAnnual <- rbindlist(list(unique(cbind(SurveyData[year==2012 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2012 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
unique(cbind(SurveyData[year==2013 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2013 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
unique(cbind(SurveyData[year==2014 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2014 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))

AddDestGen <- rbind(cbind(AddDestGenAnnual, data.table(monsoon = rep("NEM",nrow(AddDestGenAnnual)))), cbind(AddDestGenAnnual, data.table(monsoon = rep("SWM",nrow(AddDestGenAnnual)))))

save(AddDestSim,  AddDestGen, file="script_output/SurveyData/for_likelihood_functions/2022-01-20_AddDestTables.Rdata")