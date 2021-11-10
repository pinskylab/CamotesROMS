SimConn <- fread(file="~/oceanography/script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm08DayPLD.csv")[year %in% c(2012, 2013, 2014) & dest != "CAI"] #filter out CAI as a destination for now, not very well spatially defined

##each year will require a different set of survey data, so make a list of each and index by site for fast look up
#SampledTable <- SurveyData[prop_anem_samp >0, c("year", "site")]#previously named PropSampTable
#
##make sure all sampled sites are represented when joining the survey data to the sampled simulation- this chunk has the tables to add to a subsampled particle table. no need for the full
#SampTable <- rbind(SurveyData[prop_anem_samp >0 & year %in% c(2012, 2013, 2014), c("year", "site")][, .(source=site, dest=site, year=year)][, #will join to the simulated sampling table by source and dest, so make those each a column from site and preserve the year variable as a key
#     c("year", "source", "dest")][, monsoon := "NEM"], SurveyData[prop_anem_samp >0 & year %in% c(2012, 2013, 2014), c("year", "site")][, .(source=site, dest=site, year=year)][, #will join to the simulated sampling table by source and dest, so make those each a column from site and preserve the year variable as a key
#     c("year", "source", "dest")][, monsoon := "SWM"])
#
#UnqSurvey <- unique(SampTable, by=c("source", "dest", "year", "monsoon"))#add in the diff Monsoon seasons so there are complete parentage matrices later
#
#make sure all possible routes are represented
AddDestSim <- rbind(rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "NEM"],
    rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "SWM"])

AddDestAllYearsSim <- unique(AddDestSim, by=c("source","dest"))[, -c("year","monsoon")] 

AddDestAnnualSim <- unique(AddDestSim, by=c("source","dest", "year"))[, -"monsoon"] 

AddDestMonsoonSim <- unique(AddDestSim[, -"year"], by=c("source", "dest", "monsoon"))
AddDestAnnualMonsoonSim <- unique(rbind(rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2011],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "NEM"],
    rbindlist(list(unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2011],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2012],
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2013], 
        unique(SimConn[, .(source, dest)], by=c("source", "dest"))[, year := 2014]))[, monsoon := "SWM"])
    , by=c("source", "dest","year", "monsoon"))


##add in sites that aren't represented in every year by the genetics for the biophysical PARENTAGE matrix used to fit the biophysical kernel
AddDestGenAnnual <- rbindlist(list(unique(cbind(SurveyData[year==2012 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2012 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2012],  #what destinations were sampled, for use with unassigned table
unique(cbind(SurveyData[year==2013 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2013 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2013],
unique(cbind(SurveyData[year==2014 & prop_anem_samp >0 ][, .(offs_site=site)], 
                  SurveyData[year==2014 & prop_anem_samp >0 ][ , .(par_site=site)]), by=c("offs_site", "par_site"))[, year := 2014]))

AddDestGen <- cbind(PropSamp[prop_anem_samp > 0 & year==2014, .(site)][, .(offs_site=site)], PropSamp[prop_anem_samp > 0 & year==2014, .(site)][, .(par_site=site)])

AddDestGenMonsoon <- rbind(cbind(AddDestGen, data.table(monsoon = rep("NEM",nrow(AddDestGen)))), cbind(AddDestGen, data.table(monsoon = rep("SWM",nrow(AddDestGen)))))

save(AddDestSim,  AddDestAllYearsSim, AddDestAnnualSim, AddDestMonsoonSim, AddDestAnnualMonsoonSim, AddDestGenAnnual, AddDestGen, AddDestGenMonsoon, file="~/oceanography/script_output/SurveyData/for_likelihood_functions/2021-11-04_AddDestTables.Rdata")