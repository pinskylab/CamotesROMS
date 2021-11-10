#make a parentage matrix for the whole biophysical results- CAI included with Other as "unknown", use this to fit a kernel
format_biophys_parentage_matrix <- function(SimConn, AddDestGenAnnual, total_release_days, AddDestAllYearsSim){
    
	AllYearsRec <- SimConn[, .(total_particles_rec = .N), by = c("source", 
        "dest")]
    AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released) * 
        as.numeric(total_release_days)), by = c("source")], by = "source")
    AllYearsRecInt <- rbind(AddDestAllYearsSim[!AllYearsRec, 
        on = .(source, dest)][, `:=`(total_particles_rec, 0)], 
        AllYearsRec)
    AllYearsReleaseInt <- rbind(AddDestAllYearsSim[!AllYearsRelease, 
        on = .(source)][, `:=`(total_particles_released, 0)][, 
        -"dest"], AllYearsRelease)
    AllYearsRec <- unique(AllYearsRecInt[source == "Other" | 
        source == "CAI", `:=`(source, "unknown")][, `:=`(total_particles_rec, 
        sum(total_particles_rec)), by = c("source", "dest")], 
        by = c("source", "dest"))
    AllYearsRelease <- unique(AllYearsReleaseInt[source == "Other" | 
        source == "CAI", `:=`(source, "unknown")][, `:=`(total_particles_released, 
        sum(total_particles_released)), by = "source"], by = "source")
    AllYearsNormConn <- AllYearsRelease[AllYearsRec, on = "source"][, 
        `:=`(source_norm_rec, total_particles_rec/total_particles_released)]
    FullBiophysMat <- as.matrix(rbind(dcast(AllYearsNormConn[source %in% 
        AddDestGenAnnual[, par_site] & dest %in% AddDestGenAnnual[, 
        offs_site]][source != "unknown", .(source, dest, total_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "total_particles_rec")[, 
        -"source"], dcast(AllYearsNormConn[dest %in% AddDestGenAnnual[, 
        offs_site]][source == "unknown", .(source, dest, total_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "total_particles_rec")[, 
        -"source"]))
    FullBiophysMat[is.na(FullBiophysMat)] <- 0
   # dim(FullBiophysMat)
    AnnualRec <- SimConn[, .(annual_particles_rec = .N), by = c("source", 
        "dest", "year")]
    AnnualRelease <- unique(SimConn[, .(annual_particles_released = as.numeric(daily_particles_released) * 
        as.numeric(num_release_days_annual)), by = c("source", 
        "year")], by = c("source", "year"))
    AnnualRecInt <- rbind(AddDestAnnualSim[!AnnualRec, on = .(source, 
        dest, year)][, `:=`(annual_particles_rec, 0)], AnnualRec)
    AnnualReleaseInt <- rbind(AddDestAnnualSim[!AnnualRelease, 
        on = .(source, year)][, `:=`(annual_particles_released, 
        0)][, -"dest"], AnnualRelease)
    AnnualRec <- unique(AnnualRecInt[source == "Other" | source == 
        "CAI", `:=`(source, "unknown")][, `:=`(annual_particles_rec, 
        sum(annual_particles_rec)), by = c("source", "dest", 
        "year")], by = c("source", "dest", "year"))
    AnnualRelease <- unique(AnnualReleaseInt[source == "Other" | 
        source == "CAI", `:=`(source, "unknown")][, `:=`(annual_particles_released, 
        sum(annual_particles_released)), by = c("source", "year")], 
        by = c("source", "year"))
    AnnualNormConn <- AnnualRelease[AnnualRec, on = c("source", 
        "year")][, `:=`(source_norm_rec, annual_particles_rec/annual_particles_released)]
    #summary(AnnualNormConn[, .(sum = sum(source_norm_rec)), by = c("year", 
    #    "source")][, sum])
    AnnualBiophysMat2012 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
        2012][source %in% AddDestGenAnnual[year == 2012, par_site] & 
        dest %in% AddDestGenAnnual[year == 2012, offs_site]][source != 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"], dcast(AnnualNormConn[year == 2012][dest %in% 
        AddDestGenAnnual[year == 2012, offs_site]][source == 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"]))
    AnnualBiophysMat2012[is.na(AnnualBiophysMat2012)] <- 0
    AnnualBiophysMat2013 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
        2013][source %in% AddDestGenAnnual[year == 2013, par_site] & 
        dest %in% AddDestGenAnnual[year == 2013, offs_site]][source != 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"], dcast(AnnualNormConn[year == 2013][dest %in% 
        AddDestGenAnnual[year == 2013, offs_site]][source == 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"]))
    AnnualBiophysMat2013[is.na(AnnualBiophysMat2013)] <- 0
    #dim(AnnualBiophysMat2013)
    AnnualBiophysMat2014 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
        2014][source %in% AddDestGenAnnual[year == 2014, par_site] & 
        dest %in% AddDestGenAnnual[year == 2014, offs_site]][source != 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"], dcast(AnnualNormConn[year == 2014][dest %in% 
        AddDestGenAnnual[year == 2014, offs_site]][source == 
        "unknown", .(source, dest, annual_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "annual_particles_rec")[, 
        -"source"]))
    AnnualBiophysMat2014[is.na(AnnualBiophysMat2014)] <- 0
    #dim(AnnualBiophysMat2014)
    
 	MonsoonRec <- SimConn[, .(monsoon_particles_rec = .N), by = c("source", 
        "dest", "monsoon")]
    MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released) * 
        as.numeric(num_release_days_seasonal)), by = c("source", 
        "monsoon")], by = c("source", "monsoon"))
    MonsoonRec <- SimConn[, .(monsoon_particles_rec = .N), by = c("source", 
        "dest", "monsoon")]
    MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released) * 
        as.numeric(num_release_days_seasonal)), by = c("source", 
        "monsoon")], by = c("source", "monsoon"))
    MonsoonRecInt <- rbind(AddDestMonsoonSim[!MonsoonRec, on = .(source, 
        dest, monsoon)][, `:=`(monsoon_particles_rec, 0)], MonsoonRec)
    MonsoonReleaseInt <- rbind(AddDestMonsoonSim[!MonsoonRelease, 
        on = .(source, monsoon)][, `:=`(monsoon_particles_released, 
        0)][, -"dest"], MonsoonRelease)
    MonsoonRec <- unique(MonsoonRecInt[source == "Other" | source == 
        "CAI", `:=`(source, "unknown")][, `:=`(monsoon_particles_rec, 
        sum(monsoon_particles_rec)), by = c("source", "dest", 
        "monsoon")], by = c("source", "dest", "monsoon"))
    MonsoonRelease <- unique(MonsoonReleaseInt[source == "Other" | 
        source == "CAI", `:=`(source, "unknown")][, `:=`(monsoon_particles_released, 
        sum(monsoon_particles_released)), by = c("source", "monsoon")], 
        by = c("source", "monsoon"))
    MonsoonNormConn <- MonsoonRelease[MonsoonRec, on = c("source", 
        "monsoon")][, `:=`(source_norm_rec, monsoon_particles_rec/monsoon_particles_released)]
    MonsoonBiophysMatNEM <- as.matrix(rbind(dcast(MonsoonNormConn[source %in% 
        AddDestGenAnnual[year == 2014, par_site] & dest %in% 
        AddDestGenAnnual[year == 2014, offs_site]][monsoon == 
        "NEM" & source != "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
        -"source"], dcast(MonsoonNormConn[monsoon == "NEM"][dest %in% 
        AddDestGenAnnual[year == 2014, offs_site]][source == 
        "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
        -"source"]))
    MonsoonBiophysMatNEM[is.na(MonsoonBiophysMatNEM)] <- 0
    #dim(MonsoonBiophysMatNEM)
    MonsoonBiophysMatSWM <- as.matrix(rbind(dcast(MonsoonNormConn[source %in% 
        AddDestGenAnnual[year == 2014, par_site] & dest %in% 
        AddDestGenAnnual[year == 2014, offs_site]][monsoon == 
        "SWM" & source != "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
        -"source"], dcast(MonsoonNormConn[monsoon == "SWM"][dest %in% 
        AddDestGenAnnual[year == 2014, offs_site]][source == 
        "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
        dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
        -"source"]))
    MonsoonBiophysMatSWM[is.na(MonsoonBiophysMatSWM)] <- 0
    ##dim(MonsoonBiophysMatSWM)
    #matrices <- list(AnnualBiophysMat2012 = AnnualBiophysMat2012, 
    #    AnnualBiophysMat2013 = AnnualBiophysMat2013, AnnualBiophysMat2014 = AnnualBiophysMat2014, 
    #    FullBiophysMat = FullBiophysMat, MonsoonBiophysMatNEM = MonsoonBiophysMatNEM, 
    #    MonsoonBiophysMatSWM = MonsoonBiophysMatSWM)
#AllYearsRec <- SimConn[ , .(total_particles_rec = .N), by= c("source","dest")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
#	AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released)*as.numeric(total_release_days)), by= c("source")], by="source") #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs
#
#	AllYearsRecInt <- rbind(AddDestAllYearsSim[!AllYearsRec, on =.(source, dest)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
#	    , total_particles_rec:=0 ], AllYearsRec) 
#
#	AllYearsReleaseInt <- rbind(AddDestAllYearsSim[!AllYearsRelease, on =.(source)][ #what combos are not appearing because we didn't sample particles, but the route is possible based on our survey sampling
#	    , total_particles_released:=0 ][, -"dest"], AllYearsRelease) 
#
#	#collapse CAI and "Other" together into "unknown" before making the normalized matrices
#	AllYearsRec <- unique(AllYearsRecInt[source== "Other" | source == "CAI", source := "unknown"][
#	    , total_particles_rec := sum(total_particles_rec),  by=c("source", "dest")], by=c("source", "dest"))
#	AllYearsRelease <- unique(AllYearsReleaseInt[source== "Other" | source == "CAI", source := "unknown"][
#	    , total_particles_released := sum(total_particles_released), by="source"], by="source")
#
#	#join recruited and released tables together and make a column for the normalized values
#	AllYearsNormConn <-  AllYearsRelease[AllYearsRec, on="source"][
#	    , source_norm_rec := total_particles_rec/total_particles_released]
#
#
#FullBiophysMat <- as.matrix(rbind(dcast(AllYearsNormConn[source %in% AddDestGenAnnual[, par_site] & dest %in% AddDestGenAnnual[, offs_site]][source !="unknown", .(source, dest, total_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="total_particles_rec")[#use sum to count the matches for each id variable combo, that populated the cells of the matrix
#    ,-"source"], #remove the source column after casting
#      dcast(AllYearsNormConn[dest %in% AddDestGenAnnual[, offs_site]][source  == "unknown", .(source, dest, total_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#        order(source, dest)], source ~ dest, value.var="total_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#FullBiophysMat[is.na(FullBiophysMat)] <- 0 #change NAs to zerosdim(FullBiophysMat)
#dim(FullBiophysMat)
##parentage matrix by years and monsoons
#AnnualBiophysMat2012 <- as.matrix(rbind(dcast(AnnualNormConn[year==2012][source %in% AddDestGenAnnual[year==2012, par_site] & dest %in% AddDestGenAnnual[year==2012, offs_site]][source !="unknown", .(source, dest, annual_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="annual_particles_rec")[#use sum to count the matches for each id variable combo, that populated the cells of the matrix
#    ,-"source"], #remove the source column after casting
#dcast(AnnualNormConn[year==2012][dest %in% AddDestGenAnnual[year==2012, offs_site]][source  == "unknown", .(source, dest, annual_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#        order(source, dest)], source ~ dest, value.var="annual_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#AnnualBiophysMat2012[is.na(AnnualBiophysMat2012)] <- 0 #change NAs to zerosdim(AnnualBiophysMat2012)
#
#AnnualBiophysMat2013 <- as.matrix(rbind(dcast(AnnualNormConn[year==2013][source %in% AddDestGenAnnual[year==2013, par_site] & dest %in% AddDestGenAnnual[year==2013, offs_site]][source !="unknown", .(source, dest, annual_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="annual_particles_rec")[#use sum to count the matches for each id variable combo, that populated the cells of the matrix
#    ,-"source"], #remove the source column after casting
#dcast(AnnualNormConn[year==2013][dest %in% AddDestGenAnnual[year==2013, offs_site]][source  == "unknown", .(source, dest, annual_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#        order(source, dest)], source ~ dest, value.var="annual_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#AnnualBiophysMat2013[is.na(AnnualBiophysMat2013)] <- 0 #change NAs to zeros
#dim(AnnualBiophysMat2013)
#
#AnnualBiophysMat2014 <- as.matrix(rbind(dcast(AnnualNormConn[year==2014][source %in% AddDestGenAnnual[year==2014, par_site] & dest %in% AddDestGenAnnual[year==2014, offs_site]][source !="unknown", .(source, dest, annual_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="annual_particles_rec")[#use sum to count the matches for each id variable combo, that populated the cells of the matrix
#    ,-"source"], #remove the source column after casting
#dcast(AnnualNormConn[year==2014][dest %in% AddDestGenAnnual[year==2014, offs_site]][source  == "unknown", .(source, dest, annual_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#        order(source, dest)], source ~ dest, value.var="annual_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#AnnualBiophysMat2014[is.na(AnnualBiophysMat2014)] <- 0 #change NAs to zeros
#dim(AnnualBiophysMat2014)
#
#MonsoonBiophysMatNEM <- as.matrix(rbind(dcast(MonsoonNormConn[source %in% AddDestGenAnnual[year==2014, par_site] & dest %in% AddDestGenAnnual[year==2014, offs_site]][monsoon=="NEM" & source != "unknown", .(source, dest, monsoon_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="monsoon_particles_rec")[
#    ,-"source"], #remove the source column after casting
#      dcast(MonsoonNormConn[monsoon=="NEM"][dest %in% AddDestGenAnnual[year==2014, offs_site]][source == "unknown", .(source, dest, monsoon_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#          order(source, dest)], source ~ dest, value.var="monsoon_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#MonsoonBiophysMatNEM[is.na(MonsoonBiophysMatNEM)] <- 0 #change NAs to zeros
#dim(MonsoonBiophysMatNEM)
#
#MonsoonBiophysMatSWM <- as.matrix(rbind(dcast(MonsoonNormConn[source %in% AddDestGenAnnual[year==2014, par_site] & dest %in% AddDestGenAnnual[year==2014, offs_site]][monsoon=="SWM" & source != "unknown", .(source, dest, monsoon_particles_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
#    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
#        , source ~ dest, value.var="monsoon_particles_rec")[
#    ,-"source"], #remove the source column after casting
#      dcast(MonsoonNormConn[monsoon=="SWM"][dest %in% AddDestGenAnnual[year==2014, offs_site]][source == "unknown", .(source, dest, monsoon_particles_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
#          order(source, dest)], source ~ dest, value.var="monsoon_particles_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#MonsoonBiophysMatSWM[is.na(MonsoonBiophysMatSWM)] <- 0 #change NAs to zeros
#dim(MonsoonBiophysMatSWM)

matrices <- list(AnnualBiophysMat2012=AnnualBiophysMat2012, AnnualBiophysMat2013=AnnualBiophysMat2013, AnnualBiophysMat2014=AnnualBiophysMat2014,FullBiophysMat=FullBiophysMat, MonsoonBiophysMatNEM=MonsoonBiophysMatNEM, MonsoonBiophysMatSWM=MonsoonBiophysMatSWM)

return(matrices)
}
#fwrite(FullBiophysMat, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2012-14ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2012, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2012ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2013, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2013ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2014, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2014ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNEM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrixNEM2012-14ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatSWM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrixSWM2012-14ForROMSComp15DayPLD.csv")