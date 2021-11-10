#make a normalized version that is num assigned from a source to a destination/ num released from that source
format_biophys_normalized_matrix <- function(SimConn, total_release_days,AddDestAllYearsSim,AddDestAnnualSim, AddDestMonsoonSim){

AllYearsRec <- SimConn[ , .(total_particles_rec = .N), by= c("source","dest")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released)*as.numeric(total_release_days)), by= c("source")], by="source") #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

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
#dim(FullBiophysMatNorm )
FullBiophysMatNorm[is.na(FullBiophysMatNorm)] <- 0 #change NAs to zeros

#make annual matrices with all of the particle data
AnnualRec <- SimConn[ , .(annual_particles_rec = .N), by= c("source","dest", "year")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
AnnualRelease <- unique(SimConn[, .(annual_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_annual)), by= c("source", "year")], by= c("source", "year")) #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

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
summary(AnnualNormConn[, .(sum=sum(source_norm_rec)), by=c("year", "source")][, sum])#nothing should exceed 1, it doesn't- great

#cast into wide format for each year

AnnualBiophysMatNorm2012 <- as.matrix(rbind(dcast(AnnualNormConn[year==2012 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2012][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#dim(AnnualBiophysMatNorm2012)
AnnualBiophysMatNorm2012[is.na(AnnualBiophysMatNorm2012)] <- 0 #change NAs to zeros

AnnualBiophysMatNorm2013 <- as.matrix(rbind(dcast(AnnualNormConn[year==2013 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2013][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#dim(AnnualBiophysMatNorm2013)
AnnualBiophysMatNorm2013[is.na(AnnualBiophysMatNorm2013)] <- 0 #change NAs to zeros

AnnualBiophysMatNorm2014 <- as.matrix(rbind(dcast(AnnualNormConn[year==2014 & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(AnnualNormConn[year==2014][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#dim(AnnualBiophysMatNorm2014)
AnnualBiophysMatNorm2014[is.na(AnnualBiophysMatNorm2014)] <- 0 #change NAs to zeros

#make monsoon matrices with all of the particle data
MonsoonRec <- SimConn[ , .(monsoon_particles_rec = .N), by= c("source","dest", "monsoon")] #all particles recruiting along each route FILTER HERE FOR TIME PERIOD***
MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_seasonal)), by= c("source", "monsoon")], by= c("source", "monsoon")) #calculate the number of particles released over the time frame by multiplyig the release days by the number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

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
#dim(MonsoonBiophysMatNormNEM)
MonsoonBiophysMatNormNEM[is.na(MonsoonBiophysMatNormNEM)] <- 0 #change NAs to zeros

MonsoonBiophysMatNormSWM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon=="SWM" & source != "unknown", .(source, dest, source_norm_rec)][ #for assigned particles (not from "Other") keep the source/dest columns that will be expanded into wide form to become the connectivity matrix. Filtering for time period etc can be done in i here.
    order(source, dest)] #keep sites in alphabetical order so the matrix is correctly formatted!
        , source ~ dest, value.var="source_norm_rec")[
    ,-"source"], #remove the source column after casting
      dcast(MonsoonNormConn[monsoon=="SWM"][source == "unknown", .(source, dest, source_norm_rec)][ #this is to cast the "unassigned row for the model parentage, which is anyting from "Other"
          order(source, dest)], source ~ dest, value.var="source_norm_rec")[,-"source"]))#bind these two cast wide form data tables (assigned and unassigned particles) and then turn into a matrix to be used in the likelihood functions
#dim(MonsoonBiophysMatNormSWM)
MonsoonBiophysMatNormSWM[is.na(MonsoonBiophysMatNormSWM)] <- 0 #change NAs to zeros

matrices <- list(AnnualBiophysMatNorm2012=AnnualBiophysMatNorm2012, AnnualBiophysMatNorm2013=AnnualBiophysMatNorm2013, AnnualBiophysMatNorm2014=AnnualBiophysMatNorm2014,FullBiophysMatNorm=FullBiophysMatNorm, MonsoonBiophysMatNormNEM=MonsoonBiophysMatNormNEM, MonsoonBiophysMatNormSWM=MonsoonBiophysMatNormSWM)

return(matrices)
}
#fwrite(FullBiophysMatNorm, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2012-14ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2012, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2012ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2013, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2013ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2014, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2014ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNormNEM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrixNEM2012-14ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNormSWM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysNormConnMatrixSWM2012-14ForROMSComp15DayPLD.csv")