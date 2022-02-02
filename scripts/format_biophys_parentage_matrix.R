#make a parentage matrix for the whole biophysical results- CAI included with Other as "unknown", use this to fit a kernel
format_biophys_parentage_matrix <- function(SimConn, total_release_days, AddDestSim){
    
	#double checked, this is correct***
	AllYearsRec <- SimConn[, .(total_particles_rec = .N), by = c("source", 
	        "dest")]
	AllYearsRec <- rbind(unique(AddDestSim, by=c("source", "dest"))[, -c("year", "monsoon")][!AllYearsRec, on = .(source, 
	        dest)][, total_particles_rec := 0], AllYearsRec)

	AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released) * 
	        as.numeric(total_release_days)), by = c("source")], by = "source")
	AllYearsNormConn <- AllYearsRelease[AllYearsRec, on = "source"][ 
	        , source_norm_rec:= total_particles_rec/total_particles_released]
	FullBiophysMat <- as.matrix(rbind(dcast(AllYearsNormConn[source != "unknown", .(source, dest, total_particles_rec)][order(source, 
	        dest)], source ~ dest, value.var = "total_particles_rec")[, 
	        -"source"], dcast(AllYearsNormConn[source == "unknown", .(source, dest, total_particles_rec)][order(source, 
	        dest)], source ~ dest, value.var = "total_particles_rec")[, 
	        -"source"]))
	FullBiophysMat[is.na(FullBiophysMat)] <- 0
	#sum(FullBiophysMat)
	
	#make annual matrices with all of the particle data
	AnnualRec <- SimConn[ , .(annual_particles_rec = .N), by= c("source","dest", "year")] #all particles recruiting 	along each route FILTER HERE FOR TIME 	PERIOD***
	AnnualRelease <- unique(SimConn[, .(annual_particles_released = 	as.numeric(daily_particles_released)*as.numeric(num_release_days_annual)), by= 	c("source", "year")], by= 	c("source", "year")) #calculate the number of particles released over the time frame by multiplyig the release days 	by the 	number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid 	integer overflow NAs
	
	AnnualRec <- rbind(unique(AddDestSim, by=c("source", "dest", "year"))[, -"monsoon"][!AnnualRec, on = .(source, 
	        dest, year)][, annual_particles_rec := 0], AnnualRec)
	AnnualRelease <- rbind(unique(AddDestSim, by=c("source", "year"))[, -c("dest","monsoon")][!AnnualRelease, 
	        on = .(source, year)][, annual_particles_released :=  0], AnnualRelease)
	AnnualNormConn <- AnnualRelease[AnnualRec, on = c("source", 
	        "year")][, source_norm_rec := annual_particles_rec/annual_particles_released]
	#sum(AnnualNormConn$annual_particles_rec)
	AnnualBiophysMat2012 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2012][source != "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2012][source == 
	    "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"]))
	AnnualBiophysMat2012[is.na(AnnualBiophysMat2012)] <- 0
	AnnualBiophysMat2013 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2013][source != "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2013][source == 
	    "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"]))
	AnnualBiophysMat2013[is.na(AnnualBiophysMat2013)] <- 0
	AnnualBiophysMat2014 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2014][source != "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2014][source == 
	    "unknown", .(source, dest, annual_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "annual_particles_rec")[, 
	    -"source"]))
	AnnualBiophysMat2014[is.na(AnnualBiophysMat2014)] <- 0

	MonsoonRec <- SimConn[, .(monsoon_particles_rec = .N), by = c("source", 
	    "dest", "monsoon")]
	MonsoonRelease <- unique(SimConn[, .(monsoon_particles_released = as.numeric(daily_particles_released) * 
	    as.numeric(num_release_days_seasonal)), by = c("source", 
	    "monsoon")], by = c("source", "monsoon"))
	MonsoonRec <- rbind(unique(AddDestSim, by=c("source", "dest", "monsoon"))[, -"year"][!MonsoonRec, on = .(source, 
	    dest, monsoon)][, monsoon_particles_rec :=0], MonsoonRec)
	MonsoonRelease <- rbind(unique(AddDestSim, by=c("source", "monsoon"))[, -c("dest", "year")][!MonsoonRelease, 
	    on = .(source, monsoon)][, monsoon_particles_released :=0], MonsoonRelease)
	MonsoonNormConn <- MonsoonRelease[MonsoonRec, on = c("source", 
	    "monsoon")][, source_norm_rec := monsoon_particles_rec/monsoon_particles_released]
	MonsoonBiophysMatNEM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon == 
	    "NEM"][source != "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
	    -"source"], dcast(MonsoonNormConn[monsoon == "NEM"][source == 
	    "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
	    -"source"]))
	MonsoonBiophysMatNEM[is.na(MonsoonBiophysMatNEM)] <- 0
	#dim(MonsoonBiophysMatNEM)
	  MonsoonBiophysMatSWM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon == 
	    "SWM"][source != "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
	    -"source"], dcast(MonsoonNormConn[monsoon == "SWM"][source == 
	    "unknown", .(source, dest, monsoon_particles_rec)][order(source, 
	    dest)], source ~ dest, value.var = "monsoon_particles_rec")[, 
	    -"source"]))
	MonsoonBiophysMatSWM[is.na(MonsoonBiophysMatSWM)] <- 0

matrices <- list(AnnualBiophysMat2012=AnnualBiophysMat2012, AnnualBiophysMat2013=AnnualBiophysMat2013, AnnualBiophysMat2014=AnnualBiophysMat2014,FullBiophysMat=FullBiophysMat, MonsoonBiophysMatNEM=MonsoonBiophysMatNEM, MonsoonBiophysMatSWM=MonsoonBiophysMatSWM)

return(matrices)
}
#fwrite(FullBiophysMat, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2012-14ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2012, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2012ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2013, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2013ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMat2014, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrix2014ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNEM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrixNEM2012-14ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatSWM, file="~/oceanography/script_output/ROMSDataTables/20210917_BioPhysParentageMatrixSWM2012-14ForROMSComp15DayPLD.csv")