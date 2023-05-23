#make a normalized version that is num assigned from a source to a destination/ num released from that source
format_biophys_normalized_matrix <- function(SimConn, total_release_days,AddDestSim){

	AllYearsRec <- SimConn[, .(total_particles_rec = .N), by = c("source", 
	        "dest")]
	AllYearsRec <- rbind(unique(AddDestSim, by=c("source", "dest"))[, -c("year", "monsoon")][!AllYearsRec, on = .(source, 
		        dest)][, total_particles_rec := 0], AllYearsRec)
	AllYearsRelease <- unique(SimConn[, .(total_particles_released = as.numeric(daily_particles_released) * 
	        as.numeric(total_release_days)), by = c("source")], by = "source")
	AllYearsNormConn <- AllYearsRelease[AllYearsRec, on = "source"][ 
	        , source_norm_rec:= total_particles_rec/total_particles_released]
	FullBiophysMatNorm <- as.matrix(rbind(dcast(AllYearsNormConn[source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	        dest)], source ~ dest, value.var = "source_norm_rec")[, 
	        -"source"], dcast(AllYearsNormConn[source == "unknown", .(source, dest, source_norm_rec)][order(source, 
	        dest)], source ~ dest, value.var = "source_norm_rec")[, 
	        -"source"]))
	FullBiophysMatNorm[is.na(FullBiophysMatNorm)] <- 0
	#sum(FullBiophysMatNorm)
	
	#make annual matrices with all of the particle data
	AnnualRec <- SimConn[ , .(annual_particles_rec = .N), by= c("source","dest", "year")] #all particles recruiting along each route FILTER HERE FOR TIME 	PERIOD***
	AnnualRelease <- unique(SimConn[, .(annual_particles_released = as.numeric(daily_particles_released)*as.numeric(num_release_days_annual)), by= 	c("source", "year")], by= c("source", "year")) #calculate the number of particles released over the time frame by multiplyig the release days by the 	number of particles released daily. fread() converts big numbers to integers so specify as numeric to avoid integer overflow NAs

	AnnualRec <- rbind(unique(AddDestSim, by=c("source", "dest", "year"))[, -"monsoon"][!AnnualRec, on = .(source, 
	        dest, year)][, annual_particles_rec := 0], AnnualRec)
	AnnualRelease <- rbind(unique(AddDestSim, by=c("source", "year"))[, -c("dest","monsoon")][!AnnualRelease, 
	        on = .(source, year)][, annual_particles_released :=  0], AnnualRelease)
	AnnualNormConn <- AnnualRelease[AnnualRec, on = c("source", 
	        "year")][, source_norm_rec := annual_particles_rec/annual_particles_released]
	#sum(AnnualNormConn$annual_particles_rec)
	AnnualBiophysMatNorm2012 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2012][source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2012][source == 
	    "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"]))
	AnnualBiophysMatNorm2012[is.na(AnnualBiophysMatNorm2012)] <- 0
	AnnualBiophysMatNorm2013 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2013][source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2013][source == 
	    "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"]))
	AnnualBiophysMatNorm2013[is.na(AnnualBiophysMatNorm2013)] <- 0
	AnnualBiophysMatNorm2014 <- as.matrix(rbind(dcast(AnnualNormConn[year == 
	    2014][source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"], dcast(AnnualNormConn[year == 2014][source == 
	    "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"]))
	AnnualBiophysMatNorm2014[is.na(AnnualBiophysMatNorm2014)] <- 0

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
	MonsoonBiophysMatNormNEM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon == 
	    "NEM"][source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"], dcast(MonsoonNormConn[monsoon == "NEM"][source == 
	    "unknown", .(source, dest,source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"]))
	MonsoonBiophysMatNormNEM[is.na(MonsoonBiophysMatNormNEM)] <- 0
	#dim(MonsoonBiophysMatNormNEM)
	  MonsoonBiophysMatNormSWM <- as.matrix(rbind(dcast(MonsoonNormConn[monsoon == 
	    "SWM"][source != "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"], dcast(MonsoonNormConn[monsoon == "SWM"][source == 
	    "unknown", .(source, dest, source_norm_rec)][order(source, 
	    dest)], source ~ dest, value.var = "source_norm_rec")[, 
	    -"source"]))
	MonsoonBiophysMatNormSWM[is.na(MonsoonBiophysMatNormSWM)] <- 0
	
matrices <- list(AnnualBiophysMatNorm2012=AnnualBiophysMatNorm2012, AnnualBiophysMatNorm2013=AnnualBiophysMatNorm2013, AnnualBiophysMatNorm2014=AnnualBiophysMatNorm2014,FullBiophysMatNorm=FullBiophysMatNorm, MonsoonBiophysMatNormNEM=MonsoonBiophysMatNormNEM, MonsoonBiophysMatNormSWM=MonsoonBiophysMatNormSWM)

return(matrices)
}
#fwrite(FullBiophysMatNorm, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2012-14ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2012, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2012ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2013, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2013ForROMSComp15DayPLD.csv")
#fwrite(AnnualBiophysMatNorm2014, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrix2014ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNormNEM, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrixNEM2012-14ForROMSComp15DayPLD.csv")
#fwrite(MonsoonBiophysMatNormSWM, file="script_output/ROMSDataTables/20210917_BioPhysNormConnMatrixSWM2012-14ForROMSComp15DayPLD.csv")