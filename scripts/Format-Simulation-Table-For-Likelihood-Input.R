Packages <- c("dplyr",  "nleqslv", "broom","cubature", "geosphere", "data.table",  "ggplot2", "bbmle", "stringr",  "lubridate", "RColorBrewer", "viridis")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#read in the necessary input data for likelihood functions
load(file= "script_output/SurveyData/2021-11-04_InputSurveyData.Rdata")
load(file= "script_output/ROMSDataTables/2021-11-04_SeededParticleInfo.Rdata")

##prep biophysical connectivity matrix
SimConn <- fread(file="script_output/ROMSDataTables/SimConnectivityTableWithMetaLongForm08DayPLD.csv")[dest != "CAI"] #filter out CAI as a destination for now, not very well spatially defined
#outside of the loop, trim this to only be the destinations we sampled
SourceJoin <- SurveyData[SimConn, on = .(site = source, year=year_sampled)]
setnames(SourceJoin, skip_absent=TRUE, c("site", "n_offs_gen", "prop_anem_samp", "total_anems", "num_females"), c("source", "source_num_rec_sampled_annual",  "source_prop_anem_samp", "source_total_anems", "source_num_females"))
DestJoin <- SurveyData[SourceJoin, on = .(site = dest, year)]
setnames(DestJoin, skip_absent=TRUE, c("site", "n_offs_gen", "prop_anem_samp", "total_anems", "num_females"), c("dest", "dest_num_rec_sampled_annual",  "dest_prop_samp", "dest_total_anems", "dest_num_females"))

SimConn <- DestJoin[source %!in% unrealistic_sources & dest %!in% unrealistic_sources ][#& year %in% c(2012, 2013, 2014)#sand flats and Pangasugan are not realistic source or destination sites because there's almost no habitat. Safe to drop, but keep the rest of the possibilities so we can subsample iteratively all possibilities.
    , daily_particles_released := as.numeric(daily_particles_released)] #change from integer to numeric
SimConn <- ReleaseDays[SimConn, on=.(year_sampled=year, sim_monsoon)]#join in the info for number of release days in the time frame
SimConn <- kernels[Year %in% c("2012", "2013", "2014")][, year:=as.integer(Year)][,c("year", "NumParentageMatches")][SimConn, on=.(year=year_sampled)]#add in a column for the observed number of parentage matches
#rename the monsoon column in the full table for consistency
setnames(SimConn, c("sim_monsoon", "NumParentageMatches"), c("monsoon", "num_route_parentage_matches")) #get rid of upper case and inconsistent naming
setcolorder(SimConn, c("particle_id", "source", "dest", "year", "monsoon", "date"))

#at this point, we can make the raw number assignment matrix, but we want to make a normalized version that is num assigned from a source to a destination/ num released from that source
#fwrite(SimConn, file="script_output/ROMSDataTables/SimConnectivityTableCompleteMetaLongForm08DayPLD.csv")