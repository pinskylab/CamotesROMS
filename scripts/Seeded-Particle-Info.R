Packages <- c("dplyr",  "nleqslv", "broom","cubature", "geosphere", "data.table",  "ggplot2", "bbmle", "stringr",  "lubridate", "RColorBrewer", "viridis")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#make a table of the dates of release for simulations, for calculating the number of particles released in each time frame
season1 <- data.table(date=seq(as.Date("2010-10-01"), as.Date("2011-05-31"), by="days"))

season2 <- data.table(date=seq(as.Date("2011-10-01"), as.Date("2012-05-31"), by="days"))

season3 <- data.table(date=seq(as.Date("2012-10-01"), as.Date("2013-05-31"), by="days"))

season4 <- data.table(date=seq(as.Date("2013-10-01"), as.Date("2014-04-18"), by="days"))

AllDates <- rbind(season1, season2, season3, season4)

#mark the monsoon seasons, based on the same criteria I used for the parentage indirectly through the growth estimates
NEM <- c(11, 12, 1, 2, 3, 4, 5, 6)
SWM <- c(7, 8, 9, 10)

AllDates[,date := ymd(date)][, #format as ymd
             sim_monsoon := ifelse(month(date) %in% NEM, "NEM", "SWM")][,#mark monsoon season based on month
             sim_year:=year(date)][,#add year column
            year_sampled:= ifelse(date %in% season1$date, 2011, ifelse(date %in% season2$date, 2012, ifelse(date %in% season3$date, 2013, 2014)))]#and then add a year_sampled for the empircal sampling season of that particle

ReleaseDays <- AllDates[, .(num_release_days_seasonal=.N), by=c("year_sampled", "sim_monsoon")][, num_release_days_annual:= sum(num_release_days_seasonal), by=year_sampled]

total_release_days <- AllDates[year_sampled %in% c(2012, 2013, 2014), .N]#for the all year kernel- how many days of the simulation conincide with our particle sampling?
total_release_days #should be 687

save(AllDates, ReleaseDays, total_release_days, file= "script_output/ROMSDataTables/2021-11-04_SeededParticleInfo.Rdata")