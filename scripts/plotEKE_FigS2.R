# R script to make Fig. S2 (plot of EKE by site)

# load libraries and read in EKE calculations
require(data.table)
require(ggplot2)

VelFields <- fread("ROMS/data/EKE.csv")

# averages across sites (for Table S3)
means <- colMeans(VelFields[,2:ncol(VelFields)])
means # averages across sites
apply(VelFields[,2:ncol(VelFields)], MARGIN=2, FUN=sd) # standard deviations across sites


mean(means[c('Jun2011-May2012', 'Jun2012-May2013', 'Jun2013-May2014')]) # average across years
sd(VelFields[,c(`Jun2011-May2012`, `Jun2012-May2013`, `Jun2013-May2014`)]) # std dev across years

#change to long format for plotting. rename some varibales.
VelFieldsLong <- melt(VelFields, id.vars = "Location")[ 
	, Location := gsub("_", " ", Location)][
	variable == "NEM_Amihan", variable := 'Amihan'][
	variable == "SWM_Habagat", variable := "Habagat"][
    variable == "Jun2011-May2012", variable := "2012"][
    variable == "Jun2012-May2013", variable := "2013"][
    variable == "Jun2013-May2014", variable := "2014"][
    , time_scale := ifelse(variable %in% c("2012", "2013", "2014"), "annual", "monsoonal")][
    Location == "Tamakin Dacot", Location := "Tomakin Dako"]

#order sites N-S
VelFieldsLong$Location <- factor(VelFieldsLong$Location, levels = rev(c("Palanas", "Wangag", "Magbangon", "Cabatoan", "Caridad Cemetery", "Caridad Proper", "Hicgop", "Hicgop South", "Sitio Tugas", "Elementary School", "Sitio Lonas", "San Agustin", "Poroc San Flower", "Poroc Rose",  "Pangasugan", "Visca", "Gabas", "Tomakin Dako", "Haina", "Sitio Baybayon")))
VelFieldsLong$variable <- factor(VelFieldsLong$variable, levels = c('2012', '2013', '2014', 'Amihan', 'Habagat'))

EKE_avg <- ggplot()+
    geom_point(data=VelFieldsLong, aes(x=value, y=Location, color=variable, shape=variable), size=2)+
    labs(y="Site", x = expression ("Average eddy kinetic energy"~cm^2/s^2))+
    theme_bw()+
	facet_wrap(~ time_scale, ncol=2)
EKE_avg

ggsave(filename="Camotes_EKE_Avg.pdf",  plot=EKE_avg, path="script_output/plots/", units="cm", width=15, height=15)