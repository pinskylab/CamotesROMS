#https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Script_outputs/cumulative_prop_hab_sampled_by_site.RData)
load("~/parentage/r_data/cumulative_prop_hab_sampled_by_site.RData")

#NB I checked that this data.table version matched what I made with the data.frame in all previous code 4/23/2021
PropSamp <- as.data.table(cumulative_prop_hab_sampled_by_site)[
        site=="Caridad Proper" | site=="Sitio Lonas" | site=="Sitio Tugas", total_possible_sample_anems:=max(total_anems_sampled), by=site][#total possible anems comes from the metal tag count, so when we didn't use metal tags (2012), then use the max total sampled (our best guess for the number of anems at the site) as the max possible  
        ,PropAnemSamp := total_prop_hab_sampled_anems_tidied][#recalculate the prop anems sampled, and rename
        is.nan(PropAnemSamp), PropAnemSamp := total_anems_sampled/total_possible_sample_anems][
        is.nan(PropAnemSamp), PropAnemSamp :=0][#now the Nans are true 0s, meaning we didn't sample there, rather than just we sampled 4 anems but have an artifical 0 total possible sampled anems because there were no metal tags
        ,.(site, time_frame, end_year, PropAnemSamp)][#drop extranneous columns
        ,site := gsub(". ", ".", site, fixed=T), by=site]#fix space in Magbangon names, fixed=TRUE for non-regex

#now sum together the Magbangon values- they are one site in the ROMS model but two in the genetics
PropSamp <- unique(PropSamp[site %in% c("S.Magbangon", "N.Magbangon"), Combine :="yes"][#mark the rows with a Magbangon site
    Combine=="yes", PropAnemSamp:=sum(PropAnemSamp), by=end_year][#sort them by their year, and then add the values for PropAnemSamp
    site %in% c("S.Magbangon", "N.Magbangon"), site := "Magbangon"][#rename the sites to be just "Magbangon"
    ,-"Combine"])#filter out the repeat Magbangon values (they are now duplicates of one another) and drop the combine column because we don't need it anymore
#check to make sure that this is correct/consistent with unaltered data
sum(cumulative_prop_hab_sampled_by_site[site  %in% c("S. Magbangon", "N. Magbangon") & time_frame==2012, total_prop_hab_sampled_anems_tidied])==PropSamp[site=="Magbangon" & time_frame==2012, PropAnemSamp] #should be TRUE

summary(PropSamp[,PropAnemSamp]) #nothing should be more than 1

#fwrite(PropSamp, file="~/oceanography/script_output/SurveyData/ProportionHabitatSampled.csv")
