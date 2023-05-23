#clownfish metadata
load("empirical_data/genetics/total_sampling_across_years.RData")
load("empirical_data/genetics/sampled_area_each_year.RData")

#download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "empirical_data/genetics/fish-obs.RData")
fish_obs <- readRDS("empirical_data/genetics/fish-obs.RData") 
load("empirical_data/genetics/site_dist_info.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/anem_db.RData?raw=true", destfile = "empirical_data/genetics/anem_db.RData")
load("empirical_data/genetics/anem_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/dives_db.RData?raw=true", destfile = "empirical_data/genetics/dives_db.RData")
load("empirical_data/genetics/dives_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/fish_db.RData?raw=true", destfile = "empirical_data/genetics/dives_db.RData")
load("empirical_data/genetics/fish_db.RData")
load("empirical_data/genetics/gps_db.RData")

##read in the sites that we sampled each year
N_gen_par <- read.table(file="empirical_data/genetics/all_parents_corrected.txt", header = TRUE, stringsAsFactors = F) %>%#not sure that I need the parents here
    mutate(fish_indiv=as.character(fish_indiv))
N_gen_offs <- read.table(file="empirical_data/genetics/all_offspring_corrected.txt", header=T, stringsAsFactors = F) %>%
    mutate(fish_indiv=as.character(fish_indiv))



#gather the summary of total offspring sampled
#from Allison, putting all the meta data together (Constants_database_common_functions.R)
##### Match up other relevant info (site, date, fish_indiv, etc.) to fish in the clownfish table
# Pull out year and month into a separate column in dives_db
dives_db_processed <- dives_db %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(month = as.integer(substring(date,6,7))) %>%
  mutate(dive_date = date(date))

# Pull all APCL caught or otherwise in the clownfish table
allfish_fish <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, anem_table_id, recap, tag_id, color, sex, size, fish_obs_time, fish_notes) %>%
  filter(fish_spp == 'APCL') %>%
  mutate(size = as.numeric(size))  # make the size numeric (rather than chr) so can do means and such

# and their corresponding anemones
allfish_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_notes) %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, year, month, site, gps, dive_notes) %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# join together
allfish_caught <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish_caught <- left_join(allfish_caught, allfish_dives, by="dive_table_id")

# add in the gen_ids and fish_indiv (now in a separate gen_id table) - gen_id only comes in the time the fish was sequenced, not at all captures
allfish_caught <- left_join(allfish_caught, (fish_obs %>% select(fish_table_id, gen_id, fish_indiv)), by = "fish_table_id") %>%
    select(fish_indiv, sample_id, site) %>%
    mutate()

N_gen_offs_annual  <- left_join(N_gen_offs, allfish_caught, by=c("fish_indiv", "sample_id")) %>% 
    group_by(year, site) %>%
    summarise(n_offs_gen=n()) %>%
    ungroup()

N_gen_offs_annual$site <- gsub(". ", ".", N_gen_offs_annual$site, fixed=TRUE)


#combine N/S Magbangon in the genetic sampling data

AnnualRecsSamp <- bind_rows(N_gen_offs_annual %>%
                        filter(site !="N.Magbangon" & site!="S.Magbangon"),
                    N_gen_offs_annual %>%
                        mutate(Magbangon=ifelse(site=="N.Magbangon" | site=="S.Magbangon", "Magbangon", "no")) %>%
                        filter(Magbangon=="Magbangon") %>%
                        group_by(year, Magbangon) %>%
                        summarise(sum_offs=sum(n_offs_gen)) %>%
                        rename(site="Magbangon", n_offs_gen="sum_offs")) %>%
                        mutate(year=as.character(year)) %>%
                        ungroup()
sum(AnnualRecsSamp$n_offs_gen) #should be 791

#fwrite(AnnualRecsSamp, file="script_output/SurveyData/AnnualRecruitsSampled.csv")
