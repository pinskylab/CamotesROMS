#03/08/2018 code to identify sites in ROMS model vertices list from Liz Drenkard

####SET UP#######
library(dplyr)
library(rgdal)
library(tidyr)

###read txt file with ROMS vertices
vertices <- read.table("ROMS/input/Camotes_Sea_ROMS_Grid_Vertices.txt", header=FALSE, sep=" ")

#label rows and columns so I can find which match site domains.
dim(vertices)
rows <- c(seq(from =1, to=76104, by=1))
cols <- c("one_lon", "two_lon", "three_lon", "four_lon", "one_lat", "two_lat", "three_lat", "four_lat", "id_number")
vertices$id_number <- NA
vertices$id_number <- rows
colnames(vertices) <- cols

#now trim the vertices to get rid of points outside of the sites
test <- vertices %>% filter(one_lon <=124.81000 & one_lon >=124.70000) %>% filter(one_lat <=10.8900 & one_lat >=10.6200)
test2 <- test %>% filter(two_lon <=124.81000 & two_lon >=124.70000) %>% filter(two_lat <=10.8900 & two_lat >=10.6200)
test3 <- test2 %>% filter(three_lon <=124.81000 & three_lon >=124.70000) %>% filter(three_lat <=10.8900 & three_lat >=10.6200)
test4 <- test3 %>% filter(four_lon <=124.81000 & four_lon >=124.70000) %>% filter(four_lat <=10.8900 & four_lat >=10.6200)

tidy_vertices <- gather(test4, "id_lon", "lon", 1:4)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 1:4)

#arrange by id number, ID of grid cell in text file and take the 1st, 6th, 11th, and 16th rows to get the matching coordinates. should be four lat/lon rows (four corners) for each row number
tidy_vertices2$id_number <- as.factor(tidy_vertices2$id_number)
target <- c(1, 6, 11, 16)
no_rep <- tidy_vertices2 %>% group_by(id_number) %>% arrange(id_number) %>% filter(row_number(id_number) %in% target)

write.csv(no_rep, "ROMS/input/camotes_vertices_sites.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)

#read in row id numbers from QGIS visual survey to get vertice points
results <- read.table(file="ROMS/input/camotes_site_vertices.txt", header=TRUE)

site_vertices <- inner_join(results, vertices, by="id_number")
write.table(site_vertices, "ROMS/input/camotes_vertices_sites_results.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)

#make a text file with near north and near south site points 
test <- vertices %>% filter(one_lon <=124.81000 & one_lon >=124.60000) %>% filter(one_lat <=10.6330 | one_lat >=10.8700)
test2 <- test %>% filter(two_lon <=124.81000 & two_lon >=124.60000) %>% filter(two_lat <=10.6330  | two_lat >=10.8700)
test3 <- test2 %>% filter(three_lon <=124.81000 & three_lon >=124.60000) %>% filter(three_lat <=10.6330 | three_lat >=10.8700)
test4 <- test3 %>% filter(four_lon <=124.81000 & four_lon >=124.60000) %>% filter(four_lat <=10.6330 | four_lat >=10.8700)

tidy_vertices <- gather(test4, "id_lon", "lon", 1:4)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 1:4)

#arrange by id number, ID of grid cell in text file and take the 1st, 6th, 11th, and 16th rows to get the matching coordinates. should be four lat/lon rows (four corners) for each row number
tidy_vertices <- gather(test4, "id_lon", "lon", 1:4)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 1:4)

tidy_vertices2$id_number <- as.factor(tidy_vertices2$id_number)
target <- c(1, 6, 11, 16)
no_rep <- tidy_vertices2 %>% group_by(id_number) %>% arrange(id_number) %>% filter(row_number(id_number) %in% target)

write.csv(no_rep, "ROMS/input/camotes_vertices_nearby_sites.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)

#just all vertices
tidy_vertices <- gather(vertices, "id_lon", "lon", 1:4)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 1:4)

tidy_vertices2$id_number <- as.factor(tidy_vertices2$id_number)
target <- c(1, 6, 11, 16)
no_rep <- tidy_vertices2 %>% group_by(id_number) %>% arrange(id_number) %>% filter(row_number(id_number) %in% target)

write.csv(no_rep, "ROMS/input/camotes_vertices_all.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)

####make a csv for camotes or cuatro islas

#read in row id numbers from QGIS visual survey to get vertice points
camotes <- read.table(file="cuatro_islas.txt", header=TRUE)

site_vertices_camotes <- inner_join(camotes, vertices, by="id_number")
site_vertices_camotes<- distinct(site_vertices_camotes)
write.table(site_vertices_camotes, "ROMS/input/cuatro_islas_results.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


#make csv to read in and check the vertices

tidy_vertices <- gather(site_vertices_camotes, "id_lon", "lon", 2:5)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 2:5)

tidy_vertices2$id_number <- as.factor(tidy_vertices2$id_number)
target <- c(1, 6, 11, 16)
no_rep_islas <- tidy_vertices2 %>% group_by(id_number) %>% arrange(id_number) %>% filter(row_number(id_number) %in% target)

write.csv(no_rep_islas, "ROMS/input/cuatro_islas_results.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)

####make a csv for all 

#read in row id numbers from QGIS visual survey to get vertice points
results2 <- read.table(file="ROMS/input/camotes_nearby_sites.txt", header=TRUE)

site_vertices_all <- inner_join(results2, vertices, by="id_number")
write.table(site_vertices_all, "ROMS/input/camotes_nearby_sites_results.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


#make csv to read in and check the vertices
all_results <- bind_rows(site_vertices, site_vertices_all)

tidy_vertices <- gather(all_results, "id_lon", "lon", 3:6)
tidy_vertices2 <- gather(tidy_vertices, "id_lat", "lat", 3:6)

tidy_vertices2$id_number <- as.factor(tidy_vertices2$id_number)
target <- c(1, 6, 11, 16)
no_rep_all <- tidy_vertices2 %>% group_by(id_number) %>% arrange(id_number) %>% filter(row_number(id_number) %in% target)

write.csv(no_rep_all, "ROMS/input/camotes_vertices_all.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)


####UNNECESSARY######
#now load in site maps
map_sites <- readOGR("site_hulls", "site_hulls")
map_names <- readOGR("site_names", "Sitenames")

par(new=TRUE)
plot(map_sites)
par(new=TRUE)
plot(map_names)
points(vertices)
text(edits$lon, edits$lat, labels=edits$site, col="black", cex=.7, offset=1.8, pos=4)


