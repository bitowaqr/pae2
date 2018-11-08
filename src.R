install_n_load <- function(package){
  for(i in 1:length(package)){
    if(eval(parse(text=paste("require(",package[i],")")))==0) {
      install.packages(package)
    }
  }
  return (eval(parse(text=paste("require(",package,")"))))
}
required_packages<-c("ggplot2","cowplot","shiny","leaflet","sp","raster","rgeos","geosphere")
install_n_load(required_packages)

  # load data
    poly_df = raster::shapefile("./input/polygons")
    parkrun_marker = raster::shapefile("./input/marker_england")
    greenspaces = raster::shapefile("./input/greenspaces")

  # select subset of parks
    greenspace_area = area(greenspaces)
    greenspace_area_km2 = greenspace_area/1000^2
    enough_space_to_run = greenspace_area_km2 >= 0.05
    greenspaces = subset(greenspaces,enough_space_to_run)

  # coordinates of objects
    poly.coord = coordinates(poly_df)
    colnames(poly.coord) = c("lon","lat")
    parkrun_marker.coord = coordinates(parkrun_marker)
    colnames(parkrun_marker.coord) = c("lon","lat")
    greenspace.coord = coordinates(greenspaces)
    colnames(greenspace.coord) = c("lon","lat")  


# function to find a new park to minimize the objective function
distanc0r = function(candidates, # new park coorindates matrix
                      radius= 10000,  # exclude distance in meters
                      population = poly_df$pop, # population per lsoa
                      ref_sys = poly.coord, # lsoa coordinates matrix
                      objective = "x", # can be set to x^2 or x < 2
                      events.coord = parkrun_marker.coord #parkrunevent coord matrix
                     ){ 
  # transform distance (not a real objectve function)
  objective_function = function(x) {eval(parse(text=objective))}
  
  # baseline distances and distance sum
  dist_0 = distm(ref_sys,events.coord)
  # dist in meters
  dist_0.min_raw = apply(dist_0,1,function(x) min(x))
  # dist transformed (of applicable)
  dist_0.min = objective_function(dist_0.min_raw)
  dist_0.min.pop = dist_0.min * population
  sum.0 = sum(dist_0.min.pop)
  
  # the algorithm only selects lsoa to reevaluate around the new park
  # however, not to miss some outlier lsoa, lsoa with a distance > radius
  # are also always evaluated
  complement = dist_0.min_raw > radius
  
  # used to define the area from where lsoa are evaluated
  degree = seq(0,360,length.out = 9)
  
  # loop over all candiddate park coordinates:
  len = length(candidates[,1])
  sum.i = rep(-1,times=len)
  timecheck = 0
  
  for(i in 1:len ){
    # just a progress report
    progress = round(i/len,4)*100
    if(abs(progress - timecheck)>1){
      timecheck = progress
      cat("\n",progress,"%")
    }
      
    # define area around new park candidiate
      limits = destPoint(candidates[i,],degree,radius)
      min.lon = min(limits[,"lon"])
      max.lon = max(limits[,"lon"])
      min.lat = min(limits[,"lat"])
      max.lat = max(limits[,"lat"])
      
      # exclude areas outside area
      subset_ref.index =  
        ref_sys[,"lon"] < max.lon &
        ref_sys[,"lon"] > min.lon &
        ref_sys[,"lat"] < max.lat &
        ref_sys[,"lat"] > min.lat
      # but include those that have base distance > radius
      subset_ref.complement = subset_ref.index | complement
      subset_ref.M = ref_sys[subset_ref.complement,]
      population_ref = population[subset_ref.complement]
      
      comparators_dist = dist_0.min_raw[subset_ref.complement]
      comparators_dist_pop = dist_0.min.pop[subset_ref.complement]
      
      dist_i_raw = distm(subset_ref.M,candidates[i,])
      dist_i = objective_function(dist_i_raw)
      # select lsoa with decreased travel distance
      gainers = comparators_dist > dist_i_raw
      # compute transformed dist * pop
      dist_i.gainers.pop = dist_i[gainers] * population_ref[gainers]
      # use baseline for those that did not decrease
      dist_i.nongainer.pop = comparators_dist_pop[!gainers]
      # use baseline for those that were not re-evalauted because not near to park
      dist_0.min.pop_unref = dist_0.min.pop[!subset_ref.complement]
    # new sum of distances * population
      sum.i[i] = sum(dist_i.gainers.pop,dist_i.nongainer.pop,dist_0.min.pop_unref)
  }
  
  winner_index = which(sum.i == min(sum.i))
  
  return(list("sum.0" = sum.0,
              "dist_0.min" = dist_0.min,
              "dist_0.min.pop" = dist_0.min.pop,
              "winner_index" = winner_index,
              "winner" = candidates[winner_index,],
              "sum.i" = sum.i))
         }



# select candidate parks coordinates
candidates = greenspace.coord

# copy parkrun events coordinates to be appended with new events
parkrun_marker.coord_append = parkrun_marker.coord

# give rownames to later identify which parks were used
full_len = dim(candidates)[1]
full_names = 1:full_len
rownames(candidates) = full_names

# store results
winner_i = data.frame(index = NA, sum = NA,lon = NA, lat =NA)
ingame = rep(T,times=full_len) 

for(k in 11:25){
  # find best new park
  test = distanc0r(candidates = candidates,
                   radius= 10000, 
                   population = poly_df$pop,
                   ref_sys = poly.coord,
                   events.coord = parkrun_marker.coord_append)
  
  # store results
  winner_i[k,] = c(rownames(candidates)[test$winner_index],
                   min(test$sum.i),
                   test$winner[1],
                   test$winner[2])
  
  # implement new park in the list of parkrun events
  parkrun_marker.coord_append = rbind(test$winner,parkrun_marker.coord_append)
  
  # remove the park from candidate parks 
  candidates = candidates[-test$winner_index,]
  
  # backup results
  save(winner_i,file=paste("savegame",Sys.time(),".rdata",sep=""))
}


# select the park polygons
selected = full_names %in% winner_i$index
selected_parks = subset(greenspaces,
                        selected)

# create SpatialPoints for top new park coordinates
winner_i$lon = as.numeric(winner_i$lon)
winner_i$lat= as.numeric(winner_i$lat)
g2 = data.frame(winner_i[,c("lon","lat")])
coordinates(g2 ) = ~lon + lat
proj4string(g2)=proj4string(poly_df)

# compute new distances from lsoa to new + old parkrun events
dist_new = distm(poly.coord,rbind(parkrun_marker.coord,coordinates(g2)))
dist_new.min = apply(dist_new,1,function(x) min(x))
# in km
poly_df$new_mn_dstn = round(dist_new.min/1000,1)


# select lsoa fill colors for old and new distances 
  q_dists = as.numeric(quantile(poly_df$mn_dstn,probs = c(seq(0,1,by=0.2))))
  pal_dist <- colorBin("RdYlGn", domain = poly_df$mn_dstn, bins = q_dists,reverse =T)
  pal_dist_new <- colorBin("RdYlGn", domain = poly_df$new_mn_dstn, bins = q_dists,reverse =T)
  
  # LEAFLET MAP
  map = leaflet() %>%
    # provider base tiles
    addProviderTiles(providers$Stamen.Toner, group = "Toner Map"
                     ) %>%
    addTiles(group = "OSM Map"
             ) %>%
    addProviderTiles("CartoDB.Positron",group= "Carto Map", options= providerTileOptions(opacity = 0.99)
                     ) %>%
    # set view to sheffield
    setView(lng = -1.43, lat = 53.36, zoom = 7
    ) %>%
    # layer control
    addLayersControl(
      baseGroups = c("Carto Map","Toner Map","OSM Map"),
      overlayGroups = c("parkrun Events","Distance","new Events","new Parks","new Distance"),
      options = layersControlOptions(collapsed = T,autoZIndex=T)
    ) %>%
    
    # add established parkrun events
    addCircleMarkers(
      group = "parkrun Events",
      data = parkrun_marker,
      radius = 5,
      fillColor = "blue",
      stroke = FALSE, fillOpacity = 0.9,
      popup = paste("Course:",parkrun_marker$Club,"<br>",
                    "Established:", parkrun_marker$Estblsh,"<br>",
                    "Age in years:", parkrun_marker$Age_yrs,"<br>",
                    "Mean participants:", round(parkrun_marker$Mn_prtc),"<br>",
                    "Mean volunteers:", round(parkrun_marker$Mn_vlnt),"<br>")
      ) %>% 
    # add top candidates
    addCircleMarkers(
      group = "new Events",data = g2,
      radius = 6,fillColor = "cyan",stroke = FALSE, fillOpacity = 0.9
      ) %>% 
    # add distances from lsoa to old events
    addPolygons(data = poly_df, group = "Distance",
                color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                weight = 0.1, fillOpacity = 0.5,fillColor = ~pal_dist(mn_dstn),
                highlight = highlightOptions(
                  weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE),
                popup = paste(
                  poly_df$name,"<br>","Nearest event:", poly_df$nrst_vn,"<br>",
                  "Distance: ",poly_df$mn_dstn," km <br>","new Distance: ",poly_df$new_mn_dstn," km <br>",
                              "SIMD score:", poly_df$a,"<br>","Pop density:", poly_df$pp_dnst,"<br>","Population:", poly_df$pop)
                ) %>%
    # add distances from lsoa to new and old events
    addPolygons(data = poly_df, group = "new Distance",
                color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                weight = 0.1, fillOpacity = 0.5,fillColor = ~pal_dist_new(new_mn_dstn),
                highlight = highlightOptions(
                  weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE),
                popup = paste(
                  poly_df$name,"<br>","Nearest event:", poly_df$nrst_vn,"<br>",
                  "Distance: ",poly_df$mn_dstn," km <br>","new Distance: ",poly_df$new_mn_dstn," km <br>",
                  "SIMD score:", poly_df$a,"<br>","Pop density:", poly_df$pp_dnst,"<br>","Population:", poly_df$pop)
                ) %>%
    # add candidate parks as parks
    addPolygons(data = selected_parks, group = "new Parks",
                color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                weight = 0.1, fillOpacity = 0.5,fillColor = "cyan",
                highlight = highlightOptions(
                  weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE)
                ) %>%
    # hide all layers except old parkrun events
    hideGroup("Distance") %>%
    hideGroup("new Distance") %>%
    hideGroup("new Parks") %>%
    hideGroup("new Events")


  # save the map as a stand-aline html file
  htmlwidgets::saveWidget(map,file="index.html",title = "parkrun, access & equity",selfcontained = F)
  
# map

  
  ### ANALYSIS