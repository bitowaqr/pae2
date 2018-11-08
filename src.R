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





### find best parks for new events
  # find a list of single best parks
  first.picks = distanc0r(candidates = greenspace.coord,
                          radius= 10000, 
                          population = poly_df$pop,
                          ref_sys = poly.coord,
                          events.coord = parkrun_marker.coord)
  
  first_picks = order(first.picks$sum.i)
  top_s = 100
  first_picks_selection = 1:length(first_picks) %in% first_picks[1:top_s] 
  first_selection_parks = subset(greenspaces,
                                 first_picks_selection)
  reduction_perc = round((first.picks$sum.i - first.picks$sum.0)/first.picks$sum.0,5)*100
  first_selection_parks$dist_reduction = reduction_perc[order(first.picks$sum.i,decreasing = F)][1:top_s] 
  first_selection_parks.coord = coordinates(first_selection_parks)

  ### find the consecutive top 25 
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
    
    for(k in 1:25){
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


    
    ## assess results
    # select the park polygons
    selected = full_names %in% winner_i$index
    selected_parks = subset(greenspaces,
                            selected)
    selected_parks.coord = coordinates(selected_parks)
    
    # compute new distances from lsoa to new + old parkrun events
    dist_new = distm(poly.coord,rbind(parkrun_marker.coord,selected_parks.coord))
    dist_new.min = apply(dist_new,1,function(x) min(x))
    poly_df$new_mn_dstn = round(dist_new.min/1000,1) # in km
    
    # assess resiults per parkrun event (new or old)
    m.temp = matrix(nrow = dim(dist_new)[1],
                    ncol= dim(dist_new)[2],
                    data =F)
    for(row in 1:dim(dist_new)[1]){
      m.temp[row,]   = ifelse(dist_new[row,] == min(dist_new[row,]), poly_df$pop[row]*min(dist_new[row,]),0)
    }
    event_pop = colSums(m.temp)
    event_lsoa_n = apply(m.temp,2,function(x) sum(x != 0))
    
    parkrun_marker$lsoa_served = event_lsoa_n[1:dim(parkrun_marker.coord)[1]]
    parkrun_marker$pop_served = event_pop[1:dim(parkrun_marker.coord)[1]]
    
    selected_parks$lsoa_served = event_lsoa_n[(length(event_lsoa_n)-dim(winner_i)[1]+1):length(event_lsoa_n)]
    selected_parks$pop_served = event_pop[(length(event_lsoa_n)-dim(winner_i)[1]+1):length(event_lsoa_n)]
    
    selected_parks$dist_reduction = round((as.numeric(winner_i$sum) - first.picks$sum.0)/first.picks$sum.0,5) *100

    for(i in length(selected_parks$dist_reduction):2){
      selected_parks$dist_reduction[i] = selected_parks$dist_reduction[i] - selected_parks$dist_reduction[i-1]
    }
    
    # parks first or consecutive
    selected_parks_any = subset(greenspaces,
                            selected | first_picks_selection)
    
    
### MAP RESULTS
      # select lsoa fill colors for old and new distances 
        q_dists = as.numeric(quantile(poly_df$mn_dstn,probs = c(seq(0,1,by=0.2))))
        pal_dist <- colorBin("RdYlGn", domain = poly_df$mn_dstn, bins = q_dists,reverse =T)
        pal_dist_new <- colorBin("RdYlGn", domain = poly_df$new_mn_dstn, bins = q_dists,reverse =T)
        
        # tag
        rr <- tags$div(
          'Schneider et al.',tags$em('parkrun, access, equity.'), '2018. Data and code available at:',tags$a(href="https://github.com/bitowaqr/pae2", "https://github.com/bitowaqr/pae2",target="_blank")  
        )  
        
        

        
        
        
        
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
          # add tag
          addControl(rr, position = "bottomright") %>%
        
          # layer control
          addLayersControl(
            baseGroups = c("Carto Map","Toner Map","OSM Map"),
            overlayGroups = c("parkrun Events","Distance","new Events (first)","new Events (consecutive)","candidate parks","new Distance"),
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
                          "Mean volunteers:", round(parkrun_marker$Mn_vlnt),"<br>",
                          "Served population (new):", parkrun_marker$pop_served,"<br>",
                          "Served LSOA (new):", parkrun_marker$lsoa_served,"<br>")
            ) %>% 
          # add top candidates
          addCircleMarkers(
            group = "new Events (consecutive)", 
            lng = selected_parks.coord[,1],
            lat = selected_parks.coord[,2],
            radius = 6,fillColor = "purple",stroke = FALSE, fillOpacity = 0.9,
            popup = paste("Park:", selected_parks$distName1,"<br>",
                          "Pick no:", 1:length(selected_parks$id),"<br>",
                          "Reduced pop*dist:",selected_parks$dist_reduction,"% <br>",
                          "Served population:", selected_parks$pop_served,"<br>",
                          "Served LSOA:", selected_parks$lsoa_served,"<br>")
            ) %>% 
          addCircleMarkers(
            group = "new Events (first)",
            lng = first_selection_parks.coord[,1],
            lat = first_selection_parks.coord[,2],
            radius = 6,fillColor = "violet",stroke = FALSE, fillOpacity = 0.9,
            popup = paste("Park:", first_selection_parks$distName1,"<br>",
                          "Pick no:", 1:length(first_selection_parks$distName1),"<br>",
                          "Reduced pop*dist:",first_selection_parks$dist_reduction,"%")
          ) %>% 
          # add distances from lsoa to old events
          addPolygons(data = poly_df, group = "Distance",
                      color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                      weight = 0.1, fillOpacity = 0.5,fillColor = ~pal_dist(mn_dstn),
                      highlight = highlightOptions(
                        weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE),
                      popup = paste(
                        poly_df$name,"<br>",
                        "Nearest event:", poly_df$nrst_vn,"<br>",
                        "Distance: ",poly_df$mn_dstn," km <br>",
                        "new Distance: ",poly_df$new_mn_dstn," km <br>",
                        "SIMD score:", poly_df$a,"<br>",
                        "Pop density:", poly_df$pp_dnst,"<br>",
                        "Population:", poly_df$pop)
                      ) %>%
          # add distances from lsoa to new and old events
          addPolygons(data = poly_df, group = "new Distance",
                      color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                      weight = 0.1, fillOpacity = 0.5,fillColor = ~pal_dist_new(new_mn_dstn),
                      highlight = highlightOptions(
                        weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE),
                      popup = paste(
                        poly_df$name,"<br>",
                        "Nearest event:", poly_df$nrst_vn,"<br>",
                        "Distance: ",poly_df$mn_dstn," km <br>",
                        "new Distance: ",poly_df$new_mn_dstn," km <br>",
                        "SIMD score:", poly_df$a,"<br>",
                        "Pop density:", poly_df$pp_dnst,"<br>",
                        "Population:", poly_df$pop)
                      ) %>%
          # add candidate parks as parks
          addPolygons(data = selected_parks_any, group = "candidate parks",
                      color = "gray", smoothFactor = 0,stroke = T,opacity = 0.5,
                      weight = 0.1, fillOpacity = 0.5,fillColor = "cyan",
                      highlight = highlightOptions(weight = 1,color = "cyan",opacity = 1,
                                                   bringToFront = FALSE,sendToBack = TRUE),
                      popup = paste("Park:", selected_parks_any$distName1,"<br>")
                      ) %>%
          # add legend
          addLegend("bottomleft", pal = pal_dist, values = poly_df$mn_dstn,
                    opacity = 0.7, title="Distance quintiles",group = c("Distance"),
                    labFormat = labelFormat(
                      suffix  = "km") 
                    )%>%
          addLegend("bottomleft", pal = pal_dist_new, values = poly_df$new_mn_dstn,
                    opacity = 0.7, title="Distance quintiles",group = c("new Distance"),
                    labFormat = labelFormat(
                      suffix  = "km") 
          )%>%
          # hide all layers except old parkrun events
          hideGroup("Distance") %>%
          hideGroup("new Distance") %>%
          hideGroup("candidate parks") %>%
          hideGroup("new Events (consecutive)") %>%
          hideGroup("new Events (first)")
      
        
        # map
        
        # save the map as a stand-aline html file
        htmlwidgets::saveWidget(map,file="index.html",title = "parkrun, access & equity",selfcontained = F)
        
        

  
  ### ANALYSIS
        # tbc