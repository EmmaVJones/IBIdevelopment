# pretty close
pal <- colorNumeric(c("#ff0000", #23
                      "#fde9aa",#71?
                      "#e8d1d1",#21?
                      "#4f7e48",#41
                      "#e29e8c",#22 # good
                      
                      "#dcca8f",#52
                      "#38814e",#43?
                      
                      "#FFFFCC",
                      "#fbf65d"#81
), 
values(landcover2016),
na.color = "transparent")

sort(unique(values(landcover2016)))

# missing red
pal <- colorNumeric(c("#e29e8c",#22 # good"
                      "#fde9aa",#71?
                      
                      "#e8d1d1",#21?
                      "#4f7e48",#41
                      "#FFFFCC",
                      "#e29e8c",#22 # good"
                      "#dcca8f",#52
                      "#38814e",#43?
                      
                      "#fbf65d"#81
), 
sort(unique(values(landcover2016))),
#values(landcover2016),
na.color = "transparent")




pal <- colorFactor(c('#cc807a', 
                      '#d6541c',
                      '#ed0707',
                      '#0bb32c',
                      '#0b8005',
                      '#97e657',
                      '#d9ad29',
                      '#f7d368',
                      '#ede505'), 
                    c('Developed, Open Space', 'Developed, Low Intensity', 
                      'Developed, Medium Intensity', 'Deciduous Forest',
                      'Evergreen Forest', 'Mixed Forest', 'Shrub/Scrub',
                      'Grassland/Herbaceous','Pasture/Hay'),
                    na.color = 'transparent')

qpal <- colorQuantile(c('#cc807a', 
                        '#d6541c',
                        '#ed0707',
                        '#0bb32c',
                        '#0b8005',
                        '#97e657',
                        '#d9ad29',
                        '#f7d368',
                        '#ede505'), unique(landcover2016), n = 9)


CreateWebMap(maps = c("Topo","Imagery","Hydrography"), collapsed = TRUE) %>%
  setView(-80.02, 37.41, zoom=13) %>%
  addRasterImage(landcover2016, colors = qpal, opacity = 0.8, group = 'NLCD 2016') %>%
  addPolygons(data= HAM,  color = 'black', weight = 1,
              fillColor='blue', fillOpacity = 0.3,stroke=0.1,
              group="Watershed") %>%
  addPolylines(data = HAM_NHD,  color = 'blue', weight =3,
               group="1:100k NHD", label = '1:100k NHD',
               popup=leafpop::popupTable(HAM_NHD, zcol=c('Strahler Order (1:100k NHD):'))) %>%
  addCircleMarkers(data = HAMsite, color='orange', fillColor='black', radius = 5,
                   fillOpacity = 0.5,opacity=0.5,weight = 1,stroke=T, group="Station",
                   label = ~StationID,
                   popup=paste('StationID:', HAMsite$StationID)) %>%
  addLegend(pal = pal, values = c('Developed, Open Space', 'Developed, Low Intensity', 
                                   'Developed, Medium Intensity', 'Deciduous Forest',
                                   'Evergreen Forest', 'Mixed Forest', 'Shrub/Scrub',
                                   'Grassland/Herbaceous','Pasture/Hay'),
            title = "NLCD 2016 Classification", group = 'NLCD 2016') %>% 
  hideGroup('NLCD 2016') %>%
  addLayersControl(
    baseGroups=c("Topo","Imagery","Hydrography"),
    overlayGroups = c('Station',"1:100k NHD", 'Watershed','NLCD 2016'),   
    options=layersControlOptions(collapsed=T),
    position='topleft')





qpal <- colorFactor(c('#cc807a', 
                      '#d6541c',
                      '#ed0707',
                      '#0bb32c',
                      '#0b8005',
                      '#97e657',
                      '#d9ad29',
                      '#f7d368',
                      '#ede505'), unique(landcover2016),na.color = 'transparent')

qpal <- colorBin(c('#cc807a', 
                   '#d6541c',
                   '#ed0707',
                   '#0bb32c',
                   '#0b8005',
                   '#97e657',
                   '#d9ad29',
                   '#f7d368',
                   '#ede505'), unique(landcover2016),bins=9,na.color = 'transparent')


