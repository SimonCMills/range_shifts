# Extract elevation and treecover for each mountain range in selection (currently 
# running for northern Andes)

## Packages ----
library(rgee); library(sf); library(dplyr)

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")$select("treecover2000")

# rectangle designating neotropics
geometry <- ee$Geometry$Rectangle(
    coords = c(-113,-23.5,-34,23.5),
    proj = "EPSG:4326",
    geodesic = FALSE
)

# Peruvian oriental cordillera is disjunct with a northern and southern half- 
# which needs splitting into two mountain ranges. Polygons for intersecting by:
north_split <- ee$Geometry$Rectangle(
    coords = c(-79, -8.6, -75, -4),
    proj = "EPSG:4326",
    geodesic = FALSE
)

south_split <- ee$Geometry$Rectangle(
    coords = c(-79, -25, -60, -9),
    proj = "EPSG:4326",
    geodesic = FALSE
)

# Map$addLayer(north_split) +
#     Map$addLayer(mountains) +
#     Map$addLayer(south_split)

ele <- ALOS$reproject(tc$projection())
mountains <- ee$FeatureCollection("users/scmills/GMBAMountainInventory_v1_2-World")

# get Andes
mountains_sub <- mountains$filterBounds(geometry);

# mountain subset, excluding Cordillera Oriental de Peru
mountains_sub2 <- mountains_sub$filter(ee$Filter$neq("Name", "Cordillera Oriental Peru Bolivia"))

# extract just the CO Peru and create feature
CO_Peru <- mountains_sub$filter(ee$Filter$eq("Name", "Cordillera Oriental Peru Bolivia"))
test <- ee$Feature(CO_Peru$geometry())

# intersect with north south geometries, and rejoin back into feature collection
CO_Peru_north <- test$intersection(north_split, ee$ErrorMargin(1))
CO_Peru_north2 <- CO_Peru_north$set("Name", "Cordillera Oriental Peru Bolivia- North")
CO_Peru_south <- test$intersection(south_split, ee$ErrorMargin(1))
CO_Peru_south2 <- CO_Peru_north$set("Name", "Cordillera Oriental Peru Bolivia- South")
new_msub <- mountains_sub2$merge(CO_Peru_north2)
new_msub2 <- new_msub$merge(CO_Peru_south2)

# buffer
buffer <- function(feature) feature$buffer(10000)
mountains_buffered <- new_msub2$map(buffer)

## Extract elevations ----
# note: will take a lot longer to run than the loop speed implies; this just 
# sets off the exports which will run in the background
featlist <- mountains_buffered$getInfo()["features"]
save_prefix_tc <- "tc2000_60m_10kbuffer"
save_prefix_ele <- "ele_jaxa_60m_10kbuffer"

full_name <- lapply(featlist$features, function(x) x$properties$Name) %>%
    unlist

simple_name <- full_name %>%
    gsub("�", "", .) %>%
    gsub(" ", "_", .) %>%
    gsub("\\(|\\)", "", .) %>%
    sub("ü", "u", .) %>%
    gsub("\\/", "--", .)

already_downloaded <- list.files("data/mountain_ranges/") %>%
    gsub(".*buffer(.*)_2020.*", "\\1", .) %>% 
    unique

to_download <- which(!(simple_name %in% already_downloaded))
if(length(to_download) == 0) stop("no additional ranges to download")

# loop through undownloaded ranges and export to drive
for (i in to_download) {
    # get single feature from the feature list
    feat_i <- featlist$features[[i]]
    # clip elevation to this single mountain range
    ele_i <- ele$clip(ee$Feature(feat_i))
    tc_i <- tc$clip(ee$Feature(feat_i))
    
    # export
    save_name <- simple_name[i]
        
    # ..elevation
    task_i <- ee_image_to_drive(image = ele_i, 
                                description = paste0(save_prefix_ele, save_name), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 1e9)
    
    task_i$start()
    # ..treecover
    task_i <- ee_image_to_drive(image = tc_i, 
                                description = paste0(save_prefix_tc, save_name), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 1e9)
    
    task_i$start()
    # counter
    cat(paste0("\n", save_name))
}

# Save range polygons ----
# calculate proportion overlap with tropics, and save polygons
mountains_sf <- ee_as_sf(new_msub2)
geom_sf <- ee_as_sf(geometry)

# note only approximate intersection with lat-lon
summ <- mountains_sf %>% 
    mutate(total_area = st_area(.)) %>%
    st_intersection(., geom_sf) %>% 
    mutate(intersect_area = st_area(.), 
           prop_in_tropics = intersect_area/total_area)

mountains_final <- summ %>% 
    as_tibble() %>%
    dplyr::select(Name, prop_in_tropics) %>%
    mutate(prop_in_tropics = as.numeric(prop_in_tropics)) %>%
    left_join(mountains_sf, .) %>%
    left_join(., tibble(Name = full_name, range = simple_name))

# ggplot(mountains_final) +
#     geom_sf(aes(fill=prop_in_tropics)) +
#     scale_fill_gradient2(midpoint=.5)

saveRDS(mountains_final, "data/mountain_polygons.rds")
