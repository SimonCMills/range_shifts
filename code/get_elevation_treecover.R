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

ele <- ALOS$reproject(tc$projection())
# countries <- ee$FeatureCollection("USDOS/LSIB/2013")
# countries_subset <- countries$filter(ee$Filter$inList(opt_leftField = 'cc', 
#                                                       opt_rightValue = list('CO', 'EC')))
# major_watersheds <- ee$FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_3")
mountains <- ee$FeatureCollection("users/scmills/GMBAMountainInventory_v1_2-World")
# get Andes
mountains_sub <- mountains$filterBounds(geometry);

buffer <- function(feature) feature$buffer(10000)
mountains_buffered <- mountains_sub$map(buffer)

## Extract elevations ----
# note: will take a lot longer to run than the loop speed implies; this just 
# sets off the exports which will run in the background
featlist <- mountains_buffered$getInfo()["features"]
save_prefix_tc <- "tc2000_60m_10kbuffer"
save_prefix_ele <- "ele_jaxa_60m_10kbuffer"

# In total, 123 mountain ranges intersect with tropics.
for (i in 1:length(featlist$features)) {
    # get single feature from the feature list
    feat_i <- featlist$features[[i]]
    # clip elevation to this single mountain range
    ele_i <- ele$clip(ee$Feature(feat_i))
    tc_i <- tc$clip(ee$Feature(feat_i))
    
    # extract name (need to remove unrecognised character)
    name_i <- gsub(" ", "_", gsub("�", "",feat_i$properties$Name))
    # export
    save_name <- gsub("\\(|\\)", "", name_i)
    save_name <- gsub("ü", "u", save_name)
    save_name <- gsub("\\/", "--", save_name)
        
    # ..elevation
    task_i <- ee_image_to_drive(image = ele_i, 
                                description = paste0(save_prefix_ele, save_name), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 3e8)
    
    task_i$start()
    # ..treecover
    task_i <- ee_image_to_drive(image = tc_i, 
                                description = paste0(save_prefix_tc, save_name), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 3e8)
    
    task_i$start()
    # counter
    cat(paste0("\n", name_i))
}

# calculate proportion overlap with tropics, and save polygons
mountains_sf <- ee_as_sf(mountains_sub)
geom_sf <- ee_as_sf(geometry)

summ <- mountains_sf %>% 
    mutate(total_area = st_area(.)) %>%
    st_intersection(., geom_sf) %>% 
    mutate(intersect_area = st_area(.), 
           prop_in_tropics = intersect_area/total_area)

mountains_final <- summ %>% 
    as_tibble() %>%
    dplyr::select(Name, prop_in_tropics) %>%
    mutate(prop_in_tropics = as.numeric(prop_in_tropics)) %>%
    left_join(mountains_sf, .)

saveRDS(mountains_final, "data/mountain_polygons.rds")
