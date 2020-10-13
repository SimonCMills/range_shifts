# Extract treecover metrics for buffer around points 

## Packages ----
library(rgee)

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")
countries <- ee$FeatureCollection("USDOS/LSIB/2013")
countries_subset <- countries$filter(ee$Filter$inList(opt_leftField = 'cc', 
                                                      opt_rightValue = list('CO', 'EC')))
mountains <- ee$FeatureCollection("users/scmills/GMBAMountainInventory_v1_2-World")

# get Andes
Andes <- mountains$filterBounds(countries_subset);

## Extract elevations ----
featlist <- Andes$getInfo()["features"]

for (i in 1:length(featlist$features)) {
    # get single feature from the feature list
    feat_i <- featlist$features[[i]]
    # clip elevation to this single mountain range
    ele_i <- ALOS$clip(ee$Feature(feat_i))
    
    # extract name
    name_i <-   feat_i$properties$Name
    # export
    task_i <- ee_image_to_drive(image = ele_i, 
                      description = paste0(save_prefix, gsub(" ", "_", name_i)), 
                      folder = "rgee_exports", 
                      scale = 60)
    task_i$start()
    cat(paste0("\n\n\n", name_i))
    ee_monitoring(task_i)
}
