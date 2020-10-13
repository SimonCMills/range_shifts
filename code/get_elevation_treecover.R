# Extract elevation and treecover for each mountain range in selection (currently 
# running for northern Andes)

## Packages ----
library(rgee)

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")$select("treecover2000")
countries <- ee$FeatureCollection("USDOS/LSIB/2013")
countries_subset <- countries$filter(ee$Filter$inList(opt_leftField = 'cc', 
                                                      opt_rightValue = list('CO', 'EC')))
mountains <- ee$FeatureCollection("users/scmills/GMBAMountainInventory_v1_2-World")

# get Andes
Andes <- mountains$filterBounds(countries_subset);

## Extract elevations ----
# note: will take a lot longer to run than the loop speed implies; this just 
# sets off the exports which will run in the background
featlist <- Andes$getInfo()["features"]
save_prefix_tc <- "tc2000_60m_"
save_prefix_ele <- "ele_jaxa_60m_"

for (i in 1:length(featlist$features)) {
    # get single feature from the feature list
    feat_i <- featlist$features[[i]]
    # clip elevation to this single mountain range
    ele_i <- ALOS$clip(ee$Feature(feat_i))
    tc_i <- tc$clip(ee$Feature(feat_i))
    
    # extract name (need to remove unrecognised character)
    name_i <- gsub(" ", "_", gsub("�", "",feat_i$properties$Name))
    # export
    # ..elevation
    task_i <- ee_image_to_drive(image = ele_i, 
                                description = paste0(save_prefix_ele, name_i), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 3e8)
    
    task_i$start()
    # ..treecover
    task_i <- ee_image_to_drive(image = tc_i, 
                                description = paste0(save_prefix_tc, name_i), 
                                folder = "rgee_exports", 
                                scale = 60, maxPixels = 3e8)
    
    task_i$start()
    # counter
    cat(paste0("\n", name_i))
}
