# create dataframe of range metadata:
# min&max elevation, spatial extent (of raster), proportion in neotropics
# 
# NOTE: run copy_to_RD.sh first (which will transfer rgee files, and then update
# research drive)

# housekeeping
library(raster); library(dplyr); library(data.table)

## get raster info ----
fnames_ele <- list.files("data/mountain_ranges/", full.names=T, pattern = "ele_jaxa")
fnames_tc <- list.files("data/mountain_ranges/", full.names=T, pattern = "tc2000")
range_names <- gsub(".*60m_10kbuffer(.*)_2020.*", "\\1", fnames_ele)

df_rangeinfo <- data.table(range = range_names)
for(i in 1:length(range_names)) {
    print(i)
    ele <- raster(fnames_ele[i])
    tc <- raster(fnames_tc[i])
    
    # sanity check
    if(any(dim(ele) != dim(tc))) stop("forest and elevation have different dimensions")
    
    # elevation range
    max_val_i <- max(ele[]) 
    min_val_i <- min(ele[]) 
    
    # classify tc (50% threshold)
    tc_classified <- tc
    tc_classified[tc_classified < 50] <- 0
    tc_classified[tc_classified >= 50] <- 1
    
    
    amt_area_i <- length(which(ele[] > 300))
    amt_forest_i <- length(which(tc_classified[ele[] > 300] == 1))
    # df_rangeinfo[i, max_for_ele := max(ele[tc_classified[] == 1])]
    
    df_rangeinfo[i, `:=`(max_val = max_val_i, 
                         min_val = min_val_i, 
                         nrow = nrow(ele), 
                         ncol=ncol(ele), 
                         ncell=ncell(ele), 
                         area_forest = amt_forest_i, 
                         area_gt300 = amt_area_i)]
}

saveRDS(df_rangeinfo, "data/metadata_intermediate_version.rds")
