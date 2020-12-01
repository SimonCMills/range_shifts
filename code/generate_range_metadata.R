# create dataframe of range metadata (min&max elevation, proportion in neotropics)

# housekeeping
library(raster); library(dplyr); library(sf)

## (1) copy all unique ranges ----
# some duplicated downloads which need removing
# get elevations and treecover for mountain ranges
fnames <- list.files("D:/Drive/rgee_exports/", full.names = T, pattern="60m_10kbuffer")
print(paste0("n_fnames: ", length(fnames)))

range <- gsub("60m_10kbuffer", "", fnames) %>%
    gsub(".*_exports/(.*)_2020.*", "\\1", .)

# 119 mountain ranges in ranges in total
length(range[!duplicated(range)])/2

# copy
file.copy(fnames[!duplicated(range)], "data/mountain_ranges/")

## get raster info ----
## (max ele and dimensions)
fnames_2 <- list.files("data/mountain_ranges/", full.names=T, pattern = "ele_jaxa")
range_2 <- gsub(".*60m_10kbuffer(.*)_2020.*", "\\1", fnames_2)

df_rangeinfo <- data.table(range = range_2)
for(i in 1:length(range_2)) {
    r <- raster(fnames_2[i])
    max_val_i <- max(r[]) 
    min_val_i <- min(r[]) 
    df_rangeinfo[i, `:=`(max_val = max_val_i, 
                         min_val = min_val_i, 
                         nrow = nrow(r), 
                         ncol=ncol(r), 
                         ncell=ncell(r))]
}

# 95 ranges have max elevation >= 1500
nrow(df_rangeinfo[max_val >= 1500])

## get raster info (2) ----
## is there treecover? (rules out a couple of ranges in Chile)
fnames_3 <- list.files("data/mountain_ranges/", full.names=T, pattern = "tc2000")

for(i in 1:length(range_2)) {
    ele <- raster(fnames_2[i])
    tc <- raster(fnames_3[i])
    # sanity check: layers must have same dimensions
    if(any(dim(ele) != dim(tc))) stop("forest and elevation have different dimensions")
    tc_classified <- tc
    tc_classified[tc_classified < 50] <- 0
    tc_classified[tc_classified >= 50] <- 1
    
    amt_area_i <- length(which(ele[] > 300))
    amt_forest_i <- length(which(tc_classified[ele[] > 300] == 1))
    df_rangeinfo[i, amt_forest := amt_forest_i]
    df_rangeinfo[i, amt_gt300 := amt_area_i]
}

# 84 ranges have max elevation >= 1500 & have >=10% forest  
df_rangeinfo[, pct_forest_300 := round(amt_forest/amt_gt300, 3)]
df_rangeinfo[max_val >=1500 & pct_forest_300 > .1]

mountain_ranges <- readRDS("data/mountain_polygons.rds") %>%
    mutate(range = gsub("\\(|\\)", "", Name) %>%
               gsub("ü", "u", .) %>% 
               gsub("\\/", "--", .) %>%
               gsub(" ", "_", .))

mountains2 <- left_join(mountain_ranges, df_rangeinfo)
saveRDS(mountains2, "data/range_metadata.rds")

## 119 in tropics
mountains2 %>%
    filter(prop_in_tropics > .5) %>%
    nrow

## 93 in tropics & max ele >= 1500
mountains2 %>%
    filter(prop_in_tropics > .5, max_val >= 1500) %>%
    nrow

## 85 in tropics & max ele >= 1500
mountains2 %>%
    filter(prop_in_tropics > .5, max_val >= 1500, pct_forest_300 >= .1) %>%
    nrow


