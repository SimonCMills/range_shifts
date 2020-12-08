# Read in metada generated from raster files, and merge with mountain range 
# polygons to generate final metadata file. 
# note: had to split this code up: in the end needed running on the cluster, and 
# installing sf there was proving a faff. 

#
file.copy("Y:/edwards_lab1/User/bo1scm/range_shifts/data/metadata_intermediate_version.rds", 
          "data/")

df_rangeinfo <- readRDS("data/metadata_intermediate_version.rds")
df_rangeinfo
# join with range polygons
mountain_ranges <- readRDS("data/mountain_polygons.rds") %>%
    mutate(range = gsub("\\(|\\)", "", Name) %>%
               gsub("ü", "u", .) %>% 
               gsub("\\/", "--", .) %>%
               gsub(" ", "_", .))

mountains2 <- left_join(mountain_ranges, df_rangeinfo) %>%
    mutate(pct_forest_300 = round(area_forest/area_gt300, 3)) %>%
    dplyr::select(Country:area_gt300, pct_forest_300, everything())

saveRDS(mountains2, "data/range_metadata.rds")

## 124 in tropics
mountains2 %>%
    filter(prop_in_tropics > .5) %>%
    nrow

## 100 in tropics & max ele >= 1500
mountains2 %>%
    filter(prop_in_tropics > .5, max_val >= 1500) %>%
    nrow

## 94 in tropics & max ele >= 1500
mountains2 %>%
    filter(prop_in_tropics > .5, max_val >= 1500, pct_forest_300 >= .05) %>%
    nrow
