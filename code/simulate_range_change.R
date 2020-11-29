# read graphs and calculate range change
library(data.table); library(igraph); library(raster)
source("code/functions_igraph.R")

cpu_info <- as.numeric(commandArgs(T))
print(cpu_info)
id_number <- cpu_info[1]

if(id_number == 1) ids <- 1:4
if(id_number == 2) ids <- 5:8
if(id_number == 3) ids <- 9:12
if(id_number == 4) ids <- 13:16
fnames <- list.files("outputs/", full.names = T, pattern="graph_undir")
range_names <- unique(gsub(".*_undir_(.*)", "\\1", fnames))

# get highest forested elevation
ele_gain <- 365
n_steps <- 10
increment <- ele_gain/n_steps

for(i in ids) {
    # read in the ith graph
    graph_i <- readRDS(fnames[i])
    # get highest forested elevation
    max_ele <- max(graph_i$df_vertices$ele)
    upr_vec <- seq(1100, max_ele - ele_gain, len=20)
    
    # iterate across starting ranges
    for(j in 1:length(upr_vec)) {
        # get range accessed through movement, and range that is movement free
        new_range <- get_new_range(graph_i, upr_limit = upr_vec[j], range_width = 800, increment, n_steps)
        fm_range <- get_fm_range(graph_i, upr_limit = upr_vec[j], range_width = 800, increment, n_steps)
        
        # save
        save_name <- paste0("outputs/rangeshifts/", range_names[i], "_upr", upr_vec[j], "_rw800.rds")
        saveRDS(list(range_move = new_range, range_fm = fm_range), save_name)
        
        # summarise range changes
        range_change <-  new_range %>% 
            mutate(abs_t = abs(t), 
                   t0 = sum(sign(t) == 0) + sum(sign(t) == -1)) %>% 
            group_by(abs_t) %>%
            summarise(t0 = unique(t0), 
                      gain = sum(sign(t) == 1), 
                      loss = sum(sign(t) == -1), 
                      diff = gain - loss) 
        
        range_change_fm <-  fm_range %>% 
            mutate(abs_t = abs(t), 
                   t0 = sum(sign(t) == 0) + sum(sign(t) == -1)) %>% 
            group_by(abs_t) %>%
            summarise(t0 = unique(t0), 
                      gain = sum(sign(t) == 1), 
                      loss = sum(sign(t) == -1), 
                      diff = gain - loss) 
        range_change_summary <- full_join(range_change, range_change_fm, 
                                          by=c("abs_t", "t0"), 
                                          suffix=c("_move", "_fm")) %>%
            mutate(upr = upr_vec[j], 
                   range_width = 800, 
                   range = range_names[i])
        
        save_name <- paste0("outputs/rangeshift_summaries/", range_names[i], "_upr", 
                            upr_vec[j], "_rw800.rds")
        
        saveRDS(range_change_summary, save_name)
    }    
}