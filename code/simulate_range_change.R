# read graphs and calculate range change

# housekeeping
library(data.table); library(igraph); library(raster)
source("code/functions_igraph.R")

# pass in range width
cpu_info <- as.numeric(commandArgs(T))
rwidth <- cpu_info[1]

# get fnames (and order by file size- so smallest are run first)
fnames <- list.files("outputs/", full.names = T, pattern="graph_undir")
fnames <- fnames[order(file.info(fnames)$size)]

# get range names
range_names <- unique(gsub(".*_undir_(.*)", "\\1", fnames))

# get highest forested elevation
lapse_rate <- 0.0055 # per metre
ele_gain <- round(2/lapse_rate, 0)
n_steps <- 20
increment <- ele_gain/n_steps

for(i in 1:length(range_names)) {
    # read in the ith graph
    graph_i <- readRDS(fnames[i])
    print(range_names[i])
    # get highest forested elevation
    max_ele <- max(graph_i$df_vertices$ele)
    # each species' starting upper limit
    upr_vec <- seq(1100, max_ele - ele_gain, len=20)
    
    # iterate across starting ranges
    for(j in 1:length(upr_vec)) {
        print(paste0("species: ", j))
        # get range accessed through movement, and range that is movement free
        new_range <- get_new_range(graph_i, upr_limit = upr_vec[j], range_width = rwidth, increment, n_steps)
        fm_range <- get_fm_range(graph_i, upr_limit = upr_vec[j], range_width = rwidth, increment, n_steps)
        
        # save
        save_name <- paste0("outputs/rangeshifts/", range_names[i], "_upr", upr_vec[j], "_rw", rwidth,".rds")
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
                   range_width = rwidth, 
                   range = range_names[i])
        
        save_name <- paste0("outputs/rangeshift_summaries/", range_names[i], "_upr", 
                            upr_vec[j], "_rw", rwidth, ".rds")
        
        saveRDS(range_change_summary, save_name)
    }    
}
