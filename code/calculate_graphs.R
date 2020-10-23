# script to generate graphs & rangeshifts

# get elevations and treecover for mountain ranges
fnames <- list.files("data/mountain_ranges/", full.names = T, pattern="60m")
print(length(fnames))
# housekeeping
library(raster); library(dplyr); library(data.table); library(igraph)
source("code/functions_igraph.R")

# colour palette
# cols <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", 
#                            "#FEE090", "#FDAE61", "#F46D43", "#D73027"))

range_names <- unique(gsub(".*60m_10kbuffer(.*)_2020.*", "\\1", fnames))

# should just run overnight- biggest mountain ranges will take 15 minutes or so
# 40GB RAM and see how it gets on. If expand to tropics in general then maybe 
# worth looping in parallel, but still talking about what, 100? 150? ranges, so 
# still only 25 hours run time (and actually most mountains will be substantially 
# faster)
for (i in 1:length(range_names)) {
    range_i <- range_names[i]
    fnames_i <- fnames[grepl(range_i, fnames)]
    ele <-  raster(fnames_i[1])
    tc  <-  raster(fnames_i[2])
    
    # sanity check: layers must have same dimensions
    if(any(dim(ele) != dim(tc))) stop("forest and elevation have different dimensions")
    
    # classify forest
    tc_classified <- tc
    tc_classified[tc_classified < 50] <- 0
    tc_classified[tc_classified >= 50] <- 1
    
    # get permitted cells (adding 0-valued margin)
    permitted_tc <- as.matrix(tc_classified)#[,1:1000]
    permitted_tc[get_margins(permitted_tc)] <- 0
    # permitted_cells <- which(permitted_tc[] == 1)
    ele_vec <- as.numeric(as.matrix(ele))
    n_row <- dim(permitted_tc)[1]
    n_col <- dim(permitted_tc)[2]
    
    # remove everything below 200m 
    to_remove <- which(ele_vec < 200)
    permitted_tc[to_remove] <- 0
    permitted_cells <- which(permitted_tc[] == 1)
    

    time_graph <- system.time(
        graph_out <- generate_graph_undir(ele_vec, permitted_cells, n_row, n_col, rule=1)
    )
    
    # save graph (really slow to write but will save memory overhead in future 
    # and is fast to read)
    graph_out$time <- time_graph
    saveRDS(graph_out, paste0("outputs/graph_undir_", range_i))
    rm(graph_out)
    
    # now generate directed graph (for climate connectivity metrics)
    time_graph <- system.time(
        graph_out <- generate_graph(ele_vec, permitted_cells, n_row, n_col, rule=1)
    )
    
    # save graph (really slow to write but will save memory overhead in future 
    # and is fast to read)
    graph_out$time <- time_graph
    saveRDS(graph_out, paste0("outputs/graph_dir_", range_i))
    rm(graph_out)
}
