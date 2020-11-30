# script to generate graphs (both directed- for simulating range changes- and 
# undirected- for calculating climate connectivity)

# get elevations and treecover for mountain ranges
fnames <- list.files("data/mountain_ranges/", full.names = T, pattern="60m_10kbuffer")
print(paste0("n_fnames: ", length(fnames)))

# housekeeping
library(raster); library(dplyr); library(data.table); library(igraph)
source("code/functions_igraph.R")

range_names <- unique(gsub(".*60m_10kbuffer(.*)_2020.*", "\\1", fnames))

# which graphs are already calculated? (avoid re-calculating)
graphs_calculated <- list.files("./outputs/", "graph_undir") %>%
    gsub("graph_undir_", "", .)

# range metadata
# arrange to calculate smallest graphs first, and remove those already calculated
# Note: record maximum RAM overhead after run?
range_mtd <- readRDS("data/range_metadata.rds") %>% 
    dplyr::select(-geometry) %>%
    as_tibble %>%
    filter(max_val >= 1500, prop_in_tropics >= .5, !(range %in% graphs_calculated))  %>%
    arrange(desc(ncell))

# loop across ranges and calculate graphs
for (i in 1:length(range_mtd$range)) {
    range_i <- range_mtd$range[i]
    print(range_i)
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
    permitted_tc <- as.matrix(tc_classified)
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
}
