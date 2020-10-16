# Growing set of functions for modelling spread through a network (built on igraph)
# Needs consolidation, and finalising of testing (though so far everything seems
# to be performing as expected). 
# Currently: 
# (1) generate_graph: create network from permitted cells, elevations, and movement
# rule. 
#   TO DO: think the most efficient way forward is to have a single graph 
#   generated and saved once, and then anything subsequent to simply refer to this
#   Currently most functions will rebuild graph each time. 
#
# (2) get_range: from graph, get all accessible cells and append their elevations. 
# Will be used to estimate the #cells in range through time
#
# (3) get_max_ele:get the maximum elevation accessible from all cells
#
# - a bunch of possibly redundant stuff. 
#
# (4) get_adjacent: get adjacencies under varying gap-crossing rules- currently
# written manually from 0-2 cell gap-crossing; possibly needs to be able to 
# generate adjacencies for an arbitrarily large gap-crossing number. 
# 
#
# (5) get_margins: creates impassable border at edge of raster (to prevent 
# adjacencies wrapping around to other side of matrix)
#
# Note: bizarre bug in igraph function that will generate duplicate vertices for
# *some* values of xe+05. So far have seen x values of 5, 7, 9, and 10. In all
# cases observed so far, the e+05 version only has incoming edges, while the
# x00000 version only has outgoing edges. The workaround is to identify these
# failure cases, reroute the incoming edges to the x00000 version (creating a
# correctly connected vertex in the graph), and then remove the e+05 version.
# 
# It apparently only affects a tiny number of cases (e.g. 4 duplicated nodes in
# a network of ~2.7e06 vertices), so is going to have an ignorable impact on
# results, but fixing regardless. Need to write a reproducible example and submit
# as an error.
# 
# The check to catch further cases that do not follow the observed rule is to
# compare the number of vertices listed in dat2 (from and to columns) to the
# length of the names in the graph
#

# dependencies
library(dplyr); library(data.table); library(igraph)

# generate the graph of all cells (travelling up- or along-slope, i.e from smaller
# values to larger) 
generate_graph <- function(ele_vec, permitted_cells, n_row, n_col, rule) {
  # generate adjacencies
  adj <- data.table(from = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
                    to = get_adjacent(permitted_cells, n_row, n_col, rule))
  # remove non-permitted adjacencies
  dat <- adj[to %in% permitted_cells,] 
  # this drops cells without *any* permitted adjacencies, but these are needed 
  # for calculating occupied area (but don't need including in graph)
  disc_cells <- data.table(
    id_cell = permitted_cells[!(permitted_cells %in% unique(dat$from))]
    )
  disc_cells[, id_vertex := NA]
  disc_cells[, ele := ele_vec[id_cell]]
  disc_cells[, disc := TRUE]
  
  # get elevations for from and to vertices
  dat[,c("from_ele", "to_ele") := list(ele_vec[from], ele_vec[to])]
  
  # if adjacent nodes don't have same elevation, drop the lower elevation 
  # retention index; if they are the same elevation retain both (as path needs
  # to be bidirectional)
  dat2 <- rbind(dat[from_ele > to_ele,], dat[from_ele == to_ele,])
  class(dat2$from)
  # generate graph from adjacencies
  g <- graph_from_data_frame(dat2[,list(from, to)], directed = T)
  
  # generate lookup table
  # note: This won't ever matter for my use, but treating as an integer will mean 
  # indexing tops out somewhere around 2e+09 (as 32-bit)
  vertex_names <- names(V(g))
  graph_vertices <- data.table(id_cell = as.integer(vertex_names[]), 
                               id_vertex = 1:length(vertex_names)) 
  graph_vertices[,ele := ele_vec[id_cell]]
  graph_vertices[,disc := FALSE]
  
  ##############################################################################
  # workaround for igraph bug: get all incoming and outgoing edges from vertices
  # that have been duplicated (incoming always associated with xe+ version and
  # outgoing always associated with x0 version). Identify the vertices associated
  # with incoming edges, create incoming edges from these vertices to the x0 
  # version of the duplicated nodes, and then delete the xe+ version from the 
  # graph. See more extensive discussion of the issue in unit testing markdown file. 
  
  # get duplicated vertices
  duplicated_vertices <- 
    graph_vertices[id_cell %in% graph_vertices[duplicated(id_cell),id_cell]]
  setorder(duplicated_vertices, id_cell, id_vertex)
  
  # x0.. version has the outgoing vertices
  incoming <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="in")
  # xe+.. version has incoming vertices
  outgoing <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="out")
  
  # retain only the x0 version: get the vertices reached from the xe+ version, 
  # add these edges to the x0 version, and then delete the xe+ vertices from the 
  # graph
  in_xe_ver <- incoming[lapply(incoming, length) != 0]
  
  # create vector of edge starts (as.integer coerces e+05 form to 00 form)
  in_vert_from <- rep(names(in_xe_ver), unlist(lapply(in_xe_ver, length))) %>%
    as.integer
  
  # vector of edge ends
  in_vert_to <- lapply(in_xe_ver, function(x) as.integer(names(x))) %>% unlist %>%
    as.integer
  
  # interleave these two vectors
  edges_to_add <- as.vector(rbind(in_vert_from, in_vert_to))
  
  # add new edges to graph
  g <- add_edges(g, edges_to_add)
  
  # remove xe+ vertices from graph
  g <- delete_vertices(g, 
                       duplicated_vertices$id_vertex[grepl("\\+", names(outgoing))])
  
  ## rerun lookup table code block
  # generate lookup table
  # note: This won't ever matter for my use, but treating as an integer will mean 
  # indexing tops out somewhere around 2e+09 (as 32-bit)
  vertex_names <- names(V(g))
  graph_vertices <- data.table(id_cell = as.integer(vertex_names[]), 
                               id_vertex = 1:length(vertex_names)) 
  graph_vertices[,ele := ele_vec[id_cell]]
  graph_vertices[,disc := FALSE]
  ##############################################################################
  
  # bind the dataframe of all the cells in the graph to the dataframe of the 
  # disconnected cells not in the graph
  for_return <- rbind(graph_vertices[!duplicated(id_cell),], disc_cells)
  
  # sanity check: this final dataframe containing the cells that are disconnected 
  # from the graph *and* all the cells in the graph, should be equal to the total
  # number of forested cells provided at the outset
  if(nrow(for_return) != length(permitted_cells)) stop("length mismatch")
  
  # return graph, data.table of vertices, and the two lists of ingoing and
  # outgoing edges from duplicated cells (errors should be caught by the above
  # check: but can visually check these lists after the event to confirm there
  # aren't additional failure cases)
  list(graph = g, 
       df_vertices = for_return[], 
       dup_in_edges = incoming, 
       dup_out_edges = outgoing)
}

# given the graph, associated vector of elevations, and cells in range at outset
# calculate all accessible future cells (and append associated elevations)
get_range <- function(graph_output, ele_vec, start_cells) {
  # unpack graph outputs
  graph_vertices <- graph_output$df_vertices
  
  # identify start vertices
  start_vertices <- graph_vertices[id_cell %in% start_cells & disc==FALSE, id_vertex]
  
  # run bfs
  # fix igraph behaviour to correctly store father vertices
  igraph_options(add.vertex.names=F)
  bfs_out <- bfs(graph_output$graph, root = start_vertices, unreachable = F, neimode="in")
  
  # now generate ordering for doing the membership calculation: NAs aren't reachable
  graph_vertices[as.numeric(bfs_out$order), ordering := 1:.N]
  graph_vertices[, start_cell := ifelse(id_cell %in% start_cells, T, F)]
  setorder(graph_vertices, ordering)
  
  # return vertices, dropping unreached cells (but retain start vertices)
  return(graph_vertices[!is.na(ordering) & disc == FALSE,])
}

# get maximum elevation accessible from each cell
get_max_ele <- function(ele_vec, permitted_cells, n_row, n_col, rule) {
    # generate adjacencies, then remove adjacencies that aren't permitted
  dat <- data.table(from = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
                    to = get_adjacent(permitted_cells, n_row, n_col, rule)) %>%
    .[to %in% permitted_cells] 

  # get elevations for from and to vertices
  dat[,c("from_ele", "to_ele") := list(ele_vec[from], ele_vec[to])]
  
  # if adjacent nodes don't have same elevation, drop the lower elevation 
  # retention index; if they are the same elevation retain both (as path needs
  # to be bidirectional)
  dat2 <- rbind(dat[from_ele > to_ele,], dat[from_ele == to_ele,])
  # order to create necessary ordering for dfs: dfs algorithm will jump to 
  # new roots in order of vertices: by giving edge df in order of elevation, 
  # this ensures it will always start from the next highest unreached vertex
  setorder(dat2, -from_ele, -to_ele)
  
  # generate graph from adjacencies
  g <- graph_from_data_frame(dat2, directed = T)
  
  # generate lookup table
  graph_vertices <- V(g)$name %>%
    data.table(id_cell = as.integer(.), id_vertex = 1:length(.), 
               max_ele=as.double(rep(NA, length(.)))) 
  graph_vertices[,ele := ele_vec[id_cell]]
  
  # run dfs
  # housekeeping: fix igraph behaviour to correctly store father vertices
  igraph_options(add.vertex.names=F)
  dfs_out <- dfs(g, root = 1, unreachable = T, neimode="out", dist=T, father=T)

  # add distances to df with original ordering
  graph_vertices[,father := dfs_out$father]
  # now generate ordering for doing the membership calculation
  graph_vertices[as.numeric(dfs_out$order), ordering := 1:.N]
  setorder(graph_vertices, ordering)
  # get membership and get maximum accessible elevation under upslope movement
  graph_vertices[,membership := cumsum(is.na(father))]
  graph_vertices[,max_ele := max(ele), by="membership"]
  
  return(graph_vertices[])
}

sim_new_range_igraph <- function(n_iter, ele_vec, lwr_vec, upr_vec, permitted_cells, 
                                    start_cells, n_row, n_col, rule) {
    ## part 1: get full graph of all cells that are connected to starting range
    ## by habitat- any cell that is not connected by habitat is never accessible
    dat <- data.frame(x1 = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
                      x2 = get_adjacent(permitted_cells, n_row, n_col, rule)) %>%
        filter(x2 %in% permitted_cells)
    graph <- graph_from_data_frame(dat)
    
    # identify membership of all cells 
    members <- components(graph)$membership %>%
        data.frame(membership = ., id_cell = as.integer(names(.)))
    # which cells are not connected to the starting cells (i.e range t[0])
    disconnected <- which(!(members$membership %in% 
                                members$membership[members$id_cell %in% start_cells]))
    
    # calculate the graph of cells connected to the starting range cells
    connected <- delete.vertices(graph, disconnected)
    # create df of all connected cells, appending elevations to this. 
    # these are the ids of cells that are connected to start cells: only search 
    # these for potential connections
    members_cnct <- components(connected)$membership %>%
        tibble(membership = ., id_cell = as.integer(names(.)), 
               eles = ele_vec[id_cell])
    
    ## part 2: iterate through time, updating range & recalculating graph
    # list to catch range[t]
    catch <- vector("list", n_iter)
    # initialise range[0]
    start_t <- start_cells
    
    for (t in 1:n_iter) {
        # get cell ids that fall outside elevational limts
        outside_range <- members_cnct$id_cell[members_cnct$eles <= lwr_vec[t] | 
                                                  members_cnct$eles> upr_vec[t+1]]
        # get vertex ids for these cells
        to_drop <- which(members_cnct$id_cell %in% outside_range)
        # drop these 
        graph_t <- delete.vertices(connected, to_drop)
        
        # recalculate membership
        members_t <- components(graph_t)$membership %>%
            data.frame(membership = ., id_cell = as.integer(names(.)))
        # get cell ids for all cells not isolated to starting range[t]
        linked_t <- members_t$id_cell[members_t$membership %in% 
                                          members_t$membership[members_t$id_cell %in% start_t]]
        # update range and store (dropping cells that are no longer in range)
        start_t <- linked_t[ele_vec[linked_t] > lwr_vec[t+1]]
        catch[[t]] <- start_t
    }
    return(catch)
} 

# contagious spread: not called in main function, but isolates the mechanism for 
# movement in a single time-slice 
spread_igraph <- function(permitted_vec, start_cells, n_row, n_col, rule=1) {
    # permitted_cells <- which(permitted_vec == 1)
    dat <- data.frame(x1 = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
                      x2 = get_adjacent(permitted_cells, n_row, n_col, rule)) %>%
        filter(x2 %in% permitted_cells)
    graph <- graph_from_data_frame(dat)
    members <- components(graph)$membership %>%
        data.frame(membership = ., id_cell = as.integer(names(.)))
    members$id_cell[members$membership %in% members$membership[members$id_cell %in% start_cells]]
}

# for each cell, get all adjacent cells
# this is a bit unwieldy, but works and is fast
get_adjacent <- function(cells, n_row, n_col, rule=1) {
    if (rule == 1) {
        adjacencies_i <- c(cells-n_row - 1, 
                           cells-n_row, 
                           cells-n_row+1, 
                           cells-1, 
                           cells+1, 
                           cells+n_row-1, 
                           cells+n_row, 
                           cells+n_row+1)
        return(adjacencies_i)
    }
    if(rule ==2) {
        adjacencies_i <- c(cells-n_row - 1, 
                           cells-n_row, 
                           cells-n_row+1, 
                           cells-1, 
                           cells+1, 
                           cells+n_row-1, 
                           cells+n_row, 
                           cells+n_row+1,
                           cells - 2*n_row - 2, 
                           cells - 2*n_row - 1, 
                           cells - 2*n_row, 
                           cells - 2*n_row + 1,
                           cells - 2*n_row + 2, 
                           cells - n_row - 2, 
                           cells - n_row + 2, 
                           cells - 2, 
                           cells + 2, 
                           cells + n_row - 2, 
                           cells + n_row + 2,
                           cells + 2*n_row - 2, 
                           cells + 2*n_row - 1, 
                           cells + 2*n_row, 
                           cells + 2*n_row + 1,
                           cells + 2*n_row + 2)
        return(adjacencies_i)
    }
    if(rule == 3) {
        adjacencies_i <- c(cells-n_row - 1, 
                           cells-n_row, 
                           cells-n_row+1, 
                           cells-1, 
                           cells+1, 
                           cells+n_row-1, 
                           cells+n_row, 
                           cells+n_row+1,
                           cells - 2*n_row - 2, 
                           cells - 2*n_row - 1, 
                           cells - 2*n_row, 
                           cells - 2*n_row + 1,
                           cells - 2*n_row + 2, 
                           cells - n_row - 2, 
                           cells - n_row + 2, 
                           cells - 2, 
                           cells + 2, 
                           cells + n_row - 2, 
                           cells + n_row + 2,
                           cells + 2*n_row - 2, 
                           cells + 2*n_row - 1, 
                           cells + 2*n_row, 
                           cells + 2*n_row + 1,
                           cells + 2*n_row + 2,
                           cells - 3*n_row - 3, 
                           cells - 3*n_row - 2, 
                           cells - 3*n_row - 1, 
                           cells - 3*n_row, 
                           cells - 3*n_row + 1,
                           cells - 3*n_row + 2, 
                           cells - 3*n_row + 3, 
                           cells - 2*n_row - 3,
                           cells - 2*n_row + 3, 
                           cells - n_row - 3, 
                           cells - n_row + 3, 
                           cells - 3, 
                           cells + 3, 
                           cells + n_row - 3, 
                           cells + n_row + 3,
                           cells + 2*n_row - 3, 
                           cells + 2*n_row + 3,
                           cells + 3*n_row - 3, 
                           cells + 3*n_row - 2, 
                           cells + 3*n_row - 1, 
                           cells + 3*n_row, 
                           cells + 3*n_row + 1,
                           cells + 3*n_row + 2, 
                           cells + 3*n_row + 3)
        return(adjacencies_i)
    }
}

# get matrix margins
get_margins <- function(matrix) {
    dims <- dim(matrix)
    bottom_right <- prod(dims)
    top_right <- (bottom_right - dims[1])
    c(1:dims[1], # first column 
      top_right:bottom_right, # last column
      seq(1, top_right, dims[1]), # top row
      seq(dims[1], bottom_right, dims[1])) # bottom row
}

## DEFUNCT ----
# for each cell, get all adjacent cells
# get_adjacent <- function(start_cells, n_row, n_col, rule=1) {
#     # if (rule == 1) {
#     adjacencies_i <- c(start_cells-n_row - 1, 
#                        start_cells-n_row, 
#                        start_cells-n_row+1, 
#                        start_cells-1, 
#                        start_cells+1, 
#                        start_cells+n_row-1, 
#                        start_cells+n_row, 
#                        start_cells+n_row+1)
#     # return(adjacencies_i)
#     # }
#     if(rule %in% c(2,3)) {
#         adjacencies_i <- c(adjacencies_i, 
#                            c(start_cells - 2*n_row - 2, 
#                              start_cells - 2*n_row - 1, 
#                              start_cells - 2*n_row, 
#                              start_cells - 2*n_row + 1,
#                              start_cells - 2*n_row + 2, 
#                              start_cells - n_row - 2, 
#                              start_cells - n_row + 2, 
#                              start_cells - 2, 
#                              start_cells + 2, 
#                              start_cells + n_row - 2, 
#                              start_cells + n_row + 2,
#                              start_cells + 2*n_row - 2, 
#                              start_cells + 2*n_row - 1, 
#                              start_cells + 2*n_row, 
#                              start_cells + 2*n_row + 1,
#                              start_cells + 2*n_row + 2))
#     }
#     if(rule == 3) {
#         adjacencies_i <- c(adjacencies_i, 
#                            c(start_cells - 3*n_row - 3, 
#                              start_cells - 3*n_row - 2, 
#                              start_cells - 3*n_row - 1, 
#                              start_cells - 3*n_row, 
#                              start_cells - 3*n_row + 1,
#                              start_cells - 3*n_row + 2, 
#                              start_cells - 3*n_row + 3, 
#                              start_cells - 2*n_row - 3,
#                              start_cells - 2*n_row + 3, 
#                              start_cells - n_row - 3, 
#                              start_cells - n_row + 3, 
#                              start_cells - 3, 
#                              start_cells + 3, 
#                              start_cells + n_row - 3, 
#                              start_cells + n_row + 3,
#                              start_cells + 2*n_row - 3, 
#                              start_cells + 2*n_row + 3,
#                              start_cells + 3*n_row - 3, 
#                              start_cells + 3*n_row - 2, 
#                              start_cells + 3*n_row - 1, 
#                              start_cells + 3*n_row, 
#                              start_cells + 3*n_row + 1,
#                              start_cells + 3*n_row + 2, 
#                              start_cells + 3*n_row + 3))
#     }
#     adjacencies_i
# }
# 
# 
# sim_new_range <- function(n_iter, ele_vec, lwr_vec, upr_vec, permitted_vec, 
#                           range_init, n_row, n_col, rule=1) {
#     # catch ids of range[t]
#     range_catch <- as.list(rep(NA, n_iter))
#     # initialise range 
#     range_t <- range_init
#     
#     for(t in 1:n_iter) {
#         # in elevational range and in permitted cell (i.e. habitat)
#         permitted_in_range <- as.numeric(ele_vec > lwr_vec[t] & 
#                                              ele_vec <= upr_vec[t+1] & 
#                                              permitted_vec ==1) 
#         # get range t1 (prior to removing lost range)
#         range_t1 <- spread(permitted_in_range, range_t, n_row, n_col, rule = 1)
#         # remove cells that have dropped out
#         to_drop <- which(df_gradient$xy < lwr_vec[t+1])
#         range_t1 <- range_t1[!(range_t1 %in% to_drop)]
#         # store range id and update for next iteration
#         range_catch[[t]] <- range_t1
#         range_t <- range_t1
#     }
#     return(range_catch)
# }
#
#
# spread function
# spread <- function(permitted_vec, start_cells, n_row, n_col, rule) {
#     # initial conditions
#     old_cells <- start_cells
#     adjacent_permitted <- 1
#     
#     while(length(adjacent_permitted) != 0) {
#         # for each cell, get all adjacent cells
#         adjacencies_i <- c(start_cells-n_row - 1, 
#                            start_cells-n_row, 
#                            start_cells-n_row+1, 
#                            start_cells-1, 
#                            start_cells+1, 
#                            start_cells+n_row-1, 
#                            start_cells+n_row, 
#                            start_cells+n_row+1)
#         if(rule %in% c(2,3)) {
#             adjacencies_i <- c(adjacencies_i, 
#                                c(start_cells - 2*n_row - 2, 
#                                  start_cells - 2*n_row - 1, 
#                                  start_cells - 2*n_row, 
#                                  start_cells - 2*n_row + 1,
#                                  start_cells - 2*n_row + 2, 
#                                  start_cells - n_row - 2, 
#                                  start_cells - n_row + 2, 
#                                  start_cells - 2, 
#                                  start_cells + 2, 
#                                  start_cells + n_row - 2, 
#                                  start_cells + n_row + 2,
#                                  start_cells + 2*n_row - 2, 
#                                  start_cells + 2*n_row - 1, 
#                                  start_cells + 2*n_row, 
#                                  start_cells + 2*n_row + 1,
#                                  start_cells + 2*n_row + 2))
#         }
#         if(rule == 3) {
#             adjacencies_i <- c(adjacencies_i, 
#                                c(start_cells - 3*n_row - 3, 
#                                  start_cells - 3*n_row - 2, 
#                                  start_cells - 3*n_row - 1, 
#                                  start_cells - 3*n_row, 
#                                  start_cells - 3*n_row + 1,
#                                  start_cells - 3*n_row + 2, 
#                                  start_cells - 3*n_row + 3, 
#                                  start_cells - 2*n_row - 3,
#                                  start_cells - 2*n_row + 3, 
#                                  start_cells - n_row - 3, 
#                                  start_cells - n_row + 3, 
#                                  start_cells - 3, 
#                                  start_cells + 3, 
#                                  start_cells + n_row - 3, 
#                                  start_cells + n_row + 3,
#                                  start_cells + 2*n_row - 3, 
#                                  start_cells + 2*n_row + 3,
#                                  start_cells + 3*n_row - 3, 
#                                  start_cells + 3*n_row - 2, 
#                                  start_cells + 3*n_row - 1, 
#                                  start_cells + 3*n_row, 
#                                  start_cells + 3*n_row + 1,
#                                  start_cells + 3*n_row + 2, 
#                                  start_cells + 3*n_row + 3))
#         }
#         # drop adjacencies that are start cells or outside the matrix
#         adjacency_unique <- unique(adjacencies_i)
#         adjacency_clean <- adjacency_unique[!(adjacency_unique %in% old_cells) & 
#                                                 !(adjacency_unique < 1) & 
#                                                 !(adjacency_unique > n_row*n_col)]
#         # drop non-start-cell adjacencies that aren't permitted
#         adjacent_permitted <- adjacency_clean[permitted_vec[adjacency_clean] == 1]
#         # update start_cells with new starts
#         start_cells <- adjacent_permitted
#         # update 'seen' cells to exclude in future
#         old_cells <- c(old_cells, adjacent_permitted)
#     }
#     return(old_cells)
# }


