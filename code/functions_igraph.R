# functions for contagious spread through connected habitat & contagious spread
# while tracking niche
# notes: 
#   - niche-tracking through time doesn't depend on the contagious spread function
#   - the order is take range[t], do all movement into cells that are accessible 
#   in the t:t+1 interval, and *then* remove cells that have dropped out of the 
#   lower range edge 
plot_mat <- function(x, ...) plot(raster(x), ...)

# this function has been so *incredibly* headache inducing, but sorted now
get_max_ele <- function(ele_vec, permitted_cells, n_row, n_col, rule=2) {
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
  # order to create necessary ordering for dfs (expand on this..)
  setorder(dat2, -from_ele, -to_ele)
  
  # generate graph from adjacencies
  g <- graph_from_data_frame(dat2, directed = T)
  
  # generate lookup table
  graph_vertices <- V(g)$name %>%
    data.table(id_cell = as.integer(.), id_vertex = 1:length(.), 
               max_ele=as.double(rep(NA, length(.)))) 
  graph_vertices[,ele := ele_vec[id_cell]]
  
  # run dfs
  dfs_out <- dfs(g, root = 1, unreachable = T, neimode="out", dist=T)

  # add distances to df with original ordering
  graph_vertices[,dist := dfs_out$dist]
  # now generate ordering for doing the membership calculation
  graph_vertices[as.numeric(dfs_out$order), ordering := 1:.N]
  setorder(graph_vertices, ordering)
  # get membership and get maximum accessible elevation under upslope movement
  graph_vertices[,membership := cumsum(dist==0)]
  graph_vertices[,max_ele := max(ele), by="membership"]
  
  return(graph_vertices[])
}

sim_new_range_igraph <- function(n_iter, ele_vec, lwr_vec, upr_vec, permitted_cells, 
                                    start_cells, n_row, n_col, rule=1) {
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


