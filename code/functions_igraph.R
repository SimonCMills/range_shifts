# functions for contagious spread through connected habitat & contagious spread
# while tracking niche
# note: functions don't depend on each other (have retained the former in case 
# it ends up being useful, but probably won't need it now)

sim_new_range_igraph_v3 <- function(n_iter, ele_vec, lwr_vec, upr_vec, permitted_cells, 
                                    start_cells, n_row, n_col) {
    ## part 1: get full graph of all cells that are connected to starting range
    ## by habitat- any cell that is not connected by habitat is never accessible
    # to add: functionality that allows gap-crossing 
    # - just a question of expanding the initial dat to include dist=2 and 3
    dat <- data.frame(x1 = rep(permitted_cells, 8), 
                      x2 = c(permitted_cells-n_row - 1, 
                             permitted_cells-n_row, 
                             permitted_cells-n_row+1, 
                             permitted_cells-1, 
                             permitted_cells+1, 
                             permitted_cells+n_row-1, 
                             permitted_cells+n_row, 
                             permitted_cells+n_row+1)) %>%
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
        # recalculate isolated cells
        linked_t <- members_t$id_cell[members_t$membership %in% 
                                          members_t$membership[members_t$id_cell %in% start_t]]
        # update range and store
        start_t <- linked_t[ele_vec[linked_t] > lwr_vec[t+1]]
        catch[[t]] <- start_t
    }
    return(catch)
} 

# contagious spread
spread_igraph2 <- function(permitted_vec, start_cells) {
    # permitted_cells <- which(permitted_vec == 1)
    dat <- data.frame(x1 = rep(permitted_cells, 8), 
                      x2 = c(permitted_cells-n_row - 1, 
                             permitted_cells-n_row, 
                             permitted_cells-n_row+1, 
                             permitted_cells-1, 
                             permitted_cells+1, 
                             permitted_cells+n_row-1, 
                             permitted_cells+n_row, 
                             permitted_cells+n_row+1)) %>%
        filter(x2 %in% permitted_cells)
    graph <- graph_from_data_frame(dat)
    graph.data.frame(dat)
    members <- components(graph)$membership %>%
        data.frame(membership = ., id_cell = as.integer(names(.)))
    members$id_cell[members$membership %in% members$membership[members$id_cell %in% start_cells]]
}

