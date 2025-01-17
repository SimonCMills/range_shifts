# Growing set of functions for modelling spread through a network (built on igraph)
#
# Currently: 
# (1) generate_graph_undir: create undirected graph for simulating range movement
#
# (2) get_new_range: from graph, get all accessible cells and append their elevations. 
# Will be used to estimate the #cells in range through time
#
# (3) get_fm_range: from graph, get all cells at different elevations, i.e. 
# 'free-movement' version. 
#
# (4) get_adjacent: get adjacencies under varying gap-crossing rules- currently
# written manually from 0-2 cell gap-crossing; possibly needs to be able to 
# generate adjacencies for an arbitrarily large gap-crossing number. 
# 
# (5) get_margins: creates impassable border at edge of raster (to prevent 
# adjacencies wrapping around to other side of matrix)

# prevent R from using e notation (bug in R/igraph): see bug_reprex for more info
# Have preserved original solution in DEFUNCT code at end
options(scipen=99)

# dependencies
library(dplyr); library(data.table); library(igraph)

# (1) generate the graph of all cells accessible through forest
generate_graph_undir <- function(ele_vec, permitted_cells, n_row, n_col, rule) {
  # generate adjacencies
  adj <- data.table(from = rep(permitted_cells, (rule*2 + 1)^2 - 1),
                    to = get_adjacent(permitted_cells, n_row, n_col, rule))
  # remove non-permitted adjacencies
  dat2 <- adj[to %in% permitted_cells,]

  # get elevations for from and to vertices
  dat2[,c("from_ele", "to_ele") := list(ele_vec[from], ele_vec[to])]

  # remove pathways that travel upwards or downwards by more than 200m within
  # 60m (i.e. a gradient of 75 degrees, i.e. a cliff)
  dat2 <- dat2[abs(from_ele - to_ele) < 200,]

  # remove redundant object (frees up several GB in the case of Santa Marta)
  rm(adj)

  # order to create necessary ordering for dfs: dfs algorithm will jump to
  # new roots in order of vertices: by giving edge df in order of elevation,
  # this ensures it will always start from the next highest unreached vertex
  setorder(dat2, -from_ele, -to_ele)

  # this drops cells without *any* permitted adjacencies, but need to retain for
  # for calculating occupied area (keep them in graph because it simplifies
  # process downstream).
  dropped_cells <- permitted_cells[!(permitted_cells %in% unique(dat2$from))]

  # generate graph from adjacencies
  g <- graph_from_data_frame(dat2[,list(from, to)], directed = FALSE)

  # remove redundant object (after graph is generated)
  # Not sure exactly what the total saving is due to reference in data.table, but
  # adj is [216e+06, 2] and dat2 is [4, 201e+06], all cols integer-valued in the
  # case of Cordillera Central de Ecuador (large mountain range), corresponding
  # to an object size of 4GB for dat2.
  rm(dat2)

  # add dropped vertices (those without neighbours)
  g <- g + vertices(dropped_cells)

  # generate lookup table
  # note: This won't ever matter for my use, but treating as an integer will mean
  # indexing tops out somewhere around 2e+09 (as 32-bit)
  vertex_names <- names(V(g))
  graph_vertices <- data.table(id_cell = as.integer(vertex_names[]),
                               id_vertex = 1:length(vertex_names))
  graph_vertices[, ele := ele_vec[id_cell]]
  graph_vertices[, isolated := FALSE]
  graph_vertices[id_cell %in% dropped_cells, isolated := TRUE]

  # sanity check: this final dataframe containing the cells that are disconnected
  # from the graph *and* all the cells in the graph, should be equal to the total
  # number of forested cells provided at the outset
  if(nrow(graph_vertices) != length(permitted_cells)) stop("length mismatch")

  # return graph, data.table of vertices
  list(graph = g,
       df_vertices = graph_vertices)
}


# (2) given the graph, associated lookup table with elevations, and a start range
# calculate all accessible future cells (and append associated elevations)
get_new_range <- function(graph_out, upr_limit, range_width, increment, n_steps) {
  # copy vertices
  df_vertices <- copy(graph_out$df_vertices)
  # identify cells in range at t==0
  df_vertices[ele >= upr_limit - range_width  & ele < upr_limit, t := 0] 
  
  # initialise start cells at t == 0 (only searching upper perimeter)
  start_cells <- df_vertices[t == 0 & ele >= upr_limit - 200, id_cell]
  
  # No cell is more than 200m asl different from neighbour (by design: >200m
  # adjacencies are removed). Logically therefore, only need to search within
  # 200m of the upper edge. Either a cell is within 200m of an edge and therefore
  # will be able to move up if a pathway exists once the boundary is relaxed, or
  # it is more than 200m from the elevational limit, in which case movement is not
  # allowed anyway
  for(i in 1:20) {
    # catch case where there are no start vertices
    # first trim to available subgraph in timestep i
    # note: lower range never exceeds the upper limit of the range at t == 0,
    # which are already stored at the outset. Can therefore set the lower bound
    # as a static x["upr"] - 200
    retain_vertices <- df_vertices[ele >= upr_limit - 200 + increment*(i-1) &
                                     ele < upr_limit + increment*i, id_vertex]
    
    # catch case where there are no start vertices that can access new cells 
    # note: there might still be a start cell, but it isn't within 200 m of an 
    # upper elevational limit, and by definition therefore can't access new cells. 
    if(length(retain_vertices) != 0) {
    subgraph_i <- induced_subgraph(graph_out$graph, retain_vertices)

    # generate new lookup table
    vertex_names <- names(V(subgraph_i))
    subgraph_vertices <- data.table(id_cell = as.integer(vertex_names[]),
                                    id_vertex = 1:length(vertex_names))

    # generate start vertices and run bfs
    start_vertices <- subgraph_vertices[id_cell %in% start_cells, id_vertex]
    bfs_i <- bfs(subgraph_i, start_vertices, unreachable = FALSE)

    # get cells that are reached by algo (this includes start cells)
    reached_cells <- subgraph_vertices[id_vertex %in% as.integer(bfs_i$order), id_cell]
    } else {
      reached_cells <- 0
    }

    # update range info: all cells that are not in the start cells AND
    # are in the reached cells get the t-column updated (i.e. reached at time t)
    df_vertices[!(id_cell %in% start_cells) & id_cell %in% reached_cells, t := i]

    # update start cells for next iteration
    # note: reached_cells include both start_cells and new cells (as algorithm
    # will expand downwards, and by definition is linked to all start cells)
    start_cells <- subgraph_vertices[id_cell %in% reached_cells, id_cell]
    
    # record cells that have dropped out of range
    df_vertices[t == 0 & ele < upr_limit - range_width  + increment*i, t := -i]
  }

  # return dataframe of reached cells
  return(df_vertices[!is.na(t),])
}

# (3) range under 'free movement' (fm)
get_fm_range <- function(graph_out, upr_limit, range_width, increment, n_steps) {
  # copy vertices
  df_vertices <- copy(graph_out$df_vertices)
  # identify cells in range at t==0
  df_vertices[ele >= upr_limit - range_width  & ele < upr_limit, t := 0] 

  for(i in 1:n_steps) {
    # new cells
    df_vertices[!(id_cell %in% df_vertices[!is.na(t), id_cell]) & 
                  ele > upr_limit - range_width + increment*(i-1) &
                  ele < upr_limit + increment*i, t := i]
    
    # record cells that have dropped out of range
    df_vertices[t == 0 & ele < upr_limit - range_width  + increment*i, t := -i]
  }
  
  # return dataframe of reached cells
  return(df_vertices[!is.na(t),])
}

# (4) for each cell, get all adjacent cells
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

# (5) get matrix margins
get_margins <- function(matrix) {
  dims <- dim(matrix)
  bottom_right <- prod(dims)
  top_right <- (bottom_right - dims[1])
  c(1:dims[1], # first column 
    top_right:bottom_right, # last column
    seq(1, top_right, dims[1]), # top row
    seq(dims[1], bottom_right, dims[1])) # bottom row
}

# (6) contagious spread: not called in main function, but isolates the mechanism for 
# movement in a single time-slice 
spread_igraph <- function(permitted_vec, start_cells, n_row, n_col, rule=1) {
  dat <- data.frame(x1 = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
                    x2 = get_adjacent(permitted_cells, n_row, n_col, rule)) %>%
    filter(x2 %in% permitted_cells)
  graph <- graph_from_data_frame(dat)
  members <- components(graph)$membership %>%
    data.frame(membership = ., id_cell = as.integer(names(.)))
  members$id_cell[members$membership %in% members$membership[members$id_cell %in% start_cells]]
}

## DEFUNCT ----
# functions that are not used in MS: all relate to calculating climate 
# connectivity metrics. In comparison to published versions, this will consider
# *all* pathways by which a colder cell can be reached by only moving in a single
# direction along a gradient.

# generate the graph of all cells (travelling up- or along-slope, i.e from smaller
# values to larger)
generate_graph <- function(ele_vec, permitted_cells, n_row, n_col, rule) {
  
  # prevent R from using e notation (bug in R/igraph): see bug_reprex for more info
  options(scipen=99)
  # generate adjacencies
  adj <- data.table(from = rep(permitted_cells, (rule*2 + 1)^2 - 1),
                    to = get_adjacent(permitted_cells, n_row, n_col, rule))
  # remove non-permitted adjacencies
  dat <- adj[to %in% permitted_cells,]

  # get elevations for from and to vertices
  dat[,c("from_ele", "to_ele") := list(ele_vec[from], ele_vec[to])]

  # remove redundant object (frees up several GB in the case of Santa Marta)
  rm(adj)

  dat2 <- rbind(dat[from_ele > to_ele,],
                dat[from_ele == to_ele,])

  # remove pathways that travel upwards or downwards by more than 200m within
  # 60m (i.e. a gradient of 75 degrees, i.e. a cliff)
  dat2 <- dat2[abs(from_ele - to_ele) < 200,]

  # order to create necessary ordering for dfs: dfs algorithm will jump to
  # new roots in order of vertices: by giving edge df in order of elevation,
  # this ensures it will always start from the next highest unreached vertex
  setorder(dat2, -from_ele, -to_ele)

  # this drops cells without *any* permitted adjacencies, but need to retain for
  # for calculating occupied area (keep them in graph because it simplifies
  # process downstream).
  dropped_cells <- permitted_cells[!(permitted_cells %in% unique(dat2$from))]

  # generate graph from adjacencies
  g <- graph_from_data_frame(dat2[,list(from, to)], directed = T)

  # remove redundant object (after graph is generated)
  # Not sure exactly what the total saving is due to reference in data.table, but
  # adj is [216e+06, 2] and dat2 is [4, 201e+06], all cols integer-valued in the
  # case of Cordillera Central de Ecuador (large mountain range), corresponding
  # to an object size of 4GB for dat2.
  rm(dat2)

  # add dropped vertices (those without neighbours)
  g <- g + vertices(dropped_cells)

  # generate lookup table
  # note: This won't ever matter for my use, but treating as an integer will mean
  # indexing tops out somewhere around 2e+09 (as 32-bit)
  vertex_names <- names(V(g))
  graph_vertices <- data.table(id_cell = as.integer(vertex_names[]),
                               id_vertex = 1:length(vertex_names))
  graph_vertices[, ele := ele_vec[id_cell]]
  graph_vertices[, isolated := FALSE]
  graph_vertices[id_cell %in% dropped_cells, isolated := TRUE]

  # sanity check: this final dataframe containing the cells that are disconnected
  # from the graph *and* all the cells in the graph, should be equal to the total
  # number of forested cells provided at the outset
  if(nrow(graph_vertices) != length(permitted_cells)) stop("length mismatch")

  # return graph, data.table of vertices, and the two lists of ingoing and
  # outgoing edges from duplicated cells (errors should be caught by the above
  # check: but can visually check these lists after the event to confirm there
  # aren't additional failure cases)
  list(graph = g,
       df_vertices = graph_vertices)
}

# get maximum elevation accessible from each cell
get_max_ele <- function(graph_output) {
  graph_vertices <- graph_output$df_vertices
  # run dfs
  # housekeeping: fix igraph behaviour to correctly store father vertices
  igraph_options(add.vertex.names=F)
  dfs_out <- dfs(graph_output$graph, root = 1, unreachable = T, neimode="out", dist=T, father=T)
  
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


# old version of generate_graph that contains workaround for igraph/R e-notation 
# issue- identify duplicated vertices, reroute, and then remove dupes. The simpler 
# scipen solution is much better.
#
# generate_graph <- function(ele_vec, permitted_cells, n_row, n_col, rule) {
#   
#   # generate adjacencies
#   adj <- data.table(from = rep(permitted_cells, (rule*2 + 1)^2 - 1), 
#                     to = get_adjacent(permitted_cells, n_row, n_col, rule))
#   # remove non-permitted adjacencies
#   dat <- adj[to %in% permitted_cells,] 
#   
#   # get elevations for from and to vertices
#   dat[,c("from_ele", "to_ele") := list(ele_vec[from], ele_vec[to])]
#   
#   # this drops cells without *any* permitted adjacencies, but need to retain for
#   # for calculating occupied area (keep them in graph because it simplifies 
#   # process downstream). 
#   dropped_cells <- permitted_cells[!(permitted_cells %in% unique(dat$from))]
#   
#   dat2 <- rbind(dat[from_ele > to_ele,],
#                 dat[from_ele == to_ele,])
#   
#   # order to create necessary ordering for dfs: dfs algorithm will jump to 
#   # new roots in order of vertices: by giving edge df in order of elevation, 
#   # this ensures it will always start from the next highest unreached vertex
#   setorder(dat2, -from_ele, -to_ele)
#   
#   # generate graph from adjacencies
#   g <- graph_from_data_frame(dat2[,list(from, to)], directed = T)
#   g <- g + vertices(dropped_cells)
#   
#   # generate lookup table
#   # note: This won't ever matter for my use, but treating as an integer will mean 
#   # indexing tops out somewhere around 2e+09 (as 32-bit)
#   vertex_names <- names(V(g))
#   graph_vertices <- data.table(id_cell = as.integer(vertex_names[]), 
#                                id_vertex = 1:length(vertex_names)) 
#   graph_vertices[, ele := ele_vec[id_cell]]
#   graph_vertices[, isolated := FALSE]
#   graph_vertices[id_cell %in% dropped_cells, isolated := TRUE]
#   
#   ##############################################################################
#   # workaround for igraph bug: get all incoming and outgoing edges from vertices
#   # that have been duplicated (incoming always associated with xe+ version and
#   # outgoing always associated with x0 version). Identify the vertices associated
#   # with incoming edges, create incoming edges from these vertices to the x0 
#   # version of the duplicated nodes, and then delete the xe+ version from the 
#   # graph. See more extensive discussion of the issue in unit testing markdown file. 
#   
#   # logical check: in principle a small mountain range might not have this issue
#   igraph_duplication <- FALSE
#   
#   if(length(vertex_names) != length(permitted_cells)) {
#     igraph_duplication <- TRUE
#     # get duplicated vertices
#     duplicated_vertices <- 
#       graph_vertices[id_cell %in% graph_vertices[duplicated(id_cell),id_cell]]
#     setorder(duplicated_vertices, id_cell, id_vertex)
#     
#     # x0.. version has the outgoing vertices
#     incoming <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="in")
#     # xe+.. version has incoming vertices
#     outgoing <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="out")
#     
#     # retain only the x0 version: get the vertices reached from the xe+ version, 
#     # add these edges to the x0 version, and then delete the xe+ vertices from the 
#     # graph
#     in_xe_ver <- incoming[lapply(incoming, length) != 0]
#     
#     # create vector of edge starts (as.integer coerces e+05 form to 00 form)
#     in_vert_from <- rep(names(in_xe_ver), unlist(lapply(in_xe_ver, length))) %>%
#       as.integer
#     
#     # vector of edge ends
#     in_vert_to <- lapply(in_xe_ver, function(x) as.integer(names(x))) %>% unlist %>%
#       as.integer
#     
#     # interleave these two vectors
#     edges_to_add <- as.vector(rbind(in_vert_from, in_vert_to))
#     
#     # add new edges to graph
#     g <- add_edges(g, graph_vertices[match(edges_to_add, id_cell), id_vertex])
#     
#     # remove xe+ vertices from graph
#     g <- delete_vertices(g, 
#                          duplicated_vertices$id_vertex[grepl("\\+", names(outgoing))])
#     
#     ## rerun lookup table code block
#     # generate lookup table
#     # note: This won't ever matter for my use, but treating as an integer will mean 
#     # indexing tops out somewhere around 2e+09 (as 32-bit)
#     vertex_names <- names(V(g))
#     graph_vertices <- data.table(id_cell = as.integer(vertex_names[]), 
#                                  id_vertex = 1:length(vertex_names)) 
#     graph_vertices[,ele := ele_vec[id_cell]]
#     graph_vertices[, ele := ele_vec[id_cell]]
#     graph_vertices[, isolated := FALSE]
#     graph_vertices[id_cell %in% dropped_cells, isolated := TRUE]
#   }
#   ##############################################################################
#   
#   # # bind the dataframe of all the cells in the graph to the dataframe of the 
#   # # disconnected cells not in the graph
#   # for_return <- rbind(graph_vertices[!duplicated(id_cell),], disc_cells)
#   
#   # sanity check: this final dataframe containing the cells that are disconnected 
#   # from the graph *and* all the cells in the graph, should be equal to the total
#   # number of forested cells provided at the outset
#   if(nrow(graph_vertices) != length(permitted_cells)) stop("length mismatch")
#   
#   # return graph, data.table of vertices, and the two lists of ingoing and
#   # outgoing edges from duplicated cells (errors should be caught by the above
#   # check: but can visually check these lists after the event to confirm there
#   # aren't additional failure cases)
#   list(graph = g, 
#        df_vertices = graph_vertices, 
#        dup_in_edges = ifelse(igraph_duplication == FALSE, NA, incoming), 
#        dup_out_edges = ifelse(igraph_duplication == FALSE, NA, outgoing))
# }
