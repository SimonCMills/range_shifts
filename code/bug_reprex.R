# Reprex for bizarre bug in igraph function that will generate duplicate vertices 
# for with e-notation (e.g. 1e+05). In all cases observed so far, the e-notation 
# version only has incoming edges, while the x00000 version only has outgoing edges. 
# Two workarounds (now using the first): 
# (1) options(scipen=99) prevents R from using e-notation and apparently fixes 
# the issue
# (2) identify failure cases, reroute the incoming edges to the x00000 version 
# (creating a correctly connected vertex in the graph), and then remove the e+05 version.
# 
# It only affects a tiny number of cases (e.g. 4 duplicated nodes in a network of 
# ~2.7e06 vertices), so is going to have an ignorable impact on results, but fixing 
# regardless. 
# 
# The check to catch further cases that do not follow the observed rule is to
# compare the number of vertices listed in dat2 (from and to columns) to the
# length of the names in the graph

# code example that mirrors the pipeline that discovered this issue
library(data.table); library(igraph)

# function to get adjacent cells in a matrix
get_adjacent <- function(cells, n_row, n_col) {
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

# function to get the margins of a matrix (i.e. 1-deep outer margin of cells)    
get_margins <- function(matrix) {
    dims <- dim(matrix)
    bottom_right <- prod(dims)
    top_right <- (bottom_right - dims[1])
    c(1:dims[1], # first column 
      top_right:bottom_right, # last column
      seq(1, top_right, dims[1]), # top row
      seq(dims[1], bottom_right, dims[1])) # bottom row
}


# (1) Before creating the failure case, produce a much smaller graph that 
# has the correct behaviour

# generate a matrix of 1-valued cells
test_mat <- matrix(1, ncol=100, nrow=100)

# remove edge cells to prevent the adjacencies wrapping around the edges
test_mat[get_margins(test_mat)] <- 0

# plot: all black cells are those that should be represented in the graph, and 
# each of these cells should each be linked to their immediately adjacent neighbours 
# (including diagonals - see get_adjacent function)
image(test_mat, asp=1, col=c("red", "black"))

# calculate the adjacency dataframe to calculate a graph from
permitted_cells <- which(test_mat[] == 1)
n_row <- dim(test_mat)[1]
n_col <- dim(test_mat)[2]

# full set of adjacencies
adj <- data.table(from = rep(permitted_cells, (1*2 + 1)^2 - 1), 
                  to = get_adjacent(permitted_cells, n_row, n_col))
# remove those that are 0-valued
adj_permitted <- adj[to %in% permitted_cells,] 

# calculate graph
g <- graph_from_data_frame(adj_permitted[,list(from, to)], directed = T)

# get vertex names 
vertex_names <- names(V(g))
graph_vertices <- data.table(name_vertex = vertex_names, 
                             id_cell = as.integer(vertex_names), 
                             id_vertex = 1:length(vertex_names)) 
setorder(graph_vertices, id_cell)

# looks good: same number of vertices in graph as there are 1-valued cells in the 
# original matrix
print(paste0("n_vertices: ", nrow(graph_vertices)))
print(paste0("n_cells: ", sum(test_mat)))


## (2) failure case. Code is identical to the above, save for the dimensions of 
## the matrix being much larger (1000 rather than 100), and the image() function 
## is commented out. 

# generate a matrix of 1-valued cells
test_mat <- matrix(1, ncol=1200, nrow=1200)

# remove edge cells to prevent the adjacencies wrapping around the edges
test_mat[get_margins(test_mat)] <- 0

# plot: all black cells are those that should be represented in the graph, and 
# each of these cells should each be linked to their immediately adjacent neighbours 
# (including diagonals - see get_adjacent function)
# image(test_mat, asp=1, col=c("red", "black"))

# calculate the adjacency dataframe to calculate a graph from
permitted_cells <- which(test_mat[] == 1)
n_row <- dim(test_mat)[1]
n_col <- dim(test_mat)[2]

# full set of adjacencies
adj <- data.table(from = rep(permitted_cells, (1*2 + 1)^2 - 1), 
                  to = get_adjacent(permitted_cells, n_row, n_col))
# remove those that are 0-valued
adj_permitted <- adj[to %in% permitted_cells,] 

# calculate graph
g <- graph_from_data_frame(adj_permitted[,list(from, to)], directed = T)

# get vertex names 
vertex_names <- names(V(g))
graph_vertices <- data.table(name_vertex = vertex_names, 
                             id_cell = as.integer(vertex_names), 
                             id_vertex = 1:length(vertex_names)) 
setorder(graph_vertices, id_cell)

# looks good: same number of vertices in graph as there are 1-valued cells in the 
# original matrix
print(paste0("n_vertices: ", nrow(graph_vertices)))
print(paste0("n_cells: ", sum(test_mat)))

print(paste0("n_extra_vertices: ", nrow(graph_vertices) - sum(test_mat)))

# (3) What are these extra vertices?
# get duplicated vertices
duplicated_vertices <- 
    graph_vertices[id_cell %in% graph_vertices[duplicated(id_cell),id_cell]]
setorder(duplicated_vertices, id_cell, id_vertex)

# the 7 additional vertices arise through duplication 
nrow(duplicated_vertices)
print(duplicated_vertices)

# xe+.. version has the incoming edges
incoming <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="in")
incoming[unlist(lapply(incoming, function(x) length(x) != 0))]
# x0.. version has outgoing edges
outgoing <- adjacent_vertices(g, duplicated_vertices$id_vertex, mode="out")
outgoing[unlist(lapply(outgoing, function(x) length(x) != 0))]
