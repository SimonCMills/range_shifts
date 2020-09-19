## check functions perform as expected
# note: the full raster has to have a margin of "not-permitted" (i.e. 0-valued)
# cells, otherwise it will spread across edges. Probably best to hard code this 
# into function?

library(raster); library(igraph); library(dplyr)
source("code/functions_igraph.R")

## plot matrices
# rasters index cells horizontally rather than vertically, which will mess up code
# written for matrices, but plotting them without rotation and safely (by reducing
# resolution for large matrices is useful), so use this. 
plot_mat <- function(x, ...) plot(raster(x), ...)

## test 1 ----
# does the algorithm spread contagiously?
# read 65*50 matrix that has "permitted movement" cells spaced contiguously, 
# 1-cell removed, and 2-cell removed
test_mat <- as.matrix(raster("data/pikachu.png"))
test_mat <- as.matrix(raster::disaggregate(raster(test_mat), fact=1))
test_mat[test_mat[] == 5] <- 1
test_mat[test_mat[] > 1] <-0
permitted_movement <- test_mat

# for spread function
start_cells <- min(which(test_mat[] == 1))
n_row <- dim(permitted_movement)[1]
n_col <- dim(permitted_movement)[2]
permitted_cells <- which(permitted_movement == 1)

# run
fill1 <- spread_igraph(permitted_cells, start_cells, n_row, n_col, rule=1)
fill2 <- spread_igraph(permitted_cells, start_cells, n_row, n_col, rule=2)
fill3 <- spread_igraph(permitted_cells, start_cells, n_row, n_col, rule=3)

# fill points reached by algorithm
mat1 <- test_mat; mat1[fill1] <- 20
mat2 <- test_mat; mat2[fill2] <- 20
mat3 <- test_mat; mat3[fill3] <- 20

# save results
png("figures/test_spread_algo_1.png")
par(mfrow=c(2,2))
permitted_plot <- permitted_movement
permitted_plot[start_cells] <- 20
plot(raster(permitted_plot), asp=1, main = "Permitted movement & start")
plot(raster(mat1), asp=1, main = "Spread, 1-cell search radius")
plot(raster(mat2), asp=1, main = "Spread, 2-cell search radius")
plot(raster(mat3), asp=1, main = "Spread, 3-cell search radius")
dev.off()

## test 2----
# Does function correctly track niche through time?
# create an elevational gradient
df_gradient <- expand.grid(x = 1:n_col, y = 1:n_row) %>%
    mutate(xy = (x*y), xy_scale = (xy-min(xy))/max(xy-min(xy)) * 1000) %>%
    arrange(x, y)

# for passing to simulation
ele_vec <- df_gradient$xy_scale
start_cells <- min(which(test_mat[] == 1))
permitted_cells <- which(permitted_movement == 1)
n_iter <- 300
lwr_vec <- 0 + 0:(n_iter+1)
upr_vec <- 100 + 0:(n_iter+1)

# simulate
ranges_sim <- sim_new_range_igraph(n_iter, ele_vec, lwr_vec, upr_vec, 
                                   permitted_cells, start_cells, n_row, n_col, rule=2)

# save 
png("figures/test_rangesim_rule2.png", height=200, width=110, units="mm", res=100)
par(mfrow = c(4,2))
plot_mat(matrix(ele_vec, ncol=n_col, nrow=n_row), main="Elevation gradient")
mat_2 <- test_mat; mat_2[mat_2 == 1] <- ele_vec[mat_2==1]
plot_mat(mat_2, main="Permitted movement & gradient")
# plot results
mat1 <- test_mat; mat1[ranges_sim[[1]]] <- 20
plot_mat(mat1, main="t[1]")
mat1 <- test_mat; mat1[ranges_sim[[50]]] <- 20
plot_mat(mat1, main="t[50]")
mat1 <- test_mat; mat1[ranges_sim[[100]]] <- 20
plot_mat(mat1, main="t[100]")
mat1 <- test_mat; mat1[ranges_sim[[150]]] <- 20
plot_mat(mat1, main="t[150]")
mat1 <- test_mat; mat1[ranges_sim[[200]]] <- 20
plot_mat(mat1, main="t[200]")
mat1 <- test_mat; mat1[ranges_sim[[250]]] <- 20
plot_mat(mat1, main="t[250]")
dev.off()

# generate gif
# library(animation)
# saveHTML({
#     for(i in seq(1, 100, 5)) {
#         mat_t <- test_mat; mat_t[ranges_sim[[i]]] <- 10
#         plot_mat(mat_t)
#     }
# },title = "Simulated spread, rule=1", description = desc, verbose = FALSE)
