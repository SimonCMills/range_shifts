# # copy earth engine outputs to local (note: need to run this locally before 
# # transfer to cluster)
# fnames_drive <- list.files("D:/Drive/rgee_exports/", full.names = T, pattern="60m_10kbuffer")
# fnames_drive
# file.remove(list.files("data/mountain_ranges/", full.names=T))
# file.copy(fnames_drive, "data/mountain_ranges/", overwrite = T)

library(data.table); library(igraph); library(raster)
source("code/functions_igraph.R")
plot_mat <- function(x) plot(raster(x), asp=dim(x)[1]/dim(x)[2])
i <- 16
fnames <- list.files("data/mountain_ranges/", full.names = T, pattern="60m_10kbuffer")
print(paste0("n_fnames: ", length(fnames)))
range_names <- unique(gsub(".*60m_(.*)_2020.*", "\\1", fnames))
range_i <- range_names[i]
fnames_i <- fnames[grepl(range_i, fnames)]
ele <-  raster(fnames_i[1])
tc <-  raster(fnames_i[2])

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
ele_vec <- as.integer(as.matrix(ele))
n_row <- dim(permitted_tc)[1]
n_col <- dim(permitted_tc)[2]

# remove everything below 300m (removes 4% of cells in SM)
to_remove <- which(ele_vec < 300)
permitted_tc[to_remove] <- 0
permitted_cells <- which(permitted_tc[] == 1)

# remove all non-essential objects with non-trivial memory requirements
# rm(permitted_tc, tc_classified)

# calculate graph
time_graph_2 <- system.time(
    graph_out_2 <- generate_graph_undir(ele_vec, permitted_cells, n_row, n_col, rule=1)
)

x <- c(lwr = 300, upr=1100)
# define increment of change
n_steps <- 10
increment <- 400/n_steps

new_range <- get_new_range(graph_out_2, x, increment, 10)
# remap lost range values for plotting
new_range[t < 0, t := -(11 + t)]
# get expansion into stars
temp3 <- matrix(NA, n_row, n_col);
temp3[new_range[!is.na(t), id_cell]] <- new_range[!is.na(t), t]
rtemp3 <- raster(temp3)
crs(rtemp3) <- crs(tc)
extent(rtemp3) <- extent(tc)
library(stars)
temp3_stars <- st_as_stars(rtemp3)
tc_class_stars <- st_as_stars(tc_classified)
names(tc_class_stars) <- "fill2"

downsampling <- 1
# plot (ggplot masks `:=` so only loading here)
library(ggplot2)
ggplot() + 
    geom_stars(data = tc_class_stars, downsample=downsampling)+
    scale_fill_gradientn(colours = c("white", "grey70")) +
    ggnewscale::new_scale_fill() +
    geom_stars(data = temp3_stars, downsample=downsampling) +
    scale_fill_viridis_b(option = "inferno", na.value = NA, breaks=seq(-10, 10, 1)) +
    coord_equal() +
    theme_void() +
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0)) 
# save (post-process in AI)
ggsave("figures/Upslope_range_shifts_Santa_Marta.pdf", height=200, width=200, units="mm")
