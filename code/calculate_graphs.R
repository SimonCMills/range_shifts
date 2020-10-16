# script to generate graphs & rangeshifts

# get elevations and treecover for mountain ranges
fnames <- list.files("data/mountain_ranges/", full.names = T)
# housekeeping
library(raster); library(dplyr); library(data.table); library(igraph)
source("code/functions_igraph.R")

# colour palette
# cols <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", 
#                            "#FEE090", "#FDAE61", "#F46D43", "#D73027"))

range_names <- unique(gsub(".*60m_(.*)_2020.*", "\\1", fnames))

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
    permitted_cells <- which(permitted_tc[] == 1)
    ele_vec <- as.numeric(as.matrix(ele))
    n_row <- dim(permitted_tc)[1]
    n_col <- dim(permitted_tc)[2]

    time_graph <- system.time(
        graph_out <- generate_graph(ele_vec, permitted_cells, n_row, n_col, rule=1)
    )
    
    # save graph (really slow to write but will save memory overhead in future 
    # and is fast to read)
    graph_out$time <- time_graph
    saveRDS(graph_out, paste0("outputs/graph_", range_i))
    
    
    # png(paste0("figures/permitted_and_ele_", range_i, ".png"), res=400, width=200, height=150, units="mm")
    # par(mfrow=c(1,2))
    # plot(tc_classified)
    # ele_mat_plot <- permitted_tc; ele_mat_plot[] <- NA
    # ele_mat_plot[permitted_cells] <- ele_vec[permitted_cells]
    # plot_mat(ele_mat_plot, asp=n_row/n_col)
    # dev.off()
    
    # elevational ranges to simulate movement from (ranging between lower of 300 
    # and upper of highest forested elevation in a given mountain range)
    # 322 m is a 2 degree shift with adibatic lapse rate of 6.2km^-1
    max_forested_ele <- max(ele_vec[permitted_cells])
    ele_lims <- data.table(lwr = seq(300, max_forested_ele-800-322, len = 40), 
                           upr = seq(1100, max_forested_ele-322, len = 40))
    
    # get upslope range from varying start ranges
    diffs <- apply(ele_lims, 1, 
                   function(x) {
                       in_ele_range <- which(ele_vec >= x["lwr"] & ele_vec < x["upr"])
                       start_i <- in_ele_range[in_ele_range %in% permitted_cells]
                       range_i <- get_range(graph_out, ele_vec, start_i)
                       range_i[ele >= x["lwr"] + 322 & ele < x["upr"] + 322]
                       })
    
    ######
    full <- rbindlist(diffs, idcol="id")   
    full[, lwr := ele_lims$lwr[id]]
    full[, upr := ele_lims$upr[id]]
    
    ele_lims[, n_start := sum(ele_vec[permitted_cells] >= lwr & 
                                  ele_vec[permitted_cells] < upr), by="lwr"]
    ele_lims[, n_pot := sum(ele_vec[permitted_cells] >= lwr + 322 & 
                                ele_vec[permitted_cells] < upr + 322), by="lwr"]
    
    actual_ranges <- full[ ,.(n_act = .N), by="lwr"]
    
    range_extents <- merge(ele_lims, actual_ranges, all.x=T)
    range_extents[is.na(n_act), n_act := 0]
    saveRDS(list(df_rangeshift = full, 
                 df_summary = range_extents), paste0("outputs/rangeshift_", range_i))
}
