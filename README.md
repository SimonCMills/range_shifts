# Code associated with MS, "Upslope range shifts in neotropics". 
Data accessed through earth engine; derived products available at ...

## Main code files:
[functions_igraph](code/functions_igraph.R): contains the main functions for generating graphs forest connectivity, simulating upslope movement, etc.

[get_elevation_treecover](code/get_elevation_treecover.R): extracts and saves rasters of elevation and treecover for all neotropical mountain ranges, and saves the shapefile of the mountain polygons (it also calculates the proportion of each moutain range that falls in the tropics). 

[generate_range_metadata](code/generate_range_metadata.R): calculates some additional metadata for each mountain range (max&min elevation, proportion forested, etc.).

[calculate_graphs](code/calculate_graphs.R): calculates graphs for all neotropical mountain ranges.

[simulate_range_change](code/simulate_range_change.R): simulates upslope range shifts (depends on graphs)

## Ancillary code:
[example_range_shift_Santa_Marta](code/example_range_shift_Santa_Marta.R): example script that runs through the entire pipeline for a single smallish mountain range- Santa Marta (CO). 

[test_movement_functions.Rmd](code/test_movement_functions.Rmd): unit tests for graph creation and simulating movement to confirm expected behaviour. 