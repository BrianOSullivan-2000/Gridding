
################################################################################

## Make a spatial mesh for various methods - I most often    ##
## use this approach for spatially embedded Random Forests   ##

################################################################################

make_mesh <- function(island_outline, grid, a, q, 
                      plot_mesh = F, buffer = 5000, coast_limit = 5000){
  
  # island_outline        -outline of the region of interest (probably Ireland)
  #                       this should just be points with east/north coords
  # grid                  -grid for finding minimum distances
  # a                     -maximum area of each triangle in the mesh
  # q                     -minimum angled allowed in the mesh
  
  # plot_mesh             -this will plot the mesh if needed
  # buffer                -how much padding away from the border is given
  # coast_limit           -minimum distance allowed for mesh from the border
  
  ## RTriangle package makes the mesh
  ## concaveman simplifies spatial objects
  require(RTriangle); require(concaveman); require(sf); require(sp)
  
  ## Convert the island outline to a spatial object with lat/lon coordinates
  island_points <- st_as_sf(island_outline, coords = c("x", "y"), crs = 29903)
  island_points <- st_transform(island_points, crs = 4326)
  
  ## Simplify the mesh and add a buffer
  boundary <- concaveman(island_points, concavity = 10)
  boundary <- st_buffer(boundary, dist = buffer)
  boundary <- st_simplify(boundary, dTolerance = 5000)
  
  ## Get the vertices and edges of the boundary
  boundary_nodes <- st_cast(x = boundary, "POINT", crs = 4326)
  boundary_nodes <- st_coordinates(x = boundary_nodes)
  boundary_nodes <- data.frame(lon = boundary_nodes[,1], lat = boundary_nodes[,2])
  boundary_nodes <- boundary_nodes[-nrow(boundary_nodes), ]
  boundary_segments <- cbind(1:nrow(boundary_nodes), c(2:nrow(boundary_nodes), 1))
  
  ## Make the mesh
  boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)
  mesh_regular <- triangulate(p = boundary_pslg, a=a, q=q)
  
  ## Pull the nodes from the mesh
  nodes <- as.data.frame(mesh_regular$P); colnames(nodes) <- c("x", "y")
  nodes <- st_as_sf(nodes, coords = c("x", "y"), crs = 4326)
  nodes <- st_transform(nodes, crs = 29903)
  nodes <- st_coordinates(nodes)
  
  ## We'll need minimum distances to remove nodes too far from boundary
  min_dists <- apply(spDists(as.matrix(nodes), 
                             as.matrix(grid[c("east", "north")])), 1, min)
  
  ## Plotting option
  if(plot_mesh){
    require(ggplot2)
    
    mesh_regular$E <- 
      mesh_regular$E[apply(matrix(mesh_regular$E %in% 
                                    which(min_dists < coast_limit), ncol = 2), 1, all), ]
    edges_df <- data.frame(
      x1 = nodes[mesh_regular$E[, 1], 1], y1 = nodes[mesh_regular$E[, 1], 2],
      x2 = nodes[mesh_regular$E[, 2], 1], y2 = nodes[mesh_regular$E[, 2], 2]
    )
    
    plot_nodes <- as.data.frame(nodes)
    plot_nodes$min_dists <- min_dists
    
    mesh_plot <- ggplot() +
      theme_void() + geom_point(data = island_outline, aes(x, y), 
                                colour = "gray50", size = 0.5) +
      geom_segment(data = edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), 
                   color = "red", linewidth = 0.5, alpha = 0.2) +
      geom_point(data = plot_nodes[min_dists < coast_limit, ], aes(x = X, y = Y), size = 2) +
      
      theme(aspect.ratio = 1, 
            axis.title = element_blank(), 
            axis.ticks = element_blank(), axis.text = element_blank(),
            panel.background = element_rect(fill = "white"))
    print(mesh_plot)
  }
  
  ## Remove nodes too far from boundary
  nodes <- nodes[min_dists < coast_limit, ]
  
  return(nodes)
}