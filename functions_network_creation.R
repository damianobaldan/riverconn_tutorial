#' Function to identify river network outlet
#'
#' @param river_network a polyline shapefile containing the river network
#' (must be of class sf)
#' @param catchment_mask a polygon shapefile with the extent of the catchment
#' (must be of class sf)
#'
#' @return the index of the river network that is corresponding to the outlet
#' (index points to the row in the river_network_shapefile)
#' @export
#'
#' @examples
#'
get_outlet <- function(river_network, catchment_mask, distance) {
  catchment_mask %>%
    st_union() %>%
    st_cast("LINESTRING") %>%
    st_is_within_distance(river_network, .,
                          dist = distance, sparse = FALSE) %>%
    match(TRUE, .)
}

#' Function to create an object of class igraph with the river network attributes and
#' the right directionality
#'
#' @param network_links a point shapefile containing the river joins and dams
#' (must be of class sf)
#' @param river_net_simplified a polyline shapefile containing the river
#' network already sliced (must be of class sf). Must contain an ID column named 'NodeID'.
#' @param outlet integer indicating the reach in river_net_simplified that is
#' the outlet of the catchment
#'
#' @return 'igraph' object of the river network
#' @export
#'
#' @examples
#'
create_network <- function(network_links, river_net_simplified, outlet) {
  
  # - use spatial join to get information about the links
  # - create a distance matrics in long format
  # - create temporaty network: it includes triangular cliques corresponding
  # to the joints
  # - loop over the triangles to remove the edges that are not good for having a
  # nicely directed network
  # - create a new undirected graph without loops
  # - define the outlet
  # - create a directed graph by setting directions based on the distance to the outlet
  # Check where river_net_simplified is overlapping with network links,
  # retain the information and save to network_links
  
  network_links_df <- st_is_within_distance(network_links, river_net_simplified,
                                            dist = 0.01) %>%
    lapply(.,
           FUN = function(x){data.frame("from" = x[1], "to" = x[2], "to2" = x[3]) }) %>%
    do.call(rbind,.) %>%
    cbind(network_links %>% st_drop_geometry())
  
  # Get a full distance matrix
  # - note: directions are messed up!
  # - note: for the joints I am creating triangles (i.e. there are 3 nodes and 3 links):
  # this is to be corrected afterwards when I decide the directionality
  full_net_links_df <- rbind(
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to2) %>%
      rename(from = from, to = to),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to) %>%
      rename(from = from, to = to2),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-from) %>%
      rename(from = to, to = to2) ,
    network_links_df %>%
      filter(is.na(to2)) %>%
      dplyr::select(-to2)
  ) %>%
    #dplyr::select(from, to, type, id_dam, pass_u, pass_d) %>%
    mutate(id_links = 1:nrow(.))
  
  # Create vertices vector
  vertices <- river_net_simplified %>% 
    st_drop_geometry %>%
    mutate(name = NodeID) %>%
    mutate(length = as.numeric(length))
  
  # Create a graph based on this distance matrix and adjust directions
  # note: add vertex name so it's easier to identify them in the next steps
  river_temp_2 <- igraph::graph_from_data_frame(
    d = full_net_links_df,
    v = vertices,
    directed = FALSE) %>%
    igraph::simplify(remove.loops = FALSE, remove.multiple = TRUE, edge.attr.comb="first")
  
  # Remove triangles from graph: this way the joints are keeping the right directionality
  cliq3 <- cliques(river_temp_2, 3, 3)
  links_to_remove <- list()
  for (i in 1:length(cliq3)){
    full_dist <- distances(river_temp_2, v =cliq3[[i]], to = outlet)
    links_to_remove[[i]] <- data.frame(
      "from" = rownames(full_dist)[full_dist != min(full_dist)],
      "to" = rownames(full_dist)[full_dist != min(full_dist)] %>% rev)
  }
  links_to_remove <- do.call(rbind, links_to_remove)
  
  # Create a graph with the edges to remove and subtract it from the original graph
  graph_to_remove <- graph_from_data_frame(d = links_to_remove,
                                           directed = FALSE)
  river_temp <- difference(river_temp_2, graph_to_remove )
  
  # Make the graph directional
  # the algorithms checks edges distances to outlet and assigns direction
  river_temp_df <- igraph::as_data_frame(river_temp, what = "edges") %>%
    mutate(from = as.numeric(from), to = as.numeric(to))
  out <- list()
  for (iii in 1:nrow(river_temp_df)) {
    df_iter <- river_temp_df[iii, ]
    delta <- distances(river_temp, v = df_iter$from, to = outlet) -
      distances(river_temp, v = df_iter$to, to = outlet)
    if (delta > 0) {from = df_iter$from; to = df_iter$to}
    if (delta < 0) {from = df_iter$to; to = df_iter$from}
    out[[iii]] <- data.frame("from" = from, "to" = to,
                             "type" = df_iter$type, 
                             "id_links" = df_iter$id_links, 
                             "id_barrier"=df_iter$id_barrier,
                             "pass_u" = df_iter$pass_u,
                             "pass_d" = df_iter$pass_d)
  }
  
  river_graph_df <- do.call(rbind, out)
  river_graph <- graph_from_data_frame(
    d = river_graph_df,
    directed = TRUE,
    v = vertices)
  
  # Return result
  return(river_graph)
}


#' Function to snap points to a dense polyline
#'
#' @param x a point shapefile that contains the points to be snapped
#' (must be of class sf)
#' @param y a polyline shapefile that contains the river network
#' (must be of class sf)
#' @param max_dist maximum distance for snapping points
#'
#' @return a data frame that mirrors x, but with the geometries snapped to y
#' (the result of st_distance(x, y) is a vector of 0s)
#' @export
#' 
#' @details this function needs a dense polyline shapefile
# (i.e. a shapefile that has small -short- lines).
# This shapefile can be built in ArcGIS using the functions to generate
# equally spaced points along a polyline and then to
# split line at points (points should be sampled with 100 - 500 m intervals )
#'
#' @examples
#' 
snap_to_river <- function(x, y, max_dist = 1000){
  
  ## NOTE: This function is needed because sf is not good at snapping
  
  # Simplify y to check distance
  y_simple <- y %>% st_union()
  # Make the river shapefile dense: extract points at each vertex
  y_points <- y %>%
    st_as_sf %>%
    st_cast("POINT") %>%
    mutate(id = 1:nrow(.))
  # Make the river shapefile dense: add splits at each vertex
  y_complex <- lwgeom::st_split(y, y_points) %>%
    st_collection_extract(.,"LINESTRING") %>%
    st_as_sf() %>%
    mutate(length = st_length(.))
  ## Check the distribution length of the simplified reaches
  #hist(y_complex$length, nclass = 300)
  # Initialize output list
  out <- list()
  # Loop over x rows
  for (i in 1:nrow(x)){
    # Select the point that will be snapped
    x_i.snap <- x[i,]
    # Extract the nearest feature
    nf <- y_complex[st_nearest_feature(x[i,], y_complex),]
    # Check the distance between the point to snap and the nearest feature
    # if yes, then snap, otherwise do not snap and keep the original point
    if ( as.numeric(st_distance(x[i,], nf)) < max_dist) {
      # # Snap the point to the nearest feature NOTE: st_snap is not working properly
      # x_i.snap <- st_snap(x[i,], nf %>% st_cast("POINT"), tolerance = max_dist)
      # Transform the nearest polyline to points
      nf.points <- st_cast(nf, "POINT")
      # Select the point that has the minimum distance to the point to snap
      nf.points.min <- nf.points[ which.min(st_distance(x[i,],nf.points)),]
      # Substitute the geometry of the point to snap with
      # the geometry of the closest point in the nearest feature (nf) object
      x_i.snap$geometry <- nf.points.min$geometry
      # # Check the distance
      # st_distance(x_i.snap, y_simple)
    }
    # Create output data frame
    out[[i]] <- x_i.snap
    #cat(paste0(i, " "))
  }
  # out_x should contain the same metadata of dams
  out_x <- do.call(rbind, out)
  return(out_x)
}


#' Check if dams are located on the headwaters of the network
#'
#' @param dams a sf point shapefile that is snapped to shape
#' @param shape a sf polyline shapefile
#'
#' @return a data.frame that contains all the attributes in the dams shapefile plus two new columns. 
#' The new columns are: 'n_intersections', reporting how many segments of the polyline are touching the dam, 
#' and dam_headwater that flags dams that are located in the headwaters (if value is TRUE). In this column, 
#' a TRUE value identifies dams that should be removed from the analysis.
#' @export
#' 
#' @description this functions identifies dams that are located in the headwaters. The workflow is as follows:
#' 1) a circular buffer is created around each dam; 2) sf::st_intersect is used to check how many sections of the 
#' river polyline intersect the buffer; 3) if the dam does not have any section of the polyline upstream (i.e. the 
#' buffer intersects the polyline only once), then it is located in the headwaters. These dams are flagged in the 
#' output so the user can remove them
#'
#' @examples
headwaters_dam <- function(dams, shape){
  
  # Create circular buffer around dams
  dams_buffer <- st_buffer(dams, 1)
  
  # Break shapefile at points
  shape_points <- shape %>%
    st_as_sf %>%
    st_cast("POINT") %>%
    mutate(id = 1:nrow(.))
  
  shape_complex <- lwgeom::st_split(shape, shape_points) %>%
    st_collection_extract(.,"LINESTRING") %>%
    st_as_sf()
  
  # Intersect with polyline
  intersections_mat <- st_intersects(dams_buffer, shape_complex)
  
  # Check if some element in intersection_mat is smaller than 2
  # This means that the buffer is crossed only one time by the river network,
  # meaning that the dam does not have an upstream segment
  
  dams_df <- dams %>% 
    st_drop_geometry() %>% 
    mutate(
      n_intersections = do.call(rbind, lapply(intersections_mat, FUN = length)),
      flag_headwater = ifelse(n_intersections == 1, TRUE, FALSE))
  
  return(dams_df)
  
}


#' Check if there are issues with confluences in hydrosheds shapefile 
#'
#' @param shape a sf polyline
#'
#' @return a sf point object with the confluences and two additional fields: n_confluences, containing the number
#' of reaches that enter and exit the confluence, and flag_confluences 
#' that flags the 'suspect' confluences with a TRUE value
#' @export
#'
#' @examples
multiple_confluences <- function(shape){
  
  # Break shapefile at points
  river_to_points <- shape %>%
    mutate(NodeID = 1:nrow(.)) %>%
    st_as_sf %>%
    st_cast("POINT") %>%
    mutate(id = 1:nrow(.))
  
  joins_selection <- river_to_points %>%
    st_equals() %>%
    Filter(function(x){length(x) > 2}, .) %>%
    lapply(., FUN = min) %>%
    unlist() %>%
    unique()
  
  river_joins <- river_to_points %>%
    dplyr::filter(id %in% joins_selection) %>%
    dplyr::select(NodeID) %>% 
    st_as_sf
  
  # Create buffer around shape_points
  river_joins_buffer <- river_joins %>%
    st_buffer(1)
  
  # Intersect with polyline
  intersections_mat <- st_intersects(river_joins_buffer, shape)
  
  # Add information to the joints shapefile
  river_joins$n_confluences <-  as.vector(do.call(rbind, lapply(intersections_mat, FUN = length)) )
  river_joins$flag_confluences <- ifelse(river_joins$n_confluences <=3, FALSE, TRUE)
  
  return(river_joins)
  
}


#' Detects gaps in the river network polyline that might cause the function network_creation to fail
#'
#' @param shape the river network shapefile to be used as input to the network_creation function
#'
#' @return a point shapefile with the position of the gaps in the river network
#'
gaps_detection <- function(shape){
  
  # Break shapefile at points
  river_to_points <- shape %>%
    st_as_sf %>%
    st_cast("POINT") %>%
    mutate(id = 1:nrow(.))
  
  # Retain all confluences, both good and 'broken' ones
  joins_selection <- river_to_points %>%
    st_is_within_distance(dist = 0.0001) %>%
    Filter(function(x){length(x) > 2 }, .) %>%
    lapply(., FUN = min) %>%
    unlist() %>%
    unique()
  
  all_river_joins <- river_to_points %>%
    filter(id %in% joins_selection)
  
  # Select casted points that are next to each confluence
  close_points <- st_is_within_distance(all_river_joins, river_to_points, dist = 1)
  
  # Loop over each points triplet and calculate the distance
  out_dist <- list()
  for (i in 1:length(close_points)){
    out_dist[[i]] <- river_to_points %>%
      filter(id %in% close_points[[i]]) %>%
      st_distance() %>%
      max()
  }
  
  # Problematic points are extracted
  problematic_points <- close_points[which(out_dist > 0)] %>%
    unlist() %>%
    unique()
  
  points_out <- river_to_points %>%
    filter(id %in% problematic_points)
  
  return(points_out)
  
}


#' Creates a river network shapefile with the component attribute that 
#' can be useful for identifying disconnected chunks
#'
#' @param network_links a point shapefile containing the river joins and dams
#' (must be of class sf)
#' @param river_net_simplified a polyline shapefile containing the river
#' network already sliced (must be of class sf).  Must contain an ID column named 'NodeID'.
#' @return sf object that copies river_net_simplified with an additional 
#' 'component' field that can be plotted to spot isolated components
#' @export
#'
#' @examples
#'
check_components <- function(network_links, river_net_simplified) {
  
  # - use spatial join to get information about the links
  # - create a distance matrics in long format
  # - create temporaty network: it includes triangular cliques corresponding
  # to the joints
  # - loop over the triangles to remove the edges that are not good for having a
  # nicely directed network
  # - create a new undirected graph without loops
  # - extract the graph components and join with original shapefile
  
  network_links_df <- st_is_within_distance(network_links, 
                                            river_net_simplified,
                                            dist = 0.01) %>%
    lapply(.,
           FUN = function(x){data.frame("from" = x[1], "to" = x[2], "to2" = x[3]) }) %>%
    do.call(rbind,.) %>%
    cbind(network_links %>% st_drop_geometry())
  
  # Get a full distance matrix
  # - note: directions are messed up!
  # - note: for the joints I am creating triangles (i.e. there are 3 nodes and 3 links):
  # this is to be corrected afterwards when I decide the directionality
  full_net_links_df <- rbind(
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to2) %>%
      rename(from = from, to = to),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to) %>%
      rename(from = from, to = to2),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-from) %>%
      rename(from = to, to = to2) ,
    network_links_df %>%
      filter(is.na(to2)) %>%
      dplyr::select(-to2)
  ) %>%
    #dplyr::select(from, to, type, id_dam, pass_u, pass_d) %>%
    mutate(id_links = 1:nrow(.))
  
  # Create vertices vector
  vertices <- river_net_simplified %>% 
    st_drop_geometry %>%
    mutate(name = 1:nrow(.)) %>%
    mutate(length = as.numeric(length))
  
  # Create a graph based on this distance matrix and adjust directions
  # note: add vertex name so it's easier to identify them in the next steps
  river_temp_2 <- igraph::graph_from_data_frame(
    d = full_net_links_df,
    v = vertices,
    directed = FALSE) %>%
    igraph::simplify(remove.loops = FALSE, remove.multiple = TRUE, edge.attr.comb="first")
  
  # Get the components of river_temp
  river_net_simplified_comp <- river_net_simplified %>%
    mutate(NodeID = 1:nrow(.)) %>%
    dplyr::select(NodeID)
  river_net_simplified_comp$component <- as.vector(components(river_temp_2)$membership)
  
  # Return output
  return(river_net_simplified_comp)
  
}


#' Creates a river network shapefile with the n_triangles field that 
#' can be useful for identifying loops for manual cleaning of the network
#'
#' @param network_links a point shapefile containing the river joins and dams
#' (must be of class sf)
#' @param river_net_simplified a polyline shapefile containing the river
#' network already sliced (must be of class sf).  Must contain an ID column named 'NodeID'.
#' @return sf object that mirrors river_net_simplified with an additional 
#' 'n_triangles' field that can be plotted to spot isolated components
#' @export
#'
#' @examples
#'
loops_identification <- function(network_links, river_net_simplified) {
  
  # - use spatial join to get information about the links
  # - create a distance matrics in long format
  # - create temporaty network: it includes triangular cliques corresponding
  # to the joints
  # - loop over the triangles to remove the edges that are not good for having a
  # nicely directed network
  # - create a new undirected graph without loops
  # - extract the graph components and join with original shapefile
  
  network_links_df <- st_is_within_distance(network_links, river_net_simplified, dist = 0.1) %>%
    lapply(.,
           FUN = function(x){data.frame("from" = x[1], "to" = x[2], "to2" = x[3]) }) %>%
    do.call(rbind,.) %>%
    cbind(network_links %>% st_drop_geometry())
  
  # Get a full distance matrix
  # - note: directions are messed up!
  # - note: for the joints I am creating triangles (i.e. there are 3 nodes and 3 links):
  # this is to be corrected afterwards when I decide the directionality
  full_net_links_df <- rbind(
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to2) %>%
      rename(from = from, to = to),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-to) %>%
      rename(from = from, to = to2),
    network_links_df %>%
      filter(!is.na(to2)) %>%
      dplyr::select(-from) %>%
      rename(from = to, to = to2) ,
    network_links_df %>%
      filter(is.na(to2)) %>%
      dplyr::select(-to2)
  ) %>%
    #dplyr::select(from, to, type, id_dam, pass_u, pass_d) %>%
    mutate(id_links = 1:nrow(.))
  
  # Create vertices vector
  vertices <- river_net_simplified %>% 
    st_drop_geometry %>%
    mutate(name = 1:nrow(.)) %>%
    mutate(length = as.numeric(length))
  
  # Create a graph based on this distance matrix and adjust directions
  # note: add vertex name so it's easier to identify them in the next steps
  river_temp_2 <- igraph::graph_from_data_frame(
    d = full_net_links_df,
    v = vertices,
    directed = FALSE) # %>% igraph::simplify(remove.loops = FALSE, remove.multiple = TRUE, edge.attr.comb="first")
  
  # Get the components of river_temp
  river_net_simplified_comp <- river_net_simplified %>%
    mutate(NodeID = 1:nrow(.)) %>%
    dplyr::select(NodeID)
  river_net_simplified_comp$n_triangles <- count_triangles(river_temp_2)
  
  # Return output
  return(river_net_simplified_comp)
  
}




