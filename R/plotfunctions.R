#' @title Print method for "gaussdiag" class
#'
#' @description Plots score diagnostics that are the output of the \code{ng.check} function.
#' It including the contribution of each index in the latent field to the overall Bayes factor (BF)
#' and the obsedved BF sensitivity compared against the reference distribution.
#'
#' @param x List. Output of the \code{ng.score} function
#' @param n Integer. Plots the index number of the \code{n} indeces with the largest BF sensitivity.
#' @param samples Optional vector. Samples of the reference distribution.
#' @return A list of ggplots:
#' \itemize{
#'    \item \code{First:}  The contribution of each index to the overallBF sensitivity: \eqn{d_i(\mathbf{y})}.
#'    \item \code{Second:} Reference distribution of the BF sensitivity vs observed BF sensitivity.
#'  }
#' @export
plot.gaussdiag <- function(x, n = rep(0,nrow(x$BF.check)), samples = NULL, ...){
  comp.names <- rownames(x$BF.check)
  n.comps    <- length(comp.names)

  data.plot  <- matrix(NA,0,5)
  for(i in 1:n.comps){
    s.i               <- x$BF.index[[i]]
    N                 <- length(s.i)
    indeces           <- order(abs(s.i), decreasing = TRUE)[seq_len(n[i])]
    labels            <- rep("", length(s.i))
    labels[indeces]   <- as.character(indeces)
    alpha             <- rep(1,length(s.i))
    alpha[indeces]    <- 0
    data.plot         <- rbind(data.plot, data.frame(index     = 1:N,
                                                     score     = s.i,
                                                     alpha     = alpha,
                                                     labels    = labels,
                                                     component = rep(comp.names[i],N) ))
  }

  plot.ranking <- ggplot(data = data.plot,  aes(x = index, y = score, color = component, label= labels)) +
    geom_point(size=2.5, aes(alpha=alpha, color = component, fill = component), stroke = 0) +
    geom_text(show.legend = F) +
    #scale_color_brewer(palette = "Dark2")+
    #scale_fill_brewer(palette = "Dark2")+
    xlab("index") +
    ylab(TeX("d_i(\\mathbf{y})")) +
    ggtitle("Score for each index") +
    scale_alpha(guide = 'none') +
    theme_bw() +
    theme(text = element_text(size = 15))


  color.vector <- hue_pal()(n.comps)
  #  color.vector <- brewer.pal(max(3,n.comps),"Dark2")



  if(!is.null(samples)){
    sds               <- unlist(lapply(samples, function(x) sd(x) ))
  }
  else{
    if(!is.null(x$BF.check$var.ref.mode)){
      sds  <- sqrt(x$BF.check$var.ref.mode)}
    else{return(plot.ranking)}
  }

  plot.ref <- list()

  for(i in 1:n.comps){

    if(!is.null(samples)){
      plot.ref[[i]] <- ggplot(data = data.frame(samples = samples[[i]])) +
        geom_histogram(aes(x=samples, y = ..density.., colour="s(\\mathbf{y}^{pred})"), fill = "grey")+
        stat_function(fun = dnorm,
                      args = list(0, sds[i]),
                      size = 1.5, aes(colour = "s(\\mathbf{y}^{pred})")) +
        geom_vline(xintercept = x$BF.check$s0.mode[i],
                   colour = color.vector[i],
                   aes(colour = "s(\\mathbf{y})"),
                   size=1.5) +
        scale_colour_manual("",
                            values = c("s(\\mathbf{y}^{pred})" = "black", "s(\\mathbf{y})" = color.vector[i]),
                            labels = c(TeX("s(\\mathbf{y}^{pred})"),TeX("s(\\mathbf{y})")))+
        xlab("score") +
        ylab("count") +
        xlim(min(-4*sds[i],x$BF.check$s0.mode[i]*1.1), max(4*sds[i],x$BF.check$s0.mode[i]*1.1)) +
        ggtitle(paste0("Ref. and obs. score for ", comp.names[i]) ) +
        theme_bw() +
        theme(legend.position = "none") +
        theme(legend.text.align = 0, text = element_text(size = 15))

    }
    else{
      plot.ref[[i]] <- ggplot() +
        stat_function(fun = dnorm,
                      args = list(0, sds[i]),
                      size = 1.5, aes(colour = "s(\\mathbf{y}^{pred})")) +
        geom_vline(xintercept = x$BF.check$s0.mode[i],
                   colour = color.vector[i],
                   aes(colour = "s(\\mathbf{y})"),
                   size=1.5) +
        scale_colour_manual("",
                            values = c("s(\\mathbf{y}^{pred})" = "black", "s(\\mathbf{y})" = color.vector[i]),
                            labels = c(TeX("s(\\mathbf{y}^{pred})"),TeX("s(\\mathbf{y})")))+
        xlab("score") +
        ylab("density") +
        xlim(min(-4*sds[i], x$BF.check$s0.mode[i]*1.1), max(4*sds[i], x$BF.check$s0.mode[i]*1.1)) +
        ggtitle(paste0("Ref. and obs. score for ", comp.names[i]) ) +
        theme_bw()   +
        theme(legend.position = "none") +
        theme(legend.text.align = 0, text = element_text(size = 15))

    }


  }


  return(list(plot.ranking = plot.ranking, plot.ref = plot.ref))

}


#' @title Plot areal data
#'
#' @description Generates a Leaflet widget to plot areal data. Requires \code{leaflet} and
#' \code{leaflet.extra} package.
#'
#' @param SpatialPolygons Spatial Polygons object.
#' @param data Vector of data to plot.
#' @param palette Character. Palette to use in Leaflet widget.
#' @param title Character. Title for the legend.
#' @return Leaflet widget
#' @export
areal.plot <- function(map, data, palette = "YlOrRd", title = "Value"){

  leaf       <- require("leaflet")
  leaf.extra <- require("leaflet.extras")
  if( (leaf == FALSE) || (leaf.extra == FALSE) ){
    stop("This function requires the packages: leaflet, leaflet.extras")
  }


  data           <- data.frame(data = data)
  rownames(data) <- sapply(slot(map, "polygons"), function(x){slot(x, "ID")})

  #rownames(data) <- paste("ID",1:N,sep="")
  #sapply(slot(map, "polygons"), function(x){slot(x, "ID")})

  N              <- length(data)
  names          <- colnames(data)
  map            <- SpatialPolygonsDataFrame(map, data)

  pal <- colorNumeric(palette = palette, domain = data)

  labels <- sprintf("<strong>%s: </strong> %s", rownames(data), round(data$data,5)) %>% lapply(htmltools::HTML)


  l <- leaflet(map) %>%
    addPolygons(
      color = "grey", weight = 1, fillColor = ~ pal(data),
      fillOpacity = 0.5,
      highlightOptions = highlightOptions(weight = 4),
      label = labels,
      labelOptions = labelOptions(
        style =
          list(
            "font-weight" = "normal",
            padding = "3px 8px"
          ),
        textsize = "15px", direction = "auto"
      )
    ) %>%
    addLegend_decreasing(
      pal = pal, values = ~data, opacity = 0.5, title = title,
      position = "bottomright", decreasing = TRUE
    )

  return(l)
}

#' @title Plot geostat data
#'
#' @description Generates a Leaflet widget to plot geostat data. Requires \code{leaflet} and
#' \code{leaflet.extra} package.
#'
#' @param data Data frame. Should countain 3 collumns with the latitude, longitude and measurements.
#' @param mesh inla.mesh object.
#' @param n Integer. Plots the index number of the \code{n} indeces with the largest score.
#' @param palette Character. Palette to use in Leaflet widget.
#' @param domain Vector. Domain to use in palette scale
#' @param tile.provider Character. Name of the Leaflet tile provider.
#' @return Leaflet widget
#' @export
geo.plot <- function(data, mesh = NULL, n = 0,
                     palette = c("blue","white","red"), domain = data[,3],
                     tile.provider = providers$Stamen.Terrain){


  leaf       <- require("leaflet")
  leaf.extra <- require("leaflet.extras")
  if( (leaf == FALSE) || (leaf.extra == FALSE) ){
    stop("This function requires the packages: leaflet, leaflet.extras")
  }


  plot_data           <- data
  V                   <- data[,3]
  N                   <- nrow(data)
  var.name            <- colnames(data)[3]
  colnames(plot_data) <- c("lon","lat","value")

  suppressWarnings(
    if(!is.null(mesh)){
      mesh_map <- inla.mesh2sp(mesh)$triangles
      proj4string(mesh_map) <- CRS("+proj=longlat +datum=WGS84")
      mesh_map_proj <- spTransform(mesh_map, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    }
  )


  sel    <- rep(0,N)
  sel[order(V, decreasing = TRUE)[seq_len(n)]] <- 1 #do 95% quantiles include the value 1?


  #Palletes
  pal.V  <- colorNumeric(palette = palette, domain = domain)

  l <- leaflet(plot_data) %>%

    addProviderTiles(tile.provider) %>%

    addPolygons(data=mesh_map_proj, weight = 0.5, fill = FALSE, color = "#0A0708")  %>%

    addCircleMarkers(data = plot_data, lng = ~lon, lat = ~lat, color="black",
                     stroke= TRUE, fill = TRUE, fillColor = ~pal.V(value),
                     radius = 6, fillOpacity = TRUE, weight = sel*4.5 + 0.5, opacity = 1) %>%

    addLegend_decreasing(pal = pal.V, values = plot_data$value, opacity = 0.8, title = var.name,
                         position = "bottomright", decreasing = TRUE) %>%
    addFullscreenControl()

  return(l)
}

#leaflet decreasing label
addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft","topleft"),
                                  pal, values, na.label = "NA", bins = 7, colors,
                                  opacity = 0.5, labels = NULL, labFormat = labelFormat(),
                                  title = NULL, className = "info legend", layerId = NULL,
                                  group = NULL, data = getMapData(map), decreasing = FALSE) {

  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors))
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula"))
      title <- deparse(values[[2]])
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] ==
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins))
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1)
        pretty(values, bins)
      else bins
      if (length(bins) > 2)
        if (!all(abs(diff(bins, differences = 2)) <=
                 sqrt(.Machine$double.eps)))
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values)))
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels))
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)),
                 na_color = na.color, na_label = na.label, opacity = opacity,
                 position = position, type = type, title = title, extra = extra,
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}

inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
  }

  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)

  list(triangles = triangles, vertices = vertices)
}


