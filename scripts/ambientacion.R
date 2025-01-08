# Paso 1: Configuración inicial
if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
if (!requireNamespace("e1071", quietly = TRUE)) install.packages("e1071")
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("smoothr", quietly = TRUE)) install.packages("smoothr")
if (!requireNamespace("units", quietly = TRUE)) install.packages("units")
if (!requireNamespace("lwgeom", quietly = TRUE)) {
  install.packages("lwgeom")
}
devtools::install_github("michaeldorman/nngeo")


# Cargar librerías
library(lwgeom)
library(terra)
library(sf)
library(dplyr)
library(psych)
library(e1071)
library(smoothr)
library(units)
library(ggplot2)
library(future)
library(future.apply)
library(leaflet)
library(htmlwidgets)

# Configuración de paralelización
plan(multicore, workers = parallel::detectCores() - 1)

setwd("C:/Users/sanicosi/Desktop/geospatial/amenabar")

# Definir CRS deseado (EPSG:32720 - UTM zona 20S WGS84)
crs_metrico <- "EPSG:32720"

log_message <- function(message) {
  cat(Sys.time(), "-", message, "\n")
}

# Paso 1: Cargar polígono y transformar al CRS EPSG:32720
log_message("Cargando y procesando polígono.")
poligono <- tryCatch({
  st_read("POLYGON.shp")
}, error = function(e) {
  stop("Error al cargar el polígono: ", e$message)
})

if (is.na(st_crs(poligono))) {
  stop("El polígono no tiene CRS definido. Verifica el archivo.")
} else if (st_crs(poligono)$epsg != 32720) {
  poligono <- st_transform(poligono, crs = crs_metrico)
}
poligono <- st_simplify(poligono, dTolerance = 5)  # Simplificar geometría
print(st_crs(poligono))
# Ploteo del polígono
plot(st_geometry(poligono), main = "Polígono Original", col = "lightblue", border = "darkblue")

# Paso 2: Cargar y procesar rasters
log_message("Cargando y procesando rasters.")
raster1 <- tryCatch({
  rast("NDVI_Image_2023.tif")
}, error = function(e) {
  stop("Error al cargar raster 1: ", e$message)
})

raster2 <- tryCatch({
  rast("NDVI_Image.tif")
}, error = function(e) {
  stop("Error al cargar raster 2: ", e$message)
})

# Reproyectar rasters
raster1 <- project(raster1, crs_metrico)
plot(raster1)
raster2 <- project(raster2, crs_metrico)
plot(raster2)
raster_stack <- c(raster1, raster2)

# Recortar rasters al polígono
raster_recortado <- crop(raster_stack, vect(poligono))
raster_recortado <- mask(raster_recortado, vect(poligono))
writeRaster(raster_recortado, "raster_recortado.tif", overwrite = TRUE)
# Ploteo de los rasters recortados
plot(raster_recortado, main = "Rasters Recortados")


# Paso 3: Generar grilla ajustada
log_message("Generando grilla ajustada.")
cellsize <- 5
grid <- st_make_grid(
  poligono,
  cellsize = cellsize,
  what = "centers"
) %>%
  st_intersection(poligono)
# Ploteo de la grilla
# Ploteo de los rasters recortados con grilla y polígono superpuestos
# Ploteo de raster, grilla y polígono usando ggplot
raster_df <- as.data.frame(raster_recortado[[1]], xy = TRUE)
raster_df <- raster_df[complete.cases(raster_df), ]  # Eliminar valores NA

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = NDVI)) +
  scale_fill_viridis_c(name = "NDVI", option = "C") +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.2) +
  geom_sf(data = poligono, fill = NA, color = "red", linewidth = 1.8) +
  labs(title = "Superposición de Raster, Grilla y Polígono", x = "Easting", y = "Northing") +
  theme_minimal()


# Paso 4: Calcular NDVI promedio y normalizar
log_message("Calculando NDVI promedio y normalizando.")
ndvi_mean <- mean(raster_recortado)
mean_value <- global(ndvi_mean, "mean", na.rm = TRUE)$mean  # Extraer la media como número
ndvi_mean_normalizado <- app(ndvi_mean, fun = function(x) x / mean_value)


# Extraer valores del raster normalizado
valores_ndvi <- terra::extract(ndvi_mean_normalizado, st_coordinates(grid))

# Crear un DataFrame con coordenadas y NDVI
matriz_datos <- data.frame(
  x = st_coordinates(grid)[, 1],
  y = st_coordinates(grid)[, 2],
  ndvi = valores_ndvi
) %>% na.omit()

# Paso 5: Desnormalizar NDVI para calcular las medias originales
log_message("Desnormalizando NDVI y calculando medias originales.")
max_raster1 <- global(raster1, "max", na.rm = TRUE)$max
max_raster2 <- global(raster2, "max", na.rm = TRUE)$max

matriz_datos <- matriz_datos %>%
  mutate(
    ndvi_desnormalizado = lyr.1 * ((max_raster1 + max_raster2) / 2)  # Promedio de valores máximos originales
  )



# Paso 6: Clustering (k-means)
log_message("Ejecutando clustering.")
silhouette_analysis <- function(data, max_k = 10) {
  library(cluster)
  avg_silhouette <- numeric(max_k - 1)
  for (k in 2:max_k) {
    kmeans_model <- kmeans(data, centers = k, nstart = 25)
    sil <- silhouette(kmeans_model$cluster, dist(data))
    avg_silhouette[k - 1] <- mean(sil[, 3])
  }
  return(avg_silhouette)
}

determine_optimal_k <- function(data, max_k = 10) {
  silhouette <- silhouette_analysis(data, max_k)
  silhouette[1] <- -Inf  # Ignorar K=2
  return(which.max(silhouette) + 1)
}

optimal_k <- determine_optimal_k(matriz_datos[, "lyr.1", drop = FALSE], max_k = 5)
clustering <- kmeans(matriz_datos[, "lyr.1", drop = FALSE], centers = 3, nstart = 25)
matriz_datos$cluster <- as.factor(clustering$cluster)

# Calcular la media desnormalizada por cluster
media_ndvi <- matriz_datos %>%
  group_by(cluster) %>%
  summarise(media_ndvi_original = mean(ndvi_desnormalizado, na.rm = TRUE))

print(media_ndvi)

# Paso 7: Crear raster de clustering
log_message("Creando raster de clustering.")
matriz_datos <- matriz_datos %>%
  group_by(x, y) %>%
  summarise(cluster = mean(as.numeric(cluster), na.rm = TRUE)) %>%
  ungroup()

raster_result <- rast(matriz_datos[, c("x", "y", "cluster")])
crs(raster_result) <- crs_metrico
writeRaster(raster_result, "raster_cluster.tif", overwrite = TRUE)
# Ploteo del raster de clustering




# Ploteo de clustering con ggplot
# Convertir cluster a factor
matriz_datos$cluster <- as.factor(matriz_datos$cluster)

# Ploteo con ggplot
ggplot(matriz_datos, aes(x = x, y = y, fill = cluster)) +
  geom_tile() +
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "Clustering Resultante", x = "Oeste", y = "Sur") +
  theme_minimal()

# Paso 7: Crear raster de clustering
log_message("Creando raster de clustering.")
matriz_datos <- matriz_datos %>%
  group_by(x, y) %>%
  summarise(cluster = mean(as.numeric(cluster), na.rm = TRUE)) %>%
  ungroup()

raster_result <- rast(matriz_datos[, c("x", "y", "cluster")])
crs(raster_result) <- crs_metrico
writeRaster(raster_result, "raster_cluster.tif", overwrite = TRUE)

# Paso 8: Suavizar raster


# Paso 8: Suavizar raster
log_message("Suavizando raster.")
smoothness <- 3

raster_suavizado <- focal(raster_result, w = matrix(1, nrow = smoothness, ncol = smoothness), fun = modal, na.rm = TRUE)
writeRaster(raster_suavizado, "raster_suavizado.tif", overwrite = TRUE)
# Ploteo del raster suavizado
plot(raster_suavizado, main = "Raster Suavizado", col = topo.colors(5))

# Recortar el raster suavizado por el polígono original
log_message("Recortando el raster suavizado por la extensión del polígono original.")
raster_suavizado_recortado <- mask(raster_suavizado, vect(poligono))

# Guardar el raster recortado
writeRaster(raster_suavizado_recortado, "raster_suavizado_recortado.tif", overwrite = TRUE)

# Ploteo del raster suavizado y recortado
plot(raster_suavizado_recortado, main = "Raster Suavizado y Recortado", col = topo.colors(5))


# Paso 9: Vectorizar y filtrar por área
log_message("Vectorizando y filtrando por área.")
vectorizado <- as.polygons(raster_suavizado_recortado)
vectorizado_sf <- st_as_sf(vectorizado)

vectorizado_sf <- vectorizado_sf %>%
  mutate(area = st_area(.))  # Recalcula el área en m²

log_message("Convirtiendo el área de m² a hectáreas.")
vectorizado_sf <- vectorizado_sf %>%
  mutate(area_ha = as.numeric(area) / 10000)  # 1 ha = 10,000 m²

vectorizado_sf <- vectorizado_sf %>%
  mutate(cluster_original = cluster) %>%
  st_cast("POLYGON")

# Crear una copia del vectorizado_sf para modificarla
vectorizado_actualizado <- vectorizado_sf

# Iterar sobre cada cluster
for (i in unique(vectorizado_sf$cluster)) {
  log_message(paste("Procesando cluster:", i))
  
  # Extraer geometría del cluster actual
  current_cluster <- vectorizado_sf %>% filter(cluster == i)
  
  # Extraer geometrías de los otros clusters
  other_clusters <- vectorizado_sf %>% filter(cluster != i)
  
  # Calcular la diferencia geométrica
  diferencia <- st_difference(st_union(current_cluster$geometry), st_union(other_clusters$geometry))
  
  # Actualizar la geometría del cluster actual
  vectorizado_actualizado <- vectorizado_actualizado %>%
    mutate(geometry = if_else(cluster == i, diferencia, geometry))
}

# Verificar el resultado final
print(vectorizado_actualizado)

# Visualizar con ggplot
ggplot() +
  geom_sf(data = vectorizado_actualizado, aes(fill = cluster), color = "black") +
  labs(title = "Clusters sin Intersecciones", x = "Easting", y = "Northing") +
  theme_minimal()


vectorizado_sf <- vectorizado_sf %>%
  mutate(area_individual = as.numeric(st_area(geometry))/10000)  # Calcular área en m²

# Filtrar áreas mayores a 1 hectárea
log_message("Filtrando áreas mayores a 1 hectárea.")
vectorizado_sf <- vectorizado_sf %>%
  filter(area_individual > 1)

vectorizado_sf <- vectorizado_sf %>%
  mutate(geometry = fill_holes(geometry, threshold = units::set_units(2e4, "m^2")))  # Rellenar huecos menores a 1 ha



# Asegurarse de que ambas columnas 'cluster' tengan el mismo tipo
media_ndvi <- media_ndvi %>%
  mutate(cluster = as.character(cluster))

vectorizado_sf <- vectorizado_sf %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(media_ndvi, by = c("cluster" = "cluster"))

print(st_crs(vectorizado_sf))         # Confirmar CRS
print(summary(vectorizado_sf$area))  # Resumen de áreas en hectáreas
print(vectorizado_sf)                 # Visualizar datos


# Suavizar bordes de la capa vectorial
log_message("Suavizando bordes de la capa vectorial.")
vectorizado_sf <- vectorizado_sf %>%
  mutate(geometry = st_simplify(geometry, dTolerance = 5))
vector_smoothness <-1

#vectorizado_sf <- smooth(vectorizado_sf, method = "ksmooth", smoothness = vector_smoothness)
log_message("Aplicando suavizado local a cada geometría.")
#vectorizado_sf <- vectorizado_sf %>%
#  group_by(cluster) %>%
#  mutate(geometry = smooth(geometry, method = "ksmooth", smoothness = 0.5)) %>%
#  ungroup()

# Corregir posibles problemas de geometría
vectorizado_sf <- st_make_valid(vectorizado_sf)

vectorizado_sf <- vectorizado_sf %>%
  mutate(cluster_original = cluster) %>%
  st_cast("POLYGON")

vectorizado_sf <- vectorizado_sf %>%
  group_by(cluster) %>%
  summarise(
    media_ndvi_original = mean(media_ndvi_original, na.rm = TRUE),
    geometry = st_union(geometry)
  ) %>%
  ungroup()


# Recortar el vectorizado por el polígono original
log_message("Recortando el vectorizado por la extensión del polígono original.")
vectorizado_recortado <- st_intersection(vectorizado_sf, poligono)

# Verificar el resultado
log_message("Verificando el resultado final.")
print(vectorizado_recortado)

# Guardar el vectorizado recortado con huecos rellenados
st_write(vectorizado_recortado, "vectorizado_recortado.shp", append = FALSE)

# Ploteo del vectorizado recortado
ggplot() +
  geom_sf(data = vectorizado_recortado, aes(fill = media_ndvi_original), color = "black") +
  geom_sf(data = poligono, fill = NA, color = "red", size = 1.2) +  # Contorno del polígono original
  scale_fill_gradient2(low = "red", mid = "blue", high = "darkgreen",
                       midpoint = mean(vectorizado_recortado$media_ndvi_original, na.rm = TRUE),
                       name = "Valor NDVI") +
  labs(title = "Vectorizado Recortado con Huecos Rellenos", x = "Easting", y = "Northing") +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))





# Ploteo del vectorizado recortado
ggplot() +
  geom_sf(data = vectorizado_recortado, aes(fill = media_ndvi_original), color = "black") +
  geom_sf(data = poligono, fill = NA, color = "black", size = 1.8) +  # Contorno del polígono original
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen",
                       midpoint = mean(vectorizado_recortado$media_ndvi_original, na.rm = TRUE),
                       name = "Valor NDVI") +
  labs(title = "Vectorizado Recortado por el Polígono", x = "Oeste", y = "Sur") +
  theme_minimal()

# Guardar el vectorizado recortado en un archivo
st_write(vectorizado_recortado, "vectorizado_recortado.shp", append = FALSE)




log_message("Proceso completado exitosamente.")




# Convertir el vectorizado recortado a un objeto comprensible para Leaflet
vectorizado_recortado_leaflet <- st_transform(vectorizado_recortado, crs = 4326)  # Transformar a WGS 84 (EPSG:4326)



# Crear el mapa Leaflet
mapa_leaflet <- leaflet() %>%
  addProviderTiles("Esri.WorldImagery", group = "Satélite") %>%
  addPolygons(
    data = vectorizado_recortado_leaflet,
    fillColor = ~colorNumeric("RdYlGn", vectorizado_recortado_leaflet$media_ndvi_original)(media_ndvi_original),
    weight = 1,
    color = "black",
    fillOpacity = 0.9,
    label = ~paste0("NDVI: ", round(media_ndvi_original, 2))
  ) %>%
  addLayersControl(
    baseGroups = c("Satélite"),
    options = layersControlOptions(collapsed = TRUE)
  ) %>%
  addLegend(
    pal = colorNumeric("RdYlGn", vectorizado_recortado_leaflet$media_ndvi_original),
    values = vectorizado_recortado_leaflet$media_ndvi_original,
    title = "NDVI",
    position = "bottomright"
  )

# Guardar el mapa como un archivo HTML
saveWidget(mapa_leaflet, file = "mapa_leaflet.html", selfcontained = TRUE)

# Mensaje de confirmación
cat("El mapa interactivo se ha guardado como 'mapa_leaflet.html'.\n")






### todo este flujo como funcion



### opcion para usar raster si terra no funciona


procesar_datos_geoespaciales_raster <- function(
    lote_name = "victoria",
    base_path = "C:/Users/sanicosi/Desktop/geospatial",
    k = 3,
    cellsize = 5,
    base_output_dir = "output",
    crs_metrico = "+proj=utm +zone=20 +south +datum=WGS84",
    simplify_tolerance = 5,
    min_area_ha = 1,
    smooth_window_size = 3,
    color_palette = "RdYlGn",
    plot_steps = TRUE,
    output_format = "shp"
) {
  library(sf)
  library(raster)
  library(dplyr)
  library(ggplot2)
  library(leaflet)
  library(htmlwidgets)
  library(rmapshaper)
  
  log_message <- function(message) {
    cat(Sys.time(), "-", message, "\n")
  }
  
  # Crear directorio de salida específico para el lote
  output_dir <- file.path(base_output_dir, lote_name)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_message(paste("Directorio de salida:", output_dir))
  
  # Construir rutas automáticamente
  poligono_path <- file.path(base_path, lote_name, "POLYGON.shp")
  # Detectar automáticamente los archivos raster en la carpeta
  raster_files <- list.files(file.path(base_path, lote_name), pattern = "\\.tif$", full.names = TRUE)
  
  # Validar que hay al menos dos archivos raster
  if (length(raster_files) < 2) {
    stop("Se requieren al menos dos archivos raster en la carpeta del lote.")
  }
  
  # Asignar los primeros dos rasters encontrados
  raster1_path <- raster_files[1]
  raster2_path <- raster_files[2]
  
  
  # Llamar a la lógica principal de la función
  log_message(paste("Procesando lote:", lote_name))
  log_message(paste("Ruta del polígono:", poligono_path))
  log_message(paste("Ruta del raster 1:", raster1_path))
  log_message(paste("Ruta del raster 2:", raster2_path))
  
  # Paso 1: Cargar polígono y transformar CRS
  log_message("Cargando y procesando polígono.")
  poligono <- tryCatch({
    st_read(poligono_path)
  }, error = function(e) {
    stop("Error al cargar el polígono: ", e$message)
  })
  
  if (is.na(st_crs(poligono))) {
    stop("El polígono no tiene CRS definido. Verifica el archivo.")
  } else {
    poligono <- st_transform(poligono, crs_metrico)
  }
  poligono <- st_simplify(poligono, dTolerance = simplify_tolerance)
  
  if (plot_steps) {
    plot(st_geometry(poligono), main = "Polígono Original", col = "lightblue", border = "darkblue")
  }
  
  # Paso 2: Cargar y procesar rasters con raster
  log_message("Cargando y procesando rasters.")
  raster1 <- raster(raster1_path)
  raster2 <- raster(raster2_path)
  
  raster1 <- projectRaster(raster1, crs = crs_metrico)
  raster2 <- projectRaster(raster2, crs = crs_metrico)
  
  raster_stack <- stack(raster1, raster2)
  raster_recortado <- crop(raster_stack, poligono)
  raster_recortado <- mask(raster_recortado, poligono)
  writeRaster(raster_recortado, file.path(output_dir, "raster_recortado.tif"), overwrite = TRUE)
  
  if (plot_steps) {
    plotRGB(raster_recortado, main = "Rasters Recortados")
  }
  
  # Paso 3: Generar grilla ajustada
  log_message("Generando grilla ajustada.")
  grid <- st_make_grid(poligono, cellsize = cellsize, what = "centers") %>%
    st_intersection(poligono)
  
  # Paso 4: Calcular NDVI promedio
  log_message("Calculando NDVI promedio.")
  ndvi_mean <- calc(raster_recortado, mean, na.rm = TRUE)
  valores_ndvi <- raster::extract(ndvi_mean, as.data.frame(st_coordinates(grid)))
  
  matriz_datos <- data.frame(
    x = st_coordinates(grid)[, 1],
    y = st_coordinates(grid)[, 2],
    ndvi = valores_ndvi
  ) %>% na.omit()
  
  # Paso 5: Clustering
  log_message(paste("Ejecutando clustering con K =", k))
  clustering <- kmeans(matriz_datos$ndvi, centers = k, nstart = 25)
  matriz_datos$cluster_id <- as.factor(clustering$cluster)
  
  # Paso 6: Crear raster de clustering
  log_message("Creando raster de clustering.")
  raster_cluster <- rasterFromXYZ(matriz_datos[, c("x", "y", "cluster_id")])
  crs(raster_cluster) <- crs_metrico
  writeRaster(raster_cluster, file.path(output_dir, "raster_cluster.tif"), overwrite = TRUE)
  
  # Paso 7: Suavizar raster
  log_message("Suavizando raster.")
  focal_matrix <- matrix(1, nrow = smooth_window_size, ncol = smooth_window_size)
  raster_suavizado <- focal(raster_cluster, w = focal_matrix, fun = modal, na.rm = TRUE)
  writeRaster(raster_suavizado, file.path(output_dir, "raster_suavizado.tif"), overwrite = TRUE)
  
  # Paso 8: Vectorizar y procesar resultados
  log_message("Vectorizando y guardando resultados.")
  vectorizado <- rasterToPolygons(raster_suavizado, dissolve = TRUE)
  vectorizado_sf <- st_as_sf(vectorizado)
  
  vectorizado_sf <- vectorizado_sf %>%
    rename(cluster_id = layer)
  
  # Calcular área y filtrar
  vectorizado_sf <- vectorizado_sf %>%
    mutate(area_ha = as.numeric(st_area(geometry)) / 10000) %>%
    filter(area_ha > min_area_ha)
  
  log_message("Incorporando valores medios de NDVI.")
  
  # Calcular media de NDVI por cluster
  media_ndvi <- matriz_datos %>%
    group_by(cluster_id) %>%
    summarise(media_ndvi_original = mean(ndvi, na.rm = TRUE))
  
  media_ndvi <- media_ndvi %>%
    mutate(cluster_id = as.character(cluster_id))
  
  vectorizado_sf <- vectorizado_sf %>%
    mutate(cluster_id = as.character(cluster_id))
  
  # Unir medias de NDVI al vectorizado
  vectorizado_sf <- vectorizado_sf %>%
    left_join(media_ndvi)
  
  log_message("Valores medios de NDVI asignados correctamente.")
  
  labels <- switch(
    as.character(k),
    "3" = c("Bajo", "Medio", "Alto"),
    "4" = c("Muy Bajo", "Bajo", "Medio", "Alto"),
    "5" = c("Muy Bajo", "Bajo", "Medio", "Alto", "Muy Alto"),
    stop("k debe ser entre 3 y 5.")
  )
  
  vectorizado_sf <- vectorizado_sf %>%
    mutate(
      ambiente = case_when(
        k == 3 ~ case_when(
          media_ndvi_original < quantile(media_ndvi_original, 1/3) ~ "Bajo",
          media_ndvi_original < quantile(media_ndvi_original, 2/3) ~ "Medio",
          TRUE ~ "Alto"
        ),
        k == 4 ~ case_when(
          media_ndvi_original < quantile(media_ndvi_original, 1/4) ~ "Muy Bajo",
          media_ndvi_original < quantile(media_ndvi_original, 2/4) ~ "Bajo",
          media_ndvi_original < quantile(media_ndvi_original, 3/4) ~ "Medio",
          TRUE ~ "Alto"
        ),
        k == 5 ~ case_when(
          media_ndvi_original < quantile(media_ndvi_original, 1/5) ~ "Muy Bajo",
          media_ndvi_original < quantile(media_ndvi_original, 2/5) ~ "Bajo",
          media_ndvi_original < quantile(media_ndvi_original, 3/5) ~ "Medio",
          media_ndvi_original < quantile(media_ndvi_original, 4/5) ~ "Alto",
          TRUE ~ "Muy Alto"
        ),
        TRUE ~ NA_character_  # Manejo de errores si k no es 3, 4 o 5
      )
    )
  
  
  log_message("Clasificación de ambiente agregada al vectorizado.")
  
  # Suavizar bordes
  log_message("Suavizando bordes de las geometrías.")
  vectorizado_sf <- vectorizado_sf %>%
    mutate(geometry = st_simplify(geometry, dTolerance = 1))
  
  # Recortar vectorizado por el polígono original
  log_message("Recortando el vectorizado por el polígono original.")
  vectorizado_recortado <- st_intersection(vectorizado_sf, poligono)
  
  log_message("Guardando archivo vectorizado recortado.")
  output_vector_file <- file.path(output_dir, paste0(lote_name, "_vectorizado_recortado.", tolower(output_format)))
  
  if (tolower(output_format) == "shp") {
    st_write(vectorizado_sf, output_vector_file, append = FALSE)
  } else if (tolower(output_format) == "geojson") {
    st_write(vectorizado_sf, output_vector_file, driver = "GeoJSON", append = FALSE)
  } else {
    stop("Formato de salida no soportado. Use 'shp' o 'geojson'.")
  }
  
  log_message("Proceso completado exitosamente.")

  
  color_palette <- c(
    "Muy Bajo" = "purple",
    "Bajo" = "red",
    "Medio" = "yellow",
    "Alto" = "darkgreen",
    "Muy Alto" = "blue"
  )
  
  ambientes_presentes <- unique(procesamiento$vectorizado_sf$ambiente)
  colores_usados <- color_palette[ambientes_presentes]
  
  plot_ndvi <- ggplot() +
    geom_sf(data = vectorizado_sf, aes(fill = media_ndvi_original), color = "black") +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = mean(vectorizado_sf$media_ndvi_original, na.rm = TRUE), name = "NDVI") +
    labs(title = "NDVI Promedio por Ambiente") +
    theme_minimal()
  
  plot_ambientes <- ggplot() +
    geom_sf(data = procesamiento$vectorizado_sf, aes(fill = ambiente), color = "black", lwd = 0.5) +
    scale_fill_manual(values = color_palette, name = "Ambiente") +
    labs(title = "Clasificación de Ambientes", x = "Oeste", y = "Sur") +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(list(
    vectorizado_sf = vectorizado_sf,
    plot_ndvi = plot_ndvi,
    plot_ambientes = plot_ambientes
  ))
  
  log_message("Proceso completado exitosamente.")
}




procesamiento = procesar_datos_geoespaciales_raster(
  lote_name = "las_balas",  # Nombre del lote definido manualmente
  base_path = "C:/Users/sanicosi/Desktop/geospatial",
  k = 5,
  cellsize = 5,
  base_output_dir = "C:/Users/sanicosi/Desktop/geospatial/output",
  simplify_tolerance = 1,
  min_area_ha = 0.5,
  smooth_window_size = 33,
  color_palette = "YlGnBu",
  plot_steps = FALSE,
  output_format = "SHP"
)

print(procesamiento)


# Definir el directorio base y nombre del lote
base_output_dir <- "C:/Users/sanicosi/Desktop/geospatial/output"
lote_name <- "las_balas"

# Crear la ruta completa al archivo vectorizado
vector_file_path <- file.path(base_output_dir, lote_name, paste0(lote_name, "_vectorizado_recortado.shp"))

# Ruta del archivo guardado
shapefile_path <-vector_file_path#"C:/Users/sanicosi/Desktop/geospatial/output/victoria/vectorizado_recortado.shp"

# Cargar el shapefile
vectorizado_recortado <- st_read(shapefile_path)

# Ver la estructura del archivo cargado
print(vectorizado_recortado)

# Plotear el shapefile usando ggplot2
ggplot() +
  geom_sf(data = vectorizado_recortado, aes(fill = md_ndv_), color = "black", lwd = 0.5) +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", 
                       midpoint = mean(vectorizado_recortado$ md_ndv_, na.rm = TRUE), 
                       name = "NDVI") +
  labs(title = "Vectorizado Recortado", x = "Easting", y = "Northing") +
  theme_minimal() +
  theme(legend.position = "right")

# Definir los colores para cada posible ambiente
color_palette <- c(
  "Muy Bajo" = "purple",
  "Bajo" = "red",
  "Medio" = "yellow",
  "Alto" = "darkgreen",
  "Muy Alto" = "blue"
)

# Filtrar los colores según los ambientes presentes en los datos
ambientes_presentes <- unique(vectorizado_recortado$ambient)
colores_usados <- color_palette[ambientes_presentes]

# Plotear el shapefile con los ambientes clasificados
ggplot() +
  geom_sf(data = vectorizado_recortado, aes(fill = ambient), color = "black", lwd = 0.5) +
  scale_fill_manual(values = colores_usados, name = "Ambiente") +
  labs(title = "Clasificación de Ambientes", x = "Oeste", y = "Sur") +
  theme_minimal() +
  theme(legend.position = "right")



#############################################################
#
#--------------Generacion de prescripciones-----------------#
#
#
#############################################################


generar_prescripcion <- function(
    vectorizado_sf, 
    tipo_prescripcion = "semilla", 
    cultivo = "maiz", 
    fraccionar_fertilizacion = FALSE,
    n_inicial_kg_ha = 0,
    n_suelo_kg_ha = 0,  # Nitrógeno aportado por el suelo
    rendimientos = c(11, 12, 13, 14, 15), 
    plantas_por_ambiente = c(45000, 50000, 55000, 60000, 70000),
    kg_n_por_ton = 20,
    fuente_nitrogeno = "Urea",
    fertilizantes = data.frame(
      producto = c("Urea", "Nitrato de Amonio", "Nitrato de Amonio Calcáreo (CAN)", 
                   "Sulfato de Amonio", "Sulfonitrato de Amonio", "Fosfato Diamónico (DAP)"),
      equivalencia_n = c(0.46, 0.32, 0.27, 0.21, 0.26, 0.18)
    )
) {
  library(dplyr)
  library(ggplot2)
  
  # Validar tipo de prescripción
  if (!tipo_prescripcion %in% c("semilla", "fertilizante")) {
    stop("El tipo de prescripción debe ser 'semilla' o 'fertilizante'.")
  }
  
  # Validar cultivo para fraccionamiento
  if (tipo_prescripcion == "fertilizante" && cultivo != "maiz" && fraccionar_fertilizacion) {
    stop("El fraccionamiento de fertilización solo está disponible para el cultivo de maíz.")
  }
  
  # Validar fuente de nitrógeno
  if (!fuente_nitrogeno %in% fertilizantes$producto) {
    stop("La fuente de nitrógeno especificada no es válida.")
  }
  
  # Detectar la cantidad de ambientes en los datos y ordenarlos por NDVI promedio
  vectorizado_sf <- vectorizado_sf %>%
    arrange(md_ndv_)  # Ordenar por NDVI promedio
  
  # Asignar rendimientos y dosis según el orden de NDVI
  n_ambientes <- nrow(vectorizado_sf)
  if (n_ambientes < 3 || n_ambientes > 5) {
    stop("La cantidad de ambientes debe estar entre 3 y 5.")
  }
  
  rendimientos_ajustados <- rendimientos[1:n_ambientes]
  plantas_ajustadas <- plantas_por_ambiente[1:n_ambientes]
  
  vectorizado_sf <- vectorizado_sf %>%
    mutate(
      rendimiento_ton_ha = rendimientos_ajustados[row_number()],
      area_ha = as.numeric(st_area(geometry)) / 10000  # Calcular superficie en hectáreas
    )
  
  if (tipo_prescripcion == "semilla") {
    vectorizado_sf <- vectorizado_sf %>%
      mutate(
        dosis_semillas_plantas_ha = plantas_ajustadas[row_number()]
      )
    fill_column <- "dosis_semillas_plantas_ha"
    fill_label <- "Dosis Semillas (Plantas/ha)"
  } else if (tipo_prescripcion == "fertilizante") {
    equivalencia_n <- fertilizantes %>%
      filter(producto == fuente_nitrogeno) %>%
      pull(equivalencia_n)
    
    vectorizado_sf <- vectorizado_sf %>%
      mutate(
        dosis_nitrogeno_kg_ha = pmax(0, rendimiento_ton_ha * kg_n_por_ton - n_suelo_kg_ha),
        dosis_fertilizante_kg_ha = dosis_nitrogeno_kg_ha / equivalencia_n,
        total_fertilizante_kg = dosis_fertilizante_kg_ha * area_ha  # Cálculo del fertilizante total
      )
    
    # Calcular refertilización si está habilitado
    if (fraccionar_fertilizacion) {
      vectorizado_sf <- vectorizado_sf %>%
        mutate(
          dosis_refertilizacion_kg_ha = pmax(0, dosis_nitrogeno_kg_ha - n_inicial_kg_ha),
          dosis_fertilizante_refertilizacion_kg_ha = dosis_refertilizacion_kg_ha / equivalencia_n,
          total_fertilizante_refertilizacion_kg = dosis_fertilizante_refertilizacion_kg_ha * area_ha
        )
      fill_column <- "dosis_fertilizante_refertilizacion_kg_ha"
      fill_label <- paste("Refertilización (", fuente_nitrogeno, ")", sep = "")
    } else {
      fill_column <- "dosis_fertilizante_kg_ha"
      fill_label <- paste("Dosis", fuente_nitrogeno, "(kg/ha)")
    }
  }
  
  # Convertir la columna de dosis a factor para colores discretos
  vectorizado_sf[[fill_column]] <- factor(round(vectorizado_sf[[fill_column]], 2), 
                                          levels = sort(unique(round(vectorizado_sf[[fill_column]], 2))))
  
  # Generar colores únicos para cada dosis
  color_palette_discreta <- scales::hue_pal()(length(levels(vectorizado_sf[[fill_column]])))
  names(color_palette_discreta) <- levels(vectorizado_sf[[fill_column]])
  
  # Plotear el mapa con colores discretos asignados a cada dosis
  plot <- ggplot() +
    geom_sf(data = vectorizado_sf, aes_string(fill = fill_column), color = "black", lwd = 0.5) +
    scale_fill_manual(values = color_palette_discreta, name = fill_label) +
    labs(title = paste("Prescripción:", tipo_prescripcion, "-", fuente_nitrogeno), 
         x = "Oeste", y = "Sur") +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(list(vectorizado_sf = vectorizado_sf, plot = plot))
}



### Prescripcion de semillas
resultado <- generar_prescripcion(
  vectorizado_sf = vectorizado_recortado,
  tipo_prescripcion = "semilla",
  cultivo = "maiz",
  fraccionar_fertilizacion = FALSE,
  n_inicial_kg_ha = 150,
  rendimientos = c( 12, 13, 14),
  plantas_por_ambiente = c(45000,  55000, 70000),
  fuente_nitrogeno = "Urea"
)

print(resultado)


### prescripcion de ferti
resultado <- generar_prescripcion(
  vectorizado_sf = vectorizado_recortado,
  tipo_prescripcion = "fertilizante",
  cultivo = "maiz",
  fraccionar_fertilizacion = FALSE,
  n_inicial_kg_ha = 50,
  n_suelo_kg_ha = 30,  # Ejemplo: 30 kg/ha de nitrógeno aportado por el suelo
  rendimientos = c(12, 13, 14),
  plantas_por_ambiente = c(45000, 55000, 70000),
  fuente_nitrogeno = "Urea"
)

print(resultado$vectorizado_sf)  # Ver tabla con áreas y dosis ajustadas
print(resultado$plot)  # Ver el mapa de prescripción

library(DT)

datatable(resultado$vectorizado_sf, 
          extensions = 'Buttons', 
          options = list(
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
            pageLength = 5
          ),
          caption = "Prescripción de Fertilizantes por Ambiente")
