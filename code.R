########################################
#        Movement Patterns, Home Range, 
# and Habitat Selection of Crocodylians in 
#  Tayrona National Natural Park, Colombia.
#
#######################################

rm(list = ls())
dev.off()

# 0. Initial processing process----
## 0.1 Libraries----
library(readr)
library(ggplot2)
library(move)
library(dplyr)
library(sp)
library(rgdal)
library(rgeos)  
library(raster)
library(adehabitatHR)
library(scico)
library(scales)

#library(terra)
library(sf)
library(amt)
library(ctmm)
#library(parallel)
#library(INLA)
#library(corrplot)
library(stringr)
library(purrr)
library(broom)
library(ggthemes)
library(reshape2)
library(corrplot)
library(ggspatial)
library(patchwork)
library(cowplot)

library(lubridate)
library(tidyr)
library(mgcv)
library(gratia)
library(marginaleffects)
library(viridis) 

## 0.2 Creating folders----
loc.data <- paste0("./DATA/")
loc.fig <- paste0("./Figures/")
loc.output <- paste0("./dBBMMoutput/")
loc.output1 <- paste0("./iSSFoutput/")

## 0.3 Calling data and defining parameters
dir()

all.data <- read_csv(file = "Data.csv")
movement_data <- read.csv(file = "movements.csv")

crs.proj <- CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
comment(crs.proj) <- NULL 

data_presence <-all.data

telemetry <- read.csv("Data.csv")
telemetry1 <- read.csv("Data_wgs84.csv")
telemetry <- telemetry %>% filter(!id %in% c("CrAc77", "CrAc4"))
telemetry1 <- telemetry1 %>% filter(!id %in% c("CrAc77", "CrAc4"))
# 1.0 Generation of Kernel density and dBBMM----
## 1.0.1 Generating variograms and shapefiles of home range----
telemetry1$datetime <- as.POSIXct(telemetry1$datetime, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
telemetry1$timestamp <- telemetry1$datetime
telemetry1$individual <- telemetry1$id

individuals <- split(telemetry1, telemetry1$individual)
individuals <- lapply(individuals, function(df) {
  df <- df[, c("timestamp", "longitude", "latitude", "individual")]
  colnames(df) <- c("timestamp", "longitude", "latitude", "individual")
  return(df)
})

data_list <- lapply(individuals, function(df) {
  as.telemetry(df, projection = "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")
})

varioDataList <- lapply(names(data_list), function(id) {
  d <- data_list[[id]]
  vg <- variogram(d, fast = FALSE, CI = "Gauss")
  guess <- ctmm.guess(d, interactive = FALSE)
  fits <- ctmm.select(d, guess, verbose = TRUE, method = "pHREML")
  
  png(paste0("Figures/variogram_", id, ".png"))
  plot(vg, CTMM = fits, col.CTMM = c("red","purple","blue"), fraction = 0.65, level = 0.5)
  dev.off()
  
  varioData <- data.frame(SemiVar = vg$SVF, DOF = vg$DOF, Lag = vg$lag)
  CI.upper <- Vectorize(function(k, Alpha) qchisq(Alpha / 2, k, lower.tail = FALSE) / k)
  CI.lower <- Vectorize(function(k, Alpha) qchisq(Alpha / 2, k, lower.tail = TRUE) / k)
  varioData$SemiVarLow95 <- varioData$SemiVar * CI.lower(varioData$DOF, 0.05)
  varioData$SemiVarUpp95 <- varioData$SemiVar * CI.upper(varioData$DOF, 0.05)
  varioData$animal <- id
  
  return(list(animal = id, varioData = varioData, fits = fits))
})
varioDataAll <- map_dfr(varioDataList, function(x) {
  df <- x$varioData
  df$animal <- x$animal
  return(df)
}) %>% 
  filter(animal != "CrAc4") %>%
  mutate(
    LagDays = Lag / 86400,
    SemiVar = SemiVar / 10000,
    SemiVarLow95 = SemiVarLow95 / 10000,
    SemiVarUpp95 = SemiVarUpp95 / 10000
  )

max_lags <- varioDataAll %>%
  group_by(animal) %>%
  summarise(max_lag = max(LagDays, na.rm = TRUE))

varioDataAll <- varioDataAll %>%
  left_join(max_lags, by = "animal")

ggplot(varioDataAll, aes(x = LagDays, y = SemiVar)) +
  geom_ribbon(aes(ymin = SemiVarLow95, ymax = SemiVarUpp95),
              alpha = 0.3, fill = "gray60") +
  geom_line(size = 0.6, color = "black") +
  ggtitle("")+
  facet_wrap(~animal, scales = "free_y", ncol = 3) +
  coord_cartesian(xlim = c(0, NA)) +
  scale_x_continuous(
    name = "Time lag (days)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(
    name = "Semi-variance (ha)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(filename = "Figures/S_Figure_4.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)

varioFits <- do.call(rbind, lapply(varioDataList, function(x){
  fitDf <- as.data.frame(summary(x$fits))
  names(fitDf)[2] <- "RMSPE_m"
  fitDf$animal <- x$animal
  return(fitDf)
}))
write.csv(varioFits, "Table/ModelFits.csv", row.names = TRUE)
View(varioFits)
topfits <- varioFits %>%
  group_by(animal) %>%
  slice_min(order_by = ΔAICc, n = 1)
sapply(varioDataList, function(x) length(x$fits))
sapply(varioDataList, function(x) is.null(x$fits) || is.null(x$fits[[1]]))

akdeAreas <- list()
akdePolyList <- list()
for (i in seq_along(data_list)) {
  id <- names(data_list)[i]
  telemetry_obj <- data_list[[id]]
  fit_obj <- tryCatch(varioDataList[[i]]$fits[[1]], error = function(e) NULL)
  
  if (is.null(fit_obj) || !inherits(fit_obj, "ctmm")) {
    message("Modelo no disponible para ", id, " - se omite.")
    next
  }
  
  ud <- tryCatch({
    akde(telemetry_obj, fit_obj, weights = TRUE, res = 3)
  }, error = function(e) {
    message("Error al calcular AKDE para ", id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(ud)) next
  
  for (level in c(0.9, 0.95, 0.99)) {
    iso <- tryCatch({
      iso_poly <- SpatialPolygonsDataFrame.UD(ud, level.UD = level)
      iso_poly@data$level <- paste0(level * 100, "%")
      writeOGR(iso_poly, dsn = "Kernel", layer = paste0(id, "_", level * 100),
               driver = "ESRI Shapefile", overwrite_layer = TRUE)
      iso_poly
    }, error = function(e) {
      message("Error al exportar shapefile para ", id, " (nivel ", level, "): ", e$message)
      return(NULL)
    })
  }
  
  area_df <- tryCatch({
    df <- do.call(rbind, lapply(c(0.9, 0.95, 0.99), function(l){
      aDf <- as.data.frame(summary(ud, level.UD = l)$CI)
      aDf$contour <- l
      return(aDf)
    }))
    df$animal <- id
    akdeAreas[[id]] <- df
  }, error = function(e) {
    message("Error al calcular áreas para ", id, ": ", e$message)
    return(NULL)
  })
}
akdeAreasDf <- do.call(rbind, lapply(akdeAreas, function(x){
  units <- gsub("[[:digit:]]|[[:punct:]]", "", row.names(x))
  units <- sub("area ", "", units)
  x$units <- units
  return(x)
}))
rownames(akdeAreasDf) <- NULL

areaTable <- akdeAreasDf %>%
  mutate(Contour = case_when(
    contour == 0.90 ~ 90,
    contour == 0.95 ~ 95,
    contour == 0.99 ~ 99,
  )) %>%
  dplyr::select(Animal = animal, Contour, aKDE_low = low, aKDE_est = est, aKDE_high = high, Units = units)

write.csv(areaTable, "Table/AreaEstimates.csv", row.names = FALSE)
## 1.0.2 Graphing kernel range
shp_files <- list.files("Kernel", pattern = "\\.shp$", full.names = TRUE)
shp_files <- shp_files[!grepl("Tac4", shp_files)]
basename(shp_files)
shp_list <- lapply(shp_files, function(f) {
  s <- st_read(f, quiet = TRUE)
  name <- basename(f)
  
  id <- str_extract(name, "T[a-z]+\\d+")
  level <- str_extract(name, "(90|95|99)")
  
  s$animal <- id
  s$level <- paste0(level, "%")
  return(s)
})

shp_all <- do.call(rbind, shp_list) %>%
  filter(level %in% c("90%", "95%", "99%")) %>%
  mutate(level = factor(level, levels = c("90%", "95%", "99%")))%>%
  mutate(geometry = st_make_valid(geometry)) %>%
  group_by(animal, level) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_as_sf()
basename(shp_files)

colnames(shp_all)
points_df <- telemetry %>%
  st_as_sf(coords = c("x", "y"), crs = 32618) 
plot_list <- list()
shp_all <- shp_all %>%
  mutate(
    animal = case_when(
      str_detect(animal, "^Tca") ~ str_replace(animal, "^Tca", "CaCr"),
      str_detect(animal, "^Tac") ~ str_replace(animal, "^Tac", "CrAc"),
      TRUE ~ animal
    )
  )
animals <- unique(shp_all$animal)
unique(shp_all$animal)
unique(shp_all$level)
unique(points_df$id)
unique(shp_all$animal)
sapply(animals, function(a) {
  sum(points_df$id == a, na.rm = TRUE)
})
for (i in seq_along(animals)) {
  animal_data <- shp_all %>% filter(animal == animals[i])
  points_data <- points_df %>% filter(id == animals[i])
  
  plot_list[[i]] <- ggplot() +
    geom_sf(data = animal_data, aes(fill = level), color = "black", alpha = 0.3) +
    geom_sf(data = points_data, aes(geometry = geometry), color = "black", size = 0.6) +
    scale_fill_manual(values = c("90%" = "dodgerblue4", 
                                 "95%" = "dodgerblue1", 
                                 "99%" = "skyblue1")) +
    annotation_scale(location = "bl", 
                     bar_cols = c("black", "black"),
                     width_hint = 0.20,  
                     text_cex = 0.5,
                     unit_category = "metric",
                     unit = "m",
                     height = unit(0.1, "cm"),
                     pad_x = unit(0, "cm"),
                     pad_y = unit(0.05, "cm"),
                     line_width = 0.4) +
    ggtitle(
      str_replace(
        str_replace(animals[i], "^Tac", "CrAc"),
        "^Tca", "CaCr"
      )
    ) +
    coord_sf(datum = NA) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}
sapply(plot_list, class)

plot_list_clean <- plot_list[!sapply(plot_list, is.null)]
which(sapply(plot_list, is.null))

wrap_plots(plot_list_clean, ncol = 3) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
ggsave(filename = "Figures/Figure_4.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)


## 1.1 Preprocessing and parameter setting----
all.fort.poly <- vector(mode = "list", length = length(unique(all.data$id)))
data.var.list <- vector(mode = "list", length = length(unique(all.data$id)))
all.area.list <- vector(mode = "list", length = length(unique(all.data$id)))
dbbmm.list <- vector(mode = "list", length = length(unique(all.data$id)))

## 1.2 dBBMM for each animal----
### 1.2.1 CrAc11 (Joselito)----
data <- all.data[all.data$id == 'CrAc11', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)
?brownian.motion.variance.dyn
dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)
?brownian.bridge.dyn
dbbv@interest[timeLag(move,"mins")>410] <- FALSE
dbbmm_Tac11 <- brownian.bridge.dyn(object=dbbv, location.error=5, dimSize=500, ext=2,
                                   time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mihttp://127.0.0.1:38519/graphics/bec05dbd-14e1-4bdc-a528-e745fdf27a40.pngns"))

cont4 <- raster2contour(dbbmm_Tac11, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_tac11_contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)

### 1.2.2 CaCr7 (Lulo)----
data <- all.data[all.data$id == 'CaCr7', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)

dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)

dbbv@interest[timeLag(move,"mins")>410] <- FALSE

dbbmm_Tca7 <- brownian.bridge.dyn(object=dbbv, location.error=7, dimSize=500, ext=2,
                                   time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

cont4 <- raster2contour(dbbmm_Tca7, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_Tca7_contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)

poly4 <- rasterToPolygons(dbbmm_Tca7, fun=function(x){x >= quantile(x, probs=c(0.90, 0.95, 0.99))}, dissolve=TRUE)
writeOGR(obj = poly4, dsn = loc.output, layer = "dbbmm_Tca7_polygons", driver = "ESRI Shapefile", overwrite_layer = TRUE)

### 1.2.3 CaCr10 (Hugo)----
data <- all.data[all.data$id == 'CaCr10', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)

dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)

dbbv@interest[timeLag(move,"mins")>410] <- FALSE

dbbmm_Tca10 <- brownian.bridge.dyn(object=dbbv, location.error=7, dimSize=500, ext=2,
                                  time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

cont4 <- raster2contour(dbbmm_Tca10, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_Tca10_contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)

### 1.2.4 CaCr11 (Rosa)----
data <- all.data[all.data$id == 'CaCr11', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)

dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)

dbbv@interest[timeLag(move,"mins")>410] <- FALSE

dbbmm_Tca11 <- brownian.bridge.dyn(object=dbbv, location.error=7, dimSize=500, ext=2,
                                  time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

cont4 <- raster2contour(dbbmm_Tca11, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_Tca11_contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)

### 1.2.5 CaCr18 (Elena)----
data <- all.data[all.data$id == 'CaCr18', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)

dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)

dbbv@interest[timeLag(move,"mins")>410] <- FALSE

dbbmm_Tca18 <- brownian.bridge.dyn(object=dbbv, location.error=7, dimSize=500, ext=2,
                                   time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

cont4 <- raster2contour(dbbmm_Tca18, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_Tca18contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)

### 1.2.6 CaCr20 (La Gorda)----
data <- all.data[all.data$id == 'CaCr20', ]

data$gps.accuracy.filled <- ifelse(is.na(data$gps.accuracy), 0, data$gps.accuracy)

move <- move(x = data$x, y = data$y, proj = crs.proj, time = as.POSIXct(data$datetime))

ind.start <- Sys.time()

proj4string(move)

dbbv <- brownian.motion.variance.dyn(move, location.error=7, window.size=7, margin=3)

dbbv@interest[timeLag(move,"mins")>410] <- FALSE

dbbmm_Tca20 <- brownian.bridge.dyn(object=dbbv, location.error=7, dimSize=500, ext=2,
                                   time.step=60, margin=3, window.size = 7, proj = crs.proj)

print(paste("------ dBBMM computation time:",
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

cont4 <- raster2contour(dbbmm_Tca20, level=c(.90,.95,.99))
writeOGR(obj = cont4, dsn = loc.output, layer = "dbbmm_Tca20contours", driver = "ESRI Shapefile", overwrite_layer = TRUE)


# 2.0 Habitat use----
## 2.1 Checking autocorrelation
list.files("./Rasters")
cov.names <- c(
  "Dist_Be",
  "Dist_MCPNS",
  "Dist_PV",
  "Dist_RFB",
  "Dist_SV",
  "Dist_WPWA"
)
covariates <- list.files("./Rasters", pattern = "\\.tif$", full.names = TRUE)
covariates <- covariates[!grepl("\\.aux\\.xml$|\\.tfw$|\\.xml$", covariates)]
base_raster <- raster(covariates[1])
resampled_list <- lapply(covariates, function(f) {
  r <- raster(f)
  resample(r, base_raster, method = "bilinear")
})

rp <- stack(resampled_list)
names(rp) <- cov.names
plot(rp)
globalCor <- cor(values(rp), use = "pairwise.complete.obs")
print(globalCor)
corrplot(globalCor, method = "color", type = "upper", tl.cex = 0.8)

highCorPairs <- which(abs(globalCor) > 0.85 & abs(globalCor) < 1, arr.ind = TRUE)
apply(highCorPairs, 1, function(idx) paste(rownames(globalCor)[idx[1]], "-", colnames(globalCor)[idx[2]]))

vars_seleccionadas <- c(
  "Dist_Be",
  "Dist_MCPNS",
  "Dist_PV",
  "Dist_RFB",
  "Dist_WPWA"
)

globalCor_sub <- globalCor[vars_seleccionadas, vars_seleccionadas]

plot.new()
dev.off()
corrplot(globalCor_sub,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.5,  # REDUCIDO
         addCoef.col = "black",
         number.cex = 0.7,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         diag = FALSE)

## 2.2 

all.ssfResults <- NULL
ind.ssfResults <- NULL
c_ids <- setdiff(unique(all.data$id), c("CrAc77", "CrAc4", "CaCr11","CaCr7"))
mods <- list()

# Single-factor models
mods[[1]]  <- case_ ~ log_sl * cos_ta + strata(step_id_)
mods[[2]]  <- case_ ~ Dist_WPWA + Dist_WPWA:log_sl + Dist_WPWA:cos_ta + log_sl * cos_ta + strata(step_id_)
mods[[3]]  <- case_ ~ Dist_RFB + Dist_RFB:log_sl + Dist_RFB:cos_ta + log_sl * cos_ta + strata(step_id_)
mods[[4]]  <- case_ ~ Dist_Be + Dist_Be:log_sl + Dist_Be:cos_ta + log_sl * cos_ta + strata(step_id_)
mods[[5]]  <- case_ ~ Dist_PV + Dist_PV:log_sl + Dist_PV:cos_ta + log_sl * cos_ta + strata(step_id_)
mods[[6]]  <- case_ ~ Dist_MCPNS + Dist_MCPNS:log_sl + Dist_MCPNS:cos_ta + log_sl * cos_ta + strata(step_id_)

# Multi-factor models
mods[[7]]  <- case_ ~ Dist_WPWA + Dist_RFB + log_sl * cos_ta + strata(step_id_)
mods[[8]]  <- case_ ~ Dist_WPWA + Dist_Be + log_sl * cos_ta + strata(step_id_)
mods[[9]]  <- case_ ~ Dist_Be + Dist_PV + Dist_MCPNS + log_sl * cos_ta + strata(step_id_)
mods[[10]] <- case_ ~ Dist_RFB + Dist_PV + Dist_MCPNS+ log_sl * cos_ta + strata(step_id_)
mods[[11]] <- case_ ~ Dist_MCPNS + Dist_WPWA + Dist_RFB + Dist_Be + Dist_PV + log_sl * cos_ta + strata(step_id_)
mods[[12]] <- case_ ~ Dist_WPWA:log_sl + Dist_RFB:log_sl + Dist_Be:log_sl + Dist_PV:log_sl + log_sl * cos_ta + strata(step_id_)

## 2.3

# SSF Loop per Individual

i <- 0
ssfResults <- vector(mode = "list", length = length(c_ids))
for(crocodylian in c_ids){
  i <- i+1
  data.1 <- dplyr::filter(data_presence, id == crocodylian)
  tr1 <- mk_track(data.1, x, y, datetime, crs = crs.proj, id = id)
  
  ssf1 <- tr1 %>%
    steps() %>%
    filter(sl_ > 5)  
  
  ssf1 <- ssf1 %>% 
    random_steps(n = 200) %>% 
    extract_covariates(rp, where = "end") 
  
  ssf1 <- ssf1 %>% 
    mutate(
      Dist_WPWA = log(Dist_WPWA + 1),
      Dist_RFB = log(Dist_RFB + 1),
      Dist_Be = log(Dist_Be + 1),
      Dist_PV = log(Dist_PV + 1),
      Dist_MCPNS = log(Dist_MCPNS + 1),
      cos_ta = cos(ta_), 
      log_sl = log(sl_ + 1)
    )
  
  ssf1$id <- crocodylian 
  
  ?fit_issf()
  
  ind.ssfResults <- NULL
  iter <- 0
  for(model in mods){
    iter <- iter + 1
    ssffit <-  fit_issf(ssf1, formula = model)
    results_ssf <- tidy(ssffit$model, conf.int = TRUE, conf.level = 0.95)
    results_ssf$id <- crocodylian
    results_ssf$model <- paste0('model',iter)
    results_ssf$AIC <- AIC(ssffit$model)
    ind.ssfResults <- rbind(ind.ssfResults, results_ssf)
  }
  ssfResults[[i]] <- ind.ssfResults
} 

all.ssfResults <- do.call(rbind, ssfResults)
write_csv(all.ssfResults, "Table/All_ssf_results.csv")

# Create AIC table for all individuals
ssf.AIC.results <- all.ssfResults %>%
  group_by(model, id) %>%
  summarise(AIC = mean(AIC)) %>%
  arrange(AIC)


ssf.AIC.results
write_csv(ssf.AIC.results, "Table/ssfAIC_results.csv")

ssf.AIC.results.t2 <- ssf.AIC.results %>% 
  group_by(id) %>% 
  mutate(AIC = ifelse(AIC < min(AIC)+2, paste0("*",
                                               round(AIC, digits = 2),
                                               "*"),
                      round(AIC, digits = 2))) %>%
  mutate(model = factor(model, levels = c("model1", "model2", "model3", "model4", "model5", "model6", "model7", "model8", "model9", "model10", "model11", "model12" ))) %>%
  arrange(model)

# ssf.AIC.results.t2$AIC <- round(ssf.AIC.results$AIC)
table.2 <- reshape2::dcast(ssf.AIC.results.t2, model ~ id)
table.2 
write_csv(table.2, "Table/AIC scores and models.csv")

table.2$model <- paste0("model", 1:nrow(table.2))
table.2$formula <- c(
  "log_sl*cos_ta + strata(step_id_)",
  "Dist_WPWA + Dist_WPWA:log_sl + Dist_WPWA:cos_ta + log_sl*cos_ta + strata(step_id_)",
  "Dist_RFB + Dist_RFB:log_sl + Dist_RFB:cos_ta + log_sl*cos_ta + strata(step_id_)",
  "Dist_Be + Dist_Be:log_sl + Dist_Be:cos_ta + log_sl*cos_ta + strata(step_id_)",
  "Dist_PV + Dist_PV:log_sl + Dist_PV:cos_ta + log_sl * cos_ta + strata(step_id_)",
  "Dist_MCPNS + Dist_MCPNS:log_sl + Dist_MCPNS:cos_ta + log_sl*cos_ta + strata(step_id_)",
  "Dist_WPWA + Dist_RFB + log_sl * cos_ta + strata(step_id_)",
  "Dist_WPWA + Dist_Be + log_sl * cos_ta + strata(step_id_)",
  "Dist_Be + Dist_PV + Dist_MCPNS+ log_sl * cos_ta + strata(step_id_)",
  "Dist_RFB + Dist_PV + Dist_MCPNS+ log_sl * cos_ta + strata(step_id_)",
  "Dist_MCPNS + Dist_WPWA + Dist_RFB + Dist_Be + Dist_PV + log_sl * cos_ta + strata(step_id_)",
  "Dist_WPWA:log_sl + Dist_RFB:log_sl + Dist_Be:log_sl + Dist_PV:log_sl + log_sl * cos_ta + strata(step_id_)"
)

col_order <- c("model", "formula", setdiff(names(table.2), c("model", "formula")))
table.2 <- table.2[, col_order]

names(table.2)[names(table.2) == "model"] <- "Model #"
names(table.2)[names(table.2) == "formula"] <- "Model formula"

write.csv(table.2, "Table/AIC_scores_and_models.csv", row.names = FALSE)

# Return the top model for each individual

top.models <- ssf.AIC.results %>%
  group_by(id) %>%
  filter(AIC == min(AIC)) %>%
  arrange(model)

top.models

topmodelsdf <- as.data.frame(top.models)

write_csv(top.models, "Table/ssf_AICresults_topmodels.csv")

## 2.4 
ssf.results <- read.csv(file = "Table/All_ssf_results.csv", stringsAsFactors = FALSE)

### dist_feature plot
ssf.results %>% 
  group_by(id) %>% 
  mutate(best = ifelse(AIC <= min(AIC)+2, "*", "")) %>% 
  filter(term %in% grep("log_sl|cos_ta", term, value = TRUE,
                        invert = TRUE)) %>% 
  group_by(id) %>% 
  arrange(term) %>% 
  mutate(plot.order = row_number()) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5) +
  geom_errorbar(aes(x = plot.order, ymin = conf.low, ymax = conf.high,
                    colour = term),
                width = 0.25) +
  geom_point(aes(x = plot.order, y = estimate, colour = term,
                 shape = best)) +
  facet_wrap(
    . ~ id,
    scales = "free_y",
    labeller = labeller(
      id = c(
        "Tac11" = "CrAc11",
        "Tca10" = "CaCr10",
        "Tca18" = "CaCr18",
        "Tca20" = "CaCr20"
      )
    )
  ) +
  scale_shape_manual(values = c(3, 16)) +
  scale_x_continuous(labels = c("Be", "MCPNS", "PV", "RFB", "WPWA"),
                     breaks = c(2.5, 6.5, 10.5, 14.5, 18.5),
                     minor_breaks = c(4.5, 8.5, 12.5,16.6)
  ) +
  labs(x = "Feature", y = "Coefficient estimate (β)", colour = "Distance to\nfeature:") +
  scale_colour_manual(values = c("#009E73", "brown1", "#E69F00", "#0072B2", "#CC79A7"),
                      labels = c("Beach",
                                 "Mosaic of crops, pastures and natural spaces",
                                 "Primary vegetation",
                                 "Recreational facilities and buildings",
                                 "Weedy pastures and wooded areas")) +
  theme_bw() +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_none()
  )+ 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = 4, hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(colour = "grey65", linetype = 1),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = 2),
        legend.text = element_text(lineheight = 1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, face = 2, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(hjust = 0.5, face = 2, margin = margin(10,0,0,0)),
        plot.title = element_text(face = 4),
        strip.text.y = element_blank()
  ) +
  guides(shape = guide_none())

ggsave(file = "./Figures/Figure6.jpg", width = 200, height = 150,
       dpi = 400, units = "mm")

# 3.0 Movements----
## 3.1 Check temporary data format----
all.data$datetime <- as.POSIXct(all.data$datetime, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
data_mov <- all.data %>%
  filter(!is.na(id), !is.na(datetime)) %>% 
  filter(!id %in% c("Tac4", "Tac77")) %>%  
  arrange(id, datetime)

## 3.2 Calculate displacement between points----
data_mov <- data_mov %>%
  group_by(id) %>%
  mutate(displacement = sqrt((x - lag(x))^2 + (y - lag(y))^2)) %>%
  ungroup()

average_displacement <- data_mov %>%
  group_by(id) %>%
  summarise(mean_displacement = mean(displacement, na.rm = TRUE),
            sd_displacement = sd(displacement, na.rm = TRUE),  
            total_displacement = sum(displacement, na.rm = TRUE),
            num_observations = n()) %>%
  ungroup()

monitoring_days <- data_mov %>%
  group_by(id) %>%
  summarise(days_monitored = n_distinct(as.Date(datetime))) %>%
  ungroup()

summary_data <- average_displacement %>%
  left_join(monitoring_days, by = "id")

print(summary_data)

## 3.3 Add and check Week and Month Columns---- 
data_mov <- data_mov %>%
  mutate(week = floor_date(datetime, "week"), 
         month = floor_date(datetime, "month"))

month_colors <- data.frame(
  month = factor(month.abb, levels = month.abb), 
  color = c("darksalmon", "darksalmon", "lightblue",  
            "lightblue", "lightblue", "darksalmon",        
            "darksalmon", "darksalmon", "lightblue",   
            "lightblue", "lightblue", "darksalmon"))

weekly_displacement <- data_mov %>%
  group_by(id, week) %>%
  summarise(
    total_displacement = sum(displacement, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(id) %>%
  complete(
    week = seq(min(week), max(week), by = "week"),
    fill = list(total_displacement = 0)
  ) %>%
  ungroup() %>%
  filter(!is.na(id)) %>%
  mutate(
    month = factor(format(week, "%b"), levels = month.abb),
    week = as.Date(week)
  )

str(weekly_displacement)

weekly_displacement <- weekly_displacement %>%
  mutate(week = as.Date(week))

## 3.4 Plot the movements of all individuals----
ggplot(weekly_displacement, aes(x = week, y = total_displacement)) +
  geom_rect(data = unique(weekly_displacement[, c("week", "month")]),
            aes(xmin = floor_date(week, "month"), 
                xmax = ceiling_date(week, "month") - days(1), 
                ymin = -Inf, ymax = Inf, fill = month),
            alpha = 0.2, inherit.aes = FALSE) +  
  geom_line(size = 1.0) +
  scale_fill_manual(values = setNames(month_colors$color, month_colors$month), guide = "none") +
  labs(x = "Date (Month/Year)",y = "Total movement (m)") +
  theme_bw() +
  theme(
    strip.text = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray50", size = 0.5), 
    panel.grid.minor = element_line(color = "gray50", size = 0.5) ) +
  scale_x_date(date_breaks = "1 month",   
    date_labels = "%b %Y",   expand = c(0, 0)) +
  facet_wrap(~ id, ncol = 1, scales = "free_y") +
  geom_text(
    data = weekly_displacement %>% 
      group_by(id) %>% 
      summarise(
        min_week = min(weekly_displacement$week, na.rm = TRUE),
        max_y = max(total_displacement, na.rm = TRUE),
        .groups = "drop"
      ),
    aes(
      x = min_week + weeks(48),
      y = max_y,
      label = id
    ),
    hjust = 1,
    vjust = 1,
    fontface = "bold",
    size = 5,
    inherit.aes = FALSE
  )

ggsave(filename = "Figures/Figure_2.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)

# 4.0 Impact of environmental variables on movement----
colSums(is.na(movement_data))
str(movement_data) 
View(movement_data)
## 4.1 CrAc11----
data_Tac11 <- movement_data %>%
  filter(id %in% c("CrAc11"))


Tac11_gamma <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
            s(average_weekly_solar_radiation) + 
            s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
            s(average_weekly_maxtemperature) + s(tourists), 
          data = data_Tac11, method = "REML", 
          family = Gamma(link = "log"), select = FALSE)

Tac11_tweedie <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                       s(average_weekly_solar_radiation) + 
                       s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                       s(average_weekly_maxtemperature) + s(tourists), 
                  data = data_Tac11, method = "REML", 
                  family = tw(), select = FALSE)

Tac11_gaussian <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                        s(average_weekly_solar_radiation) + 
                        s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                        s(average_weekly_maxtemperature) + s(tourists), 
                   data = data_Tac11, method = "REML", 
                   family = gaussian, select = FALSE)

summary(Tac11_gamma)
summary(Tac11_tweedie)
summary(Tac11_gaussian)
AIC(Tac11_gamma, Tac11_tweedie, Tac11_gaussian)

appraise(Tac11_gamma, method = 'simulate') #Gamma is the best
gam.check(Tac11_gamma)


#Gamma is the best
appraise(Tac11_gamma, method = 'simulate') #Gamma is the best

ggsave(filename = "Figures/Figure_STac11.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)

A<-plot_predictions(Tac11_gamma, condition = "average_weekly_precipitation", type = 'response', rug = T) + 
  labs(title="A.",y = "Expected response", x = "Average precipitation (mm)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
C<-plot_predictions(Tac11_gamma, condition = "average_weekly_temperature", type = 'response', rug = T) + 
  labs(title="C.",y = "", x = "Average temperature (°C)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
B<-plot_predictions(Tac11_gamma, condition = "average_weekly_mintemperature", type = 'response', rug = T) + 
  labs(title="B.",y = "", x = "Minimum temperature (°C)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()

A+B+C

## 4.4 CaCr11----
data_Tca11 <- movement_data %>%
  filter(id %in% c("CaCr11"))

Tca11_gamma <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                     s(average_weekly_solar_radiation) + 
                     s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                     s(average_weekly_maxtemperature) + s(tourists), 
                   data = data_Tca11, method = "REML", 
                   family = Gamma(link = "log"), select =  TRUE)

Tca11_tweedie <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                       s(average_weekly_solar_radiation) + 
                       s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                       s(average_weekly_maxtemperature) + s(tourists), 
                     data = data_Tca11, method = "REML", 
                     family = tw(), select =  TRUE)

Tca11_gaussian <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                        s(average_weekly_solar_radiation) + 
                        s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                        s(average_weekly_maxtemperature) + s(tourists), 
                      data = data_Tca11, method = "REML", 
                      family = gaussian, select =  TRUE)

summary(Tca11_gamma)
summary(Tca11_tweedie)
summary(Tca11_gaussian)

AIC(Tca11_gamma, Tca11_tweedie, Tca11_gaussian)
gam.check(Tca11_gamma)
appraise(Tca11_gamma, method = 'simulate') #The best
ggsave(filename = "Figures/S_Figure_Tca11.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)

Aa<-plot_predictions(Tca11_gamma, condition = "average_weekly_precipitation", type = 'response', rug = T) + 
  labs(title="D.",y = "", x = "Average precipitation (mm)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
Aa
Bb<-plot_predictions(Tca11_gamma, condition = "weekly_humidity", type = 'response', rug = T) + 
  labs(title="E.",y = "", x = "Average humidity (%)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
Bb
Cc<-plot_predictions(Tca11_gamma, condition = "average_weekly_maxtemperature", type = 'response', rug = T) + 
  labs(title="F.",y = "", x = "Maximum temperature (°C)") +
  guides(x =  guide_axis(angle = 45)) +
  theme_bw()
Cc
## 4.5 CaCr18----
data_Tca18 <- movement_data %>%
  filter(id %in% c("CaCr18"))

Tca18_gamma1 <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                      s(average_weekly_solar_radiation) + 
                      s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                      s(average_weekly_maxtemperature) + s(tourists), 
                   data = data_Tca18, method = "REML", 
                   family = Gamma(link = "log"), select =  TRUE)

Tca18_tweedie1 <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                        s(average_weekly_solar_radiation) + 
                        s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                        s(average_weekly_maxtemperature) + s(tourists), 
                     data = data_Tca18, method = "REML", 
                     family = tw(), select =  TRUE)

Tca18_gaussian1 <- gam(movementss ~ s(average_weekly_precipitation) + s(weekly_humidity) + 
                         s(average_weekly_solar_radiation) + 
                         s(average_weekly_temperature) + s(average_weekly_mintemperature) + 
                         s(average_weekly_maxtemperature) + s(tourists), 
                      data = data_Tca18, method = "REML", 
                      family = gaussian, select = TRUE)

summary(Tca18_tweedie1)
summary(Tca18_gaussian1)

AIC(Tca18_tweedie1, Tca18_gaussian1)

appraise(Tca18_tweedie1, method = 'simulate')#The best 
appraise(Tca18_gaussian1, method = 'simulate') 
ggsave(filename = "Figures/Figure_STca18.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)


## 4.7 Graphing results----
(A+B+C)/(Aa+Bb+Cc)
ggsave(filename = "Figures/Figure_3.jpg", plot = last_plot(), width = 3000, height = 2500, 
       units = "px", dpi = 300)

