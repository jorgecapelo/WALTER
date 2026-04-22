## ============================================================
##  World Macrobioclimate Map — Rivas-Martínez (2004)
##  Penultimate Mediterranean criterion:
##    Ios2 < 2  AND  Iosc4 < 2  →  Mediterranean
##    Ios2 ≥ 2  AND  Iosc4 ≥ 2  →  Temperate
##    mixed                       →  Walter tiebreaker
##
##  Data:  WorldClim v2.1 at 10 arc-min  (~45 MB download)
##  Sampling: aggregated ×9  →  240×120 = 28,800 cells
##            land pixels:  ~8,400  (well below 10,000 limit)
##  Display:  NN-resampled to 0.5° for smooth rendering,
##            then projected to Robinson
##
##  Output:  world_macrobioclimates.png  (working directory)
##
##  Runs in ~2 min after data download (single-threaded R loop).
## ============================================================


# ── 0.  Packages ──────────────────────────────────────────────
required_pkgs <- c("terra", "geodata")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# ── 1.  Configuration ─────────────────────────────────────────
WC_CACHE    <- file.path(getwd(), "worldclim_cache")
AGG_FACTOR  <- 9      # 10 arc-min × 9 = 90 arc-min ≈ 1.5°
#                       → 240×120 grid, ~8,400 land cells
DISPLAY_RES <- 0.5    # degrees — NN-resampled output for display
OUT_FILE    <- "world_macrobioclimates.png"

dir.create(WC_CACHE, recursive = TRUE, showWarnings = FALSE)
message("Cache dir  : ", WC_CACHE)
message("Agg factor : ", AGG_FACTOR, "  (classification grid ≈ 1.5°)")
message("Display res: ", DISPLAY_RES, "°  (nearest-neighbour upscale)")


# ── 2.  Download & aggregate WorldClim at 10 arc-min ──────────
message("\nLoading WorldClim rasters (10 arc-min) ...")

load_wc <- function(var) {
  r <- geodata::worldclim_global(var, res = 10, path = WC_CACHE)
  message("  ", var, " loaded — aggregating by ", AGG_FACTOR, "× ...")
  aggregate(r, fact = AGG_FACTOR, fun = "mean", na.rm = TRUE)
}

r_prec <- load_wc("prec")
r_tavg <- load_wc("tavg")
r_tmin <- load_wc("tmin")
r_tmax <- load_wc("tmax")
message("  Done. Classification grid: ",
        nrow(r_tavg), " rows × ", ncol(r_tavg), " cols")


# ── 3.  Latitude raster ───────────────────────────────────────
#  Needed to determine hemisphere and apply lat-band thresholds.
r_lat <- init(r_tavg[[1]], "y")   # cell-centre latitude


# ── 4.  Classification helpers ────────────────────────────────

## 4a.  Itc correction C  (Table 3.8)
.itc_c <- function(Ic) {
  if      (Ic <=  8) 10  * (Ic -  8)
  else if (Ic <= 18) 0
  else if (Ic <= 21) 5   * (Ic - 18)
  else if (Ic <= 28) 15  + 15 * (Ic - 21)
  else if (Ic <= 46) 120 + 25 * (Ic - 28)
  else               570 + 30 * (min(Ic, 65) - 46)
}

## 4b.  Ios2  (2 hottest summer months, T > 0 only)
.ios2 <- function(T, P, sum_m) {
  Ts   <- T[sum_m]; Ps <- P[sum_m]
  h2   <- order(Ts, decreasing = TRUE)[1:2]
  pos  <- Ts[h2] > 0
  Tps2 <- sum(Ts[h2][pos]) * 10
  Pps2 <- sum(Ps[h2][pos])
  if (Tps2 > 0) (Pps2 / Tps2) * 10 else Inf
}

## 4c.  Iosc4  (prev month + 3 summer months, T > 0 only)
.iosc4 <- function(T, P, sum_m, prev_m) {
  four <- c(prev_m, sum_m); T4 <- T[four]; P4 <- P[four]
  pos  <- T4 > 0
  Tps4 <- sum(T4[pos]) * 10; Pps4 <- sum(P4[pos])
  if (Tps4 > 0) (Pps4 / Tps4) * 10 else Inf
}

## 4d.  Mediterranean?  — penultimate two-index criterion
.is_med <- function(ios2_v, iosc4_v, T, P, sum_m) {
  if (ios2_v <  2 && iosc4_v <  2) return(TRUE)
  if (ios2_v >= 2 && iosc4_v >= 2) return(FALSE)
  # Mixed: Walter tiebreaker — ≥2 consecutive P ≤ 2T in 5-month warm window
  prev_w <- (sum_m[1] - 2L) %% 12L + 1L
  next_w <-  sum_m[3]        %% 12L + 1L
  win    <- c(prev_w, sum_m, next_w)
  dry    <- P[win] <= 2 * T[win]
  runs   <- rle(dry)
  csec   <- runs$lengths[runs$values]
  length(csec) > 0 && max(csec) >= 2
}

## 4e.  Boreal conditions  (Table 3.9)
.check_boreal <- function(Ic, T_ann, Tp) {
  if      (Ic <= 11) T_ann <= 6.0 && Tp > 0
  else if (Ic <= 21) T_ann <= 5.3 && Tp > 380 && Tp <= 720
  else if (Ic <= 28) T_ann <= 4.8 && Tp > 380 && Tp <= 740
  else if (Ic <  46) T_ann <= 3.8 && Tp > 380 && Tp <= 800
  else               T_ann <= 0.0 && Tp > 380 && Tp <= 800
}

## 4f.  Classify one pixel → integer code (or NA)
##
##  Macrobioclimate codes:
##    1 = Tropical   2 = Mediterranean   3 = Temperate
##    4 = Boreal     5 = Polar
##
.classify <- function(P, T, Tmin, Tmax, lat) {

  if (any(is.na(c(P, T, Tmin, Tmax)))) return(NA_integer_)

  abs_lat  <- abs(lat)
  NH       <- lat >= 0

  ## Hemisphere-aware season indices (1-based)
  if (NH) {
    sum_m  <- 6:8; prev_m <- 5L; win_m <- c(12L, 1L, 2L)
  } else {
    sum_m  <- c(12L, 1L, 2L); prev_m <- 11L; win_m <- 6:8
  }

  ## Basic thermal parameters
  T_ann  <- mean(T)
  cold_m <- which.min(T); hot_m <- which.max(T)
  m_val  <- Tmin[cold_m]; M_val <- Tmax[cold_m]

  ## Indices
  Ic   <- T[hot_m] - T[cold_m]
  It   <- (T_ann + m_val + M_val) * 10
  Itc  <- It + .itc_c(Ic)
  Tp   <- sum(T[T > 0]) * 10
  Pp   <- sum(P[T > 0])
  Io   <- if (Tp > 0) (Pp / Tp) * 10 else 0

  ## Summer seasonality
  ios2_v  <- .ios2(T, P, sum_m)
  iosc4_v <- .iosc4(T, P, sum_m, prev_m)
  s_med   <- .is_med(ios2_v, iosc4_v, T, P, sum_m)
  Pcm1    <- sum(P[sum_m]);  Ps_w <- sum(P[win_m])
  s_rain  <- Pcm1 > Ps_w    # summer wetter than winter → Tropical pattern

  ## Decision tree
  if (abs_lat > 66) return(5L)   # Polar

  if (abs_lat > 52)
    return(if (.check_boreal(Ic, T_ann, Tp)) 4L else 3L)

  if (abs_lat > 35)
    return(if (s_med) 2L else 3L)

  if (abs_lat > 23) {
    if (s_med) return(2L)
    warm_a <- sum(c(T_ann >= 21, M_val >= 18, Itc >= 470)) >= 2
    warm_b <- sum(c(T_ann >= 25, m_val >= 10, Itc >= 580)) >= 2
    if ((warm_a || warm_b) && s_rain) return(1L)
    return(3L)
  }

  1L   # Tropical (equatorial / eutropical belt)
}


# ── 5.  Pixel-wise classification ─────────────────────────────
message("\nClassifying pixels ...")

P_m    <- as.matrix(r_prec)
T_m    <- as.matrix(r_tavg)
Tmin_m <- as.matrix(r_tmin)
Tmax_m <- as.matrix(r_tmax)
lat_v  <- as.numeric(values(r_lat))

n      <- nrow(P_m)
code_v <- integer(n)

for (i in seq_len(n)) {
  if (any(is.na(T_m[i, ]))) { code_v[i] <- NA_integer_; next }
  code_v[i] <- .classify(P_m[i,], T_m[i,], Tmin_m[i,], Tmax_m[i,], lat_v[i])
}

## Tally
land_n <- sum(!is.na(code_v))
cat("\n  Land pixels classified:", land_n, "\n")
tbl <- table(code_v, useNA = "no")
labs <- c("1"="Tropical","2"="Mediterrânico","3"="Temperado","4"="Boreal","5"="Polar")
for (k in names(tbl))
  cat(sprintf("    %d  %-15s : %5d  (%.1f%%)\n",
              as.integer(k), labs[k], tbl[k], 100*tbl[k]/land_n))


# ── 6.  Build & display categorical raster ────────────────────

## 6a.  Coarse classified raster
r_class         <- r_tavg[[1]]
values(r_class) <- code_v

## 6b.  Smooth to DISPLAY_RES° via nearest-neighbour resampling
r_display <- rast(extent = ext(r_class),
                  resolution = DISPLAY_RES,
                  crs        = crs(r_class))
r_display  <- resample(r_class, r_display, method = "near")

## 6c.  Assign levels and colour table
lev_df <- data.frame(
  value = 1:5,
  label = c("Tropical", "Mediterrânico", "Temperado", "Boreal", "Polar")
)
levels(r_display) <- lev_df
bioclima_cols <- c(
  "#E8A020",   # 1 Tropical      — warm amber
  "#C04010",   # 2 Mediterranean  — terracotta
  "#4A7C3F",   # 3 Temperate      — forest green
  "#3B6EA5",   # 4 Boreal         — steel blue
  "#B8B8D0"    # 5 Polar          — silvery lavender
)
coltab(r_display) <- data.frame(value = 1:5, col = bioclima_cols)

## 6d.  Project to Robinson for display
message("\nReprojecting to Robinson ...")
robin_crs  <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
r_robin    <- project(r_display, robin_crs, method = "near")

## 6e.  World boundaries
message("Loading world boundaries ...")
wld        <- geodata::world(path = WC_CACHE, resolution = 1)
wld_robin  <- project(wld, robin_crs)


# ── 7.  Plot ──────────────────────────────────────────────────
message("Plotting ...")

png(OUT_FILE, width = 3200, height = 1800, res = 240, bg = "#D0E8F5")

## Layout: map | legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4.2, 1.0))

## ── Panel 1: map ──────────────────────────────────────────────
par(mar = c(2.5, 1, 3.5, 0.5), bg = "#D0E8F5")
plot(r_robin,
     legend = FALSE, axes = FALSE, box = FALSE,
     col    = bioclima_cols,
     main   = "")

## Graticule (every 30°)
for (lon in seq(-180, 180, by = 30)) {
  lats <- seq(-89, 89, length.out = 180)
  xy   <- project(cbind(lon, lats), from = "EPSG:4326", to = robin_crs)
  lines(xy, col = "white", lwd = 0.4, lty = 2)
}
for (lat in seq(-60, 60, by = 30)) {
  lons <- seq(-179, 179, length.out = 360)
  xy   <- project(cbind(lons, lat), from = "EPSG:4326", to = robin_crs)
  lines(xy, col = "white", lwd = 0.4, lty = 2)
}

## Country boundaries
lines(wld_robin, col = "gray20", lwd = 0.4)

## Ocean outline (rectangle enclosing the projection)
## Approximate Robinson bounding ellipse via clipped rectangle
bb_x <- c(-17005833, 17005833)
bb_y <- c(-8625154, 8625154)

## Titles
mtext("Macrobioclimas do Mundo  —  Rivas-Martínez (2004)",
      side = 3, line = 2.0, cex = 1.05, font = 2, col = "gray10")
mtext(paste0("WorldClim v2.1  |  10 arc-min × agg. ",
             AGG_FACTOR, "×  (~", land_n, " pontos de terra)"),
      side = 3, line = 0.8, cex = 0.65, col = "gray35")

## Source line at bottom
mtext(paste0("Fronteira Med/Temp: Ios2 < 2 e Iosc4 < 2  |  ",
             "Classificação agregada a ", AGG_FACTOR*10, " arc-min"),
      side = 1, line = 1.2, cex = 0.55, col = "gray40")


## ── Panel 2: legend ───────────────────────────────────────────
par(mar = c(2.5, 0, 3.5, 1.5), bg = "#D0E8F5")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Header
text(0.5, 0.96, "Macrobioclima", cex = 0.75, font = 2,
     col = "gray10", adj = 0.5)

entries <- data.frame(
  code  = 1:5,
  label = c("Tropical", "Mediterrânico", "Temperado", "Boreal", "Polar"),
  col   = bioclima_cols,
  note  = c(
    "< 23° (equatorial) ou subtropical\nverão mais chuvoso + quente",
    "Secura estival: Ios2 < 2 e Iosc4 < 2",
    "Sem secura estival; extratropical",
    "Frio; T_ann ≤ 4–6 °C  (45–71°N/S)",
    "Muito frio; > 66° lat."
  ),
  stringsAsFactors = FALSE
)

n_e    <- nrow(entries)
y_top  <- 0.88
y_step <- 0.14

for (k in seq_len(n_e)) {
  yc <- y_top - (k - 1) * y_step
  # Colour swatch
  rect(0.04, yc - 0.04, 0.22, yc + 0.04,
       col    = entries$col[k],
       border = "gray30",
       lwd    = 0.7)
  # Label
  text(0.27, yc + 0.018, entries$label[k],
       adj = 0, cex = 0.72, font = 2, col = "gray10")
  # Note (wrapped into 2 lines via \n in the string)
  lines_note <- strsplit(entries$note[k], "\n")[[1]]
  for (j in seq_along(lines_note))
    text(0.27, yc - 0.012 - (j - 1) * 0.024,
         lines_note[j],
         adj = 0, cex = 0.50, col = "gray35")
}

# Pixel-count sub-legend
y_bot <- y_top - n_e * y_step - 0.03
text(0.05, y_bot,
     sprintf("Pixels classificados\n(terra): %d", land_n),
     adj = 0, cex = 0.52, col = "gray45")

text(0.05, y_bot - 0.10,
     "Fonte: WorldClim v2.1\n(Fick & Hijmans 2017)",
     adj = 0, cex = 0.50, col = "gray50")

dev.off()
message("\nMap saved to: ", normalizePath(OUT_FILE))
