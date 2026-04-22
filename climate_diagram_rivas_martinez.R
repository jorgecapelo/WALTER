## ============================================================
##  Walter-Lieth Climate Diagram + Rivas-Martínez Bioclimatic
##  Classification  —  WorldClim v2.1 / OpenStreetMap edition
##
##  * Geocoding  : OpenStreetMap / Nominatim
##  * Climate    : WorldClim v2.1 global rasters via geodata
##  * Elevation  : WorldClim elevation_global()
##  * Resolution : WC_RES (default 2.5 arc-min, ~700 MB)
##  * Classification: Rivas-Martínez (2004), as described in
##    Mesquita, S. (2005) MSc thesis, pp. 37-45
##
##  Macrobioclimate logic:
##
##    MEDITERRANEAN vs TEMPERATE — two-index rule:
##      Ios2  = ombrothermic index of the 2 hottest summer months
##      Iosc4 = compensated index (3 summer months + preceding month)
##
##      Ios2 < 2  AND  Iosc4 < 2  →  Mediterranean  (summer drought)
##      Ios2 >= 2 AND  Iosc4 >= 2 →  Temperate      (no summer drought)
##      mixed case                 →  Walter dry-month tiebreaker
##                                    (>= 2 consecutive months P <= 2T)
##
##    TROPICAL (subtropical belt 23-35°) — summer rain:
##      Summer trimester wetter than winter trimester AND
##      warm-enough thermal thresholds (Table 3.9).
##      These two are ANTITHETICAL:
##        Tropical      → rain follows the sun  (summer wet)
##        Mediterranean → rain opposes the sun  (winter wet)
##
##    TEMPERATE  — extratropical, no summer drought
##    BOREAL     — cold temperate (thermal thresholds, Table 3.9)
##    POLAR      — abs(lat) > 66°
##
##  Usage:
##    source("climate_diagram_rivas_martinez.R")
##    climate_diagram("Seville, Spain")
##    climate_diagram("Nairobi, Kenya", save = TRUE)
## ============================================================


# ── 0. Packages ───────────────────────────────────────────────
required_pkgs <- c("httr", "jsonlite", "terra", "geodata", "climatol")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# ── 1. Configuration ──────────────────────────────────────────
WC_RES   <- 2.5
WC_CACHE <- file.path(getwd(), "worldclim_cache")
dir.create(WC_CACHE, recursive = TRUE, showWarnings = FALSE)
message("WorldClim global cache : ", WC_CACHE)
message("WorldClim resolution   : ", WC_RES, " arc-min")


# ── 2. In-memory raster cache ─────────────────────────────────
.wc_global <- new.env(parent = emptyenv())


# ── 3. Global raster loaders ──────────────────────────────────
get_wc_global <- function(var) {
  key <- paste0("wc_", var)
  if (!exists(key, envir = .wc_global)) {
    message("  Loading WorldClim global '", var,
            "' (res = ", WC_RES, " arc-min) ...")
    r <- geodata::worldclim_global(var = var, res = WC_RES, path = WC_CACHE)
    assign(key, r, envir = .wc_global)
    message("    Done — ", nlyr(r), " layers.")
  }
  get(key, envir = .wc_global)
}

get_elev_global <- function() {
  if (!exists("elev", envir = .wc_global)) {
    message("  Loading global elevation raster ...")
    r <- geodata::elevation_global(res = WC_RES, path = WC_CACHE)
    assign("elev", r, envir = .wc_global)
    message("    Done.")
  }
  get("elev", envir = .wc_global)
}


# ── 4. Extraction helpers ─────────────────────────────────────
extract_monthly <- function(r, lon, lat) {
  pt <- terra::vect(matrix(c(lon, lat), ncol = 2), crs = "EPSG:4326")
  as.numeric(terra::extract(r, pt)[1, -1])
}
extract_value <- function(r, lon, lat) {
  pt <- terra::vect(matrix(c(lon, lat), ncol = 2), crs = "EPSG:4326")
  as.numeric(terra::extract(r, pt)[1, 2])
}


# ── 5. Geocoding ──────────────────────────────────────────────
geocode_osm <- function(place_name) {
  resp <- httr::GET(
    "https://nominatim.openstreetmap.org/search",
    query = list(q = place_name, format = "json", limit = 1, addressdetails = 1),
    httr::user_agent("WalterLeithDiagram/1.0 (R educational script)")
  )
  httr::stop_for_status(resp)
  res <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
  if (length(res) == 0 || nrow(res) == 0)
    stop("Location not found on OpenStreetMap: '", place_name, "'")
  list(lon     = as.numeric(res$lon[1]),
       lat     = as.numeric(res$lat[1]),
       display = res$display_name[1])
}


# ── 6. Optional pre-load ──────────────────────────────────────
load_worldclim_global <- function() {
  message("\n--- Pre-loading all WorldClim global rasters ---")
  get_wc_global("prec"); get_wc_global("tmax")
  get_wc_global("tmin"); get_wc_global("tavg")
  get_elev_global()
  message("--- All rasters ready ---\n")
  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════
#  RIVAS-MARTÍNEZ (2004) BIOCLIMATIC CLASSIFICATION
#  Source: Mesquita, S. (2005) MSc thesis, pp. 37-45
#
#  Input matrix rows (same as diagwl except row 2 = Tmean):
#    row 1 : monthly precipitation P    (mm)
#    row 2 : monthly mean temp     T    (°C)   ← tavg, NOT tmax
#    row 3 : monthly mean min temp Tmin (°C)
#    row 4 : monthly mean max temp Tmax (°C)
# ══════════════════════════════════════════════════════════════

# ── RM-1.  Itc correction C (Table 3.8) ───────────────────────
.rm_itc_correction <- function(Ic) {
  if      (Ic <=  8) 10  * (Ic -  8)        # oceanic: negative correction
  else if (Ic <= 18) 0                        # euoceanic: no correction needed
  else if (Ic <= 21) 5   * (Ic - 18)
  else if (Ic <= 28) 15  + 15 * (Ic - 21)
  else if (Ic <= 46) 120 + 25 * (Ic - 28)
  else               570 + 30 * (min(Ic, 65) - 46)
}

# ── RM-2a.  Ios2 — ombrothermic index of the 2 hottest summer months ──
#
#  Tps2 / Pps2: temperature and precipitation of the 2 hottest months
#  of the summer trimester (restricted to months where T > 0).
#
.rm_ios2 <- function(T, P, sum_m) {
  Ts   <- T[sum_m]; Ps <- P[sum_m]
  h2   <- order(Ts, decreasing = TRUE)[1:2]   # indices of the 2 hottest
  pos  <- Ts[h2] > 0
  Tps2 <- sum(Ts[h2][pos]) * 10
  Pps2 <- sum(Ps[h2][pos])
  if (Tps2 > 0) (Pps2 / Tps2) * 10 else Inf
}

# ── RM-2b.  Iosc4 — compensated summer ombrothermic index ─────────────
#
#  Covers the 3 summer months PLUS the immediately preceding month,
#  so that rainfall carried over from late spring can compensate for
#  early-summer dryness.
#
.rm_iosc4 <- function(T, P, sum_m, prev_m) {
  four_m <- c(prev_m, sum_m)
  T4     <- T[four_m]; P4 <- P[four_m]
  pos    <- T4 > 0
  Tps4   <- sum(T4[pos]) * 10
  Pps4   <- sum(P4[pos])
  if (Tps4 > 0) (Pps4 / Tps4) * 10 else Inf
}

# ── RM-2c.  Is Mediterranean? — primary two-index rule ────────────────
#
#  Ios2 < 2  AND  Iosc4 < 2  →  TRUE   (summer drought confirmed)
#  Ios2 >= 2 AND  Iosc4 >= 2 →  FALSE  (no summer drought)
#  Mixed case                 →  tiebreaker: Walter dry-month criterion
#                                 (>= 2 consecutive months with P <= 2T
#                                  in the 5-month warm-season window)
#
.rm_is_mediterranean <- function(Ios2, Iosc4, T, P, sum_m) {
  if (Ios2 <  2 && Iosc4 <  2) return(TRUE)
  if (Ios2 >= 2 && Iosc4 >= 2) return(FALSE)

  # Mixed case: fall back to Walter's consecutive dry-month criterion
  prev_m <- (sum_m[1] - 2) %% 12 + 1
  next_m <-  sum_m[3]       %% 12 + 1
  window <- c(prev_m, sum_m, next_m)
  dry    <- P[window] <= 2 * T[window]
  runs   <- rle(dry)
  consec <- runs$lengths[runs$values == TRUE]
  length(consec) > 0 && max(consec) >= 2
}

# ── RM-3.  Summer-rain pattern — Tropical criterion ────────────────────
#
#  Mediterranean and Tropical are ANTITHETICAL:
#    Tropical      → rain follows the sun → summer trimester is WET
#    Mediterranean → rain opposes the sun → winter trimester is WET
#
#  Criterion: summer trimester total > winter trimester total
#
.rm_summer_rain <- function(Pcm1, Ps_w) Pcm1 > Ps_w

# ── RM-4.  Boreal conditions (Table 3.9) ──────────────────────
.rm_check_boreal <- function(Ic, T_annual, Tp) {
  if      (Ic <= 11) T_annual <= 6.0 && Tp > 0
  else if (Ic <= 21) T_annual <= 5.3 && Tp > 380 && Tp <= 720
  else if (Ic <= 28) T_annual <= 4.8 && Tp > 380 && Tp <= 740
  else if (Ic <  46) T_annual <= 3.8 && Tp > 380 && Tp <= 800
  else               T_annual <= 0.0 && Tp > 380 && Tp <= 800
}

# ── RM-5.  Macrobioclima ──────────────────────────────────────
#
#  Decision tree (priority order top → bottom):
#
#  abs_lat > 66°          →  Polar
#  abs_lat 52–66°         →  Boreal (thermal thresholds) | Temperate
#  abs_lat 35–52°         →  Mediterranean (summer drought) | Temperate
#  abs_lat 23–35°  ──────────── subtropical three-way split:
#    summer drought present        → Mediterranean
#    summer rain + warm enough     → Tropical
#    otherwise                     → Temperate
#  abs_lat <= 23°         →  Tropical (always, equatorial belt)
#
.rm_macrobioclima <- function(T, P, lat, Ic, T_annual,
                               m_val, M_val, Itc, Tp,
                               Pcm1, Ps_w, sum_m, prev_m) {
  abs_lat <- abs(lat)

  Ios2   <- .rm_ios2(T, P, sum_m)
  Iosc4  <- .rm_iosc4(T, P, sum_m, prev_m)
  s_med  <- .rm_is_mediterranean(Ios2, Iosc4, T, P, sum_m)
  s_rain <- .rm_summer_rain(Pcm1, Ps_w)

  if (abs_lat > 66) return("Polar")

  if (abs_lat > 52) {
    if (.rm_check_boreal(Ic, T_annual, Tp)) return("Boreal")
    return("Temperado")
  }

  if (abs_lat > 35) {
    if (s_med) return("Mediterrânico")
    return("Temperado")
  }

  if (abs_lat > 23) {
    # Summer drought takes priority — it is the defining Med criterion
    if (s_med) return("Mediterrânico")

    # Thermal thresholds for subtropical Tropical (Table 3.9)
    warm_a <- sum(c(T_annual >= 21, M_val >= 18, Itc >= 470)) >= 2
    warm_b <- sum(c(T_annual >= 25, m_val >= 10, Itc >= 580)) >= 2

    if ((warm_a || warm_b) && s_rain) return("Tropical")
    return("Temperado")
  }

  "Tropical"
}

# ── RM-6.  Bioclima (Table 3.9) ───────────────────────────────
.rm_bioclima <- function(macro, Io, Iod2, Ic) {
  switch(macro,
    "Tropical" = {
      if      (Io >= 3.6) if (Iod2 > 2.5) "Pluvial" else "Pluviestacional"
      else if (Io >= 1.0) "Xérico"
      else if (Io >= 0.1) "Desértico"
      else                "Hiperdesértico"
    },
    "Mediterrânico" = {
      if      (Io >  2.0 && Ic <= 21) "Pluviestacional Oceânico"
      else if (Io >  2.2 && Ic >  21) "Pluviestacional Continental"
      else if (Io >= 1.0 && Ic <= 21) "Xérico Oceânico"
      else if (Io >= 1.0 && Ic >  21) "Xérico Continental"
      else if (Io >= 0.1 && Ic <= 21) "Desértico Oceânico"
      else if (Io >= 0.1 && Ic >  21) "Desértico Continental"
      else                            "Hiperdesértico"
    },
    "Temperado" = {
      if      (Io > 3.6 && Ic <= 11) "Hiperoceânico"
      else if (Io > 3.6 && Ic <= 21) "Oceânico"
      else if (Io > 3.6)             "Continental"
      else                           "Xérico"
    },
    "Boreal" = {
      if      (Io <= 3.6) "Xérico"
      else if (Ic <= 11)  "Hiperoceânico"
      else if (Ic <= 21)  "Oceânico"
      else if (Ic <= 28)  "Subcontinental"
      else if (Ic <= 46)  "Continental"
      else                "Hipercontinental"
    },
    "Polar" = {
      if      (Ic <= 11) "Hiperoceânico"
      else if (Ic <= 21) "Oceânico"
      else               "Continental"
    },
    NA_character_
  )
}

# ── RM-7.  Termotipo (Table 3.12) ─────────────────────────────
#  Note (Table 3.12 footnote): use Tp instead of Itc when
#  Ic >= 21 (continental) or Itc < 120 (very cold sites).
.rm_termotipo <- function(macro, Itc, Tp, use_Tp) {

  if (macro == "Tropical") {
    if      (Tp > 3350) "Infratropical inferior"
    else if (Tp > 3100) "Infratropical superior"
    else if (Tp > 2900) "Termotropical inferior"
    else if (Tp > 2700) "Termotropical superior"
    else if (Tp > 2400) "Mesotropical inferior"
    else if (Tp > 2100) "Mesotropical superior"
    else if (Tp > 1575) "Supratropical inferior"
    else if (Tp > 1050) "Supratropical superior"
    else if (Tp >  750) "Orotropical inferior"
    else if (Tp >  450) "Orotropical superior"
    else if (Tp >  150) "Criorotropical inferior"
    else if (Tp >    0) "Criorotropical superior"
    else                "Tropical gélido"

  } else if (macro == "Mediterrânico") {
    if (!use_Tp) {
      if      (Itc >= 515) "Inframediterrânico inferior"
      else if (Itc >= 450) "Inframediterrânico superior"
      else if (Itc >= 400) "Termomediterrânico inferior"
      else if (Itc >= 350) "Termomediterrânico superior"
      else if (Itc >= 280) "Mesomediterrânico inferior"
      else if (Itc >= 220) "Mesomediterrânico superior"
      else if (Itc >= 150) "Supramediterrânico inferior"
      else                 "Supramediterrânico superior"
    } else {
      if      (Tp > 2650) "Inframediterrânico inferior"
      else if (Tp > 2450) "Inframediterrânico superior"
      else if (Tp > 2300) "Termomediterrânico inferior"
      else if (Tp > 2150) "Termomediterrânico superior"
      else if (Tp > 1825) "Mesomediterrânico inferior"
      else if (Tp > 1500) "Mesomediterrânico superior"
      else if (Tp > 1200) "Supramediterrânico inferior"
      else if (Tp >  900) "Supramediterrânico superior"
      else if (Tp >  675) "Oromediterrânico inferior"
      else if (Tp >  450) "Oromediterrânico superior"
      else if (Tp >  130) "Crioromediterrânico inferior"
      else if (Tp >    0) "Crioromediterrânico superior"
      else                NA_character_
    }

  } else if (macro == "Temperado") {
    if (!use_Tp) {
      if      (Itc >  410) "Infratemperado"
      else if (Itc >= 350) "Termotemperado inferior"
      else if (Itc >= 290) "Termotemperado superior"
      else if (Itc >= 240) "Mesotemperado inferior"
      else if (Itc >= 190) "Mesotemperado superior"
      else if (Itc >= 120) "Supratemperado inferior"
      else                 "Supratemperado superior"
    } else {
      if      (Tp > 2350) "Infratemperado"
      else if (Tp > 2175) "Termotemperado inferior"
      else if (Tp > 2000) "Termotemperado superior"
      else if (Tp > 1700) "Mesotemperado inferior"
      else if (Tp > 1400) "Mesotemperado superior"
      else if (Tp > 1100) "Supratemperado inferior"
      else if (Tp >  800) "Supratemperado superior"
      else if (Tp >  590) "Orotemperado inferior"
      else if (Tp >  380) "Orotemperado superior"
      else if (Tp >  130) "Criorotemperado inferior"
      else if (Tp >    0) "Criorotemperado superior"
      else                NA_character_
    }

  } else if (macro == "Boreal") {
    if      (Tp > 750) "Termoboreal inferior"
    else if (Tp > 700) "Termoboreal superior"
    else if (Tp > 600) "Mesoboreal inferior"
    else if (Tp > 500) "Mesoboreal superior"
    else if (Tp > 440) "Supraboreal inferior"
    else if (Tp > 380) "Supraboreal superior"
    else if (Tp > 230) "Oroboreal inferior"
    else if (Tp >  80) "Oroboreal superior"
    else if (Tp >  40) "Crioroboreal inferior"
    else if (Tp >   0) "Crioroboreal superior"
    else               "Boreal gélido"

  } else if (macro == "Polar") {
    if      (Tp > 230) "Mesopolar inferior"
    else if (Tp >  80) "Mesopolar superior"
    else if (Tp >  40) "Suprapolar inferior"
    else if (Tp >   0) "Suprapolar superior"
    else               "Polar gélido"
  } else NA_character_
}

# ── RM-8.  Ombrotipo (Table 3.13) ─────────────────────────────
.rm_ombrotipo <- function(Io) {
  if      (Io <  0.1) "Ultra-hiperárido"
  else if (Io <  0.2) "Hiperárido inferior"
  else if (Io <  0.3) "Hiperárido superior"
  else if (Io <  0.6) "Árido inferior"
  else if (Io <  1.0) "Árido superior"
  else if (Io <  1.5) "Semiárido inferior"
  else if (Io <  2.0) "Semiárido superior"
  else if (Io <  2.8) "Seco inferior"
  else if (Io <  3.6) "Seco superior"
  else if (Io <  4.8) "Sub-húmido inferior"
  else if (Io <  6.0) "Sub-húmido superior"
  else if (Io <  9.0) "Húmido inferior"
  else if (Io < 12.0) "Húmido superior"
  else if (Io < 18.0) "Hiper-húmido inferior"
  else if (Io < 24.0) "Hiper-húmido superior"
  else                "Ultra-hiper-húmido"
}

# ── RM-9.  Continentality type (Table 3.14) ───────────────────
.rm_continentality <- function(Ic) {
  if      (Ic <=  4) "Hiperoceânico · Ultra-hiperoceânico"
  else if (Ic <=  8) "Hiperoceânico · Eu-hiperoceânico"
  else if (Ic <= 11) "Hiperoceânico · Sub-hiperoceânico"
  else if (Ic <= 14) "Oceânico · Semi-hiperoceânico"
  else if (Ic <= 17) "Oceânico · Euoceânico"
  else if (Ic <= 21) "Oceânico · Semicontinental"
  else if (Ic <= 28) "Continental · Subcontinental"
  else if (Ic <= 46) "Continental · Eucontinental"
  else               "Continental · Hipercontinental"
}

# ── RM-10.  Bioclimatic variants (Table 3.11) ─────────────────
.rm_variant <- function(macro, Ic, Io, Pcm1, Ps_w, T, P, sum_m, win_m) {
  v <- character(0)

  # Estépica: high continentality, summer drier than average,
  #           and at least one summer month with P < 3T
  if (Ic > 17 && Pcm1 > 1.1 * Ps_w && Io > 0.1 && Io < 4.8 &&
      any(P[sum_m] < 3 * T[sum_m]))
    v <- c(v, "Estépica")

  # Submediterrânica: Temperate only; at least one summer month P < 2.8T
  if (macro == "Temperado" && any(P[sum_m] < 2.8 * T[sum_m]))
    v <- c(v, "Submediterrânica")

  # Antitropical: Tropical (not Pluvial/Hiperdesértico), winter > summer rain
  if (macro == "Tropical" && Io < 3.6 && Io > 0.1 && Ps_w > Pcm1)
    v <- c(v, "Antitropical")

  if (length(v) == 0) "—" else paste(v, collapse = " + ")
}

# ── RM-11.  Master diagnosis function ─────────────────────────
#'
#' Compute and print the Rivas-Martínez bioclimatic classification.
#'
#' @param dat  4×12 matrix: rows = P (mm), Tmean (°C), Tmin (°C), Tmax (°C)
#' @param lat  Decimal latitude (+N / −S)
#' @param name Station name string (optional)
#'
rivas_martinez_diagnosis <- function(dat, lat, name = "") {

  P    <- as.numeric(dat[1, ])
  T    <- as.numeric(dat[2, ])
  Tmin <- as.numeric(dat[3, ])
  Tmax <- as.numeric(dat[4, ])

  # Hemisphere & seasonal month indices (1-based)
  NH <- lat >= 0
  if (NH) {
    sum_m  <- 6:8;         prev_m <- 5
    fall_m <- 9:11;        spr_m  <- 3:5
    win_m  <- c(12, 1, 2)
  } else {
    sum_m  <- c(12, 1, 2); prev_m <- 11
    fall_m <- 3:5;         spr_m  <- 9:11
    win_m  <- 6:8
  }

  # Basic thermal parameters
  T_annual <- mean(T)
  cold_m   <- which.min(T)
  hot_m    <- which.max(T)
  m_val    <- Tmin[cold_m]   # mean minimum of coldest month
  M_val    <- Tmax[cold_m]   # mean maximum of coldest month

  # ── Indices (Table 3.8) ─────────────────────────────────────
  Ic   <- T[hot_m] - T[cold_m]
  It   <- (T_annual + m_val + M_val) * 10
  C    <- .rm_itc_correction(Ic)
  Itc  <- It + C
  Tp   <- sum(T[T > 0]) * 10       # annual positive temperature (×10)
  Pp   <- sum(P[T > 0])             # positive precipitation (mm)
  Io   <- if (Tp > 0) (Pp / Tp) * 10 else 0

  # Ombrothermic index of the 2 driest months in the driest trimester
  all_tri <- list(c(12,1,2), c(3,4,5), c(6,7,8), c(9,10,11))
  tri_P   <- sapply(all_tri, function(m) sum(P[m]))
  dry_tri <- all_tri[[which.min(tri_P)]]
  d2      <- dry_tri[order(P[dry_tri])[1:2]]
  pos_d   <- T[d2] > 0
  Tpd2    <- sum(T[d2][pos_d]) * 10
  Ppd2    <- sum(P[d2][pos_d])
  Iod2    <- if (Tpd2 > 0) (Ppd2 / Tpd2) * 10 else Inf

  # Summer ombrothermic indices (Med/Temp boundary)
  Ios2  <- .rm_ios2(T, P, sum_m)
  Iosc4 <- .rm_iosc4(T, P, sum_m, prev_m)

  # Seasonal precipitation totals
  Pcm1 <- sum(P[sum_m])    # summer
  Pcm2 <- sum(P[fall_m])   # trimester following summer
  Pcm3 <- sum(P[spr_m])    # trimester preceding summer
  Ps_w <- sum(P[win_m])    # winter

  # ── Rain-seasonality diagnostics ────────────────────────────
  s_med  <- .rm_is_mediterranean(Ios2, Iosc4, T, P, sum_m)
  s_rain <- .rm_summer_rain(Pcm1, Ps_w)

  # ── Macrobioclima ───────────────────────────────────────────
  macro <- .rm_macrobioclima(T, P, lat, Ic, T_annual,
                              m_val, M_val, Itc, Tp,
                              Pcm1, Ps_w, sum_m, prev_m)

  # ── Sub-classifications ─────────────────────────────────────
  bio   <- .rm_bioclima(macro, Io, Iod2, Ic)
  termo <- .rm_termotipo(macro, Itc, Tp, use_Tp = (Ic >= 21 || Itc < 120))
  ombro <- .rm_ombrotipo(Io)
  cont  <- .rm_continentality(Ic)
  var   <- .rm_variant(macro, Ic, Io, Pcm1, Ps_w, T, P, sum_m, win_m)

  # ── Print ───────────────────────────────────────────────────
  rule  <- paste0(rep("─", 64), collapse = "")
  drule <- paste0(rep("═", 64), collapse = "")
  lbl   <- function(label, value)
    cat(sprintf("  %-32s %s\n", paste0(label, " :"), value))

  cat("\n", drule, "\n", sep = "")
  cat("  CLASSIFICAÇÃO BIOCLIMÁTICA DE RIVAS-MARTÍNEZ (2004)\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Estação: %s\n", name))
  cat(drule, "\n\n", sep = "")

  cat("  Índices bioclimáticos\n", rule, "\n", sep = "")
  cat(sprintf("  %-6s Continentalidade            %8.2f\n",  "Ic",   Ic))
  cat(sprintf("  %-6s Termicidade                 %8.1f\n",  "It",   It))
  cat(sprintf("  %-6s Correcção de It             %8.1f\n",  "C",    C))
  cat(sprintf("  %-6s Termicidade Compensado      %8.1f\n",  "Itc",  Itc))
  cat(sprintf("  %-6s Temp. Positiva Anual (×10)  %8.0f\n",  "Tp",   Tp))
  cat(sprintf("  %-6s Precipitação Positiva (mm)  %8.1f\n",  "Pp",   Pp))
  cat(sprintf("  %-6s Ombrotérmico Anual          %8.2f\n",  "Io",    Io))
  cat(sprintf("  %-6s Ombrotérmico bimestre seco  %8.2f\n",  "Iod2",  Iod2))
  cat(sprintf("  %-6s Ombrotérmico bimestre verão %8.2f\n",  "Ios2",  Ios2))
  cat(sprintf("  %-6s Ombrotérmico compensado     %8.2f\n",  "Iosc4", Iosc4))
  cat(sprintf("  %-6s Precip. trimestre verão (mm)%8.1f\n",  "Pcm1",  Pcm1))
  cat(sprintf("  %-6s Precip. trim. inverno (mm)  %8.1f\n",  "Ps_w",  Ps_w))
  cat(sprintf("  %-6s Precip. trim. seguinte (mm) %8.1f\n",  "Pcm2",  Pcm2))
  cat(sprintf("  %-6s Precip. trim. anterior (mm) %8.1f\n",  "Pcm3",  Pcm3))

  # Med/Temp boundary: makes the two-index rule transparent
  med_rule <- if (Ios2 < 2 && Iosc4 < 2) {
    "Mediterrânico  (Ios2<2 e Iosc4<2)"
  } else if (Ios2 >= 2 && Iosc4 >= 2) {
    "Temperado      (Ios2>=2 e Iosc4>=2)"
  } else {
    paste0("Misto — critério Walter: ", if (s_med) "Med (>=2 meses P<=2T)" else "Temperado")
  }

  cat("\n  Fronteira Mediterrânico / Temperado\n", rule, "\n", sep = "")
  cat(sprintf("  Ios2  = %5.2f  |  Iosc4 = %5.2f  →  %s\n", Ios2, Iosc4, med_rule))
  cat(sprintf("  Chuva estival (Pcm_verao>Pcm_inv) : %s\n",
              if (s_rain) "SIM  → critério Tropical"
              else        "NAO  → inverno mais chuvoso"))

  cat("\n  Classificação\n", rule, "\n", sep = "")
  lbl("Macrobioclima",    macro)
  lbl("Bioclima",         bio)
  lbl("Termotipo",        termo)
  lbl("Ombrotipo",        ombro)
  lbl("Continentalidade", cont)
  lbl("Variante",         var)
  cat(drule, "\n\n", sep = "")

  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════
#  MAIN FUNCTION: climate_diagram()
# ══════════════════════════════════════════════════════════════

#' Draw a Walter-Lieth diagram and print the Rivas-Martínez
#' bioclimatic classification for any named location.
#'
#' @param place  Character: place name (add country for best results).
#' @param save   Logical: also save a PNG? Default FALSE.
#' @param outdir Directory for saved PNGs (default = working directory).
#'
climate_diagram <- function(place, save = FALSE, outdir = ".") {

  ## ── Geocode ───────────────────────────────────────────────
  message("\n=== Processing: ", place, " ===")
  geo <- geocode_osm(place)
  lon <- geo$lon; lat <- geo$lat
  message("  OSM match : ", geo$display)

  ## ── Load rasters (cached after first call) ────────────────
  r_prec <- get_wc_global("prec"); r_tmax <- get_wc_global("tmax")
  r_tmin <- get_wc_global("tmin"); r_tavg <- get_wc_global("tavg")
  r_elev <- get_elev_global()

  ## ── Elevation ─────────────────────────────────────────────
  elev <- extract_value(r_elev, lon, lat)
  message(sprintf("  Lon: %.4f  Lat: %.4f  Elev: %s m",
                  lon, lat,
                  ifelse(is.na(elev), "N/A", round(elev))))

  ## ── Monthly climate series ────────────────────────────────
  prec <- extract_monthly(r_prec, lon, lat)
  tmax <- extract_monthly(r_tmax, lon, lat)
  tmin <- extract_monthly(r_tmin, lon, lat)
  tavg <- extract_monthly(r_tavg, lon, lat)

  # Absolute minimum approximation (WorldClim has no abs-min layer)
  abs_min <- tmin - pmax((tavg - tmin) / 2, 0)

  southern <- lat < 0
  if (southern) message("  Southern hemisphere detected")

  ## ── Walter-Lieth matrix (rows: P / Tmax / Tmin / Tabs) ───
  clim_mat <- rbind(prec, tmax, tmin, abs_min)

  ## ── Rivas-Martínez matrix (rows: P / Tmean / Tmin / Tmax)
  ##    Row 2 MUST be Tmean (tavg) for It / Tp / Io to be correct
  rm_mat <- rbind(prec, tavg, tmin, tmax)

  ## ── Station label ─────────────────────────────────────────
  station_label <- trimws(sub(",.*", "", geo$display))

  ## ── Walter-Lieth diagram ──────────────────────────────────
  if (save) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    safe  <- gsub("[^[:alnum:]_-]", "_", place)
    fpath <- file.path(outdir, paste0("climate_", safe, ".png"))
    png(fpath, width = 800, height = 600, res = 96)
    message("  Saving PNG to: ", fpath)
  }

  climatol::diagwl(
    dat    = clim_mat,
    cols   = NULL,
    stname = station_label,
    alt    = ifelse(is.na(elev), NA, round(elev)),
    per    = "1970-2000",
    mlab   = "en",
    shem   = southern,
    pcol   = "#4A90D9",
    tcol   = "#D94A4A",
    pfcol  = "#A8D4F5",
    sfcol  = "#FFD0A0",
    p3line = FALSE
  )

  if (save) { dev.off(); message("  Done.") }

  ## ── Rivas-Martínez bioclimatic diagnosis ──────────────────
  rivas_martinez_diagnosis(rm_mat, lat = lat, name = station_label)

  invisible(list(
    place     = place,
    display   = geo$display,
    lon       = lon,
    lat       = lat,
    elevation = elev,
    southern  = southern,
    clim_mat  = clim_mat,
    rm_mat    = rm_mat
  ))
}


# ── Ready message ─────────────────────────────────────────────
message("
==========================================================
 Walter-Lieth + Rivas-Martínez script loaded

 Cache dir  : ", WC_CACHE, "
 Resolution : ", WC_RES, " arc-min  (change WC_RES before sourcing)
 Data source: WorldClim v2.1 global rasters

 The FIRST call downloads ~700 MB (at 2.5 arc-min).
 Every subsequent call is instant — rasters stay in memory.

 To pre-load everything before your first diagram:
   load_worldclim_global()

 Usage:
   climate_diagram(\"City, Country\")
   climate_diagram(\"City, Country\", save = TRUE)

 Examples:
   climate_diagram(\"Seville, Spain\")          # Mediterranean
   climate_diagram(\"Nairobi, Kenya\")           # Tropical
   climate_diagram(\"Paris, France\")            # Temperate
   climate_diagram(\"Buenos Aires, Argentina\")  # SH Temperate
   climate_diagram(\"Barrow, Alaska, USA\")      # Polar/Boreal
   climate_diagram(\"Alice Springs, Australia\", save = TRUE)

 Batch (rasters loaded only once):
   load_worldclim_global()
   lapply(c(\"Cairo, Egypt\", \"Lima, Peru\", \"Oslo, Norway\"),
          climate_diagram)

 Resolution options (set WC_RES before sourcing):
   WC_RES <- 10   # ~45 MB  — fastest
   WC_RES <- 5    # ~175 MB
   WC_RES <- 2.5  # ~700 MB — default
   WC_RES <- 0.5  # ~11 GB  — highest detail
==========================================================")
