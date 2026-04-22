## ============================================================
##  Walter-Lieth Climate Diagram + Rivas-Martínez Bioclimatic
##  Classification  —  WorldClim v2.1 / OpenStreetMap edition
##
##  * Geocoding      : OpenStreetMap / Nominatim
##  * Climate        : WorldClim v2.1 global rasters via geodata
##  * Elevation      : WorldClim elevation_global()
##  * Resolution     : WC_RES (default 2.5 arc-min, ~700 MB)
##  * Classification : Rivas-Martínez (2004) — Mesquita (2005) MSc
##                     Köppen-Geiger     — Peel et al. (2007)
##                     Walter Zonobiome  — Walter (1985)
##                     Mucina (2019)     — global biome level
##
##  All outputs in English.
##
##  Usage:
##    source("climate_diagram_full.R")
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
#  Input matrix rows:
#    row 1 : monthly precipitation P    (mm)
#    row 2 : monthly mean temp     T    (°C)   ← tavg
#    row 3 : monthly mean min temp Tmin (°C)
#    row 4 : monthly mean max temp Tmax (°C)
# ══════════════════════════════════════════════════════════════

# ── RM-1.  Itc correction C (Table 3.8) ───────────────────────
.rm_itc_correction <- function(Ic) {
  if      (Ic <=  8) 10  * (Ic -  8)
  else if (Ic <= 18) 0
  else if (Ic <= 21) 5   * (Ic - 18)
  else if (Ic <= 28) 15  + 15 * (Ic - 21)
  else if (Ic <= 46) 120 + 25 * (Ic - 28)
  else               570 + 30 * (min(Ic, 65) - 46)
}

# ── RM-2a.  Ios2 ──────────────────────────────────────────────
.rm_ios2 <- function(T, P, sum_m) {
  Ts   <- T[sum_m]; Ps <- P[sum_m]
  h2   <- order(Ts, decreasing = TRUE)[1:2]
  pos  <- Ts[h2] > 0
  Tps2 <- sum(Ts[h2][pos]) * 10
  Pps2 <- sum(Ps[h2][pos])
  if (Tps2 > 0) (Pps2 / Tps2) * 10 else Inf
}

# ── RM-2b.  Iosc4 ─────────────────────────────────────────────
.rm_iosc4 <- function(T, P, sum_m, prev_m) {
  four_m <- c(prev_m, sum_m)
  T4     <- T[four_m]; P4 <- P[four_m]
  pos    <- T4 > 0
  Tps4   <- sum(T4[pos]) * 10
  Pps4   <- sum(P4[pos])
  if (Tps4 > 0) (Pps4 / Tps4) * 10 else Inf
}

# ── RM-2c.  Is Mediterranean? ─────────────────────────────────
.rm_is_mediterranean <- function(Ios2, Iosc4, T, P, sum_m) {
  if (Ios2 <  2 && Iosc4 <  2) return(TRUE)
  if (Ios2 >= 2 && Iosc4 >= 2) return(FALSE)
  prev_m <- (sum_m[1] - 2) %% 12 + 1
  next_m <-  sum_m[3]       %% 12 + 1
  window <- c(prev_m, sum_m, next_m)
  dry    <- P[window] <= 2 * T[window]
  runs   <- rle(dry)
  consec <- runs$lengths[runs$values == TRUE]
  length(consec) > 0 && max(consec) >= 2
}

# ── RM-3.  Summer-rain pattern ────────────────────────────────
.rm_summer_rain <- function(Pcm1, Ps_w) Pcm1 > Ps_w

# ── RM-4.  Boreal conditions (Table 3.9) ──────────────────────
.rm_check_boreal <- function(Ic, T_annual, Tp) {
  if      (Ic <= 11) T_annual <= 6.0 && Tp > 0
  else if (Ic <= 21) T_annual <= 5.3 && Tp > 380 && Tp <= 720
  else if (Ic <= 28) T_annual <= 4.8 && Tp > 380 && Tp <= 740
  else if (Ic <  46) T_annual <= 3.8 && Tp > 380 && Tp <= 800
  else               T_annual <= 0.0 && Tp > 380 && Tp <= 800
}

# ── RM-5.  Macrobioclimate ────────────────────────────────────
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
    return("Temperate")
  }
  if (abs_lat > 35) {
    if (s_med) return("Mediterranean")
    return("Temperate")
  }
  if (abs_lat > 23) {
    if (s_med) return("Mediterranean")
    warm_a <- sum(c(T_annual >= 21, M_val >= 18, Itc >= 470)) >= 2
    warm_b <- sum(c(T_annual >= 25, m_val >= 10, Itc >= 580)) >= 2
    if ((warm_a || warm_b) && s_rain) return("Tropical")
    return("Temperate")
  }
  "Tropical"
}

# ── RM-6.  Bioclimate (Table 3.9) ─────────────────────────────
.rm_bioclima <- function(macro, Io, Iod2, Ic) {
  switch(macro,
    "Tropical" = {
      if      (Io >= 3.6) if (Iod2 > 2.5) "Pluvial" else "Pluvioseasonal"
      else if (Io >= 1.0) "Xeric"
      else if (Io >= 0.1) "Desertic"
      else                "Hyperdesertic"
    },
    "Mediterranean" = {
      if      (Io >  2.0 && Ic <= 21) "Pluvioseasonal Oceanic"
      else if (Io >  2.2 && Ic >  21) "Pluvioseasonal Continental"
      else if (Io >= 1.0 && Ic <= 21) "Xeric Oceanic"
      else if (Io >= 1.0 && Ic >  21) "Xeric Continental"
      else if (Io >= 0.1 && Ic <= 21) "Desertic Oceanic"
      else if (Io >= 0.1 && Ic >  21) "Desertic Continental"
      else                            "Hyperdesertic"
    },
    "Temperate" = {
      if      (Io > 3.6 && Ic <= 11) "Hyperoceanic"
      else if (Io > 3.6 && Ic <= 21) "Oceanic"
      else if (Io > 3.6)             "Continental"
      else                           "Xeric"
    },
    "Boreal" = {
      if      (Io <= 3.6) "Xeric"
      else if (Ic <= 11)  "Hyperoceanic"
      else if (Ic <= 21)  "Oceanic"
      else if (Ic <= 28)  "Subcontinental"
      else if (Ic <= 46)  "Continental"
      else                "Hypercontinental"
    },
    "Polar" = {
      if      (Ic <= 11) "Hyperoceanic"
      else if (Ic <= 21) "Oceanic"
      else               "Continental"
    },
    NA_character_
  )
}

# ── RM-7.  Thermotype (Table 3.12) ────────────────────────────
.rm_termotipo <- function(macro, Itc, Tp, use_Tp) {

  if (macro == "Tropical") {
    if      (Tp > 3350) "Infratropical lower"
    else if (Tp > 3100) "Infratropical upper"
    else if (Tp > 2900) "Thermotropical lower"
    else if (Tp > 2700) "Thermotropical upper"
    else if (Tp > 2400) "Mesotropical lower"
    else if (Tp > 2100) "Mesotropical upper"
    else if (Tp > 1575) "Supratropical lower"
    else if (Tp > 1050) "Supratropical upper"
    else if (Tp >  750) "Orotropical lower"
    else if (Tp >  450) "Orotropical upper"
    else if (Tp >  150) "Cryorotropical lower"
    else if (Tp >    0) "Cryorotropical upper"
    else                "Gelidtropical"

  } else if (macro == "Mediterranean") {
    if (!use_Tp) {
      if      (Itc >= 515) "Inframediterranean lower"
      else if (Itc >= 450) "Inframediterranean upper"
      else if (Itc >= 400) "Thermomediterranean lower"
      else if (Itc >= 350) "Thermomediterranean upper"
      else if (Itc >= 280) "Mesomediterranean lower"
      else if (Itc >= 220) "Mesomediterranean upper"
      else if (Itc >= 150) "Supramediterranean lower"
      else                 "Supramediterranean upper"
    } else {
      if      (Tp > 2650) "Inframediterranean lower"
      else if (Tp > 2450) "Inframediterranean upper"
      else if (Tp > 2300) "Thermomediterranean lower"
      else if (Tp > 2150) "Thermomediterranean upper"
      else if (Tp > 1825) "Mesomediterranean lower"
      else if (Tp > 1500) "Mesomediterranean upper"
      else if (Tp > 1200) "Supramediterranean lower"
      else if (Tp >  900) "Supramediterranean upper"
      else if (Tp >  675) "Oromediterranean lower"
      else if (Tp >  450) "Oromediterranean upper"
      else if (Tp >  130) "Cryoromediterranean lower"
      else if (Tp >    0) "Cryoromediterranean upper"
      else                NA_character_
    }

  } else if (macro == "Temperate") {
    if (!use_Tp) {
      if      (Itc >  410) "Infratemperate"
      else if (Itc >= 350) "Thermotemperate lower"
      else if (Itc >= 290) "Thermotemperate upper"
      else if (Itc >= 240) "Mesotemperate lower"
      else if (Itc >= 190) "Mesotemperate upper"
      else if (Itc >= 120) "Supratemperate lower"
      else                 "Supratemperate upper"
    } else {
      if      (Tp > 2350) "Infratemperate"
      else if (Tp > 2175) "Thermotemperate lower"
      else if (Tp > 2000) "Thermotemperate upper"
      else if (Tp > 1700) "Mesotemperate lower"
      else if (Tp > 1400) "Mesotemperate upper"
      else if (Tp > 1100) "Supratemperate lower"
      else if (Tp >  800) "Supratemperate upper"
      else if (Tp >  590) "Orotemperate lower"
      else if (Tp >  380) "Orotemperate upper"
      else if (Tp >  130) "Cryorotemperate lower"
      else if (Tp >    0) "Cryorotemperate upper"
      else                NA_character_
    }

  } else if (macro == "Boreal") {
    if      (Tp > 750) "Thermoboreal lower"
    else if (Tp > 700) "Thermoboreal upper"
    else if (Tp > 600) "Mesoboreal lower"
    else if (Tp > 500) "Mesoboreal upper"
    else if (Tp > 440) "Supraboreal lower"
    else if (Tp > 380) "Supraboreal upper"
    else if (Tp > 230) "Oroboreal lower"
    else if (Tp >  80) "Oroboreal upper"
    else if (Tp >  40) "Cryoroboreal lower"
    else if (Tp >   0) "Cryoroboreal upper"
    else               "Gelidboreal"

  } else if (macro == "Polar") {
    if      (Tp > 230) "Mesopolar lower"
    else if (Tp >  80) "Mesopolar upper"
    else if (Tp >  40) "Suprapolar lower"
    else if (Tp >   0) "Suprapolar upper"
    else               "Gelidpolar"
  } else NA_character_
}

# ── RM-8.  Ombrotype (Table 3.13) ─────────────────────────────
.rm_ombrotipo <- function(Io) {
  if      (Io <  0.1) "Ultra-hyperarid"
  else if (Io <  0.2) "Hyperarid lower"
  else if (Io <  0.3) "Hyperarid upper"
  else if (Io <  0.6) "Arid lower"
  else if (Io <  1.0) "Arid upper"
  else if (Io <  1.5) "Semiarid lower"
  else if (Io <  2.0) "Semiarid upper"
  else if (Io <  2.8) "Dry lower"
  else if (Io <  3.6) "Dry upper"
  else if (Io <  4.8) "Subhumid lower"
  else if (Io <  6.0) "Subhumid upper"
  else if (Io <  9.0) "Humid lower"
  else if (Io < 12.0) "Humid upper"
  else if (Io < 18.0) "Hyperhumid lower"
  else if (Io < 24.0) "Hyperhumid upper"
  else                "Ultra-hyperhumid"
}

# ── RM-9.  Continentality type (Table 3.14) ───────────────────
.rm_continentality <- function(Ic) {
  if      (Ic <=  4) "Hyperoceanic – Ultra-hyperoceanic"
  else if (Ic <=  8) "Hyperoceanic – Eu-hyperoceanic"
  else if (Ic <= 11) "Hyperoceanic – Sub-hyperoceanic"
  else if (Ic <= 14) "Oceanic – Semi-hyperoceanic"
  else if (Ic <= 17) "Oceanic – Eu-oceanic"
  else if (Ic <= 21) "Oceanic – Semicontinental"
  else if (Ic <= 28) "Continental – Subcontinental"
  else if (Ic <= 46) "Continental – Eu-continental"
  else               "Continental – Hypercontinental"
}

# ── RM-10.  Bioclimatic variants (Table 3.11) ─────────────────
.rm_variant <- function(macro, Ic, Io, Pcm1, Ps_w, T, P, sum_m, win_m) {
  v <- character(0)

  if (Ic > 17 && Pcm1 > 1.1 * Ps_w && Io > 0.1 && Io < 4.8 &&
      any(P[sum_m] < 3 * T[sum_m]))
    v <- c(v, "Steppic")

  if (macro == "Temperate" && any(P[sum_m] < 2.8 * T[sum_m]))
    v <- c(v, "Submediterranean")

  if (macro == "Tropical" && Io < 3.6 && Io > 0.1 && Ps_w > Pcm1)
    v <- c(v, "Antitropical")

  if (length(v) == 0) "—" else paste(v, collapse = " + ")
}

# ── RM-11.  Master diagnosis function ─────────────────────────
rivas_martinez_diagnosis <- function(dat, lat, name = "") {

  P    <- as.numeric(dat[1, ])
  T    <- as.numeric(dat[2, ])
  Tmin <- as.numeric(dat[3, ])
  Tmax <- as.numeric(dat[4, ])

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

  T_annual <- mean(T)
  cold_m   <- which.min(T)
  hot_m    <- which.max(T)
  m_val    <- Tmin[cold_m]
  M_val    <- Tmax[cold_m]

  Ic   <- T[hot_m] - T[cold_m]
  It   <- (T_annual + m_val + M_val) * 10
  C    <- .rm_itc_correction(Ic)
  Itc  <- It + C
  Tp   <- sum(T[T > 0]) * 10
  Pp   <- sum(P[T > 0])
  Io   <- if (Tp > 0) (Pp / Tp) * 10 else 0

  all_tri <- list(c(12,1,2), c(3,4,5), c(6,7,8), c(9,10,11))
  tri_P   <- sapply(all_tri, function(m) sum(P[m]))
  dry_tri <- all_tri[[which.min(tri_P)]]
  d2      <- dry_tri[order(P[dry_tri])[1:2]]
  pos_d   <- T[d2] > 0
  Tpd2    <- sum(T[d2][pos_d]) * 10
  Ppd2    <- sum(P[d2][pos_d])
  Iod2    <- if (Tpd2 > 0) (Ppd2 / Tpd2) * 10 else Inf

  Ios2  <- .rm_ios2(T, P, sum_m)
  Iosc4 <- .rm_iosc4(T, P, sum_m, prev_m)

  Pcm1 <- sum(P[sum_m])
  Pcm2 <- sum(P[fall_m])
  Pcm3 <- sum(P[spr_m])
  Ps_w <- sum(P[win_m])

  s_med  <- .rm_is_mediterranean(Ios2, Iosc4, T, P, sum_m)
  s_rain <- .rm_summer_rain(Pcm1, Ps_w)

  macro <- .rm_macrobioclima(T, P, lat, Ic, T_annual,
                              m_val, M_val, Itc, Tp,
                              Pcm1, Ps_w, sum_m, prev_m)

  bio   <- .rm_bioclima(macro, Io, Iod2, Ic)
  termo <- .rm_termotipo(macro, Itc, Tp, use_Tp = (Ic >= 21 || Itc < 120))
  ombro <- .rm_ombrotipo(Io)
  cont  <- .rm_continentality(Ic)
  var   <- .rm_variant(macro, Ic, Io, Pcm1, Ps_w, T, P, sum_m, win_m)

  rule  <- paste0(rep("─", 64), collapse = "")
  drule <- paste0(rep("═", 64), collapse = "")
  lbl   <- function(label, value)
    cat(sprintf("  %-32s %s\n", paste0(label, " :"), value))

  cat("\n", drule, "\n", sep = "")
  cat("  RIVAS-MARTÍNEZ BIOCLIMATIC CLASSIFICATION (2004)\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Station: %s\n", name))
  cat(drule, "\n\n", sep = "")

  cat("  Bioclimatic Indices\n", rule, "\n", sep = "")
  cat(sprintf("  %-6s Continentality index         %8.2f\n",  "Ic",    Ic))
  cat(sprintf("  %-6s Thermicity index             %8.1f\n",  "It",    It))
  cat(sprintf("  %-6s It correction                %8.1f\n",  "C",     C))
  cat(sprintf("  %-6s Compensated thermicity       %8.1f\n",  "Itc",   Itc))
  cat(sprintf("  %-6s Annual positive temp (×10)   %8.0f\n",  "Tp",    Tp))
  cat(sprintf("  %-6s Positive precipitation (mm)  %8.1f\n",  "Pp",    Pp))
  cat(sprintf("  %-6s Annual ombrothermic index    %8.2f\n",  "Io",    Io))
  cat(sprintf("  %-6s Dry bimester ombrothermic   %8.2f\n",  "Iod2",  Iod2))
  cat(sprintf("  %-6s Summer bimester ombrothermic%8.2f\n",  "Ios2",  Ios2))
  cat(sprintf("  %-6s Compensated ombrothermic     %8.2f\n",  "Iosc4", Iosc4))
  cat(sprintf("  %-6s Summer quarter precip (mm)   %8.1f\n",  "Pcm1",  Pcm1))
  cat(sprintf("  %-6s Winter quarter precip (mm)   %8.1f\n",  "Ps_w",  Ps_w))
  cat(sprintf("  %-6s Following quarter precip (mm)%8.1f\n",  "Pcm2",  Pcm2))
  cat(sprintf("  %-6s Preceding quarter precip (mm)%8.1f\n",  "Pcm3",  Pcm3))

  med_rule <- if (Ios2 < 2 && Iosc4 < 2) {
    "Mediterranean  (Ios2<2 and Iosc4<2)"
  } else if (Ios2 >= 2 && Iosc4 >= 2) {
    "Temperate      (Ios2>=2 and Iosc4>=2)"
  } else {
    paste0("Mixed — Walter criterion: ",
           if (s_med) "Mediterranean (>=2 months P<=2T)" else "Temperate")
  }

  cat("\n  Mediterranean / Temperate Boundary\n", rule, "\n", sep = "")
  cat(sprintf("  Ios2  = %5.2f  |  Iosc4 = %5.2f  →  %s\n", Ios2, Iosc4, med_rule))
  cat(sprintf("  Summer rain (Pcm_summer > Pcm_winter) : %s\n",
              if (s_rain) "YES  → Tropical criterion"
              else        "NO   → winter wetter than summer"))

  cat("\n  Classification\n", rule, "\n", sep = "")
  lbl("Macrobioclimate", macro)
  lbl("Bioclimate",      bio)
  lbl("Thermotype",      termo)
  lbl("Ombrotype",       ombro)
  lbl("Continentality",  cont)
  lbl("Variant",         var)
  cat(drule, "\n\n", sep = "")

  invisible(list(macro = macro, bio = bio, termo = termo,
                 ombro = ombro, cont = cont, var = var,
                 Ic = Ic, It = It, Itc = Itc, Tp = Tp,
                 Pp = Pp, Io = Io, Iod2 = Iod2,
                 Ios2 = Ios2, Iosc4 = Iosc4))
}


# ══════════════════════════════════════════════════════════════
#  KÖPPEN-GEIGER CLIMATE CLASSIFICATION
#  Reference: Peel, M.C., Finlayson, B.L. & McMahon, T.A. (2007)
#    Updated world map of the Köppen-Geiger climate classification.
#    Hydrol. Earth Syst. Sci., 11, 1633–1644.
# ══════════════════════════════════════════════════════════════

koppen_geiger_class <- function(T, P, lat) {
  # T   : monthly mean temperature (°C)  – 12-element vector
  # P   : monthly precipitation (mm)     – 12-element vector
  # lat : decimal latitude (+N / –S)

  Tann  <- mean(T)
  Pann  <- sum(P)
  Tcold <- min(T)
  Thot  <- max(T)
  Pmin  <- min(P)
  n10   <- sum(T >= 10)        # months with T ≥ 10 °C

  NH <- lat >= 0
  # "Summer" = warm half-year (Apr–Sep NH, Oct–Mar SH)
  if (NH) {
    sum_idx <- 4:9
    win_idx <- c(10, 11, 12, 1, 2, 3)
  } else {
    sum_idx <- c(10, 11, 12, 1, 2, 3)
    win_idx <- 4:9
  }

  Ps     <- sum(P[sum_idx])
  Pw_    <- sum(P[win_idx])
  Pmin_s <- min(P[sum_idx])
  Pmax_s <- max(P[sum_idx])
  Pmin_w <- min(P[win_idx])
  Pmax_w <- max(P[win_idx])

  # B-type aridity threshold (Peel et al. 2007)
  if      (Ps / Pann >= 0.7)   Pth <- 2 * Tann + 28
  else if (Pw_ / Pann >= 0.7)  Pth <- 2 * Tann
  else                          Pth <- 2 * Tann + 14

  # Temperature subtype for C and D groups
  get_sub <- function() {
    if      (Thot  >= 22)  "a"   # hot summer
    else if (n10   >=  4)  "b"   # warm summer
    else if (Tcold < -38)  "d"   # extremely cold winter
    else                   "c"   # cold summer/winter
  }
  sub_desc_map <- c(a = "hot summer",
                    b = "warm summer",
                    c = "cold summer/winter",
                    d = "extremely cold winter")

  # ── Group A: Tropical (Tcold ≥ 18 °C) ───────────────────
  if (Tcold >= 18) {
    if (Pmin >= 60)
      return(list(code  = "Af",
                  group = "A – Tropical",
                  desc  = "Tropical rainforest – no dry season"))
    if (Pmin >= 100 - Pann / 25)
      return(list(code  = "Am",
                  group = "A – Tropical",
                  desc  = "Tropical monsoon"))
    # Aw vs As: which half-year is drier?
    if (Pmin_w <= Pmin_s)
      return(list(code  = "Aw",
                  group = "A – Tropical",
                  desc  = "Tropical savanna – dry winter"))
    return(list(code  = "As",
                group = "A – Tropical",
                desc  = "Tropical savanna – dry summer"))
  }

  # ── Group B: Arid (annual P < threshold) ─────────────────
  if (Pann < 10 * Pth) {
    hot <- Tann >= 18
    if (Pann < 5 * Pth) {
      code <- if (hot) "BWh" else "BWk"
      desc <- if (hot) "Hot desert" else "Cold desert"
    } else {
      code <- if (hot) "BSh" else "BSk"
      desc <- if (hot) "Hot semi-arid (steppe)" else "Cold semi-arid (steppe)"
    }
    return(list(code = code, group = "B – Arid", desc = desc))
  }

  sub   <- get_sub()
  sdesc <- sub_desc_map[sub]

  # ── Group C: Temperate (–3 < Tcold < 18, Thot ≥ 10) ─────
  if (Tcold > -3 && Thot >= 10) {
    if (Pmin_s < 40 && Pmax_w >= 3 * Pmin_s)
      return(list(code  = paste0("Cs", sub),
                  group = "C – Temperate",
                  desc  = paste0("Mediterranean – dry summer, ", sdesc)))
    if (Pmax_s >= 10 * Pmin_w)
      return(list(code  = paste0("Cw", sub),
                  group = "C – Temperate",
                  desc  = paste0("Temperate – dry winter, ", sdesc)))
    return(list(code  = paste0("Cf", sub),
                group = "C – Temperate",
                desc  = paste0("Temperate – no dry season, ", sdesc)))
  }

  # ── Group D: Continental (Tcold ≤ –3, Thot ≥ 10) ────────
  if (Thot >= 10 && Tcold <= -3) {
    if (Pmin_s < 40 && Pmax_w >= 3 * Pmin_s)
      return(list(code  = paste0("Ds", sub),
                  group = "D – Continental",
                  desc  = paste0("Continental – dry summer, ", sdesc)))
    if (Pmax_s >= 10 * Pmin_w)
      return(list(code  = paste0("Dw", sub),
                  group = "D – Continental",
                  desc  = paste0("Continental – dry winter, ", sdesc)))
    return(list(code  = paste0("Df", sub),
                group = "D – Continental",
                desc  = paste0("Continental – no dry season, ", sdesc)))
  }

  # ── Group E: Polar (Thot < 10 °C) ────────────────────────
  if (Thot >= 0)
    return(list(code = "ET", group = "E – Polar", desc = "Tundra"))
  return(list(code = "EF", group = "E – Polar", desc = "Ice cap (perpetual frost)"))
}


# ══════════════════════════════════════════════════════════════
#  WALTER ZONOBIOME / ZONOECOTONE CLASSIFICATION
#  Reference: Walter, H. (1985) Vegetation of the Earth.
#             Springer-Verlag, Berlin.
#             Walter, H. & Breckle, S.-W. (1985) Ökologie der
#             Erde, Vol. 1. Fischer, Stuttgart.
#
#  Nine zonobiomes (ZB I – IX) based on climate physiognomy.
#  Zonoecotones are transitional belts between adjacent ZBs,
#  flagged when the site meets borderline criteria.
#
#  Dry month criterion (Walter): P (mm) < 2 × T (°C), T > 0.
# ══════════════════════════════════════════════════════════════

walter_zonobiome_class <- function(T, P, lat) {
  # T   : monthly mean temperature (°C)  – 12-element vector
  # P   : monthly precipitation (mm)     – 12-element vector
  # lat : decimal latitude (+N / –S)

  Tann  <- mean(T)
  Pann  <- sum(P)
  Tcold <- min(T)
  Thot  <- max(T)
  Pmin  <- min(P)

  NH <- lat >= 0
  if (NH) {
    sum_m <- 6:8
    win_m <- c(12, 1, 2)
  } else {
    sum_m <- c(12, 1, 2)
    win_m <- 6:8
  }

  Ps <- sum(P[sum_m])
  Pw <- sum(P[win_m])

  # Walter dry months (only count when T > 0)
  dry     <- (T > 0) & (P < 2 * T)
  n_dry   <- sum(dry)
  dry_sum <- sum(dry[sum_m])
  dry_win <- sum(dry[win_m])

  # Aridity (Köppen B threshold, used as proxy for Walter's water balance)
  if      (Ps / Pann >= 0.7)           Pth <- 2 * Tann + 28
  else if ((Pann - Ps) / Pann >= 0.7)  Pth <- 2 * Tann
  else                                  Pth <- 2 * Tann + 14
  is_semiarid <- Pann < 10 * Pth
  is_desert   <- Pann < 5  * Pth

  # Helper ─────────────────────────────────────────────────
  zb_out <- function(zb, name, zono = NULL) {
    list(zb          = zb,
         name        = name,
         zonoecotone = if (!is.null(zono)) zono else "—")
  }

  # ── ZB IX: Tundra / Polar ─────────────────────────────────
  if (Thot < 10) {
    zono <- if (Thot >= 6) "ZB VIII / ZB IX" else NULL
    return(zb_out("ZB IX", "Tundra / Polar (Arctic–Alpine)", zono))
  }

  # ── ZB VIII: Boreal / Taiga ───────────────────────────────
  if (Tcold < -10 && Thot >= 10 && !is_semiarid) {
    zono <- if (Thot  <  15)  "ZB VIII / ZB IX" else
            if (Tcold >= -15) "ZB VI / ZB VIII"  else NULL
    return(zb_out("ZB VIII", "Boreal – Taiga (Cold Coniferous Forest)", zono))
  }

  # ── ZB VII: Temperate Arid ────────────────────────────────
  if (is_semiarid && Tcold < 5 && Tann < 20) {
    zono <- if (!is_desert && Tann >= 18) "ZB III / ZB VII" else
            if (!is_desert && Tcold < -5) "ZB VII / ZB VIII" else
            if (!is_desert)               "ZB VI / ZB VII"   else NULL
    return(zb_out("ZB VII",
                  "Temperate Arid – Continental Steppe / Cold Desert", zono))
  }

  # ── ZB III: Subtropical / Tropical Arid ──────────────────
  if (is_semiarid && (Tann >= 18 || (Tann >= 15 && Tcold >= 5))) {
    zono <- if (!is_desert && Ps > Pw) "ZB II / ZB III" else
            if (!is_desert && Pw > Ps) "ZB III / ZB IV" else NULL
    return(zb_out("ZB III",
                  "Subtropical Arid – Hot Desert / Semi-desert", zono))
  }

  # ── ZB I: Equatorial / Tropical Rainforest ────────────────
  if (Tcold >= 18 && n_dry == 0) {
    zono <- if (Pmin < 100) "ZB I / ZB II" else NULL
    return(zb_out("ZB I",
                  "Equatorial – Tropical Rainforest (Always Humid)", zono))
  }

  # ── ZB II: Tropical with Summer Rain ─────────────────────
  if (Tcold >= 18) {
    zono <- if (n_dry <= 2) "ZB I / ZB II" else
            if (n_dry >= 6) "ZB II / ZB III" else NULL
    return(zb_out("ZB II",
                  "Tropical – Summer Rain / Savanna / Deciduous Forest", zono))
  }

  # ── ZB IV: Mediterranean ──────────────────────────────────
  if (dry_sum >= 2 && Pw >= Ps && Tcold >= -1) {
    zono <- if (n_dry <= 2)                 "ZB IV / ZB V"   else
            if (n_dry >= 7 || is_semiarid)  "ZB III / ZB IV" else
            if (Tcold < 3)                  "ZB IV / ZB VI"  else NULL
    return(zb_out("ZB IV",
                  "Mediterranean – Winter Rain, Summer Drought (Sclerophyllous)",
                  zono))
  }

  # ── ZB V: Warm Temperate / Lauriphyllous ──────────────────
  if (Tann >= 13 && Tcold >= 1 && n_dry <= 1) {
    zono <- if (dry_sum == 1) "ZB IV / ZB V" else
            if (Tcold   <  3) "ZB V / ZB VI"  else NULL
    return(zb_out("ZB V",
                  "Warm Temperate – Lauriphyllous / Humid Subtropical", zono))
  }

  # ── ZB VI: Temperate Nemoral ──────────────────────────────
  if (Thot >= 10 && Tcold >= -15 && n_dry <= 2) {
    zono <- if (Tcold  < -10)               "ZB VI / ZB VIII" else
            if (Tcold  >= -1 && Tann >= 12) "ZB V / ZB VI"    else
            if (n_dry  >=  2)               "ZB VI / ZB VII"  else NULL
    return(zb_out("ZB VI",
                  "Temperate Nemoral – Deciduous / Mixed Forest", zono))
  }

  # ── Fallbacks ─────────────────────────────────────────────
  if (Thot >= 10)
    return(zb_out("ZB VIII", "Boreal – Taiga (Cold Coniferous Forest)", NULL))

  return(zb_out("ZB IX", "Tundra / Polar (Arctic–Alpine)", NULL))
}


# ══════════════════════════════════════════════════════════════
#  COMBINED PRINT: Köppen-Geiger + Walter Zonobiome
# ══════════════════════════════════════════════════════════════

print_koppen_walter <- function(T, P, lat, name = "") {

  rule  <- paste0(rep("─", 64), collapse = "")
  drule <- paste0(rep("═", 64), collapse = "")
  lbl   <- function(label, value)
    cat(sprintf("  %-32s %s\n", paste0(label, " :"), value))

  kg <- koppen_geiger_class(T, P, lat)
  wz <- walter_zonobiome_class(T, P, lat)

  # ── Köppen-Geiger ──────────────────────────────────────────
  cat("\n", drule, "\n", sep = "")
  cat("  KÖPPEN-GEIGER CLIMATE CLASSIFICATION\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Station: %s\n", name))
  cat(drule, "\n\n", sep = "")
  cat("  Peel, Finlayson & McMahon (2007)\n", rule, "\n", sep = "")
  lbl("Code",        kg$code)
  lbl("Group",       kg$group)
  lbl("Description", kg$desc)
  cat(drule, "\n\n", sep = "")

  # ── Walter Zonobiome / Zonoecotone ────────────────────────
  cat("\n", drule, "\n", sep = "")
  cat("  WALTER ZONOBIOME / ZONOECOTONE\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Station: %s\n", name))
  cat(drule, "\n\n", sep = "")
  cat("  Walter (1985) – Vegetation of the Earth\n", rule, "\n", sep = "")
  lbl("Zonobiome",                wz$zb)
  lbl("Name",                     wz$name)
  lbl("Zonoecotone (transition)", wz$zonoecotone)
  cat(drule, "\n\n", sep = "")

  invisible(list(koppen = kg, walter = wz))
}


# ══════════════════════════════════════════════════════════════
#  MUCINA (2019) GLOBAL BIOME CLASSIFICATION
#  Reference: Mucina, L. (2019). Biome: evolution of a crucial
#    ecological and biogeographical concept. New Phytologist,
#    222(1), 97–114.
#
#  12 global biomes — deterministic, rule-based from monthly T & P.
#  These correspond to Mucina's "physiognomic" biomes, whose
#  boundaries are driven primarily by macroclimate.
#
#  The sub-global "continental biome" level (44 units) is
#  floristically defined and requires biogeographic information
#  beyond what monthly climate alone can provide; it is therefore
#  not computed here.
#
#  Global biome codes and names:
#    PD    Polar Desert & Ice Sheet
#    TUN   Tundra & Polar Grassland
#    BFT   Boreal Forest & Taiga
#    CDS   Cold Desert & Semi-desert
#    SD    Subtropical Hot Desert
#    TGS   Temperate Grassland, Steppe & Shrubland
#    TRF   Tropical Rainforest
#    TSF   Tropical Seasonal Forest & Woodland
#    TSG   Tropical Savanna & Grassland
#    MSW   Mediterranean Shrubland & Woodland
#    WTR   Warm-temperate Rainforest & Sclerophyllous Thicket
#    TDMF  Temperate Deciduous & Mixed Forest
# ══════════════════════════════════════════════════════════════

mucina_global_biome_class <- function(T, P, lat) {
  # T   : monthly mean temperature (°C)  – 12-element vector
  # P   : monthly precipitation (mm)     – 12-element vector
  # lat : decimal latitude (+N / –S)

  Tann  <- mean(T)
  Pann  <- sum(P)
  Tcold <- min(T)
  Thot  <- max(T)

  NH <- lat >= 0
  if (NH) { sum_m <- 6:8; win_m <- c(12, 1, 2) }
  else    { sum_m <- c(12, 1, 2); win_m <- 6:8 }

  Ps  <- sum(P[sum_m])
  Pw  <- sum(P[win_m])
  f_s <- Ps / Pann

  # Walter dry months (P < 2T, T > 0)
  dry     <- (T > 0) & (P < 2 * T)
  n_dry   <- sum(dry)
  dry_sum <- sum(dry[sum_m])

  # Aridity threshold (Köppen-B method)
  if      (f_s >= 0.70)                Pth <- 2 * Tann + 28
  else if ((Pann - Ps) / Pann >= 0.70) Pth <- 2 * Tann
  else                                  Pth <- 2 * Tann + 14
  is_arid   <- Pann < 10 * Pth
  is_desert <- Pann < 5  * Pth

  # ── Polar ─────────────────────────────────────────────────
  if (Thot < 0)
    return(list(code = "PD",
                name = "Polar Desert & Ice Sheet"))
  if (Thot < 10)
    return(list(code = "TUN",
                name = "Tundra & Polar Grassland"))

  # ── Boreal – cold and humid ───────────────────────────────
  if (Tcold < -10 && !is_arid)
    return(list(code = "BFT",
                name = "Boreal Forest & Taiga"))

  # ── True deserts ──────────────────────────────────────────
  if (is_desert) {
    if (Tann < 18 && Tcold < 5)
      return(list(code = "CDS",
                  name = "Cold Desert & Semi-desert"))
    return(list(code = "SD",
                name = "Subtropical Hot Desert"))
  }

  # ── Semi-arid ─────────────────────────────────────────────
  if (is_arid) {
    if (Tann >= 18 || Tcold >= 5)
      return(list(code = "SD",
                  name = "Subtropical Hot Desert"))
    return(list(code = "TGS",
                name = "Temperate Grassland, Steppe & Shrubland"))
  }

  # ── Always-hot tropical (Tcold ≥ 18 °C) ──────────────────
  if (Tcold >= 18) {
    if (n_dry == 0)
      return(list(code = "TRF",
                  name = "Tropical Rainforest"))
    if (n_dry <= 4)
      return(list(code = "TSF",
                  name = "Tropical Seasonal Forest & Woodland"))
    return(list(code = "TSG",
                name = "Tropical Savanna & Grassland"))
  }

  # ── Mediterranean (summer drought, winter rain) ───────────
  if (dry_sum >= 2 && Pw > Ps && Tcold >= -1)
    return(list(code = "MSW",
                name = "Mediterranean Shrubland & Woodland"))

  # ── Warm-temperate rainforest (humid subtropical) ─────────
  if (Tann >= 14 && Tcold >= 3 && n_dry == 0)
    return(list(code = "WTR",
                name = "Warm-temperate Rainforest & Sclerophyllous Thicket"))

  # ── Temperate deciduous / mixed forest ────────────────────
  if (Thot >= 10 && n_dry <= 2)
    return(list(code = "TDMF",
                name = "Temperate Deciduous & Mixed Forest"))

  # ── Default: dry temperate ────────────────────────────────
  return(list(code = "TGS",
              name = "Temperate Grassland, Steppe & Shrubland"))
}


# ── Print Mucina global biome ─────────────────────────────────
print_mucina <- function(T, P, lat, name = "") {

  rule  <- paste0(rep("─", 64), collapse = "")
  drule <- paste0(rep("═", 64), collapse = "")
  lbl   <- function(label, value)
    cat(sprintf("  %-32s %s\n", paste0(label, " :"), value))

  gb <- mucina_global_biome_class(T, P, lat)

  cat("\n", drule, "\n", sep = "")
  cat("  MUCINA (2019) GLOBAL BIOME CLASSIFICATION\n")
  cat("  Mucina, L. (2019) New Phytologist 222(1): 97-114\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Station: %s\n", name))
  cat(drule, "\n\n", sep = "")
  cat("  Global Biome  (deterministic from macroclimate)\n",
      rule, "\n", sep = "")
  lbl("Code",         gb$code)
  lbl("Global Biome", gb$name)
  cat(drule, "\n\n", sep = "")

  invisible(gb)
}


# ══════════════════════════════════════════════════════════════
#  LOIDI, NAVARRO-SÁNCHEZ & VYNOKUROV (2022)
#  BIOCLIMATIC AND BIOME CLASSIFICATION
#  Reference: Loidi, J., Navarro-Sánchez, G. & Vynokurov, D. (2022).
#    Climatic definitions of the world's terrestrial biomes.
#    Vegetation Classification and Survey, 3: 231–271.
#    doi: 10.3897/VCS.86102
#
#  Parallel dual classification:
#
#  A. Bioclimatic classification
#     4 Macrobioclimates (I Cold · II Mesic · III Arid · IV Warm)
#    10 Bioclimates      (I.a Tundral · I.b Boreal · II.a Pluviestival ·
#                         II.b Aridiestival · II.c Steppic · III.a Cold
#                         desertic · III.b Temperate desertic · III.c Warm
#                         desertic · IV.a Pluviseasonal · IV.b Pluvial)
#     7 Ombrotypes (Table 3)  ·  6 Thermotypes (Table 4, lat-dependent)
#
#  B. Biotic unit classification
#     4 Domains → 7 Ecozones → 9 Biomes → 20 Subbiomes  (Table 7)
#
#  Indices (Table 5):
#    Tp    : positive temperature = Σ T_i × 10  (months T > 0)
#    Pp    : positive precipitation = Σ P_i      (months T > 0)
#    Io    : Pp / Tp  — ombrothermic index
#    BIO7  : annual temperature range = max(Tmax) − min(Tmin)  (°C)
#    BIO10 : mean T of warmest consecutive quarter  (°C)
#    BIO18 : precipitation of warmest consecutive quarter  (mm)
#    dws   : drought of warm season = BIO18 / BIO10
#    BIO15 : precipitation seasonality = 100 × sd(P) / (mean(P) + 1)
#    nfd   : frost days (approximated from monthly T and Tmin)
#
#  Macrobioclimate assignment (priority order):
#    1. Io < 1                     → III Arid
#    2. Tp < 1000                  → I   Cold
#    3. Tp > 2000                  → IV  Warm Tropical
#    4. otherwise (Tp 1000–2000)   → II  Mesic
# ══════════════════════════════════════════════════════════════


# ── LV-0.  Internal helpers ───────────────────────────────────

# Three consecutive months with the highest sum of T (wrapping)
.loidi_warm_qtr <- function(T) {
  best <- -Inf; s <- 1L
  for (i in 1:12) {
    idx <- (i - 1L + 0:2) %% 12L + 1L
    v   <- sum(T[idx])
    if (v > best) { best <- v; s <- i }
  }
  (s - 1L + 0:2) %% 12L + 1L
}

# Frost-day approximation from monthly T and Tmin
#   Months with T < 0  → ~30 frost days each
#   Months with T ≥ 0 but Tmin < 0 → ~15 frost days each
.loidi_nfd <- function(T, Tmin) {
  round(30 * sum(T < 0) + 15 * sum(Tmin < 0 & T >= 0))
}


# ── LV-1.  Ombrotype (Table 3) ───────────────────────────────
loidi_ombrotype <- function(Io) {
  if      (Io <  0.4) "Hyper-arid"
  else if (Io <  1.0) "Arid"
  else if (Io <  2.0) "Sub-arid"
  else if (Io <  3.6) "Dry"
  else if (Io <  6.0) "Sub-humid"
  else if (Io < 12.0) "Humid"
  else                "Hyper-humid"
}


# ── LV-2.  Thermotype (Table 4) ──────────────────────────────
#
#  For each latitude band the table gives the lower Tp bound of
#  each thermotype band (checked from warmest to coldest).
#  NA = that thermotype does not exist at that latitude.
#
loidi_thermotype <- function(Tp, lat) {
  a <- abs(lat)
  # Lower Tp bounds: [Infra, Thermo, Meso, Supra, Oro, Cryoro]
  lo <- if      (a <= 20) c(2900, 2300, 1700,  950, 450, 100)
        else if (a <= 40) c(2400, 2100, 1500,  900, 450, 100)
        else if (a <= 50) c(  NA, 2000, 1400,  800, 380, 100)
        else if (a <= 60) c(  NA,   NA,  800,  380,  NA, 100)
        else if (a <= 70) c(  NA,   NA,   NA,  380,  NA, 100)
        else if (a <= 80) c(  NA,   NA,   NA,   NA,  NA, 100)
        else              c(  NA,   NA,   NA,   NA,  NA,  NA)
  nms <- c("Infra", "Thermo", "Meso", "Supra", "Oro", "Cryoro")
  for (i in seq_along(lo))
    if (!is.na(lo[i]) && Tp >= lo[i]) return(nms[i])
  "Icy"
}


# ── LV-3.  Bioclimatic classification ────────────────────────
#
#  Returns a named list with all computed indices plus:
#    macrobioclimate, bioclimate, thermotype, ombrotype
#
loidi_bioclimate_class <- function(T, P, Tmin, Tmax, lat) {

  Tp    <- sum(T[T > 0]) * 10
  Pp    <- sum(P[T > 0])
  Io    <- if (Tp > 0) (Pp / Tp) * 10 else 0   # × 10 — same formula as Rivas-Martínez

  wq    <- .loidi_warm_qtr(T)
  BIO10 <- mean(T[wq])
  BIO18 <- sum(P[wq])
  dws   <- if (BIO10 > 0) BIO18 / BIO10 else Inf
  BIO7  <- max(Tmax) - min(Tmin)
  BIO15 <- 100 * sd(P) / (mean(P) + 1)
  nfd   <- .loidi_nfd(T, Tmin)
  Tann  <- mean(T)
  Thot  <- max(T)
  Tcold_m <- min(T)   # coldest monthly mean T

  # Mediterranean = winter rain maximum  +  expressive T seasonality
  # Tropical      = summer rain maximum  +  feeble     T seasonality
  #
  # These are OPPOSITE signatures; both must agree before a site is
  # assigned to the Thermocratic (Warm Tropical) domain.
  #
  # 1. Rain seasonality — hemisphere-aware warm vs cool semester totals
  #    (Loidi et al. 2022: rain seasonality determines bioclimate within
  #    every domain; Mediterranean is uniquely defined by WINTER maximum)
  NH    <- lat >= 0
  a_lat <- abs(lat)
  sum6  <- if (NH) 4:9 else c(10L, 11L, 12L, 1L, 2L, 3L)
  win6  <- if (NH) c(10L, 11L, 12L, 1L, 2L, 3L) else 4:9
  Ps6   <- sum(P[sum6])   # warm-semester precipitation (mm)
  Pw6   <- sum(P[win6])   # cool-semester precipitation (mm)
  # Walter dry months in the warm semester — confirms Mediterranean summer drought
  # even when dws > 1 (e.g. Lisbon: early-autumn rain inflates the warmest quarter)
  n_dry_sum <- sum((T[sum6] > 0) & (P[sum6] < 2 * T[sum6]))
  #   For equatorial sites (abs_lat < 10°) rainfall is bimodal and
  #   nearly equal in both semesters; require only Ps6 ≥ 75 % of Pw6
  #   (i.e. winter is NOT strongly dominant) rather than Ps6 ≥ Pw6.
  rain_trop <- if (a_lat < 10) Ps6 >= 0.75 * Pw6  else  Ps6 >= Pw6

  # 2. Thermal seasonality — range of monthly mean temperatures
  #    (T_seas < 25 °C distinguishes genuinely tropical thermal regimes
  #    from temperate/Mediterranean ones with clear winter/summer contrast;
  #    threshold set high enough to include continental monsoon areas
  #    such as the Indian subcontinent where T_seas can reach 20 °C)
  T_seas    <- Thot - Tcold_m
  tseas_trop <- T_seas < 25

  # ── Macrobioclimate (priority: Arid → Cold → Warm → Mesic) ──
  is_warm_trop <- Tp > 2000 && rain_trop && tseas_trop
  macro <- if      (Io  <  1)        "III – Arid (Xerocratic)"
            else if (Tp  < 1000)     "I – Cold (Cryocratic)"
            else if (is_warm_trop)   "IV – Warm Tropical (Thermocratic)"
            else                     "II – Mesic (Mesocratic)"

  therm <- loidi_thermotype(Tp, lat)
  ombro <- loidi_ombrotype(Io)

  # ── Bioclimate ───────────────────────────────────────────────
  bioclim <- if (grepl("Arid", macro)) {
    if      (Tann >= 18 && nfd < 20)          "III.c – Warm desertic"
    else if (nfd > 40 || (Tann < 12 && BIO7 > 25)) "III.a – Continental cold desertic"
    else                                       "III.b – Temperate desertic"

  } else if (grepl("Cold", macro)) {
    if (Thot < 10) "I.a – Tundral"  else  "I.b – Boreal"

  } else if (grepl("Warm", macro)) {
    if (BIO15 < 60) "IV.b – Tropical pluvial"  else  "IV.a – Tropical pluviseasonal"

  } else {
    # Mesic — priority: aridiestival first (summer drought defines Mediterranean),
    # then steppic (hypercontinentality without summer drought), then pluviestival.
    # Steppic and Mediterranean are antithetical: steppic has no summer drought;
    # Mediterranean has intense summer drought.  BIO7 can be large in Med areas
    # (hot summers + cold nights) so it must not be checked before dws.
    #
    # The primary aridiestival signal is dws ≤ 1.  A secondary criterion
    # (Pw6 > Ps6 AND n_dry_sum ≥ 2) catches sites like Lisbon where
    # early-autumn rain falls within the warmest consecutive quarter, pushing
    # dws just above 1 despite a genuine and prolonged summer drought.
    med_sign <- (dws <= 1) || (Pw6 > Ps6 && n_dry_sum >= 2)
    if      (med_sign)     "II.b – Temperate aridiestival (Mediterranean s.l.)"
    else if (BIO7 >= 30)   "II.c – Temperate steppic (hypercontinental)"
    else                   "II.a – Temperate pluviestival (summer-rainy)"
  }

  list(Tp = Tp, Pp = Pp, Io = Io,
       BIO7 = BIO7, BIO10 = BIO10, BIO18 = BIO18,
       dws = dws, BIO15 = BIO15, nfd = nfd,
       Ps6 = Ps6, Pw6 = Pw6, T_seas = T_seas, n_dry_sum = n_dry_sum,
       Tann = Tann, Thot = Thot,
       thermotype      = therm,
       ombrotype       = ombro,
       macrobioclimate = macro,
       bioclimate      = bioclim)
}


# ── LV-4.  Biome classification (Table 7) ─────────────────────
#
#  Returns domain, ecozone, biome, subbiome as character strings.
#  Subbiome assignment uses climate + coarse geographic context
#  (lat, lon) where needed.  Mountain subbiomes (1b, 1c, 2b, 4c)
#  require elevation data not available here; they are flagged
#  when the climate profile is consistent with them.
#
loidi_biome_class <- function(T, P, Tmin, Tmax, lat, lon) {

  bc    <- loidi_bioclimate_class(T, P, Tmin, Tmax, lat)
  macro <- bc$macrobioclimate
  bioclim <- bc$bioclimate
  Io    <- bc$Io;   Tp  <- bc$Tp;  BIO7  <- bc$BIO7
  BIO15 <- bc$BIO15; dws <- bc$dws; nfd  <- bc$nfd
  Tann  <- bc$Tann;  Thot <- bc$Thot
  Tcold <- min(T)
  Tmin_cold <- Tmin[which.min(T)]
  n_dry <- sum((T > 0) & (P < 2 * T))    # Walter dry months
  a_lat <- abs(lat)

  # ── Domain A – Cryocratic ────────────────────────────────────
  if (grepl("Cold", macro)) {
    if (Thot < 10) {
      sub <- if (a_lat > 55)
               "1a – Polar (Arctic/Antarctic) tundra"
             else if (Tann > 5 && a_lat < 30)
               "1c – Tundra of tropical mountains (cryoro belt)"
             else
               "1b – Tundra of temperate mountains (cryoro belt)"
      return(list(
        domain  = "A – Cryocratic (Domain of cold climates)",
        ecozone = "AA – Polar and Boreal ecozone",
        biome   = "1 – Biome of the tundra",
        subbiome = sub))
    } else {
      sub <- if (a_lat > 45) "2a – Lowland boreal taiga"
             else            "2b – Forests/shrublands of the temperate oro belt"
      return(list(
        domain  = "A – Cryocratic (Domain of cold climates)",
        ecozone = "AA – Polar and Boreal ecozone",
        biome   = "2 – Biome of the boreal forest",
        subbiome = sub))
    }
  }

  # ── Domain C – Xerocratic ────────────────────────────────────
  if (grepl("Arid", macro)) {
    sub <- if      (grepl("Warm",   bioclim)) "7c – Warm deserts and semi-deserts"
           else if (grepl("cold",   bioclim)) "7a – Cold deserts and semi-deserts"
           else                               "7b – Temperate deserts and semi-deserts"
    return(list(
      domain  = "C – Xerocratic (Domain of arid climates)",
      ecozone = "CA – Ecozone of the deserts and semi-deserts",
      biome   = "7 – Biome of the deserts and semi-deserts of arid regions",
      subbiome = sub))
  }

  # ── Domain D – Thermocratic ──────────────────────────────────
  if (grepl("Warm", macro)) {
    if (BIO15 < 60) {
      return(list(
        domain  = "D – Thermocratic (Domain of warm climates)",
        ecozone = "DB – Tropical pluvial, rainy all year ecozone",
        biome   = "9 – Biome of the tropical rain forests",
        subbiome = "9a – Tropical rain forests"))
    } else {
      # 8a vs 8b: xeric (very long dry season) vs pluviseasonal forest
      sub <- if (n_dry > 7 || (Io < 1.5 && BIO15 > 80))
               "8a – Tropical xeric shrublands and woodlands"
             else
               "8b – Tropical pluviseasonal forests and woodlands"
      return(list(
        domain  = "D – Thermocratic (Domain of warm climates)",
        ecozone = "DA – Tropical pluviseasonal, rainy and dry seasons",
        biome   = "8 – Biome of the tropical pluviseasonal forests and shrublands",
        subbiome = sub))
    }
  }

  # ── Domain B – Mesocratic ────────────────────────────────────
  # BC ecozone – Steppic (Biome 6)
  if (grepl("steppic", bioclim)) {
    sub <- if (Io >= 2.5) "6a – Forest-steppe"  else  "6b – Grass-steppe"
    return(list(
      domain  = "B – Mesocratic (Domain of temperate climates)",
      ecozone = "BC – Temperate hypercontinental steppic ecozone",
      biome   = "6 – Biome of the steppe",
      subbiome = sub))
  }

  # BB ecozone – Aridiestival / Mediterranean (Biome 5)
  if (grepl("aridiestival", bioclim)) {
    sub <- if (lat < -38 && lon > -76 && lon < -63)
             "5c – Patagonian shrubland"
           else if (BIO7 >= 20 && Tmin_cold < 0)
             "5b – Continental scrub and woodlands"
           else
             "5a – Oceanic sclerophyllous-microphyllous (Mediterranean s.str.)"
    return(list(
      domain  = "B – Mesocratic (Domain of temperate climates)",
      ecozone = "BB – Temperate aridiestival ecozone",
      biome   = "5 – Biome of the temperate aridiestival evergreen forests and shrublands",
      subbiome = sub))
  }

  # BA ecozone – Pluviestival (Biomes 3 or 4)
  # Key: Tcold ≥ 5 °C → mild winters → evergreen (Biome 4)
  #      Tcold <  5 °C → cold winters → deciduous (Biome 3)
  if (Tcold >= 5) {
    # Biome 4 – Temperate pluvial evergreen forest, shrublands and grasslands
    sub <- if (a_lat < 28 && Tp < 2100)
             # Mesic climate at tropical lat → mountain cloud forest
             "4c – Tropical/subtropical montane cloud lauroid & conifer evergreen forest"
           else if (a_lat > 38 && BIO7 <= 15 &&
                    (lon < -120 || (lat < -35 && lon < -68)))
             # Very oceanic, mid-high lat (Pacific NW or southern Chile)
             "4b – Conifer coastal forests"
           else
             "4a – Lauroid subtropical evergreen forest of the lowlands"
    return(list(
      domain  = "B – Mesocratic (Domain of temperate climates)",
      ecozone = "BA – Temperate pluviestival ecozone",
      biome   = "4 – Biome of the temperate pluvial evergreen forest, shrublands and grasslands",
      subbiome = sub))
  } else {
    # Biome 3 – Temperate deciduous forests
    return(list(
      domain  = "B – Mesocratic (Domain of temperate climates)",
      ecozone = "BA – Temperate pluviestival ecozone",
      biome   = "3 – Biome of the temperate deciduous forests",
      subbiome = "3a – Temperate deciduous forests"))
  }
}


# ── LV-5.  Print Loidi classification ────────────────────────
print_loidi <- function(T, P, Tmin, Tmax, lat, lon, name = "") {

  rule  <- paste0(rep("─", 64), collapse = "")
  drule <- paste0(rep("═", 64), collapse = "")
  lbl   <- function(label, value)
    cat(sprintf("  %-36s %s\n", paste0(label, " :"), value))

  bc <- loidi_bioclimate_class(T, P, Tmin, Tmax, lat)
  bm <- loidi_biome_class(T, P, Tmin, Tmax, lat, lon)

  # ── Header ────────────────────────────────────────────────
  cat("\n", drule, "\n", sep = "")
  cat("  LOIDI, NAVARRO-SÁNCHEZ & VYNOKUROV (2022)\n")
  cat("  Vegetation Classification and Survey 3: 231-271\n")
  if (nchar(trimws(name)) > 0)
    cat(sprintf("  Station: %s\n", name))
  cat(drule, "\n\n", sep = "")

  # ── Part 1: Bioclimatic classification ────────────────────
  cat("  Bioclimatic Classification\n", rule, "\n", sep = "")
  cat(sprintf("  %-6s Positive temperature (×10)       %8.0f\n", "Tp",    bc$Tp))
  cat(sprintf("  %-6s Positive precipitation (mm)      %8.1f\n", "Pp",    bc$Pp))
  cat(sprintf("  %-6s Ombrothermic index (Pp/Tp)       %8.2f\n", "Io",    bc$Io))
  cat(sprintf("  %-6s Annual temperature range (°C)    %8.1f\n", "BIO7",  bc$BIO7))
  cat(sprintf("  %-6s Drought of warm season (BIO18/10)%8.2f\n", "dws",   bc$dws))
  cat(sprintf("  %-6s Precip seasonality (CV)          %8.1f\n", "BIO15", bc$BIO15))
  cat(sprintf("  %-6s Frost days (approx.)             %8.0f\n", "nfd",   bc$nfd))
  cat(sprintf("  %-6s Warm-semester precip (mm)        %8.1f\n", "Ps6",   bc$Ps6))
  cat(sprintf("  %-6s Cool-semester precip (mm)        %8.1f\n", "Pw6",   bc$Pw6))
  cat(sprintf("  %-6s Monthly T range Thot-Tcold (°C)  %8.1f\n", "T_seas",bc$T_seas))
  cat(rule, "\n", sep = "")
  lbl("Macrobioclimate", bc$macrobioclimate)
  lbl("Bioclimate",      bc$bioclimate)
  lbl("Thermotype",      bc$thermotype)
  lbl("Ombrotype",       bc$ombrotype)

  # ── Part 2: Biotic unit classification ────────────────────
  cat("\n  Biotic Unit Classification\n", rule, "\n", sep = "")
  lbl("Domain",   bm$domain)
  lbl("Ecozone",  bm$ecozone)
  lbl("Biome",    bm$biome)
  lbl("Subbiome", bm$subbiome)
  cat(drule, "\n\n", sep = "")

  invisible(list(bioclimate = bc, biome = bm))
}


# ══════════════════════════════════════════════════════════════
#  MAIN FUNCTION: climate_diagram()
# ══════════════════════════════════════════════════════════════

#' Draw a Walter-Lieth diagram and print five classifications
#' for any named location:
#'   1. Rivas-Martínez (2004)
#'   2. Köppen-Geiger – Peel et al. (2007)
#'   3. Walter Zonobiome / Zonoecotone – Walter (1985)
#'   4. Mucina (2019) global biome
#'   5. Loidi, Navarro-Sánchez & Vynokurov (2022)
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

  abs_min <- tmin - pmax((tavg - tmin) / 2, 0)

  southern <- lat < 0
  if (southern) message("  Southern hemisphere detected")

  ## ── Walter-Lieth matrix (rows: P / Tmax / Tmin / Tabs) ───
  clim_mat <- rbind(prec, tmax, tmin, abs_min)

  ## ── Rivas-Martínez matrix (rows: P / Tmean / Tmin / Tmax) ─
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

  ## ── Köppen-Geiger + Walter Zonobiome ──────────────────────
  kw <- print_koppen_walter(tavg, prec, lat, name = station_label)

  ## ── Mucina (2019) global biome ────────────────────────────
  mu <- print_mucina(tavg, prec, lat, name = station_label)

  ## ── Loidi, Navarro-Sánchez & Vynokurov (2022) ─────────────
  lo <- print_loidi(tavg, prec, tmin, tmax, lat, lon,
                    name = station_label)

  invisible(list(
    place       = place,
    display     = geo$display,
    lon         = lon,
    lat         = lat,
    elevation   = elev,
    southern    = southern,
    clim_mat    = clim_mat,
    rm_mat      = rm_mat,
    koppen      = kw$koppen,
    walter      = kw$walter,
    mucina_glob = mu,
    loidi       = lo
  ))
}


# ── Ready message ─────────────────────────────────────────────
message("
==========================================================
 Walter-Lieth + Rivas-Martínez + Köppen-Geiger
 + Walter Zonobiome + Mucina (2019)
 + Loidi, Navarro-Sánchez & Vynokurov (2022) loaded

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
