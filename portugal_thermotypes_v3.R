## ============================================================
##  Thermotype Map — Continental Portugal
##  Rivas-Martínez (2004) · WorldClim v2.1 · geodata
##
##  Fix v3: Walter consecutive-drought criterion is now the
##  NECESSARY GATE for Mediterranean macrobioclimate.
##
##  Logic:
##    1. Walter gate  — requires ≥2 consecutive months P < 2T
##       (strict) in May–Sep window.  If absent → Temperate.
##    2. Ombrothermic confirmation — if Walter gate passes,
##       at least one of Ios2 < 2 or Iosc4 < 2 confirms Med.
##       If neither is < 2 → Temperate.
##
##  This correctly assigns high-elevation cells (Serra da Estrela,
##  Serra de Arga, etc.) to Temperate even when their summer P/T
##  ratios are numerically low due to cool temperatures.
##
##  Diagnostic (Section 5b) prints pixel counts and the full
##  transition table comparing v3 against the original script.
## ============================================================


# ── 0. Packages ───────────────────────────────────────────────
required_pkgs <- c("terra", "geodata")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# ── 1. Configuration ──────────────────────────────────────────
WC_RES   <- 2.5
WC_CACHE <- file.path(getwd(), "worldclim_cache")
dir.create(WC_CACHE, recursive = TRUE, showWarnings = FALSE)

PT_EXT <- ext(-9.6, -6.1, 36.9, 42.2)


# ── 2. Load & crop WorldClim rasters ──────────────────────────
message("Loading WorldClim rasters ...")
r_prec <- crop(geodata::worldclim_global("prec", WC_RES, WC_CACHE), PT_EXT)
r_tmin <- crop(geodata::worldclim_global("tmin", WC_RES, WC_CACHE), PT_EXT)
r_tmax <- crop(geodata::worldclim_global("tmax", WC_RES, WC_CACHE), PT_EXT)
r_tavg <- crop(geodata::worldclim_global("tavg", WC_RES, WC_CACHE), PT_EXT)
message("  Done.")


# ── 3. Portugal administrative boundaries ─────────────────────
message("Loading Portugal boundaries ...")
pt_nuts <- geodata::gadm("PRT", level = 1, path = WC_CACHE)
pt_cont <- pt_nuts[!pt_nuts$NAME_1 %in% c("Açores", "Madeira"), ]
pt_land <- aggregate(pt_cont)
message("  Done.")


# ── 4. Classification helpers ─────────────────────────────────

.itc_c <- function(Ic) {
  if      (Ic <=  8) 10  * (Ic -  8)
  else if (Ic <= 18) 0
  else if (Ic <= 21) 5   * (Ic - 18)
  else if (Ic <= 28) 15  + 15 * (Ic - 21)
  else if (Ic <= 46) 120 + 25 * (Ic - 28)
  else               570 + 30 * (min(Ic, 65) - 46)
}

.ios2 <- function(T, P) {
  Ts <- T[6:8]; Ps <- P[6:8]
  h2   <- order(Ts, decreasing = TRUE)[1:2]
  pos  <- Ts[h2] > 0
  Tps2 <- sum(Ts[h2][pos]) * 10
  Pps2 <- sum(Ps[h2][pos])
  if (Tps2 > 0) (Pps2 / Tps2) * 10 else Inf
}

.iosc4 <- function(T, P) {
  T4 <- T[5:8]; P4 <- P[5:8]
  pos  <- T4 > 0
  Tps4 <- sum(T4[pos]) * 10
  Pps4 <- sum(P4[pos])
  if (Tps4 > 0) (Pps4 / Tps4) * 10 else Inf
}


## ── Mediterranean check v3 — Walter is the necessary gate ────
##
##  A cell is Mediterranean ONLY IF:
##    (a) it has ≥2 consecutive months with P < 2T (strict) in
##        the May–Sep window  [Walter necessary condition], AND
##    (b) at least one of Ios2 < 2 or Iosc4 < 2.
##
##  Mountain cells that lack a real summer drought run (e.g.
##  Serra da Estrela, Serra de Arga) fail (a) → Temperate.
##
.is_med <- function(Ios2, Iosc4, T, P) {

  # ── Step 1: Walter gate ──────────────────────────────────
  win <- c(5, 6, 7, 8, 9)
  dry <- P[win] < 2 * T[win]          # strict: P = 2T is NOT dry
  runs <- rle(dry)
  csec <- runs$lengths[runs$values]
  has_dry_run <- length(csec) > 0 && max(csec) >= 2

  if (!has_dry_run) return(FALSE)     # no summer drought → Temperate

  # ── Step 2: ombrothermic confirmation ────────────────────
  # At least one index must signal summer aridity.
  if (Ios2 < 2 || Iosc4 < 2) return(TRUE)

  return(FALSE)                        # both indices ≥ 2 → Temperate
}


## ── Original version (non-strict, tiebreaker-only) ───────────
.is_med_orig <- function(Ios2, Iosc4, T, P) {
  if (Ios2 <  2 && Iosc4 <  2) return(TRUE)
  if (Ios2 >= 2 && Iosc4 >= 2) return(FALSE)
  win  <- c(5, 6, 7, 8, 9)
  dry  <- P[win] <= 2 * T[win]
  runs <- rle(dry)
  csec <- runs$lengths[runs$values]
  length(csec) > 0 && max(csec) >= 2
}


TERMO_LABELS <- c(
  "01_Lim" = "Inframediterrânico inf.",
  "02_Uim" = "Inframediterrânico sup.",
  "03_Ltm" = "Termomediterrânico inf.",
  "04_Utm" = "Termomediterrânico sup.",
  "05_Lmm" = "Mesomediterrânico inf.",
  "06_Umm" = "Mesomediterrânico sup.",
  "07_Lsm" = "Supramediterrânico inf.",
  "08_Usm" = "Supramediterrânico sup.",
  "09_Lom" = "Oromediterrânico inf.",
  "10_Uom" = "Oromediterrânico sup.",
  "11_Lcm" = "Crioromediterrânico",
  "12_Ltt" = "Termotemperado inf.",
  "13_Utt" = "Termotemperado sup.",
  "14_Lmt" = "Mesotemperado inf.",
  "15_Umt" = "Mesotemperado sup.",
  "16_Lst" = "Supratemperado inf.",
  "17_Ust" = "Supratemperado sup.",
  "18_Lot" = "Orotemperado inf.",
  "19_Uot" = "Orotemperado sup."
)

.assign_thermotype <- function(med, Itc, Tp, use_Tp) {
  if (med) {
    if (!use_Tp) {
      if      (Itc >= 515) "01_Lim"
      else if (Itc >= 450) "02_Uim"
      else if (Itc >= 400) "03_Ltm"
      else if (Itc >= 350) "04_Utm"
      else if (Itc >= 280) "05_Lmm"
      else if (Itc >= 220) "06_Umm"
      else if (Itc >= 150) "07_Lsm"
      else                 "08_Usm"
    } else {
      if      (Tp > 2650) "01_Lim"
      else if (Tp > 2450) "02_Uim"
      else if (Tp > 2300) "03_Ltm"
      else if (Tp > 2150) "04_Utm"
      else if (Tp > 1825) "05_Lmm"
      else if (Tp > 1500) "06_Umm"
      else if (Tp > 1200) "07_Lsm"
      else if (Tp >  900) "08_Usm"
      else if (Tp >  675) "09_Lom"
      else if (Tp >  450) "10_Uom"
      else                "11_Lcm"
    }
  } else {
    if (!use_Tp) {
      if      (Itc >  410) "12_Ltt"
      else if (Itc >= 350) "12_Ltt"
      else if (Itc >= 290) "13_Utt"
      else if (Itc >= 240) "14_Lmt"
      else if (Itc >= 190) "15_Umt"
      else if (Itc >= 120) "16_Lst"
      else                 "17_Ust"
    } else {
      if      (Tp > 2175) "12_Ltt"
      else if (Tp > 2000) "13_Utt"
      else if (Tp > 1700) "14_Lmt"
      else if (Tp > 1400) "15_Umt"
      else if (Tp > 1100) "16_Lst"
      else if (Tp >  800) "17_Ust"
      else if (Tp >  590) "18_Lot"
      else                "19_Uot"
    }
  }
}

.classify_pixel <- function(P, T, Tmin, Tmax, is_med_fn) {
  T_annual <- mean(T)
  cold_m   <- which.min(T)
  hot_m    <- which.max(T)
  m_val    <- Tmin[cold_m]
  M_val    <- Tmax[cold_m]

  Ic   <- T[hot_m] - T[cold_m]
  It   <- (T_annual + m_val + M_val) * 10
  Itc  <- It + .itc_c(Ic)
  Tp   <- sum(T[T > 0]) * 10
  Pp   <- sum(P[T > 0])

  Ios2  <- .ios2(T, P)
  Iosc4 <- .iosc4(T, P)
  med   <- is_med_fn(Ios2, Iosc4, T, P)

  use_Tp <- Ic >= 21 || Itc < 120
  .assign_thermotype(med, Itc, Tp, use_Tp)
}


# ── 5. Pixel-wise classification (both versions) ──────────────
message("Classifying pixels (v3 + original) ...")

P_m    <- as.matrix(r_prec)
T_m    <- as.matrix(r_tavg)
Tmin_m <- as.matrix(r_tmin)
Tmax_m <- as.matrix(r_tmax)

n_cells  <- nrow(P_m)
key_v3   <- character(n_cells)
key_orig <- character(n_cells)

for (i in seq_len(n_cells)) {
  T <- T_m[i, ]
  if (any(is.na(T))) {
    key_v3[i]   <- NA_character_
    key_orig[i] <- NA_character_
    next
  }
  P    <- P_m[i, ]
  Tmin <- Tmin_m[i, ]
  Tmax <- Tmax_m[i, ]
  key_v3[i]   <- .classify_pixel(P, T, Tmin, Tmax, .is_med)
  key_orig[i] <- .classify_pixel(P, T, Tmin, Tmax, .is_med_orig)
}
message("  Done.")


# ── 5b. DIAGNOSTIC ────────────────────────────────────────────
message("\n========== DIAGNOSTIC: v3 vs original ==========")

land_mask <- !is.na(key_v3)
n_land    <- sum(land_mask)
changed   <- land_mask & (key_v3 != key_orig)
n_changed <- sum(changed)

message(sprintf("  Land pixels total  : %d", n_land))
message(sprintf("  Pixels changed     : %d  (%.2f %% of land)",
                n_changed, 100 * n_changed / n_land))

if (n_changed > 0) {
  from_vec <- key_orig[changed]
  to_vec   <- key_v3[changed]
  trans_df <- as.data.frame(table(From = from_vec, To = to_vec),
                             stringsAsFactors = FALSE)
  trans_df <- trans_df[trans_df$Freq > 0, ]
  trans_df$From_label <- TERMO_LABELS[trans_df$From]
  trans_df$To_label   <- TERMO_LABELS[trans_df$To]
  trans_df <- trans_df[order(-trans_df$Freq), ]

  message("\n  Transition table (original → v3):")
  message(sprintf("  %-32s  →  %-32s  %s", "From", "To", "n pixels"))
  message(paste(rep("-", 80), collapse = ""))
  for (r in seq_len(nrow(trans_df))) {
    message(sprintf("  %-32s  →  %-32s  %d",
                    trans_df$From_label[r],
                    trans_df$To_label[r],
                    trans_df$Freq[r]))
  }

  med_codes  <- names(TERMO_LABELS)[1:11]
  temp_codes <- names(TERMO_LABELS)[12:19]
  unexpected <- changed &
    !(key_orig %in% med_codes & key_v3 %in% temp_codes)
  n_unexp <- sum(unexpected, na.rm = TRUE)
  if (n_unexp == 0) {
    message("\n  [OK] All changes are Mediterranean → Temperate.")
  } else {
    message(sprintf("\n  [WARN] %d pixel(s) changed in an unexpected direction!", n_unexp))
    # Print unexpected cases for inspection
    unexp_from <- key_orig[unexpected]
    unexp_to   <- key_v3[unexpected]
    for (j in seq_along(unexp_from))
      message(sprintf("         %s → %s",
                      TERMO_LABELS[unexp_from[j]],
                      TERMO_LABELS[unexp_to[j]]))
  }

  # Summary by macrobioclimate
  orig_med  <- sum(key_orig[land_mask] %in% med_codes,  na.rm = TRUE)
  orig_temp <- sum(key_orig[land_mask] %in% temp_codes, na.rm = TRUE)
  v3_med    <- sum(key_v3[land_mask]   %in% med_codes,  na.rm = TRUE)
  v3_temp   <- sum(key_v3[land_mask]   %in% temp_codes, na.rm = TRUE)

  message("\n  Macrobioclimate pixel counts:")
  message(sprintf("  %-20s  original: %5d  →  v3: %5d  (Δ %+d)",
                  "Mediterranean",  orig_med,  v3_med,  v3_med  - orig_med))
  message(sprintf("  %-20s  original: %5d  →  v3: %5d  (Δ %+d)",
                  "Temperate",      orig_temp, v3_temp, v3_temp - orig_temp))
}
message("=================================================\n")


# ── 6. Build categorical SpatRaster (v3) ──────────────────────
present_keys   <- sort(intersect(unique(key_v3[!is.na(key_v3)]),
                                 names(TERMO_LABELS)))
present_labels <- TERMO_LABELS[present_keys]

code_vec <- rep(NA_integer_, n_cells)
for (k in seq_along(present_keys))
  code_vec[key_v3 == present_keys[k]] <- k

r_template   <- r_tavg[[1]]
r_termo      <- r_template
values(r_termo) <- code_vec

lev_df <- data.frame(value = seq_along(present_labels),
                     label = unname(present_labels))
levels(r_termo) <- lev_df

all_colors <- c(
  "01_Lim" = "#990000", "02_Uim" = "#CC2200", "03_Ltm" = "#EE4400",
  "04_Utm" = "#FF6600", "05_Lmm" = "#FF9900", "06_Umm" = "#FFCC00",
  "07_Lsm" = "#99CC00", "08_Usm" = "#33AA33", "09_Lom" = "#009999",
  "10_Uom" = "#006699", "11_Lcm" = "#003399", "12_Ltt" = "#6699CC",
  "13_Utt" = "#5577AA", "14_Lmt" = "#446699", "15_Umt" = "#334466",
  "16_Lst" = "#553366", "17_Ust" = "#7733AA", "18_Lot" = "#9922CC",
  "19_Uot" = "#BB00EE"
)
present_colors <- all_colors[present_keys]
coltab(r_termo) <- data.frame(value = seq_along(present_colors),
                               col   = unname(present_colors))
r_termo <- mask(r_termo, pt_land)


# ── 7. Plot ───────────────────────────────────────────────────
message("Plotting ...")

out_file <- "portugal_thermotypes.png"
png(out_file, width = 1800, height = 2200, res = 200, bg = "#E8F4F8")

layout(matrix(c(1, 2), nrow = 1), widths = c(3.2, 1.8))

## Panel 1: map
par(mar = c(3, 3, 4, 0.5))
plot(r_termo, legend = FALSE, axes = TRUE, box = FALSE,
     col = present_colors, main = "", xlab = "", ylab = "")

axis(1, at = seq(-9.5, -6.5, by = 1),
     labels = paste0(seq(-9.5, -6.5, by = 1), "°"),
     cex.axis = 0.65, col.axis = "gray30", tck = -0.015)
axis(2, at = seq(37, 42, by = 1),
     labels = paste0(seq(37, 42, by = 1), "°"),
     cex.axis = 0.65, col.axis = "gray30", tck = -0.015, las = 1)

lines(pt_cont, col = "gray50", lwd = 0.6)
lines(pt_land, col = "gray10", lwd = 1.4)

mtext("Termotipos de Portugal Continental",
      side = 3, line = 2.2, cex = 1.15, font = 2, col = "gray10")
mtext("Rivas-Martínez (2004)  ·  WorldClim v2.1  ·  2.5 arc-min  [v3: Walter gate]",
      side = 3, line = 0.9, cex = 0.72, col = "gray35")

sb_x0 <- -9.4; sb_x1 <- -8.25; sb_y <- 37.15
rect(sb_x0, sb_y - 0.04, sb_x1, sb_y + 0.04, col = "gray10", border = NA)
text((sb_x0 + sb_x1) / 2, sb_y - 0.18, "~100 km", cex = 0.55, col = "gray20")

arrows(-6.35, 41.3, -6.35, 41.9, length = 0.12, lwd = 1.5, col = "gray20")
text(-6.35, 41.95, "N", cex = 0.8, font = 2, col = "gray20")

## Panel 2: legend
par(mar = c(3, 0, 4, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

n_lev <- length(present_labels)
y_top <- 0.96
y_stp <- 0.80 / n_lev

med_keys  <- present_keys[startsWith(present_keys, "0") | present_keys == "11_Lcm"]
temp_keys <- present_keys[!present_keys %in% med_keys]

draw_legend_section <- function(keys, title, y_start) {
  if (length(keys) == 0) return(y_start)
  text(0.05, y_start, title, adj = 0, cex = 0.68, font = 2, col = "gray15")
  y_start <- y_start - y_stp * 0.7
  for (k in keys) {
    idx <- which(present_keys == k)
    col <- present_colors[idx]
    lbl <- unname(TERMO_LABELS[k])
    rect(0.05, y_start - y_stp * 0.35, 0.22, y_start + y_stp * 0.35,
         col = col, border = "white", lwd = 0.5)
    text(0.26, y_start, lbl, adj = 0, cex = 0.60, col = "gray10")
    y_start <- y_start - y_stp
  }
  y_start
}

y <- draw_legend_section(med_keys,  "Macrobioclima Mediterrânico", y_top)
y <- draw_legend_section(temp_keys, "Macrobioclima Temperado",     y - y_stp * 0.5)

text(0.05, 0.04,
     paste0("Dados: WorldClim v2.1 (", WC_RES, " arc-min)\n",
            "Walter gate: ≥2 consec. meses P < 2T (maio–set)\n",
            "Confirmação: Ios2 < 2 ou Iosc4 < 2"),
     adj = 0, cex = 0.50, col = "gray45")

dev.off()
message("Map saved to: ", normalizePath(out_file))
