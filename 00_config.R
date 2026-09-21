# ==============================================================================
# 00_config.R  -  Configuración común a todos los scripts
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# Carga paquetes, fija la semilla, define rutas relativas, parámetros del
# análisis, un tema gráfico común y funciones auxiliares.
# Los demás scripts lo cargan con source("R/00_config.R").
# ==============================================================================

# ---- Paquetes (se instalan solos si faltan) -----------------------------------
paquetes <- c("tidyverse", "readxl", "lubridate", "forecast", "tseries",
              "lmtest", "strucchange", "patchwork", "scales", "zoo")
faltan <- setdiff(paquetes, rownames(installed.packages()))
if (length(faltan) > 0) install.packages(faltan, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(paquetes, library, character.only = TRUE)))

# ---- Semilla ------------------------------------------------------------------
set.seed(2026)

# ---- Rutas relativas (ejecutar desde la raíz del repositorio) -----------------
dir_raw  <- file.path("data", "raw")
dir_proc <- file.path("data", "processed")
dir_fig  <- file.path("output", "figures")
dir_tab  <- file.path("output", "tables")
for (d in c(dir_proc, dir_fig, dir_tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Parámetros del análisis --------------------------------------------------
FECHA_INICIO    <- as.Date("2012-01-01")  # primer año completo con metas de inflación (desde mayo-2011)
FECHA_FIN_TRAIN <- as.Date("2022-12-01")  # último mes de estimación; 2023-01 en adelante = evaluación
FECHA_METAS_INF <- as.Date("2011-05-01")  # anuncio formal del esquema de metas de inflación (BCP)
MAX_LAG_X       <- 6                      # máximo rezago de dln_tcn considerado en la búsqueda
UMBRAL_Z        <- 3.5                    # umbral del z-score robusto (mediana/MAD) para marcar atípicos
RAIZ_MIN        <- 1.05                   # módulo mínimo de las raíces AR y MA (evita soluciones casi no estacionarias)
INICIO_VALIDACION <- as.Date("2018-01-01")  # primer mes de la validación cruzada temporal dentro del entrenamiento

# ---- Paleta y tema ------------------------------------------------------------
col_ipc  <- "#1F4E79"
col_tcn  <- "#C0504D"
col_gris <- "grey55"
col_ok   <- "#2E7D32"

tema_tp <- theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey30"),
        plot.caption  = element_text(colour = "grey40", size = 8),
        panel.grid.minor = element_blank(),
        strip.text    = element_text(face = "bold"))
theme_set(tema_tp)

# ---- Funciones auxiliares -----------------------------------------------------
guardar_fig <- function(p, nombre, ancho = 9, alto = 5) {
  ggsave(file.path(dir_fig, nombre), plot = p, width = ancho, height = alto,
         dpi = 200, bg = "white")
}

guardar_tab <- function(df, nombre) {
  write_csv(df, file.path(dir_tab, nombre))
}

# Carga la base transformada y construye los rezagos de dln_tcn ANTES de recortar
# la muestra, de modo que el primer mes de la ventana ya tenga sus rezagos.
cargar_datos <- function() {
  ruta <- file.path(dir_proc, "base_transformada.csv")
  if (!file.exists(ruta)) stop("Falta ", ruta, ". Ejecutá primero R/01_import_transform.R")
  d <- read_csv(ruta, show_col_types = FALSE)
  for (k in 1:12) d[[paste0("dln_tcn_l", k)]] <- dplyr::lag(d$dln_tcn, k)
  d$dln_tcn_l0 <- d$dln_tcn
  d$mes <- factor(month(d$fecha), levels = 1:12)
  d
}

# z-score robusto (mediana y MAD) para detectar atípicos
z_robusto <- function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)

# Convierte un vector a objeto ts mensual
a_ts <- function(x, fecha_ini) {
  ts(x, start = c(year(fecha_ini), month(fecha_ini)), frequency = 12)
}
