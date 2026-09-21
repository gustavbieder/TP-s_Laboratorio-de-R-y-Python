# ==============================================================================
# 01_import_transform.R
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# Objetivo : importar la base mensual (IPC_imp y TCN Gs./USD), verificar su
#            calidad y aplicar (1) logaritmo natural y (2) primera diferencia
#            del logaritmo.
# Entrada  : data/raw/Base_Arimax.xlsx   (hoja "Hoja1": fecha, IPC, TCN, IPC_imp)
# Salida   : data/processed/base_transformada.csv
# Uso      : ejecutar con el directorio de trabajo en la raíz del repositorio
#            (por ejemplo, abriendo el proyecto en RStudio).
# Nota     : el archivo también trae la columna IPC (índice general), que este
#            trabajo no utiliza; se descarta al importar.
# ==============================================================================

library(readxl)
library(dplyr)
library(readr)
library(lubridate)

# ---- 0. Rutas relativas -------------------------------------------------------
ruta_raw  <- file.path("data", "raw", "Base_Arimax.xlsx")
ruta_proc <- file.path("data", "processed", "base_transformada.csv")

if (!file.exists(ruta_raw)) {
  stop("No se encuentra ", ruta_raw,
       ". Ejecutá el script desde la raíz del repositorio.")
}
dir.create(dirname(ruta_proc), showWarnings = FALSE, recursive = TRUE)

# ---- 1. Importación -----------------------------------------------------------
# La primera columna (fecha) se selecciona por posición para no depender de los
# acentos del encabezado original ("Años").
crudo <- read_excel(ruta_raw, sheet = "Hoja1")
stopifnot("faltan las columnas IPC_imp y/o TCN" = all(c("IPC_imp", "TCN") %in% names(crudo)))

base <- crudo |>
  select(fecha = 1, ipc_imp = IPC_imp, tcn = TCN) |>
  mutate(fecha = as.Date(fecha))

# ---- 2. Chequeos de integridad ------------------------------------------------
# Condiciones necesarias para poder aplicar log y diferencias sin distorsión.
indice_mes <- year(base$fecha) * 12 + month(base$fecha)

stopifnot(
  "faltan valores (NA)"          = !anyNA(base),
  "IPC_imp o TCN no positivos"   = all(base$ipc_imp > 0) && all(base$tcn > 0),
  "hay fechas duplicadas"        = !anyDuplicated(base$fecha),
  "fechas fuera de orden"        = !is.unsorted(base$fecha),
  "hay meses faltantes (huecos)" = all(diff(indice_mes) == 1)
)

# ---- 3. Transformaciones ------------------------------------------------------
# ln_*  : logaritmo natural del nivel.
# dln_* : primera diferencia del logaritmo, ln(x_t) - ln(x_{t-1}).
#         Aproxima la variación porcentual mensual en formato decimal
#         (0.01 ~ 1 %). La primera observación queda en NA por construcción.
base_transformada <- base |>
  mutate(
    ln_ipc_imp  = log(ipc_imp),
    ln_tcn      = log(tcn),
    dln_ipc_imp = ln_ipc_imp - dplyr::lag(ln_ipc_imp),
    dln_tcn     = ln_tcn - dplyr::lag(ln_tcn)
  )

# ---- 4. Verificación rápida ---------------------------------------------------
cat("Período  :", format(min(base_transformada$fecha)), "a",
    format(max(base_transformada$fecha)), "\n")
cat("Meses    :", nrow(base_transformada), "\n")
cat("NA en dln:", sum(is.na(base_transformada$dln_ipc_imp)),
    "(IPC_imp) y", sum(is.na(base_transformada$dln_tcn)),
    "(TCN) -> solo el primer mes\n")
print(summary(select(base_transformada, dln_ipc_imp, dln_tcn)))

# ---- 5. Guardado --------------------------------------------------------------
write_csv(base_transformada, ruta_proc)
cat("\nArchivo guardado en:", ruta_proc, "\n")
