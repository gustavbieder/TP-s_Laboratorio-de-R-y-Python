# ==============================================================================
# run_all.R  -  Reproduce todo el proyecto de punta a punta
# Ejecutar desde la raíz del repositorio:  source("run_all.R")   o   Rscript run_all.R
# Autor    : Gustavo Santiago Biedermann Giménez
# ==============================================================================
scripts <- c("R/01_import_transform.R",
             "R/02_eda.R",
             "R/03_seleccion_arimax.R",
             "R/04_evaluacion_resultados.R")
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n")
  source(s, encoding = "UTF-8", echo = FALSE)
}
cat("\nListo. Figuras en output/figures, tablas en output/tables.\n")
