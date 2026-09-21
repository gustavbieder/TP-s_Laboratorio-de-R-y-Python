# ==============================================================================
# 03_seleccion_arimax.R  -  Parte 3 (1/2): búsqueda del mejor ARIMAX
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# Se trabaja SOLO con la muestra de entrenamiento (2012-01 a 2022-12); 2023 en
# adelante queda reservado para la evaluación final (script 04).
#
# ETAPA 1 - Grilla y filtros
#   Se estiman todas las combinaciones de:
#     - errores ARMA(p,q) regulares, p y q entre 0 y 2
#     - estacionalidad: nada, SAR(1), SAR(2), SMA(1), SARMA(1,1) o 11 dummies mensuales
#     - rezagos de d_tcn: ninguno, {0..k} con k = 0..6 ("condicional") y
#       {1..k} con k = 1..6 ("ex-ante", no usa el TCN del mes que se pronostica)
#   Un modelo es ADMISIBLE si (i) sus residuos pasan Ljung-Box (rezagos 12 y 24,
#   p > 0,05) y (ii) sus raíces AR y MA están a más de 5 % de la circunferencia
#   unitaria (descarta soluciones casi no estacionarias, inestables al reestimar).
#
# ETAPA 2 - Validación temporal dentro del entrenamiento
#   Los 8 mejores modelos ex-ante y los 8 mejores condicionales por AICc (más los 3
#   mejores sin d_tcn y el mejor por BIC) se
#   comparan por su error de pronóstico a un paso con ventana expansiva sobre
#   2018-01 a 2022-12 (60 pronósticos). El "mejor ARIMAX" es el de menor RMSE
#   en esa validación; el AICc sirve para preseleccionar, no para decidir.
#
# Salida: output/tables/sel_*.csv, output/figures/sel_*.png, output/seleccion.rds
# ==============================================================================

source("R/00_config.R", encoding = "UTF-8")
source("R/utils_modelos.R", encoding = "UTF-8")

d  <- cargar_datos()
tr <- d |> filter(fecha >= FECHA_INICIO, fecha <= FECHA_FIN_TRAIN)
cat("Observaciones de entrenamiento:", nrow(tr), "(", format(min(tr$fecha)), "a", format(max(tr$fecha)), ")\n")

# ---- 1. Grilla de especificaciones -------------------------------------------
conjuntos_lags <- c(list(integer(0)),
                    lapply(0:MAX_LAG_X, function(k) 0:k),   # condicional (incluye rezago 0)
                    lapply(1:MAX_LAG_X, function(k) 1:k))   # ex-ante (rezagos >= 1)

estacionalidad <- tribble(
  ~P, ~Q, ~dummies, ~est_txt,
  0,  0,  FALSE,    "sin componente estacional",
  1,  0,  FALSE,    "SAR(1)",
  2,  0,  FALSE,    "SAR(2)",
  0,  1,  FALSE,    "SMA(1)",
  1,  1,  FALSE,    "SARMA(1,1)",
  0,  0,  TRUE,     "dummies mensuales"
)

grilla <- expand_grid(p = 0:2, q = 0:2,
                      est_id = seq_len(nrow(estacionalidad)),
                      lag_id = seq_along(conjuntos_lags))
cat("Modelos a estimar:", nrow(grilla), "\n")

evaluar <- function(i) {
  g <- grilla[i, ]
  e <- estacionalidad[g$est_id, ]
  s <- nueva_spec(g$p, g$q, e$P, e$Q, conjuntos_lags[[g$lag_id]], e$dummies)
  f <- ajustar_spec(tr, s)
  if (is.null(f) || !is.finite(f$aicc)) return(NULL)
  res <- residuals(f)
  fitdf <- s$p + s$q + s$P + s$Q
  tibble(p = s$p, q = s$q, P = s$P, Q = s$Q, dummies = s$dummies,
         estacionalidad = e$est_txt,
         lags = if (length(s$lags)) paste0(min(s$lags), "-", max(s$lags)) else "ninguno",
         rezago_max = if (length(s$lags)) max(s$lags) else 0L,
         tipo = tipo_x(s$lags),
         k_parametros = length(f$coef) + 1,
         aic = f$aic, aicc = f$aicc, bic = f$bic, sigma2 = f$sigma2,
         lb12_p = Box.test(res, 12, "Ljung-Box", fitdf = fitdf)$p.value,
         lb24_p = Box.test(res, 24, "Ljung-Box", fitdf = fitdf)$p.value,
         raiz_min = raiz_minima(f))
}

# ---- 2. Etapa 1: estimación de toda la grilla --------------------------------
t0 <- Sys.time()
nucleos <- max(1, parallel::detectCores() - 1)
resultados <- bind_rows(parallel::mclapply(seq_len(nrow(grilla)), evaluar, mc.cores = nucleos))
cat("Estimados:", nrow(resultados), "de", nrow(grilla), "en",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

resultados <- resultados |>
  mutate(residuos_ok = lb12_p > 0.05 & lb24_p > 0.05,
         raices_ok   = raiz_min >= RAIZ_MIN,
         admisible   = residuos_ok & raices_ok,
         etiqueta = paste0("ARMA(", p, ",", q, ")(", P, ",", Q, ")[12]",
                           if_else(dummies, "+dum", ""), " | X: ", lags)) |>
  arrange(aicc)
guardar_tab(resultados |> mutate(across(where(is.double), ~ round(.x, 4))), "sel_grilla_completa.csv")

resumen_filtro <- resultados |> group_by(tipo) |>
  summarise(estimados = n(), residuos_ok = sum(residuos_ok), raices_ok = sum(raices_ok),
            admisibles = sum(admisible), .groups = "drop")
guardar_tab(resumen_filtro, "sel_resumen_filtro.csv")
print(as.data.frame(resumen_filtro))

spec_desde_fila <- function(r) {
  nueva_spec(r$p, r$q, r$P, r$Q,
             if (r$lags == "ninguno") integer(0) else eval(parse(text = sub("-", ":", r$lags))),
             r$dummies)
}

# Mejores por AICc (criterio informativo)
mejor_aicc <- function(df) df |> filter(admisible) |> slice_min(aicc, n = 1, with_ties = FALSE)
sel_aicc <- list(
  arimax_aicc = mejor_aicc(filter(resultados, tipo != "sin_X")),
  ex_ante     = mejor_aicc(filter(resultados, tipo == "ex_ante")),
  condicional = mejor_aicc(filter(resultados, tipo == "condicional")),
  sin_x       = mejor_aicc(filter(resultados, tipo == "sin_X")),
  arimax_bic  = filter(resultados, admisible, tipo != "sin_X") |> slice_min(bic, n = 1, with_ties = FALSE)
)
guardar_tab(bind_rows(sel_aicc, .id = "criterio") |>
              select(criterio, etiqueta, tipo, k_parametros, aic, aicc, bic, lb12_p, lb24_p, raiz_min) |>
              mutate(across(where(is.double), ~ round(.x, 3))), "sel_mejores_por_aicc_bic.csv")

# ---- 3. Etapa 2: validación temporal de los preseleccionados -----------------
short_ex  <- resultados |> filter(admisible, tipo == "ex_ante")     |> slice_min(aicc, n = 8, with_ties = FALSE)
short_co  <- resultados |> filter(admisible, tipo == "condicional") |> slice_min(aicc, n = 8, with_ties = FALSE)
short_sin <- resultados |> filter(admisible, tipo == "sin_X")       |> slice_min(aicc, n = 3, with_ties = FALSE)
short <- bind_rows(short_ex, short_co, short_sin, sel_aicc$arimax_bic) |> distinct(etiqueta, .keep_all = TRUE)

fechas_cv <- tr |> filter(fecha >= INICIO_VALIDACION) |> pull(fecha)
cat("Validación temporal:", nrow(short), "modelos x", length(fechas_cv), "pronósticos\n")

validar <- function(i) {
  r <- short[i, ]
  pr <- pronostico_rolling(d, spec_desde_fila(r), fechas_cv) |>
    left_join(select(tr, fecha, obs = dln_ipc_imp), by = "fecha")
  e <- pr$obs - pr$pronostico
  tibble(etiqueta = r$etiqueta, tipo = r$tipo, k_parametros = r$k_parametros,
         aicc = r$aicc, bic = r$bic,
         rmse_cv = sqrt(mean(e^2, na.rm = TRUE)) * 100,
         mae_cv  = mean(abs(e), na.rm = TRUE) * 100,
         n_cv = sum(!is.na(e)))
}
t1 <- Sys.time()
cv <- bind_rows(parallel::mclapply(seq_len(nrow(short)), validar, mc.cores = nucleos)) |>
  arrange(rmse_cv)
cat("Validación temporal terminada en", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min\n")
guardar_tab(cv |> mutate(across(where(is.double), ~ round(.x, 4))), "sel_validacion_temporal.csv")
print(as.data.frame(cv |> mutate(across(where(is.double), ~ round(.x, 3)))))

# ---- 4. Modelos elegidos ------------------------------------------------------
mejor_cv <- function(tipos) {
  cv |> filter(tipo %in% tipos, n_cv >= 55) |> slice_min(rmse_cv, n = 1, with_ties = FALSE)
}
elegidos <- list(
  arimax_principal = mejor_cv(c("ex_ante", "condicional")),
  ex_ante          = mejor_cv("ex_ante"),
  condicional      = mejor_cv("condicional"),
  sin_x            = mejor_cv("sin_X")
)
fila <- function(et) resultados |> filter(etiqueta == et) |> slice(1)
specs <- c(map(elegidos, ~ spec_desde_fila(fila(.x$etiqueta))),
           list(arimax_aicc = spec_desde_fila(sel_aicc$arimax_aicc),
                arimax_bic  = spec_desde_fila(sel_aicc$arimax_bic)))

tabla_elegidos <- bind_rows(elegidos, .id = "rol") |>
  select(rol, etiqueta, tipo, k_parametros, aicc, bic, rmse_cv, mae_cv) |>
  mutate(across(where(is.double), ~ round(.x, 3)))
guardar_tab(tabla_elegidos, "sel_modelos_elegidos.csv")
print(as.data.frame(tabla_elegidos))

top15 <- resultados |> filter(admisible) |> slice_min(aicc, n = 15, with_ties = FALSE) |>
  mutate(daicc = aicc - min(aicc)) |>
  left_join(select(cv, etiqueta, rmse_cv), by = "etiqueta")
guardar_tab(top15 |> select(etiqueta, tipo, k_parametros, aicc, bic, daicc, rmse_cv) |>
              mutate(across(where(is.double), ~ round(.x, 3))), "sel_top15_aicc.csv")

saveRDS(list(specs = specs, elegidos = tabla_elegidos, cv = cv, grilla = resultados),
        file.path("output", "seleccion.rds"))

# ---- 5. Figuras ---------------------------------------------------------------
# (a) Cuánto mejora el ajuste al sumar rezagos del TCN
mejor_por_rezago <- resultados |> filter(admisible, tipo != "sin_X") |>
  group_by(tipo, rezago_max) |> slice_min(aicc, n = 1, with_ties = FALSE) |> ungroup()
aicc_sinx <- sel_aicc$sin_x$aicc

p_aicc <- ggplot(mejor_por_rezago, aes(rezago_max, aicc, colour = tipo)) +
  geom_hline(yintercept = aicc_sinx, linetype = "dashed", colour = "grey30") +
  annotate("text", x = MAX_LAG_X, y = aicc_sinx, label = "Mejor modelo sin d_tcn", vjust = -0.8, hjust = 1, size = 3) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_colour_manual(values = c(condicional = col_tcn, ex_ante = col_ipc),
                      labels = c(condicional = "Con d_tcn contemporáneo (rezagos 0 a k)",
                                 ex_ante = "Solo rezagos (1 a k)"), name = NULL) +
  scale_x_continuous(breaks = 0:MAX_LAG_X) +
  labs(title = "Ajuste del mejor modelo admisible según el rezago máximo de d_tcn",
       subtitle = "Menor AICc = mejor ajuste penalizado. Muestra de entrenamiento 2012-2022",
       x = "Rezago máximo k incluido", y = "AICc", caption = "Elaboración propia.") +
  theme(legend.position = "bottom")
guardar_fig(p_aicc, "sel_01_aicc_por_rezagos.png", 8.5, 4.8)

# (b) Error de validación temporal de los preseleccionados
p_cv <- cv |>
  mutate(etiqueta_corta = str_wrap(etiqueta, 34),
         elegido = etiqueta == elegidos$arimax_principal$etiqueta) |>
  ggplot(aes(rmse_cv, reorder(etiqueta_corta, -rmse_cv), colour = tipo, shape = elegido)) +
  geom_segment(aes(x = min(rmse_cv) * 0.98, xend = rmse_cv, yend = reorder(etiqueta_corta, -rmse_cv)), colour = "grey85") +
  geom_point(size = 3) +
  scale_colour_manual(values = c(condicional = col_tcn, ex_ante = col_ipc, sin_X = col_gris),
                      labels = c(condicional = "Con d_tcn contemporáneo", ex_ante = "Solo rezagos", sin_X = "Sin d_tcn"),
                      name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  labs(title = "Validación temporal de los modelos preseleccionados",
       subtitle = "RMSE de pronóstico a un mes, 2018-2022, reestimando cada mes (rombo = modelo elegido)",
       x = "RMSE (puntos porcentuales de inflación mensual)", y = NULL, caption = "Elaboración propia.") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
guardar_fig(p_cv, "sel_02_validacion_temporal.png", 9, 6.5)

cat("Selección terminada.\n")
