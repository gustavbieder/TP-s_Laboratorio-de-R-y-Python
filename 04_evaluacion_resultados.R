# ==============================================================================
# 04_evaluacion_resultados.R  -  Parte 3 (2/2): evaluación, modelo final e interpretación
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# 1) Pronóstico fuera de muestra a un paso (ventana expansiva, 2023-01 en adelante)
#    de los modelos elegidos en 03_seleccion_arimax.R frente a alternativas simples,
#    con test de Diebold-Mariano.
# 2) Estimación del ARIMAX principal en toda la ventana 2012-2026: coeficientes,
#    aporte de d_tcn (razón de verosimilitud) y diagnósticos de residuos.
# 3) Traspaso del tipo de cambio a precios: efecto por rezago y acumulado, con el
#    modelo elegido por AICc y con un modelo de rezagos 0 a 6 sin restricciones.
# 4) Estabilidad y robustez: cambio del traspaso desde 2020, otra muestra, atípicos.
# ==============================================================================

source("R/00_config.R", encoding = "UTF-8")
source("R/utils_modelos.R", encoding = "UTF-8")

d     <- cargar_datos()
sel   <- readRDS(file.path("output", "seleccion.rds"))
specs <- sel$specs
ventana     <- d |> filter(fecha >= FECHA_INICIO)
fechas_test <- ventana |> filter(fecha > FECHA_FIN_TRAIN) |> pull(fecha)
cat("Evaluación fuera de muestra:", format(min(fechas_test)), "a", format(max(fechas_test)),
    "(", length(fechas_test), "meses )\n")

spec_final    <- specs$arimax_principal   # mejor ARIMAX por validación temporal
spec_dinamico <- specs$arimax_aicc        # mejor ARIMAX por AICc
spec_largo    <- spec_dinamico; spec_largo$lags <- 0:MAX_LAG_X   # rezagos 0 a 6 sin restricciones
cat("ARIMAX principal      :", etiqueta_spec(spec_final), "\n")
cat("ARIMAX por AICc       :", etiqueta_spec(spec_dinamico), "\n")

# ---- 1. Pronóstico fuera de muestra a un paso ---------------------------------
candidatos <- list(
  "Mejor modelo sin d_tcn"                 = specs$sin_x,
  "ARIMAX principal"                       = spec_final,
  "ARIMAX ex-ante (solo rezagos de d_tcn)" = specs$ex_ante,
  "ARIMAX condicional (con d_tcn del mes)" = specs$condicional,
  "ARIMAX elegido por AICc"                = spec_dinamico
)
# Si dos roles comparten la misma especificación se estima una sola vez
etiquetas <- map_chr(candidatos, etiqueta_spec)
modelos_oos <- candidatos[!duplicated(etiquetas)]
tabla_modelos <- tibble(nombre = names(candidatos), especificacion = etiquetas,
                        se_evalua_como = names(modelos_oos)[match(etiquetas, etiquetas[!duplicated(etiquetas)])])
guardar_tab(tabla_modelos, "eval_modelos_evaluados.csv")
print(as.data.frame(tabla_modelos))

pron_modelos <- imap(modelos_oos, function(s, nombre) {
  cat("  pronosticando:", nombre, "-", etiqueta_spec(s), "\n")
  pronostico_rolling(d, s, fechas_test) |> mutate(modelo = nombre)
}) |> bind_rows()

ref <- ventana |>
  mutate(media_exp = map_dbl(fecha, ~ mean(dln_ipc_imp[fecha < .x])),
         naive     = dplyr::lag(dln_ipc_imp, 1),
         naive_est = dplyr::lag(dln_ipc_imp, 12)) |>
  filter(fecha %in% fechas_test) |>
  select(fecha, `Media histórica (expansiva)` = media_exp,
         `Mes anterior (naive)` = naive, `Mismo mes del año anterior` = naive_est) |>
  pivot_longer(-fecha, names_to = "modelo", values_to = "pronostico")

obs <- ventana |> filter(fecha %in% fechas_test) |> select(fecha, observado = dln_ipc_imp)
pron <- bind_rows(pron_modelos, ref) |>
  left_join(obs, by = "fecha") |>
  mutate(error = observado - pronostico)
guardar_tab(pron |> mutate(across(where(is.double), ~ round(.x, 6))), "eval_pronosticos.csv")

orden_modelos <- c("Media histórica (expansiva)", "Mes anterior (naive)", "Mismo mes del año anterior",
                   names(modelos_oos))
mse_media <- pron |> filter(modelo == "Media histórica (expansiva)") |> summarise(m = mean(error^2)) |> pull(m)

metricas <- pron |> group_by(modelo) |>
  summarise(n = sum(!is.na(error)),
            RMSE = sqrt(mean(error^2, na.rm = TRUE)) * 100,
            MAE  = mean(abs(error), na.rm = TRUE) * 100,
            R2_oos_vs_media = 1 - mean(error^2, na.rm = TRUE) / mse_media,
            .groups = "drop") |>
  mutate(modelo = factor(modelo, levels = orden_modelos)) |> arrange(modelo) |>
  mutate(across(where(is.double), ~ round(.x, 4)))
guardar_tab(metricas, "eval_metricas.csv")
print(as.data.frame(metricas))

# Test de Diebold-Mariano (H1: el segundo método es más preciso), con pérdida cuadrática
# (asociada al RMSE) y con pérdida absoluta (asociada al MAE; menos sensible a los saltos extremos)
dm <- function(m_ref, m_nuevo) {
  a <- pron |> filter(modelo == m_ref)   |> select(fecha, e1 = error)
  b <- pron |> filter(modelo == m_nuevo) |> select(fecha, e2 = error)
  j <- inner_join(a, b, by = "fecha") |> filter(!is.na(e1), !is.na(e2))
  prueba <- function(potencia) tryCatch(
    forecast::dm.test(j$e1, j$e2, alternative = "greater", h = 1, power = potencia),
    error = function(err) NULL)
  r2 <- prueba(2); r1 <- prueba(1)
  tibble(referencia = m_ref, modelo = m_nuevo, n = nrow(j),
         estadistico_cuadratica = if (is.null(r2)) NA_real_ else unname(r2$statistic),
         p_valor_cuadratica     = if (is.null(r2)) NA_real_ else r2$p.value,
         estadistico_absoluta   = if (is.null(r1)) NA_real_ else unname(r1$statistic),
         p_valor_absoluta       = if (is.null(r1)) NA_real_ else r1$p.value)
}
nombres_arimax <- setdiff(names(modelos_oos), "Mejor modelo sin d_tcn")
tests_dm <- bind_rows(
  map_dfr(names(modelos_oos), ~ dm("Media histórica (expansiva)", .x)),
  map_dfr(nombres_arimax,     ~ dm("Mejor modelo sin d_tcn", .x))
) |> mutate(across(where(is.double), ~ round(.x, 4)))
guardar_tab(tests_dm, "eval_diebold_mariano.csv")
print(as.data.frame(tests_dm))

# Figura: pronósticos vs observado (modelos seleccionados)
estilos <- tribble(
  ~modelo,                                  ~color,    ~tipo_linea,
  "Media histórica (expansiva)",            "grey60",  "dotted",
  "Mejor modelo sin d_tcn",                 "#7F7F7F", "dashed",
  "ARIMAX principal",                       col_ipc,   "solid",
  "ARIMAX ex-ante (solo rezagos de d_tcn)", col_tcn,   "twodash"
) |> filter(modelo %in% unique(pron$modelo))

p_pron <- pron |>
  filter(modelo %in% estilos$modelo) |>
  mutate(modelo = factor(modelo, levels = estilos$modelo)) |>
  ggplot(aes(fecha, pronostico * 100, colour = modelo, linetype = modelo)) +
  geom_line(data = obs, aes(fecha, observado * 100), inherit.aes = FALSE, colour = "black", linewidth = 1) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = setNames(estilos$color, estilos$modelo), name = NULL) +
  scale_linetype_manual(values = setNames(estilos$tipo_linea, estilos$modelo), name = NULL) +
  labs(title = "Pronóstico a un mes fuera de muestra vs. variación mensual observada del IPC_imp",
       subtitle = "Línea negra: d_ipc_imp observado. Modelos reestimados cada mes con datos hasta el mes previo",
       x = NULL, y = "Variación mensual del IPC_imp (%)", caption = "Elaboración propia con datos del BCP.") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
guardar_fig(p_pron, "eval_01_pronosticos_fuera_muestra.png", 9.5, 5.4)

p_met <- metricas |>
  pivot_longer(c(RMSE, MAE), names_to = "metrica", values_to = "valor") |>
  ggplot(aes(valor, fct_rev(modelo), fill = metrica)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", valor)), position = position_dodge(width = 0.7), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c(RMSE = col_ipc, MAE = col_tcn), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Error de pronóstico fuera de muestra (2023-01 en adelante)",
       subtitle = "En puntos porcentuales de variación mensual del IPC_imp; menor es mejor",
       x = NULL, y = NULL, caption = "Elaboración propia con datos del BCP.") +
  theme(legend.position = "bottom")
guardar_fig(p_met, "eval_02_metricas_fuera_muestra.png", 9.5, 4.8)

# ---- 2. Modelo principal en toda la ventana ------------------------------------
resumen_modelo <- function(fit, spec, sufijo) {
  coefs <- tibble(termino = names(coef(fit)), estimacion = coef(fit),
                  error_est = sqrt(diag(fit$var.coef))) |>
    mutate(z = estimacion / error_est, p_valor = 2 * pnorm(-abs(z)),
           ic95_inf = estimacion - 1.96 * error_est, ic95_sup = estimacion + 1.96 * error_est)
  guardar_tab(coefs |> mutate(across(where(is.double), ~ signif(.x, 5))), paste0("eval_coeficientes_", sufijo, ".csv"))
  # razón de verosimilitud contra el mismo modelo sin d_tcn
  spec0 <- spec; spec0$lags <- integer(0)
  fit0  <- ajustar_spec(ventana, spec0)
  k_x   <- length(spec$lags)
  lr    <- 2 * (as.numeric(fit$loglik) - as.numeric(fit0$loglik))
  ajuste <- tibble(n = fit$nobs, log_verosimilitud = as.numeric(fit$loglik), aic = fit$aic, aicc = fit$aicc,
                   bic = fit$bic, sigma_residual_pct = sqrt(fit$sigma2) * 100,
                   rezagos_x = k_x, estadistico_LR = lr, p_valor_LR = pchisq(lr, df = k_x, lower.tail = FALSE),
                   aicc_sin_x = fit0$aicc)
  guardar_tab(ajuste |> mutate(across(where(is.double), ~ signif(.x, 5))), paste0("eval_ajuste_", sufijo, ".csv"))
  list(coefs = coefs, ajuste = ajuste)
}

fit <- ajustar_spec(ventana, spec_final)
stopifnot(!is.null(fit))
r_final <- resumen_modelo(fit, spec_final, "principal")
print(as.data.frame(r_final$coefs |> mutate(across(where(is.double), ~ round(.x, 4)))))
print(as.data.frame(r_final$ajuste |> mutate(across(where(is.double), ~ round(.x, 4)))))

# Diagnósticos de residuos del modelo principal
fitdf <- spec_final$p + spec_final$q + spec_final$P + spec_final$Q
lb <- function(x, l, f = 0) Box.test(x, lag = l, type = "Ljung-Box", fitdf = f)$p.value
diagnosticos <- function(fit_r) {
  res_r <- as.numeric(residuals(fit_r))
  tibble(
    prueba = c("Ljung-Box de los residuos, rezago 12", "Ljung-Box de los residuos, rezago 24",
               "Ljung-Box de los residuos al cuadrado, rezago 12", "Jarque-Bera (normalidad)",
               "Shapiro-Wilk (normalidad)", "CUSUM (estabilidad)"),
    p_valor = c(lb(res_r, 12, fitdf), lb(res_r, 24, fitdf), lb(res_r^2, 12),
                tseries::jarque.bera.test(res_r)$p.value, shapiro.test(res_r)$p.value,
                strucchange::sctest(strucchange::efp(res_r ~ 1, type = "OLS-CUSUM"))$p.value)
  ) |> mutate(p_valor = round(p_valor, 4))
}
diagn <- diagnosticos(fit)
guardar_tab(diagn, "eval_diagnosticos.csv")
print(as.data.frame(diagn))

figura_diagnostico <- function(fit_r, titulo) {
  res_r <- residuals(fit_r)
  z <- as.numeric(res_r) / sd(as.numeric(res_r))
  df_res <- tibble(fecha = ventana$fecha, res = as.numeric(res_r) * 100, z = z)
  p1 <- ggplot(df_res, aes(fecha, res)) +
    geom_hline(yintercept = c(-2, 2) * sd(df_res$res), linetype = "dashed", colour = col_tcn) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(colour = col_ipc) +
    geom_point(data = filter(df_res, abs(z) > 3), size = 2) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = "Residuos en el tiempo", subtitle = "Puntos: |z| > 3. Líneas punteadas: ± 2 desvíos", x = NULL, y = "%")
  p2 <- ggAcf(as.numeric(res_r), lag.max = 24) + labs(title = "ACF de los residuos", y = NULL)
  p3 <- ggplot(tibble(z = z), aes(sample = z)) + stat_qq(colour = col_ipc, alpha = 0.7) + stat_qq_line() +
    labs(title = "Q-Q de residuos estandarizados", x = "Cuantiles teóricos", y = "Cuantiles muestrales")
  lbs <- tibble(rezago = (fitdf + 1):24) |> mutate(p = map_dbl(rezago, ~ lb(res_r, .x, fitdf)))
  p4 <- ggplot(lbs, aes(rezago, p)) + geom_point(colour = col_ipc, size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = col_tcn) + ylim(0, 1) +
    labs(title = "Ljung-Box según el rezago", subtitle = "Sobre la línea roja: sin autocorrelación",
         x = "Rezago", y = "p-valor")
  (p1 + p2) / (p3 + p4) + plot_annotation(title = titulo, caption = "Elaboración propia.")
}
res_z <- as.numeric(residuals(fit)) / sd(as.numeric(residuals(fit)))
guardar_tab(tibble(fecha = ventana$fecha, res_pct = as.numeric(residuals(fit)) * 100, z = res_z) |>
              filter(abs(z) > 3) |> mutate(across(where(is.double), ~ round(.x, 3))), "eval_residuos_atipicos.csv")
guardar_fig(figura_diagnostico(fit, "Diagnóstico de residuos del ARIMAX principal"),
            "eval_03_diagnosticos_residuos.png", 9.5, 6.5)

# ---- 3. Traspaso dinámico ------------------------------------------------------
# Se compara el modelo preferido por AICc con uno de rezagos 0 a 6 sin restricciones
# (mismos errores ARMA) para ver hasta qué mes se extiende el efecto.
fit_din   <- ajustar_spec(ventana, spec_dinamico)
fit_largo <- ajustar_spec(ventana, spec_largo)
stopifnot(!is.null(fit_din), !is.null(fit_largo))
r_din   <- resumen_modelo(fit_din, spec_dinamico, "aicc")
r_largo <- resumen_modelo(fit_largo, spec_largo, "largo")
print(as.data.frame(r_din$coefs |> mutate(across(where(is.double), ~ round(.x, 4)))))
print(as.data.frame(r_largo$coefs |> mutate(across(where(is.double), ~ round(.x, 4)))))

traspaso_acumulado <- function(fit_r) {
  nom_x <- grep("^tcn_l", names(coef(fit_r)), value = TRUE)
  V <- fit_r$var.coef[nom_x, nom_x, drop = FALSE]; b <- coef(fit_r)[nom_x]
  rez <- as.integer(sub("tcn_l", "", nom_x))
  map_dfr(seq_along(b), function(j) {
    u <- rep(0, length(b)); u[1:j] <- 1
    est <- sum(u * b); se <- sqrt(as.numeric(t(u) %*% V %*% u))
    tibble(rezago = rez[j], acumulado = est, se = se, inf = est - 1.96 * se, sup = est + 1.96 * se,
           p_valor = 2 * pnorm(-abs(est / se)))
  })
}
lab_din   <- paste0("Elegido por AICc (rezagos ", min(spec_dinamico$lags), "-", max(spec_dinamico$lags), ")")
lab_largo <- paste0("Sin restricciones (rezagos 0-", MAX_LAG_X, ")")
acum_all <- bind_rows(mutate(traspaso_acumulado(fit_din), modelo = lab_din),
                      mutate(traspaso_acumulado(fit_largo), modelo = lab_largo),
                      mutate(traspaso_acumulado(fit), modelo = "ARIMAX principal"))
guardar_tab(acum_all |> mutate(across(where(is.double), ~ signif(.x, 5))), "eval_traspaso_acumulado.csv")
print(as.data.frame(acum_all |> mutate(across(where(is.double), ~ round(.x, 4)))))

wald_de <- function(fit_r, nombre) {
  nx <- grep("^tcn_l", names(coef(fit_r)), value = TRUE)
  V <- fit_r$var.coef[nx, nx, drop = FALSE]; b <- coef(fit_r)[nx]
  w <- as.numeric(t(b) %*% solve(V) %*% b)
  tibble(modelo = nombre, estadistico_wald = w, gl = length(b), p_valor = pchisq(w, length(b), lower.tail = FALSE))
}
guardar_tab(bind_rows(wald_de(fit_din, lab_din), wald_de(fit_largo, lab_largo)) |>
              mutate(across(where(is.double), ~ signif(.x, 4))), "eval_wald_traspaso.csv")

coef_x <- bind_rows(mutate(r_din$coefs, modelo = lab_din), mutate(r_largo$coefs, modelo = lab_largo)) |>
  filter(grepl("^tcn_l", termino)) |> mutate(rezago = as.integer(sub("tcn_l", "", termino)))
col_mod <- setNames(c(col_ipc, col_tcn), c(lab_din, lab_largo))
p_t1 <- ggplot(coef_x, aes(rezago, estimacion, colour = modelo)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup), position = position_dodge(width = 0.5), linewidth = 0.7) +
  scale_colour_manual(values = col_mod, name = NULL) +
  scale_x_continuous(breaks = 0:MAX_LAG_X) + theme(legend.position = "none") +
  labs(title = "Efecto de cada rezago", subtitle = "Coeficiente de d_tcn (t-k) sobre d_ipc_imp (t)",
       x = "Rezago k (meses; 0 = mismo mes)", y = "Coeficiente")
p_t2 <- ggplot(filter(acum_all, modelo != "ARIMAX principal"), aes(rezago, acumulado, colour = modelo, fill = modelo)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = inf, ymax = sup), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = col_mod, name = NULL) + scale_fill_manual(values = col_mod, name = NULL) +
  scale_x_continuous(breaks = 0:MAX_LAG_X) + guides(fill = "none") +
  labs(title = "Efecto acumulado sobre el nivel de precios",
       subtitle = "Respuesta del IPC_imp (%) a una depreciación de 1 %",
       x = "Rezago k (meses; 0 = mismo mes)", y = "Aumento acumulado del IPC_imp (%)") +
  theme(legend.position = "bottom")
p_trasp <- (p_t1 + p_t2) +
  plot_annotation(title = "Traspaso del tipo de cambio a los precios de bienes importados",
                  caption = "Bandas e intervalos: IC 95 %. Ventana 2012-2026. Elaboración propia con datos del BCP.")
guardar_fig(p_trasp, "eval_04_traspaso.png", 10, 5.2)

# ---- 4. Estabilidad: ¿cambió el traspaso desde marzo de 2020? ------------------
prueba_cambio <- function(spec, nombre) {
  dat <- ventana |> mutate(post = as.numeric(fecha >= as.Date("2020-03-01")))
  Xb <- armar_xreg(dat, spec)
  nx <- colnames(Xb)
  Xi <- Xb * dat$post; colnames(Xi) <- paste0(nx, "_post")
  args <- list(a_ts(dat$dln_ipc_imp, min(dat$fecha)), order = c(spec$p, 0, spec$q),
               seasonal = list(order = c(spec$P, 0, spec$Q), period = 12),
               xreg = cbind(Xb, Xi), include.mean = TRUE)
  f <- do.call(Arima, args)
  b <- coef(f); V <- f$var.coef
  suma <- function(nom) { est <- sum(b[nom]); se <- sqrt(sum(V[nom, nom])); c(est, se) }
  pre <- suma(nx); cam <- suma(paste0(nx, "_post"))
  tibble(modelo = nombre, traspaso_previo = pre[1], se_previo = pre[2],
         cambio_post_2020 = cam[1], se_cambio = cam[2],
         p_valor_cambio = 2 * pnorm(-abs(cam[1] / cam[2])))
}
cambio <- bind_rows(prueba_cambio(spec_final, "ARIMAX principal"),
                    prueba_cambio(specs$ex_ante, "ARIMAX ex-ante"),
                    prueba_cambio(spec_largo, lab_largo)) |>
  mutate(across(where(is.double), ~ round(.x, 4)))
guardar_tab(cambio, "eval_test_cambio_traspaso.csv")
print(as.data.frame(cambio))

# ---- 5. Robustez --------------------------------------------------------------
traspaso_de <- function(fit_r, etiqueta, n_obs, modelo) {
  if (is.null(fit_r)) return(tibble(modelo = modelo, variante = etiqueta, n = n_obs, traspaso_acumulado = NA_real_,
                                    error_est = NA_real_, p_valor = NA_real_))
  nx <- grep("^tcn_l", names(coef(fit_r)), value = TRUE)
  b_r <- coef(fit_r)[nx]; V_r <- fit_r$var.coef[nx, nx, drop = FALSE]
  u <- rep(1, length(b_r)); est <- sum(b_r); se <- sqrt(as.numeric(t(u) %*% V_r %*% u))
  tibble(modelo = modelo, variante = etiqueta, n = n_obs, traspaso_acumulado = est, error_est = se,
         p_valor = 2 * pnorm(-abs(est / se)))
}

# fechas de los atípicos de d_ipc_imp dentro de la ventana (mismo criterio que el EDA)
fechas_atip_y <- ventana |> mutate(z = z_robusto(dln_ipc_imp)) |> filter(abs(z) > UMBRAL_Z) |> pull(fecha)
cat("Atípicos de d_ipc_imp con dummy de impulso:", paste(format(fechas_atip_y, "%Y-%m"), collapse = ", "), "\n")

# d_tcn winsorizado (percentiles 2,5 y 97,5 de la ventana)
d_w <- {
  lim <- quantile(ventana$dln_tcn, c(0.025, 0.975))
  x <- d |> mutate(dln_tcn = pmin(pmax(dln_tcn, lim[1]), lim[2]))
  for (k in 1:12) x[[paste0("dln_tcn_l", k)]] <- dplyr::lag(x$dln_tcn, k)
  x$dln_tcn_l0 <- x$dln_tcn
  filter(x, fecha >= FECHA_INICIO)
}
pre_covid  <- filter(ventana, fecha <  as.Date("2020-03-01"))
post_covid <- filter(ventana, fecha >= as.Date("2020-03-01"))

robustez <- function(spec, nombre_modelo) {
  inicio_comp <- d |> filter(if_all(all_of(paste0("dln_tcn_l", spec$lags)), ~ !is.na(.x)), !is.na(dln_ipc_imp)) |>
    summarise(m = min(fecha)) |> pull(m)
  comp <- filter(d, fecha >= inicio_comp)
  pre_2012 <- filter(comp, fecha < FECHA_INICIO)
  spec_imp <- spec; spec_imp$impulsos <- fechas_atip_y
  bind_rows(
    traspaso_de(ajustar_spec(ventana, spec), "Ventana 2012-2026 (base)", nrow(ventana), nombre_modelo),
    traspaso_de(ajustar_spec(comp, spec), "Muestra completa (2004-2026)", nrow(comp), nombre_modelo),
    traspaso_de(ajustar_spec(ventana, spec_imp), "Con dummies de impulso en los atípicos de d_ipc_imp", nrow(ventana), nombre_modelo),
    traspaso_de(ajustar_spec(d_w, spec), "d_tcn winsorizado (2,5 %-97,5 %)", nrow(d_w), nombre_modelo),
    traspaso_de(ajustar_spec(pre_2012, spec), "Submuestra 2004-2011", nrow(pre_2012), nombre_modelo),
    traspaso_de(ajustar_spec(pre_covid, spec), "Submuestra 2012 - feb-2020", nrow(pre_covid), nombre_modelo),
    traspaso_de(ajustar_spec(post_covid, spec), "Submuestra mar-2020 - 2026", nrow(post_covid), nombre_modelo)
  )
}
rob <- bind_rows(robustez(spec_final, "ARIMAX principal"),
                 robustez(specs$ex_ante, "ARIMAX ex-ante"),
                 robustez(spec_largo, lab_largo)) |>
  mutate(across(where(is.double), ~ round(.x, 4)))
guardar_tab(rob, "eval_robustez.csv")
print(as.data.frame(rob))

p_rob <- rob |>
  filter(!is.na(traspaso_acumulado)) |>
  mutate(variante = fct_rev(fct_inorder(variante)),
         inf = traspaso_acumulado - 1.96 * error_est, sup = traspaso_acumulado + 1.96 * error_est) |>
  ggplot(aes(traspaso_acumulado, variante)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = inf, xmax = sup), colour = col_ipc) +
  facet_wrap(~modelo, ncol = 1, scales = "free_x") +
  labs(title = "Robustez del traspaso acumulado estimado",
       subtitle = "Suma de los coeficientes de d_tcn según variante (IC 95 %)",
       x = "Traspaso acumulado (elasticidad)", y = NULL, caption = "Elaboración propia.")
guardar_fig(p_rob, "eval_05_robustez.png", 9.5, 8.5)

cat("Evaluación terminada.\n")
