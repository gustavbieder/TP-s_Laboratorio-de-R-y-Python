# ==============================================================================
# 02_eda.R  -  Parte 2: Análisis exploratorio de datos
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# Entrada : data/processed/base_transformada.csv  (salida de 01_import_transform.R)
# Salida  : output/figures/eda_*.png  y  output/tables/eda_*.csv
# ==============================================================================

source("R/00_config.R", encoding = "UTF-8")
d <- cargar_datos()
w <- d |> filter(fecha >= FECHA_INICIO)          # ventana de análisis (2012 en adelante)

# ---- 1. Estructura general y valores faltantes --------------------------------
vars <- c("fecha", "ipc_imp", "tcn", "ln_ipc_imp", "ln_tcn", "dln_ipc_imp", "dln_tcn")
estructura <- tibble(
  variable = vars,
  tipo     = map_chr(d[vars], ~ class(.x)[1]),
  n_obs    = nrow(d),
  n_na     = map_int(d[vars], ~ sum(is.na(.x))),
  prop_na  = round(n_na / n_obs, 4)
)
guardar_tab(estructura, "eda_estructura.csv")

cobertura <- tibble(
  concepto = c("Observaciones (muestra completa)", "Primer mes", "Último mes",
               "Observaciones de la ventana de análisis", "Inicio de la ventana",
               "Meses duplicados", "Meses faltantes (huecos)"),
  valor = c(nrow(d), format(min(d$fecha)), format(max(d$fecha)), nrow(w),
            format(FECHA_INICIO), sum(duplicated(d$fecha)),
            length(seq(min(d$fecha), max(d$fecha), by = "month")) - nrow(d))
)
guardar_tab(cobertura, "eda_cobertura.csv")

# ---- 2. Series en niveles -----------------------------------------------------
eventos <- tribble(
  ~fecha,                ~serie, ~etiqueta,
  as.Date("2008-10-01"), "TCN",  "Crisis financiera\nglobal (2008)",
  FECHA_METAS_INF,       "IPC_imp",  "Metas de\ninflación (may-2011)",
  as.Date("2020-03-01"), "IPC_imp",  "COVID-19\n(mar-2020)"
)

niveles <- d |>
  select(fecha, IPC_imp = ipc_imp, TCN = tcn) |>
  pivot_longer(-fecha, names_to = "serie", values_to = "valor") |>
  mutate(serie_lab = recode(serie,
                            IPC_imp = "IPC_imp: IPC de bienes importados (base dic-2017 = 100)",
                            TCN = "TCN (guaraníes por dólar, serie mensual)"))

etiquetas_ev <- eventos |>
  left_join(niveles |> group_by(serie, serie_lab) |>
              summarise(y = max(valor) * 0.97, .groups = "drop"), by = "serie")

p_niveles <- ggplot(niveles, aes(fecha, valor, colour = serie)) +
  annotate("rect", xmin = min(d$fecha), xmax = FECHA_INICIO, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 0.7) +
  geom_vline(data = eventos, aes(xintercept = fecha), linetype = "dotted", colour = "grey30") +
  geom_text(data = etiquetas_ev, aes(x = fecha, y = y, label = etiqueta),
            colour = "grey20", size = 2.8, hjust = -0.05, vjust = 1, inherit.aes = FALSE) +
  facet_wrap(~serie_lab, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(IPC_imp = col_ipc, TCN = col_tcn), guide = "none") +
  labs(title = "IPC_imp y tipo de cambio nominal, enero 2004 - agosto 2026",
       subtitle = "Zona gris: años previos a la ventana de análisis (antes de 2012)",
       x = NULL, y = NULL, caption = "Fuente: Banco Central del Paraguay. Elaboración propia.")
guardar_fig(p_niveles, "eda_01_series_niveles.png", 9, 6)

# ---- 3. Series en primeras diferencias del logaritmo, con atípicos ------------
dln_long <- d |>
  filter(!is.na(dln_ipc_imp)) |>
  select(fecha, IPC_imp = dln_ipc_imp, TCN = dln_tcn) |>
  pivot_longer(-fecha, names_to = "serie", values_to = "valor") |>
  group_by(serie) |>
  mutate(z = if_else(fecha >= FECHA_INICIO,
                     (valor - median(valor[fecha >= FECHA_INICIO])) / mad(valor[fecha >= FECHA_INICIO]),
                     NA_real_),
         atipico = !is.na(z) & abs(z) > UMBRAL_Z,
         serie_lab = recode(serie, IPC_imp = "d_ipc_imp: variación mensual del IPC_imp (%)",
                                   TCN = "d_tcn: variación mensual del TCN (%)")) |>
  ungroup()

atipicos <- dln_long |> filter(atipico) |>
  transmute(fecha, serie, variacion_pct = round(valor * 100, 2), z_robusto = round(z, 2))
guardar_tab(atipicos, "eda_atipicos.csv")

p_dln <- ggplot(dln_long, aes(fecha, valor * 100, colour = serie)) +
  annotate("rect", xmin = min(d$fecha), xmax = FECHA_INICIO, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.5) +
  geom_point(data = filter(dln_long, atipico), colour = "black", size = 2.2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  facet_wrap(~serie_lab, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(IPC_imp = col_ipc, TCN = col_tcn), guide = "none") +
  labs(title = "Primera diferencia del logaritmo (variación mensual aproximada, en %)",
       subtitle = paste0("Puntos negros: atípicos de la ventana de análisis (|z robusto| > ", UMBRAL_Z,
                         "). Zona gris: fuera de la ventana."),
       x = NULL, y = "%", caption = "Fuente: BCP. Elaboración propia. Las fechas de los atípicos figuran en output/tables/eda_atipicos.csv.")
guardar_fig(p_dln, "eda_02_series_dln.png", 9, 6)

# ---- 4. Por qué la ventana empieza en 2012: cambio de volatilidad -------------
vol <- d |> filter(!is.na(dln_ipc_imp)) |>
  transmute(fecha, IPC_imp = dln_ipc_imp, TCN = dln_tcn) |>
  mutate(across(c(IPC_imp, TCN), ~ zoo::rollapplyr(.x, 24, sd, fill = NA) * 100)) |>
  pivot_longer(-fecha, names_to = "serie", values_to = "sd24") |>
  filter(!is.na(sd24))

p_vol <- ggplot(vol, aes(fecha, sd24, colour = serie)) +
  annotate("rect", xmin = min(d$fecha), xmax = FECHA_INICIO, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = FECHA_METAS_INF, linetype = "dotted") +
  annotate("text", x = FECHA_METAS_INF, y = Inf, label = " Metas de inflación", hjust = 0, vjust = 1.5, size = 3) +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(IPC_imp = col_ipc, TCN = col_tcn), guide = "none") +
  labs(title = "Desvío estándar móvil de 24 meses de las variaciones mensuales (%)",
       subtitle = "La volatilidad del TCN es máxima en 2008-2012; la del IPC_imp cambia con el tiempo y es mayor desde 2021",
       x = NULL, y = "Desvío estándar (puntos porcentuales)", caption = "Fuente: BCP. Elaboración propia.")
guardar_fig(p_vol, "eda_03_volatilidad_movil.png", 9, 5.5)

subper <- d |> filter(!is.na(dln_ipc_imp)) |>
  mutate(periodo = case_when(fecha < FECHA_INICIO ~ "2004-2011",
                             fecha < as.Date("2020-03-01") ~ "2012 - feb-2020",
                             TRUE ~ "mar-2020 - 2026")) |>
  group_by(periodo) |>
  summarise(n = n(),
            media_ipc = mean(dln_ipc_imp) * 100, desvio_ipc = sd(dln_ipc_imp) * 100,
            media_tcn = mean(dln_tcn) * 100, desvio_tcn = sd(dln_tcn) * 100,
            .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
guardar_tab(subper, "eda_subperiodos.csv")

# ---- 5. Algo raro: meses con variación cero del IPC_imp (redondeo) ----------------
calidad_ipc_imp <- d |> filter(!is.na(dln_ipc_imp)) |>
  mutate(anio = year(fecha),
         un_decimal = abs(ipc_imp * 10 - round(ipc_imp * 10)) < 1e-6,
         cero = dln_ipc_imp == 0) |>
  group_by(anio) |>
  summarise(meses_cero = sum(cero), prop_un_decimal = mean(un_decimal), .groups = "drop") |>
  mutate(tramo = if_else(prop_un_decimal > 0.5, "Publicado con 1 decimal (desde 2018)",
                         "Empalmado hacia atrás (decimales largos)"))
guardar_tab(calidad_ipc_imp, "eda_calidad_ipc_imp.csv")

p_ceros <- ggplot(calidad_ipc_imp, aes(factor(anio), meses_cero, fill = tramo)) +
  geom_col() +
  scale_fill_manual(values = c("Publicado con 1 decimal (desde 2018)" = col_tcn,
                               "Empalmado hacia atrás (decimales largos)" = col_ipc), name = NULL) +
  scale_y_continuous(breaks = 0:4) +
  labs(title = "Meses con variación mensual del IPC_imp exactamente igual a cero, por año",
       subtitle = "Hay ceros en todo el período; desde 2018 el índice se publica con un solo decimal",
       x = NULL, y = "Cantidad de meses", caption = "Fuente: BCP. Elaboración propia.") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))
guardar_fig(p_ceros, "eda_04_ceros_ipc_imp.png", 9, 4.5)

# ---- 6. Distribuciones univariadas (ventana) ----------------------------------
normalidad <- tibble(
  serie = c("d_ipc_imp", "d_tcn"),
  asimetria = c(mean((w$dln_ipc_imp - mean(w$dln_ipc_imp))^3) / sd(w$dln_ipc_imp)^3,
                mean((w$dln_tcn - mean(w$dln_tcn))^3) / sd(w$dln_tcn)^3),
  curtosis_exceso = c(mean((w$dln_ipc_imp - mean(w$dln_ipc_imp))^4) / sd(w$dln_ipc_imp)^4 - 3,
                      mean((w$dln_tcn - mean(w$dln_tcn))^4) / sd(w$dln_tcn)^4 - 3),
  jarque_bera_p = c(tseries::jarque.bera.test(w$dln_ipc_imp)$p.value,
                    tseries::jarque.bera.test(w$dln_tcn)$p.value)
) |> mutate(across(where(is.numeric), ~ round(.x, 4)))
guardar_tab(normalidad, "eda_normalidad.csv")

grafico_dist <- function(x, nombre, color) {
  x <- x * 100
  h <- ggplot(tibble(x), aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 25, fill = color, alpha = 0.6, colour = "white") +
    stat_function(fun = dnorm, args = list(mean = mean(x), sd = sd(x)), colour = "black", linetype = "dashed") +
    labs(title = paste("Histograma de", nombre), x = "%", y = "Densidad")
  q <- ggplot(tibble(x), aes(sample = x)) +
    stat_qq(colour = color, alpha = 0.7) + stat_qq_line(colour = "black") +
    labs(title = paste("Gráfico Q-Q de", nombre), x = "Cuantiles teóricos", y = "Cuantiles muestrales")
  h + q
}
p_dist <- grafico_dist(w$dln_ipc_imp, "d_ipc_imp", col_ipc) / grafico_dist(w$dln_tcn, "d_tcn", col_tcn) +
  plot_annotation(title = "Distribución de las variaciones mensuales, ventana 2012-2026",
                  subtitle = "Línea punteada: densidad normal con la misma media y desvío",
                  caption = "Elaboración propia con datos del BCP.")
guardar_fig(p_dist, "eda_05_distribuciones.png", 9, 6.5)

# ---- 7. Estacionalidad --------------------------------------------------------
p_f_ipc <- anova(lm(dln_ipc_imp ~ mes, data = w))$`Pr(>F)`[1]
p_f_tcn <- anova(lm(dln_tcn ~ mes, data = w))$`Pr(>F)`[1]
guardar_tab(tibble(serie = c("d_ipc_imp", "d_tcn"), p_valor_F_dummies_mensuales = signif(c(p_f_ipc, p_f_tcn), 3)),
            "eda_estacionalidad_test.csv")

meses_es <- c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
est <- w |> select(mes, IPC_imp = dln_ipc_imp, TCN = dln_tcn) |>
  pivot_longer(-mes, names_to = "serie", values_to = "valor") |>
  mutate(serie_lab = recode(serie,
                            IPC_imp = paste0("d_ipc_imp (test F de igualdad de medias mensuales: p = ", signif(p_f_ipc, 2), ")"),
                            TCN = paste0("d_tcn (test F de igualdad de medias mensuales: p = ", signif(p_f_tcn, 2), ")")))
p_est <- ggplot(est, aes(mes, valor * 100, fill = serie)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~serie_lab, ncol = 1, scales = "free_y") +
  scale_x_discrete(labels = meses_es) +
  scale_fill_manual(values = c(IPC_imp = col_ipc, TCN = col_tcn), guide = "none") +
  labs(title = "Estacionalidad de las variaciones mensuales, ventana 2012-2026",
       x = NULL, y = "%", caption = "Elaboración propia con datos del BCP.")
guardar_fig(p_est, "eda_06_estacionalidad.png", 9, 6)

# ---- 8. Autocorrelación -------------------------------------------------------
acf_plot <- function(x, titulo, color) {
  ggAcf(x, lag.max = 24) + labs(title = paste("ACF de", titulo), y = NULL) +
    theme(plot.title = element_text(size = 10)) -> a
  ggPacf(x, lag.max = 24) + labs(title = paste("PACF de", titulo), y = NULL) +
    theme(plot.title = element_text(size = 10)) -> b
  a + b
}
p_acf <- acf_plot(w$dln_ipc_imp, "d_ipc_imp", col_ipc) / acf_plot(w$dln_tcn, "d_tcn", col_tcn) +
  plot_annotation(title = "Autocorrelación y autocorrelación parcial, ventana 2012-2026",
                  subtitle = "Bandas: intervalo de confianza del 95 % bajo ruido blanco",
                  caption = "Elaboración propia con datos del BCP.")
guardar_fig(p_acf, "eda_07_acf_pacf.png", 9, 6.5)

# ---- 9. Relación bivariada: d_ipc_imp(t) vs d_tcn(t-k) ----------------------------
cc <- tibble(rezago = 0:12) |>
  mutate(correlacion = map_dbl(rezago, ~ cor(w[[paste0("dln_tcn_l", .x)]], w$dln_ipc_imp, use = "complete.obs")))
guardar_tab(cc |> mutate(correlacion = round(correlacion, 3)), "eda_correlacion_cruzada.csv")
banda <- 1.96 / sqrt(nrow(w))

p_sc0 <- ggplot(w, aes(dln_tcn_l0 * 100, dln_ipc_imp * 100)) +
  geom_point(colour = col_ipc, alpha = 0.6) + geom_smooth(method = "lm", colour = col_tcn, fill = col_tcn, alpha = 0.15) +
  labs(title = "Sin rezago", subtitle = paste0("Correlación = ", round(cc$correlacion[1], 2)),
       x = "d_tcn (t), %", y = "d_ipc_imp (t), %")
p_sc1 <- ggplot(w, aes(dln_tcn_l1 * 100, dln_ipc_imp * 100)) +
  geom_point(colour = col_ipc, alpha = 0.6) + geom_smooth(method = "lm", colour = col_tcn, fill = col_tcn, alpha = 0.15) +
  labs(title = "Rezago de 1 mes", subtitle = paste0("Correlación = ", round(cc$correlacion[2], 2)),
       x = "d_tcn (t-1), %", y = "d_ipc_imp (t), %")
p_cc <- ggplot(cc, aes(rezago, correlacion)) +
  geom_col(fill = col_ipc, alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = c(-banda, banda), linetype = "dashed", colour = col_tcn) +
  scale_x_continuous(breaks = 0:12) +
  labs(title = "Correlación cruzada", subtitle = "d_ipc_imp (t) con d_tcn (t-k); líneas: banda del 95 %",
       x = "Rezago k (meses)", y = "Correlación")
p_bi <- (p_sc0 | p_sc1) / p_cc +
  plot_annotation(title = "Relación entre las variaciones mensuales del IPC_imp y del tipo de cambio",
                  caption = "Ventana 2012-2026. Elaboración propia con datos del BCP.")
guardar_fig(p_bi, "eda_08_relacion_bivariada.png", 9, 7)

# ---- 10. Raíces unitarias (ventana) -------------------------------------------
test_ur <- function(x, con_tendencia) {
  suppressWarnings(tibble(
    adf_p  = tseries::adf.test(x)$p.value,
    pp_p   = tseries::pp.test(x)$p.value,
    kpss_p = tseries::kpss.test(x, null = if (con_tendencia) "Trend" else "Level")$p.value
  ))
}
raices <- bind_rows(
  tibble(serie = "ln(IPC_imp)",  transformacion = "Nivel del logaritmo", test_ur(w$ln_ipc_imp,  TRUE)),
  tibble(serie = "ln(TCN)",  transformacion = "Nivel del logaritmo", test_ur(w$ln_tcn,  TRUE)),
  tibble(serie = "d_ipc_imp",    transformacion = "Primera diferencia",  test_ur(w$dln_ipc_imp, FALSE)),
  tibble(serie = "d_tcn",    transformacion = "Primera diferencia",  test_ur(w$dln_tcn, FALSE))
) |> mutate(across(ends_with("_p"), ~ round(.x, 3)))
guardar_tab(raices, "eda_raices_unitarias.csv")

cat("EDA terminado. Figuras en", dir_fig, "y tablas en", dir_tab, "\n")
