# ==============================================================================
# utils_modelos.R  -  Funciones para especificar y ajustar modelos ARIMAX
# ------------------------------------------------------------------------------
# Autor    : Gustavo Santiago Biedermann Giménez
# Un modelo ARIMAX aquí es una regresión de d_ipc_imp sobre rezagos de d_tcn con
# errores ARMA (estacionales o no):
#     d_ipc_t = c + sum_k b_k * d_tcn_{t-k} + [dummies] + eta_t ,  eta ~ SARMA
# Como d_ipc_imp ya es la primera diferencia del log del IPC, el componente "I" del
# ARIMA está en la transformación (equivale a un ARIMAX(p,1,q) sobre ln(IPC)).
#
# Una especificación ("spec") es una lista con:
#   p, q, P, Q   : órdenes ARMA regular y estacional (período 12)
#   lags         : rezagos de d_tcn incluidos (integer(0) = modelo sin X)
#   dummies      : TRUE si se incluyen 11 dummies mensuales (estacionalidad determinista)
#   impulsos     : fechas (Date) con dummies de impulso para atípicos (opcional)
# ==============================================================================

nueva_spec <- function(p = 0, q = 0, P = 0, Q = 0, lags = integer(0),
                       dummies = FALSE, impulsos = NULL) {
  list(p = p, q = q, P = P, Q = Q, lags = as.integer(lags),
       dummies = dummies, impulsos = impulsos)
}

# Matriz de regresores exógenos para un data frame con columnas dln_tcn_lK y fecha
armar_xreg <- function(datos, spec) {
  m <- NULL
  if (length(spec$lags) > 0) {
    m <- sapply(spec$lags, function(k) datos[[paste0("dln_tcn_l", k)]])
    m <- matrix(m, nrow = nrow(datos))
    colnames(m) <- paste0("tcn_l", spec$lags)
  }
  if (isTRUE(spec$dummies)) {
    dm <- matrix(sapply(2:12, function(mm) as.numeric(month(datos$fecha) == mm)), nrow = nrow(datos))
    colnames(dm) <- paste0("mes_", 2:12)
    m <- cbind(m, dm)
  }
  if (!is.null(spec$impulsos)) {
    im <- sapply(spec$impulsos, function(f) as.numeric(datos$fecha == f))
    im <- matrix(im, nrow = nrow(datos))
    colnames(im) <- paste0("imp_", format(spec$impulsos, "%Y_%m"))
    m <- cbind(m, im)
  }
  m
}

# Ajusta la especificación sobre un data frame; devuelve NULL si falla
ajustar_spec <- function(datos, spec) {
  y <- a_ts(datos$dln_ipc_imp, min(datos$fecha))
  X <- armar_xreg(datos, spec)
  # do.call guarda los regresores dentro de la llamada, lo que permite pronosticar después
  args <- list(y, order = c(spec$p, 0, spec$q),
               seasonal = list(order = c(spec$P, 0, spec$Q), period = 12),
               include.mean = TRUE)
  if (!is.null(X)) args$xreg <- X
  tryCatch(suppressWarnings(do.call(Arima, args)), error = function(e) NULL)
}

# Módulo de la raíz más cercana a la circunferencia unitaria (AR y MA, con la parte
# estacional expandida). Valores cercanos a 1 señalan soluciones casi no estacionarias
# o casi no invertibles, típicamente inestables al reestimar.
raiz_minima <- function(fit) {
  recortar <- function(x) { nz <- which(abs(x) > 1e-12); if (length(nz)) x[seq_len(max(nz))] else numeric(0) }
  phi   <- recortar(fit$model$phi)
  theta <- recortar(fit$model$theta)
  r_ar <- if (length(phi))   min(Mod(polyroot(c(1, -phi))))  else Inf
  r_ma <- if (length(theta)) min(Mod(polyroot(c(1, theta)))) else Inf
  min(r_ar, r_ma)
}

# Etiqueta legible de una especificación
etiqueta_spec <- function(spec) {
  x <- if (length(spec$lags) == 0) "sin X" else
    paste0("d_tcn rezagos ", min(spec$lags), "-", max(spec$lags))
  est <- if (isTRUE(spec$dummies)) " + dummies mensuales" else ""
  paste0("ARMA(", spec$p, ",", spec$q, ")(", spec$P, ",", spec$Q, ")[12], ", x, est)
}

# Tipo de modelo según el tratamiento del regresor
tipo_x <- function(lags) {
  if (length(lags) == 0) "sin_X" else if (min(lags) == 0) "condicional" else "ex_ante"
}

# Pronóstico a un paso con reajuste de coeficientes (ventana expansiva)
# datos_full : data frame completo (con rezagos ya construidos)
# fechas_test: meses a pronosticar; para cada uno se estima con todo lo anterior.
# Si la reestimación falla o cae en una solución degenerada (raíces demasiado cerca
# de la circunferencia unitaria), se reutilizan los parámetros del último ajuste válido.
pronostico_rolling <- function(datos_full, spec, fechas_test, inicio = FECHA_INICIO) {
  previo <- NULL
  map_dfr(fechas_test, function(f) {
    tr <- datos_full |> filter(fecha >= inicio, fecha < f)
    te <- datos_full |> filter(fecha == f)
    ajuste <- ajustar_spec(tr, spec)
    ok <- !is.null(ajuste) && raiz_minima(ajuste) >= RAIZ_MIN
    if (ok) {
      previo <<- ajuste
    } else if (!is.null(previo)) {
      y_tr <- a_ts(tr$dln_ipc_imp, min(tr$fecha))
      X_tr <- armar_xreg(tr, spec)
      args <- list(y_tr, model = previo)
      if (!is.null(X_tr)) args$xreg <- X_tr
      ajuste <- tryCatch(suppressWarnings(do.call(Arima, args)), error = function(e) NULL)
    } else {
      ajuste <- NULL
    }
    if (is.null(ajuste)) return(tibble(fecha = f, pronostico = NA_real_, reestimado = FALSE))
    pr <- forecast(ajuste, h = 1, xreg = armar_xreg(te, spec))
    tibble(fecha = f, pronostico = as.numeric(pr$mean), reestimado = ok)
  })
}
