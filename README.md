# TP Integrador – Módulo R: traspaso del tipo de cambio al IPC de bienes importados (ARIMAX)

Trabajo práctico integrador del Laboratorio de Programación en Python y R
(Maestría en Econometría, Universidad Torcuato Di Tella, 2026).

**Alumno:** Gustavo Santiago Biedermann Giménez

## Pregunta

¿Cuánto de una depreciación del guaraní frente al dólar se traslada al IPC de bienes importados (IPC_imp) y
ese vínculo mejora el pronóstico de su variación mensual a un mes?

## Dataset

- **Fuente:** Banco Central del Paraguay (bcp.gov.py). Fecha de descarga: *(completar)*.
- **Archivo:** `data/raw/Base_Arimax.xlsx` (hoja `Hoja1`).
- **Variables usadas:** `IPC_imp` (IPC de bienes importados, base diciembre 2017 = 100) y `TCN` (tipo de cambio
  nominal, guaraníes por dólar). Frecuencia mensual, enero 2004 a agosto 2026 (272 observaciones). El archivo trae
  también el IPC general (`IPC`), que este trabajo no utiliza.
- **Transformación:** logaritmo natural y primera diferencia → `dln_ipc_imp` (`d_ipc_imp`) y `dln_tcn` (`d_tcn`),
  variaciones mensuales aproximadas (0,01 ≈ 1 %).
- **Ventana de análisis:** enero 2012 a agosto 2026 (176 meses), primer año completo con el esquema de metas de
  inflación (vigente desde mayo de 2011). La muestra completa se usa como prueba de robustez.

## Técnica

**ARIMAX** (regresión con errores ARMA) de `d_ipc_imp` sobre `d_tcn` y sus rezagos. Se busca el mejor modelo con una
grilla de 756 especificaciones (órdenes ARMA, estacionalidad, rezagos del tipo de cambio), filtros de residuos y de
raíces, y validación temporal; luego se evalúa fuera de muestra (2023-01 a 2026-08) con Diebold-Mariano
(pérdida cuadrática y absoluta) y se analizan el traspaso acumulado, su estabilidad y su robustez.

## Estructura del repositorio

```
├── run_all.R                    # corre todo el proyecto de punta a punta
├── TP_Integrador_R.Rproj        # proyecto de RStudio (fija el directorio de trabajo)
├── R/
│   ├── 00_config.R              # paquetes, semilla, rutas, parámetros, tema gráfico
│   ├── utils_modelos.R          # funciones para especificar, ajustar y evaluar ARIMAX
│   ├── 01_import_transform.R    # importa el Excel, verifica calidad, log y primera diferencia
│   ├── 02_eda.R                 # Parte 2: análisis exploratorio (ggplot2)
│   ├── 03_seleccion_arimax.R    # Parte 3 (1/2): grilla, filtros y validación temporal
│   └── 04_evaluacion_resultados.R  # Parte 3 (2/2): fuera de muestra, diagnósticos, traspaso, robustez
├── data/
│   ├── raw/Base_Arimax.xlsx
│   ├── processed/base_transformada.csv
│   └── processed/Base_transformada_IPC_imp_TCN.xlsx   # misma base transformada, en Excel con fórmulas
├── output/
│   ├── figures/                 # gráficos (PNG)
│   └── tables/                  # tablas de resultados (CSV)
└── report/
    └── informe_TP_R.pdf         # informe en PDF (sin código)
```

## Cómo correr el proyecto

Requisitos: R ≥ 4.1 y conexión a internet la primera vez (los paquetes que falten se instalan solos):
`tidyverse`, `readxl`, `lubridate`, `forecast`, `tseries`, `lmtest`, `strucchange`, `patchwork`,
`scales`, `zoo`.

1. Abrir `TP_Integrador_R.Rproj` en RStudio (o ubicar el directorio de trabajo en la raíz del repositorio).
2. Ejecutar:

```r
source("run_all.R")
```

o, desde la terminal, en la raíz del repositorio: `Rscript run_all.R`.

Tiempo aproximado: 5 a 10 minutos (la mayor parte es la búsqueda en grilla y la validación temporal).
Todas las rutas son relativas y la semilla se fija en `R/00_config.R` (`set.seed(2026)`); el
procedimiento es determinístico. Los scripts también se pueden correr uno por uno, en orden numérico.

## Resultados principales

- **Traspaso parcial y rápido.** Una depreciación de 1 % eleva el IPC_imp cerca de 0,12 % en el mismo mes y 0,18 %
  acumulado al mes siguiente (IC 95 %: 0,11 a 0,26). El efecto se agota en dos meses. Es robusto a la muestra
  (0,13 a 0,23) y no hay evidencia de cambio desde 2020.
- **Pronóstico a un mes, 2023-2026.** Usar el tipo de cambio del mes reduce el error absoluto medio unos 14 % frente a
  la media histórica (significativo) y el RMSE unos 6 % (no significativo). Con solo tipos de cambio pasados, la
  mejora es marginal.
- **Limitaciones.** Residuos con colas pesadas y heterocedasticidad, asociación no causal y muestra de evaluación corta.

Detalle, gráficos, decisiones y limitaciones: `report/informe_TP_R.pdf`.
