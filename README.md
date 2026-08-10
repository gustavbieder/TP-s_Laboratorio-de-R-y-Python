# TP's_Laboratorio-de-R-y-Python
# TP Integrador — Módulo Python
**Laboratorio de Programación en Python y R — Maestría en Econometría (UTDT)**

## ¿Qué contiene este trabajo?

El notebook `TP_python_integrador.ipynb` resuelve las tres partes del trabajo práctico integrador:

- **Parte 1 — Robustez ante contaminación de muestras.** Simulaciones Monte Carlo (M=1000
  y M=500) que comparan sesgo y ECM empírico de la media vs. la mediana, y de MCO
  (implementado a mano con la fórmula matricial, sin `statsmodels`/`sklearn`) vs. LAD
  (`scipy.optimize.minimize`), bajo distintos niveles de contaminación con outliers.
- **Parte 2 — Paradoja de Simpson y sesgo por variable omitida.** Un proceso generador
  de datos diseñado para que la regresión marginal `Y~X` dé un coeficiente nulo, `Y~X+Z`
  dé negativo y `Y~X+W` dé positivo, con evidencia (M=1000 simulaciones y un ejercicio de
  n creciente hasta 100.000) de que el sesgo por variable omitida no desaparece con
  muestras grandes, y una ilustración visual de la paradoja.
- **Parte 3 — Análisis empírico con Gapminder.** Estudio de Paraguay: convergencia en
  esperanza de vida respecto al promedio mundial (1952-2007) y un episodio atípico de
  crecimiento del PBI per cápita respecto al promedio del continente, vinculado a la
  construcción de la represa de Itaipú (1972-1982).

Cada parte incluye al menos una celda de markdown con interpretación económica, y el
notebook fija una semilla (`SEED = 42`) al inicio para que todos los resultados sean
reproducibles.

## Cómo ejecutarlo
1. Entrá a [Google Colab](https://colab.research.google.com/).
2. `Archivo → Abrir notebook → GitHub`, pegá la URL de este repositorio y seleccioná
   `TP_python_integrador.ipynb`.
3. `Entorno de ejecución → Ejecutar todas`. El notebook corre de punta a punta sin
   modificaciones ni archivos externos: los datos se generan por simulación y el
   dataset de Gapminder viene incluido en la librería `plotly`.


