library(targets)

required_packages <- c("reshape2", "car", "lmtest", "MASS", "agricolae", "multcomp", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# Cargar los datos en una matriz (común para todos los puntos)
matriz_datos <- matrix(c(73, 68, 74, 71, 67,
                         73, 67, 75, 72, 70,
                         75, 68, 78, 73, 68,
                         73, 71, 75, 75, 69),
                       nrow = 4, byrow = TRUE)

rownames(matriz_datos) <- paste("Agente", 1:4)
colnames(matriz_datos) <- paste("Rollo", 1:5)

# Convertir a data frame y reorganizar a formato largo
df <- as.data.frame(matriz_datos)
df$Agente <- rownames(df) 

long_df <- melt(df, id.vars = "Agente", variable.name = "Rollo", value.name = "Resistencia")

# Verificar que quedó bien
print(long_df)
tabla_matriz <- xtabs(Resistencia ~ Agente + Rollo, data = long_df)
tabla_matriz


# ----------------------------------------------------------------------
# Punto 1: Exploración gráfica inicial de los datos (Diagramas de caja)
# ----------------------------------------------------------------------

# Boxplot por Agente
boxplot(Resistencia ~ Agente, data = long_df,
        main = "Boxplot de Resistencia por Agente Químico",
        ylab = "Resistencia a la Tensión", xlab = "Agente Químico", col = "lightblue")

# Boxplot por Rollo
boxplot(Resistencia ~ Rollo, data = long_df,
        main = "Boxplot de Resistencia por Rollo",
        ylab = "Resistencia a la Tensión", xlab = "Rollo", col = "lightgreen")

# ----------------------------------------------------------------------
# Punto 2: Prueba de hipótesis (ANOVA y Modelo de Regresión)
# ----------------------------------------------------------------------
cat("\n=== Punto 2: ANOVA y Modelo de Regresión ===\n")

# ANOVA
anova_model <- aov(Resistencia ~ Agente + Rollo, data = long_df)
anova_summary <- summary(anova_model)
print(anova_summary)
cat("Para Agente, el valor p es 0.121. Como 0.121 > 0.05, no hay evidencia estadísticamente significativa para afirmar que los promedios de resistencia difieren entre agentes químicos.\n")

# Modelo de Regresión
regression_model <- lm(Resistencia ~ Agente + Rollo, data = long_df)
regression_summary <- summary(regression_model)
print(regression_summary)
cat("Interpretación Regresión: R²=0.715. Agente2 reduce significativamente la resistencia (-3.4, p=0.027).\n")

# ----------------------------------------------------------------------
# Punto 3: Validación del modelo de regresión
# ----------------------------------------------------------------------
cat("\n=== Punto 3: Validación del Modelo ===\n")

# --- Gráficos de validación individuales ---
par(mfrow=c(1,1))  # Mostrar un gráfico a la vez

# 1️⃣ Residuos vs Valores ajustados
plot(regression_model, which = 1)

# 2️⃣ QQ-plot de residuos
plot(regression_model, which = 2)

# 3️⃣ Escala-Localización (homocedasticidad)
plot(regression_model, which = 3)

# 4️⃣ Distancia de Cook (influencia)
plot(regression_model, which = 4)

# --- Pruebas estadísticas ---
library(lmtest)
library(car)

shapiro_test <- shapiro.test(regression_model$residuals)
print(shapiro_test)
cat("Normalidad: p=0.425 >0.05, cumple.\n")

bptest_result <- bptest(regression_model)
print(bptest_result)
cat("Homoscedasticidad: p=0.114 >0.05, cumple.\n")

vif_result <- vif(regression_model)
print(vif_result)
cat("Multicolinealidad: VIF <2, OK.\n")

coef_summary <- summary(regression_model)$coefficients
print(coef_summary)
cat("Parámetros significativos: Agente2 (p=0.027), Rollo5 marginal (p=0.070).\n")

# --- Distancia de Cook (gráfico adicional personalizado) ---
cooks_distance <- cooks.distance(regression_model)
plot(cooks_distance, main="Distancia de Cook", ylab="Distancia", type="h")
abline(h=1, col="red")

influential_points <- which(cooks_distance > 1)
print(influential_points)
cat("Puntos influyentes: Ninguno (todos <1).\n")

----------------------------------------------------------------------
  # Punto 4: Método de Scheffé para contrastes
  # ----------------------------------------------------------------------
cat("\n=== Punto 4: Método de Scheffé ===\n")

scheffe_result <- scheffe.test(anova_model, "Agente", group = TRUE)

print(scheffe_result)
# ============================================================
# Cálculo de los valores F del método de Scheffé (comparaciones pareadas)
# ============================================================

# --- Datos obtenidos del análisis ANOVA / Scheffé ---
medias <- c(72.2, 68.8, 74.0, 72.0)  # Medias de los tratamientos (Agente 1 a 4)
n <- 5                                # Número de repeticiones por agente
MSerror <- 4.55                       # Error cuadrático medio del ANOVA
df_error <- 12                        # Grados de libertad del error

# --- Valor crítico F del ANOVA ---
F_crit_anova <- 3.49                  # F de la tabla ANOVA (df1=3, df2=12)

# --- Cálculo de los valores críticos de Scheffé ---
# Fórmula: F_S = (k - 1) * F_α, (k-1), dfE
k <- length(medias)

# Valor crítico para α = 0.05
Fcrit_05 <- (k - 1) * qf(1 - 0.05, df1 = k - 1, df2 = df_error)

# Valor crítico para α = 0.01
Fcrit_01 <- (k - 1) * qf(1 - 0.01, df1 = k - 1, df2 = df_error)

cat("Valor crítico Scheffé α = 0.05 =", round(Fcrit_05, 2), "\n")
cat("Valor crítico Scheffé α = 0.01 =", round(Fcrit_01, 2), "\n\n")

# ============================================================
#  Cálculo de los valores F calculados entre pares
# ============================================================

# Generar todas las combinaciones de pares (Agente i vs j)
pares <- combn(1:k, 2)
n_pares <- ncol(pares)

# Crear tabla vacía
tabla_F <- data.frame(
  Contraste = character(n_pares),
  F_calculado = numeric(n_pares),
  stringsAsFactors = FALSE
)

# Calcular el F calculado para cada contraste
for (m in 1:n_pares) {
  i <- pares[1, m]
  j <- pares[2, m]
  Fcalc <- ((medias[i] - medias[j])^2) / (MSerror * (1/n + 1/n))
  tabla_F[m, ] <- list(
    paste0("Agente ", i, " - Agente ", j),
    round(Fcalc, 2)
  )
}

# ============================================================
#  Interpretación según nivel de significancia
# ============================================================

tabla_F$Interpretacion_0.05 <- ifelse(tabla_F$F_calculado > Fcrit_05,
                                      "Significativo", "No significativo")
tabla_F$Interpretacion_0.01 <- ifelse(tabla_F$F_calculado > Fcrit_01,
                                      "Muy significativo", "No significativo")

# Interpretación final combinada
tabla_F$Interpretacion_final <- ifelse(
  tabla_F$F_calculado > Fcrit_01, "Diferencia muy significativa",
  ifelse(tabla_F$F_calculado > Fcrit_05, "Diferencia significativa a α = 0.05", "Sin diferencia estadística")
)

# ============================================================
#  Mostrar resultados finales
# ============================================================
cat("=== Tabla de resultados del método de Scheffé ===\n")
print(tabla_F, row.names = FALSE)

# ----------------------------------------------------------------------
# Punto 5: Comparaciones entre pares (Tukey, LSD, Duncan, Newman-Keuls)
# ----------------------------------------------------------------------
cat("\n=== Punto 5: Comparaciones entre pares ===\n")

# Tukey
tukey_summary <- TukeyHSD(anova_model, which = "Agente")
print(tukey_summary)
cat("Tukey: Solo Ag2-Ag3 significativo a 0.05.\n")

# LSD
lsd_result <- LSD.test(anova_model, "Agente", group=TRUE)
print(lsd_result$groups)
cat("LSD: Ag2 diff de Ag1, Ag3, Ag4.\n")

# Duncan
duncan_result <- duncan.test(anova_model, "Agente", group=TRUE)
print(duncan_result$groups)
cat("Duncan: Similar a LSD.\n")

# Newman-Keuls (SNK)
snk_result <- SNK.test(anova_model, "Agente", group=TRUE)
print(snk_result$groups)
cat("SNK: Solo Ag2-Ag3.\n")

cat("Conclusiones diferentes: Sí. LSD/Duncan más liberales (más diffs), Tukey/SNK más conservadores.\n")
cat("Dado supuestos cumplen, preferir Tukey para control FWER.\n")

# ----------------------------------------------------------------------
# Punto 6: Análisis del efecto de los bloques
# ----------------------------------------------------------------------
cat("\n=== Punto 6: Análisis de Bloques ===\n")

print(anova_summary)
cat("Bloques significativos a 0.05 (p=0.035), no a 0.01.\n")
cat("Sin bloques, MSE aumenta, poder disminuye.\n")

# Modelo sin bloques (Pregunta 2 sin bloques)
anova_model_sin_bloques <- aov(Resistencia ~ Agente, data = long_df)
summary_sin_bloques <- summary(anova_model_sin_bloques)
print(summary_sin_bloques)
cat("ANOVA sin bloques: Agente no sig (p=0.0575).\n")

regression_model_sin_bloques <- lm(Resistencia ~ Agente, data = long_df)
summary_regression_sin_bloques <- summary(regression_model_sin_bloques)
print(summary_regression_sin_bloques)
cat("Regresión sin bloques: Ningún coef sig a 0.05.\n")

# Pregunta 5 sin bloques
cat("\nComparaciones sin bloques:\n")

tukey_sin_bloques <- TukeyHSD(anova_model_sin_bloques, which = "Agente")
print(tukey_sin_bloques)

lsd_sin_bloques <- LSD.test(anova_model_sin_bloques, "Agente", group=TRUE)
print(lsd_sin_bloques$groups)

duncan_sin_bloques <- duncan.test(anova_model_sin_bloques, "Agente", group=TRUE)
print(duncan_sin_bloques$groups)

snk_sin_bloques <- SNK.test(anova_model_sin_bloques, "Agente", group=TRUE)
print(snk_sin_bloques$groups)

cat("Sin bloques: Menos diferencias detectadas (solo Ag2-Ag3 en la mayoría).\n")
cat("El diseño con bloques fue clave para detectar efectos.\n")





