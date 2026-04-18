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

# Convertir Agente y Rollo a factores (importante para que R los trate correctamente)
long_df$Agente <- as.factor(long_df$Agente)
long_df$Rollo  <- as.factor(long_df$Rollo)

# ANOVA del DBCA
anova_model <- aov(Resistencia ~ Agente + Rollo, data = long_df)
anova_summary <- summary(anova_model)
print(anova_summary)


p_agente <- anova_summary[[1]]["Agente", "Pr(>F)"]
p_rollo  <- anova_summary[[1]]["Rollo",  "Pr(>F)"]

#Se extraen P-valores para verificar hipótesis

cat(sprintf("\nP-valor Agente: %.6f\n", p_agente))
cat(sprintf("P-valor Rollo:  %.6f\n", p_rollo))

# Conclusión automática según alpha
for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\n--- Con α = %.2f ---\n", alpha))
  cat(sprintf("Agente: %s\n", ifelse(p_agente < alpha, "SIGNIFICATIVO ✓", "No significativo ✗")))
  cat(sprintf("Rollo:  %s\n", ifelse(p_rollo  < alpha, "SIGNIFICATIVO ✓", "No significativo ✗")))
}

# Modelo de Regresión lineal

# Asegurar que las variables categóricas sean factores
long_df$Agente <- as.factor(long_df$Agente)
long_df$Rollo  <- as.factor(long_df$Rollo)

# Ajustar el modelo de regresión lineal
regression_model <- lm(Resistencia ~ Agente + Rollo, data = long_df)

# Resumen del modelo
regression_summary <- summary(regression_model)
print(regression_summary)

# Extraer R² y R² ajustado
r2     <- regression_summary$r.squared
r2_adj <- regression_summary$adj.r.squared

cat(sprintf("\nR²          = %.4f\n", r2))
cat(sprintf("R² ajustado = %.4f\n", r2_adj))

# Extraer p-values de los coeficientes
p_values <- regression_summary$coefficients[, 4]

cat("\nP-values de los coeficientes:\n")
print(p_values)

# Evaluación con alfa = 0.05
cat("\nSignificancia con alfa = 0.05:\n")
print(p_values < 0.05)

# Evaluación con alfa = 0.01
cat("\nSignificancia con alfa = 0.01:\n")
print(p_values < 0.01)

# ----------------------------------------------------------------------
# Punto 3: Validación del modelo de regresión (DBCA)
# ----------------------------------------------------------------------

required_packages <- c("car", "lmtest", "MASS")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

long_df$Agente <- factor(long_df$Agente)
long_df$Rollo  <- factor(long_df$Rollo)

modelo_dbca <- lm(Resistencia ~ Agente + Rollo, data = long_df)

cat("\n====================================\n")
cat("PUNTO 3: VALIDACIÓN DEL MODELO DBCA\n")
cat("====================================\n")

# ----------------------------------------------------------------------
# 1. Resumen del modelo
# ----------------------------------------------------------------------
cat("\n--- Resumen del modelo ---\n")
print(summary(modelo_dbca))

cat("\n--- ANOVA del modelo ---\n")
print(anova(modelo_dbca))

# ----------------------------------------------------------------------
# 2. Gráficos de diagnóstico (los 4 clásicos de R)
# ----------------------------------------------------------------------
par(mfrow = c(2, 2))
plot(modelo_dbca, main = "Diagnóstico del modelo DBCA")
par(mfrow = c(1, 1))

# ----------------------------------------------------------------------
# 3. Extracción de residuos y medidas de diagnóstico
# ----------------------------------------------------------------------
residuos  <- resid(modelo_dbca)
ajustados <- fitted(modelo_dbca)
rstand    <- rstandard(modelo_dbca)
rstud     <- rstudent(modelo_dbca)
cook      <- cooks.distance(modelo_dbca)
leverage  <- hatvalues(modelo_dbca)

# ----------------------------------------------------------------------
# 4. Pruebas de supuestos
# ----------------------------------------------------------------------

# 4.1 Normalidad — Shapiro-Wilk
cat("\n--- Prueba de normalidad (Shapiro-Wilk) ---\n")
shapiro_res <- shapiro.test(residuos)
print(shapiro_res)

for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\nCon α = %.2f: %s\n", alpha,
              ifelse(shapiro_res$p.value > alpha,
                     "No se rechaza H0 → Residuos normales ✓",
                     "Se rechaza H0 → Evidencia de no normalidad ✗")))
}

# QQ-plot manual para visualizar normalidad
qqnorm(residuos, main = "QQ-Plot de Residuos")
qqline(residuos, col = "red", lwd = 2)

# 4.2 Homocedasticidad — Breusch-Pagan
cat("\n--- Prueba de homocedasticidad (Breusch-Pagan) ---\n")
bp_res <- bptest(modelo_dbca)
print(bp_res)

for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\nCon α = %.2f: %s\n", alpha,
              ifelse(bp_res$p.value > alpha,
                     "No se rechaza H0 → Varianza constante ✓",
                     "Se rechaza H0 → Heterocedasticidad ✗")))
}

# 4.3 Levene por Agente y por Rollo
cat("\n--- Prueba de Levene por Agente ---\n")
lev_agente <- leveneTest(Resistencia ~ Agente, data = long_df)
print(lev_agente)

cat("\n--- Prueba de Levene por Rollo ---\n")
lev_rollo <- leveneTest(Resistencia ~ Rollo, data = long_df)
print(lev_rollo)

# 4.4 Independencia — Durbin-Watson  ← AGREGADO: faltaba en tu código
cat("\n--- Prueba de independencia (Durbin-Watson) ---\n")
dw_res <- dwtest(modelo_dbca)
print(dw_res)

for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\nCon α = %.2f: %s\n", alpha,
              ifelse(dw_res$p.value > alpha,
                     "No se rechaza H0 → Residuos independientes ✓",
                     "Se rechaza H0 → Evidencia de autocorrelación ✗")))
}

# ----------------------------------------------------------------------
# 5. Significancia de parámetros
# ----------------------------------------------------------------------
cat("\n--- Coeficientes del modelo ---\n")
coef_tabla <- summary(modelo_dbca)$coefficients
print(coef_tabla)

cat("\n--- Coeficientes significativos al 5% ---\n")
sig_05 <- coef_tabla[coef_tabla[, 4] < 0.05, ]
print(sig_05)

cat("\n--- Coeficientes significativos al 1% ---\n")  # ← AGREGADO
sig_01 <- coef_tabla[coef_tabla[, 4] < 0.01, ]
print(sig_01)

# ----------------------------------------------------------------------
# 6. Detección de observaciones influyentes y atípicas
# ----------------------------------------------------------------------
cat("\n--- Residuos estandarizados ---\n")
print(round(rstand, 4))

cat("\n--- Residuos studentizados ---\n")
print(round(rstud, 4))

outliers <- which(abs(rstud) > 2)
cat("\nObservaciones con |residuo studentizado| > 2:\n")
if (length(outliers) == 0) cat("Ninguna\n") else print(outliers)

cat("\n--- Distancia de Cook ---\n")
print(round(cook, 4))

umbral_cook <- 4 / nrow(long_df)
cat(sprintf("\nUmbral Cook = %.4f\n", umbral_cook))
influyentes <- which(cook > umbral_cook)
cat("Observaciones potencialmente influyentes:\n")
if (length(influyentes) == 0) cat("Ninguna\n") else print(influyentes)

cat("\n--- Leverage ---\n")
print(round(leverage, 4))

umbral_lev <- 2 * length(coef(modelo_dbca)) / nrow(long_df)
cat(sprintf("\nUmbral leverage = %.4f\n", umbral_lev))
lev_altos <- which(leverage > umbral_lev)
cat("Observaciones con leverage alto:\n")
if (length(lev_altos) == 0) cat("Ninguna\n") else print(lev_altos)

cat("\n--- Prueba Bonferroni para outliers ---\n")
print(outlierTest(modelo_dbca))

# Gráfico de Cook's Distance  ← AGREGADO: útil para visualizar influyentes
plot(cook, type = "h", main = "Distancia de Cook",
     ylab = "Cook's Distance", xlab = "Observación")
abline(h = umbral_cook, col = "red", lty = 2)
text(influyentes, cook[influyentes], labels = influyentes, pos = 3, col = "red")

# Tabla resumen de diagnóstico
diagnostico <- data.frame(
  Obs         = 1:nrow(long_df),
  Agente      = long_df$Agente,
  Rollo       = long_df$Rollo,
  Resistencia = long_df$Resistencia,
  Ajustado    = round(ajustados, 4),
  Residuo     = round(residuos, 4),
  Rstand      = round(rstand, 4),
  Rstudent    = round(rstud, 4),
  Leverage    = round(leverage, 4),
  CookD       = round(cook, 4)
)

cat("\n--- Tabla de diagnóstico ---\n")
print(diagnostico)

# ----------------------------------------------------------------------
# 7. Forma final del modelo
# ----------------------------------------------------------------------
cat("\n--- Forma final del modelo ---\n")

if (shapiro_res$p.value > 0.05 &&
    bp_res$p.value > 0.05 &&
    dw_res$p.value > 0.05 &&
    length(outliers) == 0 &&
    length(influyentes) == 0) {
  
  cat("Todos los supuestos se cumplen. Se conserva el modelo aditivo del DBCA:\n")
  cat("Y_ij = mu + tau_i + beta_j + error_ij\n\n")
  cat("Donde:\n")
  cat("  mu     = media general\n")
  cat("  tau_i  = efecto del agente químico i (i = 1,2,3,4)\n")
  cat("  beta_j = efecto del rollo j (j = 1,2,3,4,5)\n")
  cat("  error  ~ N(0, sigma²)\n")
  
} else {
  cat("Se detectan posibles problemas. Ver sección de soluciones.\n")
}

# ----------------------------------------------------------------------
# 8. Posibles soluciones si no se cumplen los supuestos
# ----------------------------------------------------------------------
cat("\n--- Posibles soluciones si no se cumplen supuestos ---\n")
cat("1. Normalidad/varianza: transformación Box-Cox.\n")
cat("2. Transformaciones: log(Y), sqrt(Y) o 1/Y.\n")
cat("3. Heterocedasticidad persistente: mínimos cuadrados ponderados (WLS).\n")
cat("4. Observaciones influyentes: revisar errores de medición.\n")
cat("5. No normalidad severa: prueba no paramétrica de Friedman.\n")

# ----------------------------------------------------------------------
# 9. Gráfico Box-Cox
# ----------------------------------------------------------------------
cat("\n--- Gráfico Box-Cox ---\n")
bc <- boxcox(modelo_dbca, lambda = seq(-2, 2, by = 0.1))
lambda_optimo <- bc$x[which.max(bc$y)]
cat(sprintf("Lambda óptimo Box-Cox: %.4f\n", lambda_optimo))  # ← AGREGADO
  
#----------------------------------------------------------------------
# Punto 4: Método de Scheffé para contrastes
# ----------------------------------------------------------------------
cat("\n=================================\n")
cat("PUNTO 4: MÉTODO DE SCHEFFÉ\n")
cat("=================================\n")

library(agricolae)

# Asegurar factores
long_df$Agente <- factor(long_df$Agente)
long_df$Rollo  <- factor(long_df$Rollo)

# Modelo DBCA
anova_model <- aov(Resistencia ~ Agente + Rollo, data = long_df)

# Resultado de agricolae
scheffe_result <- scheffe.test(anova_model, "Agente", group = TRUE)
print(scheffe_result)

# ------------------------------------------------------------
# 1. Extraer valores correctos desde el modelo
# ------------------------------------------------------------
medias <- with(long_df, tapply(Resistencia, Agente, mean))
print(medias)

anova_tab <- anova(anova_model)
MSerror   <- anova_tab["Residuals", "Mean Sq"]
df_error  <- anova_tab["Residuals", "Df"]

k <- length(medias)                        # número de tratamientos
n <- length(unique(long_df$Rollo))         # repeticiones por tratamiento = bloques

cat("\nMSerror =", round(MSerror, 6), "\n")
cat("gl del error =", df_error, "\n")
cat("Número de tratamientos =", k, "\n")
cat("Número de repeticiones por tratamiento =", n, "\n")

# ------------------------------------------------------------
# 2. Valores críticos de Scheffé
# ------------------------------------------------------------
Fcrit_05 <- (k - 1) * qf(1 - 0.05, df1 = k - 1, df2 = df_error)
Fcrit_01 <- (k - 1) * qf(1 - 0.01, df1 = k - 1, df2 = df_error)

cat("\nValor crítico de Scheffé con α = 0.05:", round(Fcrit_05, 4), "\n")
cat("Valor crítico de Scheffé con α = 0.01:", round(Fcrit_01, 4), "\n")

# ------------------------------------------------------------
# 3. Comparaciones pareadas entre tratamientos
#    (subconjunto de contrastes)
# ------------------------------------------------------------
pares <- combn(names(medias), 2)

tabla_scheffe <- data.frame(
  Contraste = character(),
  Diferencia = numeric(),
  F_calculado = numeric(),
  Decision_0.05 = character(),
  Decision_0.01 = character(),
  stringsAsFactors = FALSE
)

for (m in 1:ncol(pares)) {
  i <- pares[1, m]
  j <- pares[2, m]
  
  diff_ij <- medias[i] - medias[j]
  
  # Para comparación entre dos medias en diseño balanceado:
  # F = (yi - yj)^2 / [MSerror * (1/n + 1/n)]
  Fcalc <- (diff_ij^2) / (MSerror * (1/n + 1/n))
  
  tabla_scheffe <- rbind(tabla_scheffe, data.frame(
    Contraste = paste(i, "vs", j),
    Diferencia = round(diff_ij, 4),
    F_calculado = round(Fcalc, 4),
    Decision_0.05 = ifelse(Fcalc > Fcrit_05, "Significativo", "No significativo"),
    Decision_0.01 = ifelse(Fcalc > Fcrit_01, "Significativo", "No significativo")
  ))
}

cat("\n=== Tabla de resultados de Scheffé ===\n")
print(tabla_scheffe, row.names = FALSE)

# ----------------------------------------------------------------------
# Punto 5: Comparaciones entre pares de medias de tratamientos
# Procedimientos: Tukey, LSD, Duncan y Newman-Keuls
# con alpha = 0.05 y alpha = 0.01
# ----------------------------------------------------------------------

library(agricolae)

long_df$Agente <- factor(long_df$Agente)
long_df$Rollo  <- factor(long_df$Rollo)

anova_model <- aov(Resistencia ~ Agente + Rollo, data = long_df)

cat("\n============================================\n")
cat("PUNTO 5: COMPARACIONES ENTRE PARES DE MEDIAS\n")
cat("============================================\n")

cat("\n--- Medias por tratamiento ---\n")
medias_agente <- with(long_df, tapply(Resistencia, Agente, mean))
print(medias_agente)

for (alpha in c(0.05, 0.01)) {
  
  cat("\n============================================\n")
  cat(sprintf("RESULTADOS CON alpha = %.2f\n", alpha))
  cat("============================================\n")
  
  cat("\n--- Tukey ---\n")
  tukey_res <- HSD.test(anova_model, "Agente", alpha = alpha, group = TRUE)
  print(tukey_res$groups)
  
  cat("\n--- LSD ---\n")
  lsd_res <- LSD.test(anova_model, "Agente", alpha = alpha, group = TRUE, p.adj = "none")
  print(lsd_res$groups)
  
  cat("\n--- Duncan ---\n")
  duncan_res <- duncan.test(anova_model, "Agente", alpha = alpha, group = TRUE)
  print(duncan_res$groups)
  
  cat("\n--- Newman-Keuls ---\n")
  snk_res <- SNK.test(anova_model, "Agente", alpha = alpha, group = TRUE)
  print(snk_res$groups)
}
# ----------------------------------------------------------------------
# Punto 6: Análisis del efecto de los bloques
# ----------------------------------------------------------------------
cat("\n====================================\n")
cat("PUNTO 6: ANÁLISIS DEL EFECTO DE LOS BLOQUES\n")
cat("====================================\n")

library(agricolae)

# Asegurar factores
long_df$Agente <- factor(long_df$Agente)
long_df$Rollo  <- factor(long_df$Rollo)

# ------------------------------------------------------------
# 1. Modelo con bloques (DBCA)
# ------------------------------------------------------------
modelo_con_bloques <- aov(Resistencia ~ Agente + Rollo, data = long_df)
anova_con_bloques <- summary(modelo_con_bloques)
print(anova_con_bloques)

p_agente_cb <- anova_con_bloques[[1]]["Agente", "Pr(>F)"]
p_rollo_cb  <- anova_con_bloques[[1]]["Rollo",  "Pr(>F)"]
mse_cb      <- anova_con_bloques[[1]]["Residuals", "Mean Sq"]

cat("\n--- Interpretación con bloques ---\n")
for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\nCon α = %.2f\n", alpha))
  cat(sprintf("Agente: %s\n",
              ifelse(p_agente_cb < alpha, "Significativo", "No significativo")))
  cat(sprintf("Rollo: %s\n",
              ifelse(p_rollo_cb < alpha, "Significativo", "No significativo")))
}

cat(sprintf("\nMSE con bloques = %.4f\n", mse_cb))

# ------------------------------------------------------------
# 2. Modelo sin bloques
# ------------------------------------------------------------
modelo_sin_bloques <- aov(Resistencia ~ Agente, data = long_df)
anova_sin_bloques <- summary(modelo_sin_bloques)
print(anova_sin_bloques)

p_agente_sb <- anova_sin_bloques[[1]]["Agente", "Pr(>F)"]
mse_sb      <- anova_sin_bloques[[1]]["Residuals", "Mean Sq"]

cat("\n--- Interpretación sin bloques ---\n")
for (alpha in c(0.05, 0.01)) {
  cat(sprintf("\nCon α = %.2f\n", alpha))
  cat(sprintf("Agente: %s\n",
              ifelse(p_agente_sb < alpha, "Significativo", "No significativo")))
}
cat(sprintf("\nMSE sin bloques = %.4f\n", mse_sb))

# ------------------------------------------------------------
# 3. Comparación del efecto de usar bloques
# ------------------------------------------------------------
cat("\n--- Comparación entre modelos ---\n")
cat(sprintf("MSE con bloques    = %.4f\n", mse_cb))
cat(sprintf("MSE sin bloques    = %.4f\n", mse_sb))
cat(sprintf("Reducción del MSE  = %.4f\n", mse_sb - mse_cb))
cat(sprintf("Razón MSE sin/con  = %.4f\n", mse_sb / mse_cb))

if (mse_sb > mse_cb) {
  cat("El uso de bloques redujo la variabilidad residual y mejoró la precisión del análisis.\n")
} else {
  cat("El uso de bloques no redujo la variabilidad residual.\n")
}

# ------------------------------------------------------------
# 4. Modelo de regresión sin bloques
# ------------------------------------------------------------
reg_sin_bloques <- lm(Resistencia ~ Agente, data = long_df)
summary_reg_sb <- summary(reg_sin_bloques)

cat("\n--- Regresión sin bloques ---\n")
print(summary_reg_sb)

coef_tab_sb <- summary_reg_sb$coefficients

cat("\nCoeficientes significativos al 5%:\n")
print(coef_tab_sb[,4] < 0.05)

cat("\nCoeficientes significativos al 1%:\n")
print(coef_tab_sb[,4] < 0.01)

# ------------------------------------------------------------
# 5. Comparaciones múltiples sin bloques
# ------------------------------------------------------------
cat("\n--- Comparaciones múltiples sin bloques ---\n")

cat("\nTukey:\n")
tukey_sb <- HSD.test(modelo_sin_bloques, "Agente", group = TRUE)
print(tukey_sb$groups)

cat("\nLSD:\n")
lsd_sb <- LSD.test(modelo_sin_bloques, "Agente", group = TRUE, p.adj = "none")
print(lsd_sb$groups)

cat("\nDuncan:\n")
duncan_sb <- duncan.test(modelo_sin_bloques, "Agente", group = TRUE)
print(duncan_sb$groups)

cat("\nNewman-Keuls:\n")
snk_sb <- SNK.test(modelo_sin_bloques, "Agente", group = TRUE)
print(snk_sb$groups)


