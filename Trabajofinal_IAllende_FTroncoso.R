# Carga de librerias
library(tidyverse)
library(knitr)
library(lmtest)
library(car)

install.packages("tseries")
library(tseries)
library(glmnet)

###PREGUNTA1##
eme_raw <- read_csv("EME8_base_full_original.csv")
pliegues <- read_csv("ids_particiones_pliegues.csv")




# Cruce de bases
df_modelo <- eme_raw %>%
  inner_join(pliegues, by = "id")
nrow(df_modelo)
# Verificar particiones
table(df_modelo$particion)

df_entrenamiento <- df_modelo %>% filter(particion == "entrenamiento")
df_prueba  <- df_modelo %>% filter(particion == "prueba")

# CANTIDADES NULAS
sum(df_entrenamiento$ganancia_final <= 0, na.rm = TRUE)


# Eliminar variables con fuga de información
df_entrenamiento <- df_entrenamiento %>%
  select(-gasto_mensual, -ingreso_mensual, -ganancia_estimada, -ganancia_declarada)

head(df_entrenamiento)

####################### PREGUNTA 2 ############################################
# a) Estadísticas resumen (Tabla formal)
resumen_ganancia <- df_entrenamiento %>%
  summarise(
    Media = mean(ganancia_final, na.rm = TRUE),
    Mediana = median(ganancia_final, na.rm = TRUE),
    DS = sd(ganancia_final, na.rm = TRUE),
    Minimo = min(ganancia_final, na.rm = TRUE),
    Maximo = max(ganancia_final, na.rm = TRUE),
    P25 = quantile(ganancia_final, 0.25, na.rm = TRUE),
    P75 = quantile(ganancia_final, 0.75, na.rm = TRUE)
  )
resumen_ganancia


# b) Gráficos de dispersión (Horas y Trabajadores)
# Gráfico 1
p1 <- ggplot(df_entrenamiento, aes(x = c6_3, y = ganancia_final)) +
  geom_point(alpha = 0.4, color = "darkblue") +
  scale_y_continuous(labels = scales::label_comma(big.mark = ".")) +
  labs(title = "Ganancia vs horas de trabajo", x = "Horas semanales", y = "Ganancia (CLP)") +
  theme_minimal()

# Gráfico 2
p2 <- ggplot(df_entrenamiento, aes(x = total_trabajadores, y = ganancia_final)) +
  geom_point(alpha = 0.4, color = "darkolivegreen") +
  scale_y_continuous(labels = scales::label_comma(big.mark = ".")) +
  labs(title = "Ganancia vs total trabajadores", x = "N° de Trabajadores", y = "Ganancia (CLP)") +
  theme_minimal()

p1
p2


# c) Boxplots comparativos

df_entrenamiento %>%
  mutate(
    Informalidad_Factor = ifelse(informalidad == 1, "Informal", "Formal"),
    Sexo_Factor = ifelse(sexo == 1, "Hombre", "Mujer")
  ) %>%
  drop_na(Informalidad_Factor, Sexo_Factor, ganancia_final) %>%
  ggplot(aes(x = Informalidad_Factor, y = asinh(ganancia_final), fill = Sexo_Factor)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(title = "Distribución de Ganancia por Formalidad y Sexo",
       x = "Formalidad del Negocio", y = "Ganancia (asinh)", fill = "Sexo") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")




#Preg 3



# 1. Asegurar que las variables categóricas sean factores antes de modelar
df_entrenamiento <- df_entrenamiento %>%
  mutate(
    across(c(cine_eme_red, informalidad, c1_b, sexo, tramo_etario, region), as.factor)
  )

# 2. Ajustar el modelo OLS inicial
modelo_inicial <- lm(ganancia_final ~ c6_3 + c7 + total_trabajadores + 
                       cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                     data = df_entrenamiento)

# Ver los resultados del modelo (estimaciones y p-values)
summary(modelo_inicial)

# 3. Análisis visual de residuos
par(mfrow = c(1, 2))
plot(modelo_inicial, which = 1, main = "Residuos vs Ajustados")
plot(modelo_inicial, which = 2, main = "Gráfico Q-Q Normal")
par(mfrow = c(1, 1))


# 4. Pruebas formales de supuestos
# Homocedasticidad (Breusch-Pagan)
bptest(modelo_inicial)

# Normalidad de los residuos (Jarque-Bera)
jarque.bera.test(na.omit(residuals(modelo_inicial)))

# Colinealidad (VIF)
vif(modelo_inicial)

# Especificación del modelo (Test RESET de Ramsey)
resettest(modelo_inicial, power = 2:3, type = "fitted")





###preg 4


# --- 4.b) Ajuste de dos modelos transformados ---

# Alternativa 1: asinh(Y / 100000)
df_entrenamiento <- df_entrenamiento %>%
  mutate(y_asinh = asinh(ganancia_final / 100000))

mod_asinh <- lm(y_asinh ~ c6_3 + c7 + total_trabajadores + 
                  cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                data = df_entrenamiento)

# Alternativa 2: Yeo-Johnson (estimando lambda óptimo en entrenamiento)
pt_yj <- powerTransform(df_entrenamiento$ganancia_final ~ 1, family = "yjPower")
lambda_opt <- pt_yj$lambda
lambda_opt # Ver el lambda estimado

df_entrenamiento <- df_entrenamiento %>%
  mutate(y_yj = yjPower(ganancia_final, lambda_opt))

mod_yj <- lm(y_yj ~ c6_3 + c7 + total_trabajadores + 
               cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
             data = df_entrenamiento)

# Diagnósticos comparativos de los modelos transformados
summary(mod_asinh)$r.squared
bptest(mod_asinh)
jarque.bera.test(na.omit(residuals(mod_asinh)))
resettest(mod_asinh, power = 2:3, type = "fitted")

summary(mod_yj)$r.squared
bptest(mod_yj)
jarque.bera.test(na.omit(residuals(mod_yj)))
resettest(mod_yj, power = 2:3, type = "fitted")

# --- 4.c) Retorno a CLP (Inversión Directa vs Media Condicional) y CV (5 pliegues) ---

# Función auxiliar para invertir Yeo-Johnson exactamente
inv_yj <- function(z, lambda) {
  y <- numeric(length(z))
  pos <- z >= 0
  if (abs(lambda) > 1e-6) {
    y[pos] <- (z[pos] * lambda + 1)^(1 / lambda) - 1
  } else {
    y[pos] <- exp(z[pos]) - 1
  }
  if (abs(2 - lambda) > 1e-6) {
    y[!pos] <- 1 - (-(2 - lambda) * z[!pos] + 1)^(1 / (2 - lambda))
  } else {
    y[!pos] <- 1 - exp(-z[!pos])
  }
  return(y)
}

# Validación Cruzada en los 5 pliegues predefinidos
pliegues_unicos <- sort(unique(df_entrenamiento$fold))
resultados_cv <- data.frame()

for (k in pliegues_unicos) {
  train_k <- df_entrenamiento %>% filter(fold != k)
  val_k   <- df_entrenamiento %>% filter(fold == k)
  
  # 1. Modelo asinh en pliegue k
  m_asinh_k <- lm(asinh(ganancia_final / 100000) ~ c6_3 + c7 + total_trabajadores + 
                    cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                  data = train_k)
  
  pred_z_asinh <- predict(m_asinh_k, newdata = val_k)
  res_asinh_k  <- na.omit(residuals(m_asinh_k))
  
  # Inversión directa (plug-in): 100000 * sinh(z_hat)
  pred_clp_asinh_dir <- 100000 * sinh(pred_z_asinh)
  
  # Media condicional (Duan's Smearing): E[sinh(z_hat + e)]
  # Usando identidad: sinh(a + e) = sinh(a)*cosh(e) + cosh(a)*sinh(e)
  pred_clp_asinh_smear <- 100000 * (sinh(pred_z_asinh) * mean(cosh(res_asinh_k)) + 
                                      cosh(pred_z_asinh) * mean(sinh(res_asinh_k)))
  
  # 2. Modelo Yeo-Johnson en pliegue k (estimando lambda SOLO en train_k para evitar fuga)
  lam_k <- powerTransform(train_k$ganancia_final ~ 1, family = "yjPower")$lambda
  train_k$y_yj_k <- yjPower(train_k$ganancia_final, lam_k)
  
  m_yj_k <- lm(y_yj_k ~ c6_3 + c7 + total_trabajadores + 
                 cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
               data = train_k)
  
  pred_z_yj <- predict(m_yj_k, newdata = val_k)
  res_yj_k  <- na.omit(residuals(m_yj_k))
  
  # Inversión directa (plug-in)
  pred_clp_yj_dir <- inv_yj(pred_z_yj, lam_k)
  
  # Media condicional (Smearing numérico: promedio de inv_yj(z_hat + e_i))
  pred_clp_yj_smear <- sapply(pred_z_yj, function(z_val) mean(inv_yj(z_val + res_yj_k, lam_k)))
  
  # Guardar errores en CLP (eliminando posibles NAs de predictores faltantes)
  idx <- !is.na(pred_clp_asinh_dir) & !is.na(val_k$ganancia_final)
  y_real <- val_k$ganancia_final[idx]
  
  resultados_cv <- rbind(resultados_cv, data.frame(
    Fold = k,
    MAE_asinh_Dir    = mean(abs(y_real - pred_clp_asinh_dir[idx])),
    RMSE_asinh_Dir   = sqrt(mean((y_real - pred_clp_asinh_dir[idx])^2)),
    MAE_asinh_Smear  = mean(abs(y_real - pred_clp_asinh_smear[idx])),
    RMSE_asinh_Smear = sqrt(mean((y_real - pred_clp_asinh_smear[idx])^2)),
    MAE_yj_Dir       = mean(abs(y_real - pred_clp_yj_dir[idx])),
    RMSE_yj_Dir      = sqrt(mean((y_real - pred_clp_yj_dir[idx])^2)),
    MAE_yj_Smear     = mean(abs(y_real - pred_clp_yj_smear[idx])),
    RMSE_yj_Smear    = sqrt(mean((y_real - pred_clp_yj_smear[idx])^2))
  ))
}

# Promedio de métricas CV en CLP
colMeans(resultados_cv[, -1])






##### PREGUNTA 5

# --- 5.a) Explorar log(1 + c6_3) y log(1 + total_trabajadores) ---
mod_p5_log <- lm(y_asinh ~ log1p(c6_3) + c7 + log1p(total_trabajadores) + 
                   cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                 data = df_entrenamiento)

summary(mod_p5_log)

# --- 5.b) Formas no lineales de meses de funcionamiento (c7) con jerarquía ---
# Usamos poly(c7, 2, raw = TRUE) y poly(c7, 3, raw = TRUE) para mantener c7 + c7^2 (+ c7^3)
mod_p5_poly2 <- lm(y_asinh ~ log1p(c6_3) + poly(c7, 2, raw = TRUE) + log1p(total_trabajadores) + 
                     cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                   data = df_entrenamiento)

mod_p5_poly3 <- lm(y_asinh ~ log1p(c6_3) + poly(c7, 3, raw = TRUE) + log1p(total_trabajadores) + 
                     cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
                   data = df_entrenamiento)

# Comparación anidada (Test F secuencial para jerarquía polinómica)
anova(mod_p5_log, mod_p5_poly2, mod_p5_poly3)
summary(mod_p5_poly2)$coefficients[1:5, ]

# --- 5.c) Comparación: Diagnósticos y Error fuera de muestra (CV 5 pliegues en CLP) ---

# 1. Diagnósticos en entrenamiento (R2 ajustado, RESET, Breusch-Pagan, Jarque-Bera)
diagnosticos_p5 <- data.frame(
  Modelo = c("Base (asinh Y, X original)", "Log(1+X)", "Log(1+X) + Poly2(c7)", "Log(1+X) + Poly3(c7)"),
  R2_Adj = c(summary(mod_asinh)$adj.r.squared,
             summary(mod_p5_log)$adj.r.squared,
             summary(mod_p5_poly2)$adj.r.squared,
             summary(mod_p5_poly3)$adj.r.squared),
  RESET_stat = c(resettest(mod_asinh, power = 2:3, type = "fitted")$statistic,
                 resettest(mod_p5_log, power = 2:3, type = "fitted")$statistic,
                 resettest(mod_p5_poly2, power = 2:3, type = "fitted")$statistic,
                 resettest(mod_p5_poly3, power = 2:3, type = "fitted")$statistic),
  RESET_pval = c(resettest(mod_asinh, power = 2:3, type = "fitted")$p.value,
                 resettest(mod_p5_log, power = 2:3, type = "fitted")$p.value,
                 resettest(mod_p5_poly2, power = 2:3, type = "fitted")$p.value,
                 resettest(mod_p5_poly3, power = 2:3, type = "fitted")$p.value),
  BP_pval = c(bptest(mod_asinh)$p.value,
              bptest(mod_p5_log)$p.value,
              bptest(mod_p5_poly2)$p.value,
              bptest(mod_p5_poly3)$p.value)
)
print(diagnosticos_p5)

# 2. Validación Cruzada en los 5 pliegues (métrica en CLP: Directa y Smearing)
formulas_p5 <- list(
  Base_asinh = asinh(ganancia_final / 100000) ~ c6_3 + c7 + total_trabajadores + cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region,
  Log_X      = asinh(ganancia_final / 100000) ~ log1p(c6_3) + c7 + log1p(total_trabajadores) + cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region,
  Log_Poly2  = asinh(ganancia_final / 100000) ~ log1p(c6_3) + poly(c7, 2, raw = TRUE) + log1p(total_trabajadores) + cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region,
  Log_Poly3  = asinh(ganancia_final / 100000) ~ log1p(c6_3) + poly(c7, 3, raw = TRUE) + log1p(total_trabajadores) + cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region
)

cv_p5 <- data.frame()

for (nombre in names(formulas_p5)) {
  form_actual <- formulas_p5[[nombre]]
  mae_dir <- numeric(length(pliegues_unicos))
  rmse_dir <- numeric(length(pliegues_unicos))
  mae_smear <- numeric(length(pliegues_unicos))
  rmse_smear <- numeric(length(pliegues_unicos))
  
  for (i in seq_along(pliegues_unicos)) {
    k <- pliegues_unicos[i]
    train_k <- df_entrenamiento %>% filter(fold != k)
    val_k   <- df_entrenamiento %>% filter(fold == k)
    
    mod_k <- lm(form_actual, data = train_k)
    pred_z <- predict(mod_k, newdata = val_k)
    res_k  <- na.omit(residuals(mod_k))
    
    pred_dir   <- 100000 * sinh(pred_z)
    pred_smear <- 100000 * (sinh(pred_z) * mean(cosh(res_k)) + cosh(pred_z) * mean(sinh(res_k)))
    
    idx <- !is.na(pred_dir) & !is.na(val_k$ganancia_final)
    y_real <- val_k$ganancia_final[idx]
    
    mae_dir[i]    <- mean(abs(y_real - pred_dir[idx]))
    rmse_dir[i]   <- sqrt(mean((y_real - pred_dir[idx])^2))
    mae_smear[i]  <- mean(abs(y_real - pred_smear[idx]))
    rmse_smear[i] <- sqrt(mean((y_real - pred_smear[idx])^2))
  }
  
  cv_p5 <- rbind(cv_p5, data.frame(
    Modelo = nombre,
    MAE_Directo_CLP   = mean(mae_dir),
    RMSE_Directo_CLP  = mean(rmse_dir),
    MAE_Smearing_CLP  = mean(mae_smear),
    RMSE_Smearing_CLP = mean(rmse_smear)
  ))
}

print(cv_p5)






# ==============================================================================
# PREGUNTA 6: Observaciones atípicas e influyentes
# ==============================================================================

# Trabajamos sobre el mejor modelo de la P5 (mod_p5_log)
# Asegurar filas completas usadas por el modelo para alinear índices con df_entrenamiento
vars_mod <- all.vars(formula(mod_p5_log))
df_p6 <- df_entrenamiento %>%
  filter(if_all(all_of(vars_mod), ~ !is.na(.)))

# Reajustar exactamente sobre df_p6 para garantizar correspondencia 1 a 1 con los IDs
mod_p6 <- lm(formula(mod_p5_log), data = df_p6)

n_obs <- nobs(mod_p6)
p_par <- length(coef(mod_p6))

# --- 6.a) Cálculo de métricas y umbrales ---
df_p6 <- df_p6 %>%
  mutate(
    r_stud = rstudent(mod_p6),
    leverage = hatvalues(mod_p6),
    cooks_d = cooks.distance(mod_p6)
  )

umbral_rstud <- 3
umbral_lev   <- 2 * p_par / n_obs
umbral_cook  <- 4 / n_obs

cat("Umbral |t_i| >", umbral_rstud, ":", sum(abs(df_p6$r_stud) > umbral_rstud), "casos\n")
cat("Umbral Leverage >", round(umbral_lev, 5), ":", sum(df_p6$leverage > umbral_lev), "casos\n")
cat("Umbral Cook >", round(umbral_cook, 6), ":", sum(df_p6$cooks_d > umbral_cook), "casos\n")

# Gráfico diagnóstico de influencia (Residuos Estudentizados vs Leverage, tamaño = Cook)
plot(df_p6$leverage, df_p6$r_stud,
     xlab = "Leverage (h_ii)", ylab = "Residuos Estudentizados (t_i)",
     main = "Diagnóstico de Influencia (proporcional a Cook)",
     pch = 21, bg = rgb(0.2, 0.4, 0.6, 0.3), cex = sqrt(df_p6$cooks_d) * 15)
abline(h = c(-3, 3), col = "red", lty = 2)
abline(v = umbral_lev, col = "blue", lty = 2)

# --- 6.b) Documentar los 5 casos más influyentes por ID (sin datos personales) ---
top5_influyentes <- df_p6 %>%
  arrange(desc(cooks_d)) %>%
  slice(1:5) %>%
  select(id, ganancia_final, y_asinh, c6_3, total_trabajadores, r_stud, leverage, cooks_d)

print(top5_influyentes)

# --- 6.c) Análisis de sensibilidad (Con vs Sin los 5 casos más influyentes) ---
ids_top5 <- top5_influyentes$id
df_p6_sin_top5 <- df_p6 %>% filter(!id %in% ids_top5)

mod_p6_sens <- lm(formula(mod_p5_log), data = df_p6_sin_top5)

# Comparar coeficientes clave y bondad de ajuste
comparacion_sens <- data.frame(
  Termino = names(coef(mod_p6)),
  Coef_Con_Todos = coef(mod_p6),
  Coef_Sin_Top5  = coef(mod_p6_sens),
  Var_Porcentual = round(((coef(mod_p6_sens) - coef(mod_p6)) / abs(coef(mod_p6))) * 100, 2),
  Pval_Con_Todos = summary(mod_p6)$coefficients[, 4],
  Pval_Sin_Top5  = summary(mod_p6_sens)$coefficients[, 4]
)

# Mostrar los principales predictores (primeros 16 términos: intercepto, continuas, educación, informalidad, sector, sexo)
print(comparacion_sens[1:16, ], row.names = FALSE)

cat("R2 Ajustado (Todos):", summary(mod_p6)$adj.r.squared, "\n")
cat("R2 Ajustado (Sin Top 5):", summary(mod_p6_sens)$adj.r.squared, "\n")





# ==============================================================================
# PREGUNTA 7: Selección clásica (Forward, Backward, Stepwise) y Criterios
# ==============================================================================

# 1. Asegurar exactamente la misma muestra (sin NAs en los predictores candidatos)
vars_p7 <- all.vars(formula(mod_p5_log))
df_p7 <- df_entrenamiento %>%
  filter(if_all(all_of(vars_p7), ~ !is.na(.)))

n_p7 <- nrow(df_p7)

# Definir modelo nulo (solo intercepto) y modelo completo (todos los predictores de P5)
mod_nulo <- lm(y_asinh ~ 1, data = df_p7)
mod_full <- lm(y_asinh ~ log1p(c6_3) + c7 + log1p(total_trabajadores) + 
                 cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region, 
               data = df_p7)

# --- 7.a) Selección por bloques usando AIC (k = 2) y BIC (k = log(n)) ---

# Con criterio AIC (k = 2)
step_forward_aic  <- step(mod_nulo, scope = list(lower = mod_nulo, upper = mod_full), 
                          direction = "forward", k = 2, trace = 0)
step_backward_aic <- step(mod_full, direction = "backward", k = 2, trace = 0)
step_both_aic     <- step(mod_nulo, scope = list(lower = mod_nulo, upper = mod_full), 
                          direction = "both", k = 2, trace = 0)

# Con criterio BIC (k = log(n)) para contrastar penalización por complejidad
step_forward_bic  <- step(mod_nulo, scope = list(lower = mod_nulo, upper = mod_full), 
                          direction = "forward", k = log(n_p7), trace = 0)
step_backward_bic <- step(mod_full, direction = "backward", k = log(n_p7), trace = 0)
step_both_bic     <- step(mod_nulo, scope = list(lower = mod_nulo, upper = mod_full), 
                          direction = "both", k = log(n_p7), trace = 0)

# Revisar rutas (tablas ANOVA de pasos) y términos seleccionados

print(step_forward_aic$anova)


print(step_backward_aic$anova)

print(step_both_aic$anova)


print(step_both_bic$anova)

# --- 7.b) Comparación formal de AIC, BIC y Cp de Mallows ---
# Varianza residual de referencia (del modelo completo)
sigma2_ref <- summary(mod_full)$sigma^2
cat("\nVarianza residual de referencia (sigma^2_full):", round(sigma2_ref, 5), "\n")
cat("Parámetros totales modelo completo (p_full):", length(coef(mod_full)), "\n")

# Función para extraer métricas comparables (p = número total de parámetros incl. intercepto)
evaluar_criterios <- function(modelo, nombre, sigma2_full, n) {
  rss <- sum(residuals(modelo)^2)
  p   <- length(coef(modelo)) # Incluye intercepto + betas
  cp  <- (rss / sigma2_full) - (n - 2 * p)
  data.frame(
    Modelo = nombre,
    p_parametros = p,
    RSS = round(rss, 2),
    R2_Adj = round(summary(modelo)$adj.r.squared, 4),
    AIC = round(AIC(modelo), 2),
    BIC = round(BIC(modelo), 2),
    Cp_Mallows = round(cp, 2)
  )
}

tabla_p7 <- rbind(
  evaluar_criterios(mod_nulo,          "0. Modelo Nulo",            sigma2_ref, n_p7),
  evaluar_criterios(mod_full,          "1. Modelo Completo",        sigma2_ref, n_p7),
  evaluar_criterios(step_forward_aic,  "2. Forward (AIC)",          sigma2_ref, n_p7),
  evaluar_criterios(step_backward_aic, "3. Backward (AIC)",         sigma2_ref, n_p7),
  evaluar_criterios(step_both_aic,     "4. Stepwise Both (AIC)",    sigma2_ref, n_p7),
  evaluar_criterios(step_forward_bic,  "5. Forward (BIC)",          sigma2_ref, n_p7),
  evaluar_criterios(step_backward_bic, "6. Backward (BIC)",         sigma2_ref, n_p7),
  evaluar_criterios(step_both_bic,     "7. Stepwise Both (BIC)",    sigma2_ref, n_p7)
)

print(tabla_p7, row.names = FALSE)




# ==============================================================================
# PREGUNTA 8: Validación cruzada detallada (5 pliegues predefinidos)
# ==============================================================================

# --- 8.a) Verificar estructura de los 5 pliegues en entrenamiento ---
table(df_entrenamiento$fold)

# --- 8.b) Reportar MAE y RMSE por pliegue y promedios en CLP ---
pliegues_unicos <- sort(unique(df_entrenamiento$fold))
tabla_cv_pliegues <- data.frame()

formula_ganadora <- asinh(ganancia_final / 100000) ~ log1p(c6_3) + c7 + log1p(total_trabajadores) + 
  cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region

for (k in pliegues_unicos) {
  train_k <- df_entrenamiento %>% filter(fold != k)
  val_k   <- df_entrenamiento %>% filter(fold == k)
  
  # Ajustar modelo SOLO con los 4 pliegues de entrenamiento
  mod_k <- lm(formula_ganadora, data = train_k)
  
  # Predecir en el pliegue de validación k
  pred_z <- predict(mod_k, newdata = val_k)
  res_k  <- na.omit(residuals(mod_k))
  
  # Retorno a CLP: Inversión directa y Smearing de Duan (estimado SOLO en train_k)
  pred_dir   <- 100000 * sinh(pred_z)
  pred_smear <- 100000 * (sinh(pred_z) * mean(cosh(res_k)) + cosh(pred_z) * mean(sinh(res_k)))
  
  idx <- !is.na(pred_dir) & !is.na(val_k$ganancia_final)
  y_real <- val_k$ganancia_final[idx]
  
  tabla_cv_pliegues <- rbind(tabla_cv_pliegues, data.frame(
    Pliegue = as.character(k),
    n_val   = sum(idx),
    MAE_Directo_CLP   = round(mean(abs(y_real - pred_dir[idx])), 0),
    RMSE_Directo_CLP  = round(sqrt(mean((y_real - pred_dir[idx])^2)), 0),
    MAE_Smearing_CLP  = round(mean(abs(y_real - pred_smear[idx])), 0),
    RMSE_Smearing_CLP = round(sqrt(mean((y_real - pred_smear[idx])^2)), 0)
  ))
}

# Agregar fila de promedios finales
fila_promedio <- data.frame(
  Pliegue = "Promedio CV",
  n_val   = round(mean(tabla_cv_pliegues$n_val), 0),
  MAE_Directo_CLP   = round(mean(tabla_cv_pliegues$MAE_Directo_CLP), 0),
  RMSE_Directo_CLP  = round(mean(tabla_cv_pliegues$RMSE_Directo_CLP), 0),
  MAE_Smearing_CLP  = round(mean(tabla_cv_pliegues$MAE_Smearing_CLP), 0),
  RMSE_Smearing_CLP = round(mean(tabla_cv_pliegues$RMSE_Smearing_CLP), 0)
)

tabla_cv_pliegues <- rbind(tabla_cv_pliegues, fila_promedio)
print(tabla_cv_pliegues, row.names = FALSE)





# ==============================================================================
# PREGUNTA 9: Regularización (Ridge, LASSO y Elastic Net)
# ==============================================================================


# --- 9.a) Codificación (matriz de diseño) y preparación sin NAs ---
formula_reg <- y_asinh ~ log1p(c6_3) + c7 + log1p(total_trabajadores) + 
  cine_eme_red + informalidad + c1_b + sexo + tramo_etario + region

vars_reg <- all.vars(formula_reg)
df_p9 <- df_entrenamiento %>%
  filter(if_all(all_of(vars_reg), ~ !is.na(.)))

# Matriz X (sin columna de intercepto [,-1], glmnet incluye intercepto y estandariza)
X_train <- model.matrix(formula_reg, data = df_p9)[, -1]
y_train <- df_p9$y_asinh
y_clp   <- df_p9$ganancia_final
folds_p9 <- df_p9$fold

# --- 9.b) Ajuste y selección de hiperparámetros usando los 5 pliegues predefinidos ---
set.seed(2026) # Reproducibilidad

# 1. Ridge (alpha = 0)
cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, foldid = folds_p9, standardize = TRUE)

# 2. LASSO (alpha = 1)
cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1, foldid = folds_p9, standardize = TRUE)

# 3. Elastic Net (grilla de alpha para elegir el mejor alpha y lambda por CV)
alphas_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9)
mejor_cvm <- Inf
cv_enet <- NULL
alpha_opt <- NA

for (a in alphas_grid) {
  fit_a <- cv.glmnet(X_train, y_train, alpha = a, foldid = folds_p9, standardize = TRUE)
  if (min(fit_a$cvm) < mejor_cvm) {
    mejor_cvm <- min(fit_a$cvm)
    cv_enet   <- fit_a
    alpha_opt <- a
  }
}

# Resumen de hiperparámetros escogidos
hiperparametros_p9 <- data.frame(
  Modelo = c("Ridge (alpha = 0)", "LASSO (alpha = 1)", paste0("Elastic Net (alpha = ", alpha_opt, ")")),
  Alpha  = c(0, 1, alpha_opt),
  Lambda_min = round(c(cv_ridge$lambda.min, cv_lasso$lambda.min, cv_enet$lambda.min), 6),
  MSE_CV_asinh = round(c(min(cv_ridge$cvm), min(cv_lasso$cvm), min(cv_enet$cvm)), 5),
  Lambda_1se = round(c(cv_ridge$lambda.1se, cv_lasso$lambda.1se, cv_enet$lambda.1se), 6),
  Vars_No_Cero_min = c(
    sum(coef(cv_ridge, s = "lambda.min")[-1] != 0),
    sum(coef(cv_lasso, s = "lambda.min")[-1] != 0),
    sum(coef(cv_enet,  s = "lambda.min")[-1] != 0)
  )
)
print(hiperparametros_p9, row.names = FALSE)

# --- 9.c) Evaluación en CLP (con sinh) y comparación de estabilidad de coeficientes ---

# Función para calcular MAE y RMSE en CLP mediante CV de 5 pliegues sin fuga
evaluar_glmnet_clp <- function(X, y_z, y_real_clp, folds, alpha_val, lambda_val) {
  pliegues <- sort(unique(folds))
  mae_dir <- numeric(length(pliegues))
  rmse_dir <- numeric(length(pliegues))
  mae_smear <- numeric(length(pliegues))
  rmse_smear <- numeric(length(pliegues))
  
  for (i in seq_along(pliegues)) {
    k <- pliegues[i]
    idx_tr <- folds != k
    idx_va <- folds == k
    
    # Ajustar glmnet solo en los 4 pliegues de entrenamiento (standardize=TRUE escala solo con train)
    mod_k <- glmnet(X[idx_tr, ], y_z[idx_tr], alpha = alpha_val, lambda = lambda_val, standardize = TRUE)
    
    pred_tr <- as.numeric(predict(mod_k, newx = X[idx_tr, ]))
    res_tr  <- y_z[idx_tr] - pred_tr
    
    pred_va <- as.numeric(predict(mod_k, newx = X[idx_va, ]))
    
    # Retorno a CLP usando seno hiperbólico
    p_dir   <- 100000 * sinh(pred_va)
    p_smear <- 100000 * (sinh(pred_va) * mean(cosh(res_tr)) + cosh(pred_va) * mean(sinh(res_tr)))
    
    mae_dir[i]    <- mean(abs(y_real_clp[idx_va] - p_dir))
    rmse_dir[i]   <- sqrt(mean((y_real_clp[idx_va] - p_dir)^2))
    mae_smear[i]  <- mean(abs(y_real_clp[idx_va] - p_smear))
    rmse_smear[i] <- sqrt(mean((y_real_clp[idx_va] - p_smear)^2))
  }
  
  c(MAE_Directo = round(mean(mae_dir), 0),
    RMSE_Directo = round(mean(rmse_dir), 0),
    MAE_Smearing = round(mean(mae_smear), 0),
    RMSE_Smearing = round(mean(rmse_smear), 0))
}

res_ridge <- evaluar_glmnet_clp(X_train, y_train, y_clp, folds_p9, 0, cv_ridge$lambda.min)
res_lasso <- evaluar_glmnet_clp(X_train, y_train, y_clp, folds_p9, 1, cv_lasso$lambda.min)
res_enet  <- evaluar_glmnet_clp(X_train, y_train, y_clp, folds_p9, alpha_opt, cv_enet$lambda.min)

tabla_clp_p9 <- rbind(
  OLS_Sin_Penalizar = c(524108, 2714044, 563726, 2679716),
  Ridge             = res_ridge,
  LASSO             = res_lasso,
  Elastic_Net       = res_enet
)
print(tabla_clp_p9)

# Comparar estabilidad de los coeficientes principales entre OLS, Ridge, LASSO y Elastic Net
coefs_comp <- data.frame(
  Variable = rownames(coef(cv_ridge, s = "lambda.min")),
  OLS         = round(as.numeric(coef(mod_p5_log)), 4),
  Ridge       = round(as.numeric(coef(cv_ridge, s = "lambda.min")), 4),
  LASSO       = round(as.numeric(coef(cv_lasso, s = "lambda.min")), 4),
  Elastic_Net = round(as.numeric(coef(cv_enet,  s = "lambda.min")), 4)
)
print(coefs_comp[1:16, ], row.names = FALSE)




# ==============================================================================
# PREGUNTA 10: Evaluación final en Conjunto de Prueba y Decisión
# ==============================================================================

# 1. Asegurar que df_prueba tenga exactamente los mismos niveles de factores que entrenamiento
df_prueba <- df_prueba %>%
  mutate(
    across(c(cine_eme_red, informalidad, c1_b, sexo, tramo_etario, region), as.factor),
    y_asinh = asinh(ganancia_final / 100000)
  )

# Verificar NAs en df_prueba en las variables del modelo
vars_final <- all.vars(formula_ganadora)
df_prueba_eval <- df_prueba %>%
  filter(if_all(all_of(vars_final), ~ !is.na(.)))

cat("Observaciones evaluadas en prueba:", nrow(df_prueba_eval), "\n")
y_test_clp <- df_prueba_eval$ganancia_final

# Función auxiliar para evaluar modelos lm transformados en prueba (CLP)
eval_lm_test <- function(modelo, newdata, y_real) {
  pred_z <- predict(modelo, newdata = newdata)
  res_tr <- na.omit(residuals(modelo))
  
  p_dir   <- 100000 * sinh(pred_z)
  p_smear <- 100000 * (sinh(pred_z) * mean(cosh(res_tr)) + cosh(pred_z) * mean(sinh(res_tr)))
  
  c(MAE_Dir = round(mean(abs(y_real - p_dir)), 0),
    RMSE_Dir = round(sqrt(mean((y_real - p_dir)^2)), 0),
    MAE_Smear = round(mean(abs(y_real - p_smear)), 0),
    RMSE_Smear = round(sqrt(mean((y_real - p_smear)^2)), 0))
}

# 2. Evaluar los modelos candidatos en Prueba

# Candidato 1: OLS Inicial (escala original CLP - Pregunta 3)
pred_ols_orig <- predict(modelo_inicial, newdata = df_prueba_eval)
res_test_c1 <- c(
  MAE_Dir   = round(mean(abs(y_test_clp - pred_ols_orig)), 0),
  RMSE_Dir  = round(sqrt(mean((y_test_clp - pred_ols_orig)^2)), 0),
  MAE_Smear = NA,
  RMSE_Smear = NA
)

# Candidato 2: OLS con asinh(Y/100000) y X original (Pregunta 4)
res_test_c2 <- eval_lm_test(mod_asinh, df_prueba_eval, y_test_clp)

# Candidato 3: OLS con asinh(Y/100000) + log1p(X) (Modelo Ganador P5/P7)
res_test_c3 <- eval_lm_test(mod_p5_log, df_prueba_eval, y_test_clp)

# Candidato 4 y 5: Ridge y Elastic Net (Pregunta 9)
X_test <- model.matrix(formula_reg, data = df_prueba_eval)[, -1]

eval_glmnet_test <- function(cv_fit, X_tr, y_tr_z, X_te, y_real) {
  pred_tr <- as.numeric(predict(cv_fit, newx = X_tr, s = "lambda.min"))
  res_tr  <- y_tr_z - pred_tr
  pred_te <- as.numeric(predict(cv_fit, newx = X_te, s = "lambda.min"))
  
  p_dir   <- 100000 * sinh(pred_te)
  p_smear <- 100000 * (sinh(pred_te) * mean(cosh(res_tr)) + cosh(pred_te) * mean(sinh(res_tr)))
  
  c(MAE_Dir = round(mean(abs(y_real - p_dir)), 0),
    RMSE_Dir = round(sqrt(mean((y_real - p_dir)^2)), 0),
    MAE_Smear = round(mean(abs(y_real - p_smear)), 0),
    RMSE_Smear = round(sqrt(mean((y_real - p_smear)^2)), 0))
}

res_test_c4 <- eval_glmnet_test(cv_ridge, X_train, y_train, X_test, y_test_clp)
res_test_c5 <- eval_glmnet_test(cv_enet,  X_train, y_train, X_test, y_test_clp)

# --- 10.a) Tabla Final de Candidatos (CV vs Prueba en CLP) ---
tabla_final_p10 <- data.frame(
  Modelo = c(
    "1. OLS Inicial (Y original, X original)",
    "2. OLS Transformado Y (asinh Y, X original)",
    "3. OLS Transformado Y + X (asinh Y + log1p X)",
    "4. Ridge (asinh Y + log1p X, lambda.min)",
    "5. Elastic Net (asinh Y + log1p X, alpha=0.1)"
  ),
  CV_MAE_Dir    = c(704215, 530757, 524108, 524521, 524126),
  CV_RMSE_Smear = c(2691100, 2653844, 2679716, 2682389, 2679900),
  Test_MAE_Dir  = c(res_test_c1["MAE_Dir"],  res_test_c2["MAE_Dir"],  res_test_c3["MAE_Dir"],  res_test_c4["MAE_Dir"],  res_test_c5["MAE_Dir"]),
  Test_RMSE_Dir = c(res_test_c1["RMSE_Dir"], res_test_c2["RMSE_Dir"], res_test_c3["RMSE_Dir"], res_test_c4["RMSE_Dir"], res_test_c5["RMSE_Dir"]),
  Test_MAE_Smear  = c(res_test_c1["MAE_Smear"],  res_test_c2["MAE_Smear"],  res_test_c3["MAE_Smear"],  res_test_c4["MAE_Smear"],  res_test_c5["MAE_Smear"]),
  Test_RMSE_Smear = c(res_test_c1["RMSE_Smear"], res_test_c2["RMSE_Smear"], res_test_c3["RMSE_Smear"], res_test_c4["RMSE_Smear"], res_test_c5["RMSE_Smear"])
)

print(tabla_final_p10, row.names = FALSE)







# Cálculo exacto CV Modelo 1 (fila 1 de tabla_final_p10)
cv_c1 <- sapply(pliegues_unicos, function(k) {
  m <- lm(ganancia_final ~ c6_3 + c7 + total_trabajadores + cine_eme_red + 
            informalidad + c1_b + sexo + tramo_etario + region, 
          data = df_entrenamiento %>% filter(fold != k))
  v <- df_entrenamiento %>% filter(fold == k)
  p <- predict(m, newdata = v)
  c(mae = mean(abs(v$ganancia_final - p)), rmse = sqrt(mean((v$ganancia_final - p)^2)))
})
tabla_final_p10[1, c("CV_MAE_Dir", "CV_RMSE_Smear")] <- round(rowMeans(cv_c1), 0)
print(tabla_final_p10, row.names = FALSE)

# Archivo de dependencias obligatorio
writeLines(capture.output(sessionInfo()), "dependencias.txt")
