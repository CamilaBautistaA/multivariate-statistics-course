
# Modelo poisson con varianza robusta

# Cargar librerías necesarias
library(sandwich)
library(lmtest)
library(dplyr)
library(ggplot2)

# Simular datos
set.seed(123)

n <- 1000

producto <- factor(sample(c("Lácteo", "Cárnico", "Vegetal"), n, replace = TRUE))
proveedor <- factor(sample(c("A", "B"), n, replace = TRUE))
temperatura <- round(rnorm(n, mean = 72, sd = 5), 1)  # temperatura de cocción

# Coeficientes simulados (lacteo base, proveedor A base)
beta0 <- -2.0
beta_prod <- c(Cárnico = 0.6, Vegetal = 0.9)  # efecto sobre RR
beta_provB <- 0.4
beta_temp <- 0.03

# Calcular eta y probabilidad
eta <- beta0 +
  ifelse(producto == "Cárnico", beta_prod["Cárnico"],
         ifelse(producto == "Vegetal", beta_prod["Vegetal"], 0)) +
  ifelse(proveedor == "B", beta_provB, 0) +
  beta_temp * temperatura

p <- exp(eta)  # log link para riesgo relativo
p[p > 0.95] <- 0.95  # evitar probabilidades mayores a 1

rechazado <- rbinom(n, size = 1, prob = p)

# Dataset final
datos <- data.frame(rechazado, producto, proveedor, temperatura)
head(datos, 10)

# ---- Modelo Poisson con varianza robusta ----
modelo <- glm(rechazado ~ producto + proveedor + temperatura,
              family = poisson(link = "log"), data = datos)

library(performance)
check_model(modelo)

# Visualisar resultados del modelo
library(sjPlot)
library(ggeffects)
fancy_plot = ggeffect(modelo) |> plot() |> plot_grid()
fancy_plot

# Errores estándar robustos
se_robust <- sqrt(diag(sandwich(modelo)))

# Estimaciones
est <- coef(modelo)
IC_log <- est + outer(c(-1, 1) * 1.96, se_robust)
RR <- exp(est)
IC <- exp(IC_log)

# Resultados
tabla_resultados <- data.frame(
  Variable = names(RR),
  RR = round(RR, 3),
  IC_inf = round(IC[1, ], 3),
  IC_sup = round(IC[2, ], 3)
)

print(tabla_resultados)

# Estimadores robustos
coeftest(modelo, vcov = sandwich)

# Riesgo relativo aproximado
exp(coef(modelo))

library(emmeans)
emmeans(modelo, pairwise ~ producto,
        type = 'response', infer = T)

emmeans(modelo, ~ producto, type = 'response') %>%
  pairs(reverse = T, infer = T)

library(gtsummary)
nice_table = tbl_regression(modelo,
                            exponentiate = TRUE,
                            add_pairwise_contrasts = TRUE,
                            weights= 'prop',
                            pairwise_reverse = FALSE) |>
  
  # pimp your table
  add_significance_stars(hide_p = FALSE, hide_se = TRUE, hide_ci = FALSE) |>
  modify_header(label = '**Predictor**') |>
  modify_caption('Table 1. Nice looking table!') |>
  modify_footnote(abbreviation = TRUE) |>
  bold_labels() |>
  italicize_levels() |>
  bold_p()

#tidy_fun = broom.helpers::tidy_parameters
nice_table


# El modelo Log Binomial tiene como objetivo estimar Riesgos Relativos, pero 
# a veces presenta problemas de convergencia

library(logbin)
modelo_logbin = logbin(rechazado ~ producto + proveedor + temperatura,
                       data = datos)

library(performance)
check_model(modelo_logbin)

fancy_plot1 = ggeffect(modelo_logbin) |> plot() |> plot_grid()
fancy_plot1

# Riesgo relativo
exp(coef(modelo_logbin))

emmeans(modelo_logbin, pairwise ~ producto,
        type = 'response', infer = T)

emmeans(modelo_logbin, ~ producto, type = 'response') %>%
  pairs(reverse = T, infer = T)

# comparación de estimadores RR
RR_logbin <- exp(coef(modelo_logbin))
RR_poisson <- exp(coef(modelo))

data.frame(
  Coeficiente = names(RR_logbin),
  RR_logbin = round(RR_logbin, 3),
  RR_poisson = round(RR_poisson, 3)
)

nice_table1 = tbl_regression(modelo_logbin,
                            exponentiate = TRUE,
                            add_pairwise_contrasts = TRUE,
                            weights= 'prop',
                            pairwise_reverse = FALSE) |>
  
  # pimp your table
  add_significance_stars(hide_p = FALSE, hide_se = TRUE, hide_ci = FALSE) |>
  modify_header(label = '**Predictor**') |>
  modify_caption('Table 1. Nice looking table!') |>
  modify_footnote(abbreviation = TRUE) |>
  bold_labels() |>
  italicize_levels() |>
  bold_p()

#tidy_fun = broom.helpers::tidy_parameters
nice_table

# RR aproximado con Poisson vs RR de Log Binomial

Poisson_model = tbl_regression(
  modelo,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  weights = 'prop',
  pairwise_reverse = FALSE)

Poisson_model

Logbinomial_model = tbl_regression(
  modelo_logbin,
  exponentiate = TRUE,
  add_pairwise_contrasts = TRUE,
  weights = 'prop',
  pairwise_reverse = FALSE)

logbinomial_model

final_table = tbl_merge(
  tbls = list(Poisson_model, Logbinomial_model),
  tab_spanner = c('Poisson', 'Logbinomial')
)

final_table


