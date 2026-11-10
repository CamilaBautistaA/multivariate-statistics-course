
#####################################################
#                                                  ##
#    Modelo de Regresiòn Logìstica Multinomial     ##
# Profesor: Clodomiro Fernando Miranda Villagòmez  ##   
#         cfmiranda@lamolina.edu.pe                ##
#                                                  ##
#####################################################

library(tinytable)
library(foreign)
library(gtExtras)
library(MASS) 
library(brant)
library(tidyverse)
library(car)
library(ggstatsplot)
library(emmeans)
library(gtsummary)
library(patchwork)
library(sjPlot)
library(jtools)
library(forestmodel)
library(ggeffects)
library(effectsize)
library(DHARMa)
library(cardx)
library(survey)
library(pscl)
library(rsample)
library(janitor)
library(performance)
library(caret)
library(randomForest)
library(vip)
library(stargazer)

# Mejor que la regresión logística binaria? La Regresión multinomial. La 
# alternativa secreta.
# En la regresión logística binaria la pregunta es si una categoría es más
# probable, ahora con la RLM se pregunta si una categoría es más probable que 
# todas las demás.

# Ejemplo: Se quiere estudiar sobre los factores que determinan la especie de 
# pinguino (Adelie, Gentoo y Chinstrap).

theme_set(theme_bw())

# importar los datos y darle el nombre d
tt(head(d))
d = d |>
  fastDummies::dummy_cols(select_columns = 'species',
                          remove_first_dummy = FALSE)
tt(head(d))
is.data.frame(d)
str(d)

d <- d %>%
  mutate(across(c(species, island, sex), as.factor)) %>%
  glimpse()

# Hemos Convertido cada categoría (Adelie, Gentoo y Chinstrap) en su propia 
# columna, con la finalidad de verificar que la RL Multinomial es un conjunto 
# de RL binarias que trabajan juntas. Usaremos un predictor (island)

# Visualizar la dummyfication
# devtools::install_github("strengejacke/strengejacke", force = TRUE)

d |>
  select(starts_with('species')) |>
  mutate_all(factor) |>
  view_df(show.frq = TRUE, show.prc = TRUE)

gt_plt_summary(d, title = 'Factores Sobre la Especie')

# La RL Multinomial es un conjunto de RL binarias que trabajan juntas (ejemplo:
# species_Adelie = 1 para Adelie y 0 para otros(Gentoo o Chinstrap))

d %>% 
  plot_frq(species)

d %>% 
  group_by(species) %>%
  plot_frq(island) %>%
  plot_grid()

tbl_summary(d,
            by=species,
            statistic = all_continuous() ~ '{mean} ({sd})')

tbl_summary(d,
            by=species)

# Visualizar que tres RL binarias y una RL multinomial (con un predictor categòrico)
# producen resultados iguales

plot_species = ggbarstats(
  d |> mutate(species = factor(species, levels = c('Adelie', 'Chinstrap', 'Gentoo'))),
  x = species, y = island, label = 'both',
  results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_Adelie = ggbarstats(d, species_Adelie, island, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_Chinstrap = ggbarstats(d, species_Chinstrap, island, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_Gentoo = ggbarstats(d, species_Gentoo, island, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_species + plot_Adelie + plot_Chinstrap + plot_Gentoo

# Visualizar resultados de RL binaria y RL multinomial con un predictor 
# categórico

sA = glm(species_Adelie ~ island, data = d, family = binomial())
A_plot = plot_model(sA, type = 'eff', terms = 'island', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

sC = glm(species_Chinstrap ~ island, data = d, family = binomial())
C_plot = plot_model(sC, type = 'eff', terms = 'island', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

sG = glm(species_Gentoo ~ island, data = d, family = binomial())
G_plot = plot_model(sG, type = 'eff', terms = 'island', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Modelo Multinomial
library(nnet)
m = multinom(species ~ island, data = d)

species_plot = plot_model(m, type = 'eff', terms = 'island', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Gráficos idénticos
(A_plot + C_plot + G_plot) / ( species_plot)

# guardar el gráfico
ggsave('model_plots.png', device = png, dpi = 999, width = 9, height = 5)

# Modelo Multinomial con un predictor categòrico y uno cuantitativo

m1 = multinom(species ~ island + body_mass_g, data = d)

# Como la variable dependiente es categórica nominal usamos una regresión 
# logística multinomial

# Supuestos

# 1. Independencia de las observaciones

# 2. Ausencia de multicolinealidad severa (VIF)
vif(m1)

# 3. Muestra suficientemente grande. Si converge no hay problemas
m2 = multinom(species ~ island + body_mass_g, data = d, trace = TRUE)

# 6. No pueden haber predictores que separen completamente una categorìa del 
#    resto (no se podrìan estimar los coeficientes o las estimaciones tienden
#   al infinito)

# 7. Las variables categòricas deben tener frecuencias suficientes para cada
#    categorìa


ggeffect(m1)
ggeffect(m1) |> plot() %>% plot_grid

#Contraste de significatividad global.
# H0: Beta2=Beta3=...=Betak=0
modelo_ls0 <- multinom(species ~ 1, data = d)

ji=modelo_ls0$deviance-m1$deviance
gl.ji=m1$edf-modelo_ls0$edf
pvalor=1-pchisq(ji,gl.ji)
print(cbind(ji,gl.ji,pvalor))

rl_table = tbl_regression(m1)
rl_table

rl_table1 = tbl_regression(m1, exponentiate = T)
rl_table1

forest_model(m1, exponentiate = TRUE)

plot_summs(m1, confint = TRUE, exp = TRUE) +
  scale_x_log10()

# Veamos que se obtienen probabilidades y Odds Ratios exactos
# Trabajando con el modelo m

EMM = emmeans(m, ~ island|species) # con la RL multinomial
REMM = regrid(EMM, transform = 'logit')
contREMM = contrast(REMM, method = 'pairwise', type = 'response', df = Inf)

EMMsA = emmeans(sA, pairwise ~ island, type = 'response') # RL binaria especie Adelie = 1

EMM
EMMsA$emmeans

contREMM
EMMsA$contrasts

# Con el siguiente còdigo se tienen gl→Inf, los resultados no son iguales a lo
# anterior
emmeans(m, ~ island|species)

contrast(REMM, method = 'pairwise', type = 'response')

# Volviendo al modelo m1

# importancia de las variables. Criterio de Anova y de randomForest
car::Anova(m1) # body_mass_g es màs importante

rf = randomForest(species ~ island + body_mass_g, data = d)

vip(rf)

# Se recomienda hacer varias RL binarias porque la interpretación es más fácil. 
# Pero si la variable dependiente tiene muchas categorías, hacer muchas RL 
# binarias puede se abrumador, en este caso es mejor usar la RL multinomial, 
# ya que tiene tres ventajas: En primer lugar da una clara descripción general 
# de los datos, en particular si se tiene múltiples predictores. En segundo 
# lugar, permite comparar las categorías directamente unas contra otras. 
# En tercer lugar tiene flexibilidad para realizar gráficos.

# Visualización de los datos de diferentes maneras. Con RLog Múltiple

emmip(m1, ~ species | island, CIs = TRUE, dodge = 0.5, type = 'response')

emmip(m1, species ~ island, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top') + # + para {patchwork}

emmip(m1, island ~ species, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top')

emtrends(m1, ~ body_mass_g | species, 'body_mass_g', weights = 'prop', infer = TRUE)


EMM = emmeans(m1, ~ body_mass_g | species,
              at = list(body_mass_g = c(3600, 6400)),
              weights = 'prop')

EMM

REMM = regrid(EMM, transform = 'logit')
contrast(REMM, method = 'revpairwise', type = 'response', df = Inf)

# Gráfico
emmip(m1, species ~ body_mass_g, CIs = TRUE, dodge = 2,
      at = list(body_mass_g = c(3500, 6500)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, species ~ body_mass_g | island, CIs = TRUE, dodge = 2,
      at = list(body_mass_g = c(3500, 6500)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, species ~ island | body_mass_g, CIs = TRUE, dodge = 2,
      at = list(body_mass_g = c(3500, 6500)), type = 'response') +
  theme(legend.position = 'top')

# Predicciones
summary(d)
summary(m1)
nuevo1 = data.frame(island = 'Biscoe', body_mass_g = 4200)
b11 = -1.89628; b21 = 1.108721; b31 = -0.2854992; b41 = 0.000213944
b12 = -18.12970; b22 = -2.285602; b32 = -2.2770720; b42 = 0.00435460

e1 = exp(b11 + b21*0 + b31*0 + b41*4200)
e2 = exp(b12 + b22*0 + b32*0 + b42*4200)
(p1 = e1/(1+e1 +e2)) # probabilidad de especie Chinstrap
(p2 = e2/(1+e1 +e2)) # probabilidad de especie Gentoo
(p3 = 1/(1+e1 +e2))  # probabilidad de especie Adelie
p1 +p2 +p3

# de otra manera
p=predict(m1,newdata=nuevo1,type="probs")
p

# Tabla de confusion

clase.pred <- predict(m1, d)
head(clase.pred)
head(d$species)
is.factor(clase.pred)
is.factor(d$species)
d1 = as.factor(d$species)
cm <- caret::confusionMatrix(clase.pred, d1)
print(cm)

# Modelo multinomial con interacciones
tt(head(d))
m5 = multinom(species ~ island * sex, data = d)

emmip(m5, sex ~ island | species, CIs = TRUE, dodge = 0.5,
      type = 'response')

emmip(m5, island ~ sex | species, CIs = TRUE, dodge = 0.5,
      type = 'response')

EMM = emmeans(m5, ~ sex | species | island)
REMM = regrid(EMM, transform = 'logit')
contrast(REMM, interaction = 'revpairwise', type = 'response', df = Inf)

# Modelo multinomial multivariado (con varios predictores)
library(glmulti)
glmulti(
  species ~ bill_length_mm + island + bill_depth_mm + flipper_length_mm + body_mass_g +
    sex,
  data = d,
  level = 1,
  method = 'd',
  fitfunction = multinom)

multinom_glmulti = glmulti(
  species ~ bill_length_mm + island + bill_depth_mm + flipper_length_mm + body_mass_g +
    sex,
  data = d,
  level = 1,
  method = 'h', # fuerza bruta
  fitfunction = multinom)

plot(multinom_glmulti, type = 's')

print(multinom_glmulti)

library(flextable)
weightable(multinom_glmulti)[1:3,] |>
  regulartable() |>
  autofit()

# Historial de los datos (tomemos el  primer modelo recomendado) y grafiquemos
m6 = multinom(species ~ 1 + bill_length_mm + island + bill_depth_mm + 
                flipper_length_mm + body_mass_g + sex,
             data = d)

ggeffect(m6) |> plot() |> plot_grid()

check_model(m6)
model_performance(m6)
tbl_regression(m6)
tbl_regression(m6, exponentiate = T)


