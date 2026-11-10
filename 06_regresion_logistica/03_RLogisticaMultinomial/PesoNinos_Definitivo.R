
#####################################################
#                                                  ##
#    Modelo de Regresiòn Logìstica Multinomial     ##
# Profesor: Clodomiro Fernando Miranda Villagòmez  ##   
#         cfmiranda@lamolina.edu.pe                ##
#                                                  ##
#####################################################

library(foreign)
library(MASS) 
library(brant)
library(tidyverse)
library(forcats)
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
library(tictoc)
library(caret)
library(randomForest)
library(vip)
library(stargazer)
theme_set(theme_bw())

# Mejor que la regresión logísticabinaria? La Regresión multinomial. La 
# alternativa secreta.
# En la regresión logística binaria la pregunta es si una categoría es más
# probable, ahora con la RLM se pregunta si una categoría es más probable que 
# todas las demás.

d <-read.spss("PesoNinos1.sav",
              use.value.labels=TRUE, 
              to.data.frame=TRUE)
head(d)
d = d[, -1]
head(d)
summary(d)

d = d %>% rename(peso = Peso_RN_1,
                 edad = Edad_Madre,
                 partos = Num_Partos,
                 tipo = Tipo_de_parto,
                 sexo = Sexo_RN,
                 talla = Talla_RN,
                 ed_gest = Edad_Gestacional)
summary(d)

d$peso = fct_recode(d$peso,
                    bajo = 'Bajo Peso',
                    normal = 'Peso Normal',
                    sobre = 'Sobre Peso')

d$tipo = fct_recode(d$tipo,
                    cesarea = 'Cesarea',
                    t_normal = 'Parto  Normal')

d$sexo = fct_recode(d$sexo,
                    mujer = 'Femenino',
                    hombre = 'Masculino')

summary(d)

d = d |>
  fastDummies::dummy_cols(select_columns = 'peso',
                          remove_first_dummy = FALSE)
head(d)

# Hemos Convertido cada categoría (bajo, normal y sobre) en su propia columna, con
# la finalidad de verificar que la RL Multinomial es un conjunto de RL binarias 
# que trabajan juntas. Usaremos un predictor (sexo = sexo del recien nacido)

d |>
  select(starts_with('peso')) |>
  mutate_all(factor) |>
  view_df(show.frq = TRUE, show.prc = TRUE)

gt_plt_summary(d, title = 'Factores Sobre el Peso de Niños')

# La RL Multinomial es un conjunto de RL binarias que trabajan juntas (ejemplo:
# peso_bajo = 1 para bajo y 0 para otros(normal o sobre))

d %>% 
  plot_frq(peso)

d %>% 
  group_by(peso) %>%
  plot_frq(sexo) %>%
  plot_grid()

tbl_summary(d,
            by=peso,
            statistic = all_continuous() ~ '{mean} ({sd})')

tbl_summary(d,
            by=peso)

# Visualizar que tres RL binarias y una RL multinomial producen resultados iguales

plot_peso = ggbarstats(
  d |> mutate(peso = factor(peso, levels = c('bajo', 'normal', 'sobre'))),
  x = peso, y = sexo, label = 'both',
  results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_bajo = ggbarstats(d, peso_bajo, sexo, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_normal = ggbarstats(d, peso_normal, sexo, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_sobre = ggbarstats(d, peso_sobre, sexo, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_peso + plot_bajo + plot_normal + plot_sobre

# Visualizar resultados de RL binaria y RL multinomial con un predictor 
# categórico

sb = glm(peso_bajo ~ sexo, data = d, family = binomial())
bajo_plot = plot_model(sb, type = 'eff', terms = 'sexo', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

sn = glm(peso_normal ~ sexo, data = d, family = binomial())
normal_plot = plot_model(sn, type = 'eff', terms = 'sexo', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

ss = glm(peso_sobre ~ sexo, data = d, family = binomial())
sobre_plot = plot_model(ss, type = 'eff', terms = 'sexo', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Modelo Multinomial
library(nnet)
m = multinom(peso ~ sexo, data = d)

peso_plot = plot_model(m, type = 'eff', terms = 'sexo', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Gráficos idénticos
(bajo_plot + normal_plot + sobre_plot) / ( peso_plot)

# guardar el gráfico
ggsave('model_plots.png', device = png, dpi = 999, width = 9, height = 5)

# Modelo Multinomial con un predictor categòrico y uno cuantitativo

m1 = multinom(peso ~ sexo + edad, data = d)

# Como la variable dependiente es categórica ordinal, veamos si también se puede
# usar una regresión logística ordinal

# Modelo de regresion logistica ordinal
mord <- polr(peso ~ sexo, data = d)

# H0: Las regresiones son paralelas (los Odds son proporcionales).
# Como el pvalor de "Omnibus" es 0.11 es mejor hacer una regresión ordinal
# pero haremos la Multinomial
brant(mord)
gofcat::brant.test(mord)

# Mètodo gràfico. Se grafican las medias de los predictores por nivel de 
# categorìa en una lìnea continua y los valores esperados asumiendo que las
# suposiciones de odds proporcionales son verdaderas en una lìnea discontinua.
# Si las lìneas continua y discontinua siguen una estrecha correlaciòn, el 
# supuesto de odds proporcionales de nuestro modelo es sòlido, pero si empiezan 
# a distanciarse el mètodo gràfico es mejor que la prueba de Brandt porque 
# muestra en què medida el supuesto es erròneo. Si el supuesto se viola 
# gravemente se puede usar como alternativa el modelo logìstico multinomial.
# En el gràfico se observa que la lìnea continua y discontinua se distancian. Por
# lo tanto, es mejor usar el modelo multinomial.
rms::plot.xmean.ordinaly(peso ~ sexo * edad, data = d)

# Supuestos

# 1. Independencia de las observaciones

# 2. Ausencia de multicolinealidad severa (VIF)
vif(m1)

# 3. Muestra suficientemente grande. Si converge no hay problemas
m2 = multinom(peso ~ sexo + edad, data = d, trace = TRUE)

# 4. Relaciòn lineal entre predictores continuos y el logit (Z). math_log resulta 
#    no significativo → la relaciòn es lineal
d$edad_log = d$edad * log(d$edad)
m3 = multinom(peso ~ sexo + edad + edad_log, data = d)
s = summary(m3)
z = s$coefficients / s$standard.errors
(pvalor = 2*(1-pnorm(abs(z))))

# 5. La razòn de probabilidades entre dos categorìas no cambia si se agrega o 
#    elimina otra categorìa (Prueba de: Hausman-McFadden o de Small-Hsiao)
library(mlogit)
d_mlogit = mlogit.data(d, choice = 'peso', shape = 'wide')
m4 = mlogit(peso ~ sexo + edad, data = d_mlogit)

# 6. No pueden haber predictores que separen completamente una categorìa del 
#    resto (no se podrìan estimar los coeficientes o las estimaciones tienden
#   al infinito)

# 7. Las variables categòricas deben tener frecuencias suficientes para cada
#    categorìa


ggeffect(m1)
ggeffect(m1) |> plot() %>% plot_grid

#Contraste de significatividad global.
# H0: Beta2=Beta3=...=Betak=0
modelo_ls0 <- multinom(peso ~ 1, data = d)

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

EMM = emmeans(m, ~ sexo|peso) # con la RL multinomial
REMM = regrid(EMM, transform = 'logit')
contREMM = contrast(REMM, method = 'pairwise', type = 'response', df = Inf)

EMMps = emmeans(ss, pairwise ~ sexo, type = 'response') # RL binaria nivel alto = 1

EMM
EMMps$emmeans

contREMM
EMMps$contrasts

# No cometer las siguientes fallas, porque los resultados no son iguales
emmeans(m, ~ sexo|peso)

contrast(REMM, method = 'pairwise', type = 'response')

# Se recomienda hacer varias RL binarias porque la interpretación es más fácil. 
# Pero si la variable dependiente tiene muchas categorías, hacer muchas RL 
# binarias puede se abrumador, en este caso es mejor usar la RL multinomial, 
# ya que tiene tres ventajas: En primer lugar da una clara descripción general 
# de los datos, en particular si se tiene múltiples predictores. En segundo 
# lugar, permite comparar las categorías directamente unas contra otras. 
# En tercer lugar tiene flexibilidad para realizar gráficos.

# Visualización de los datos de diferentes maneras. Con RLog Múltiple

emmip(m, ~ peso | sexo, CIs = TRUE, dodge = 0.5, type = 'response')

emmip(m, peso ~ sexo, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top') + # + para {patchwork}

emmip(m, sexo ~ peso, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top')

emtrends(m1, ~ edad | peso, 'edad', weights = 'prop', infer = TRUE)

EMM = emmeans(m1, ~ edad | peso,
              at = list(edad = c(23, 35)),
              weights = 'prop')

EMM

REMM = regrid(EMM, transform = 'logit')
contrast(REMM, method = 'revpairwise', type = 'response', df = Inf)

# Gráfico
emmip(m1, peso ~ edad, CIs = TRUE, dodge = 2,
      at = list(edad = c(23, 35)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, peso ~ sexo | edad, CIs = TRUE, dodge = 2,
      at = list(edad = c(23, 35)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, peso ~ edad | sexo, CIs = TRUE, dodge = 2,
      at = list(edad = c(23, 35)), type = 'response') +
  theme(legend.position = 'top')

# Predicciones
summary(d)
summary(m1)
nuevo1 = data.frame(sexo = 'hombre', edad = 25)
b11 = -0.2345847; b21 = 2.061734; b31 = 0.06638464
b12 = -9.9009485; b22 = 2.702140; b32 = 0.2565848

e1 = exp(b11 + b21*1 + b31*25)
e2 = exp(b12 + b22*1 + b32*25)
(p1 = e1/(1+e1 +e2)) # probabilidad de ses middle
(p2 = e2/(1+e1 +e2)) # probabilidad de ses high
(p3 = 1/(1+e1 +e2))  # probabilidad de ses low
p1 +p2 +p3

# de otra manera
p=predict(m1,newdata=nuevo1,type="probs")
p

# Modelo multinomial con interacciones
head(d)
m1 = multinom(peso ~ sexo * tipo, data = d)

emmip(m1, tipo ~ sexo| peso, CIs = TRUE, dodge = 0.5,
      type = 'response')

emmip(m1, sexo ~ tipo | peso, CIs = TRUE, dodge = 0.5,
      type = 'response')

EMM = emmeans(m1, ~ tipo | peso | sexo)
REMM = regrid(EMM, transform = 'logit')
contrast(REMM, interaction = 'revpairwise', type = 'response', df = Inf)

# Modelo multinomial multivariado (con varias predictoras)
library(glmulti)

glmulti(
  peso ~ edad + partos + tipo + sexo + ed_gest, data = d,
  level = 1,    # considera sòlo efectos principales
  method = 'd', # resumen simple de los candidatos
  fitfunction = multinom)

#glmulti(
#  peso ~ edad + partos + tipo + sexo + ed_gest, data = d,
#  level = 1,  # considera sòlo efectos principales
#  method = 'g', # algoritmo genètico
#  fitfunction = multinom)

tic()
multinom_glmulti = glmulti(
  peso ~ edad + partos + tipo + sexo + ed_gest, data = d,
  level = 2,  # considera hasta interacciones dobles
  method = 'h', # fuerza bruta
  fitfunction = multinom)
toc()

plot(multinom_glmulti, type = 's')

print(multinom_glmulti)

# viendo el peso de los 3 mejores modelos
library(flextable)
weightable(multinom_glmulti)[1:3,] |>
  regulartable() |>
  autofit()

# Historial de los datos (tomemos el  segundo modelo recomendado) y grafiquemos
m = multinom(peso ~ 1 + edad + partos + ed_gest + partos:edad + ed_gest:partos + tipo:partos + sexo:partos,
             data = d)

ggeffect(m) |> plot() |> plot_grid()

check_model(m)
model_performance(m)
tbl_regression(m)
tbl_regression(m, exponentiate = T)

# importancia de las variables.
rf = randomForest(peso ~ 1 + edad + partos + ed_gest + partos:edad + ed_gest:partos + tipo:partos + sexo:partos,
                  data = d)

vip(rf)
