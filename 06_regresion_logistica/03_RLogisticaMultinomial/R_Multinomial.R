
#####################################################
#                                                  ##
#    Modelo de Regresiòn Logìstica Multinomial     ##
# Profesor: Clodomiro Fernando Miranda Villagòmez  ##   
#         cfmiranda@lamolina.edu.pe                ##
#                                                  ##
#####################################################

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

# Mejor que la regresión logísticabinaria? La Regresión multinomial. La 
# alternativa secreta.
# En la regresión logística binaria la pregunta es si una categoría es más
# probable, ahora con la RLM se pregunta si una categoría es más probable que 
# todas las demás.

# Ejemplo: Se quiere estudiar sobre los factores que dan forma al estatus 
# socioeconómico de las personas en California (bajo, medio y alto).

theme_set(theme_bw())

# los datos
d = read.dta('https://stats.idre.ucla.edu/stat/data/hsbdemo.dta') |>
  fastDummies::dummy_cols(select_columns = 'ses',
                          remove_first_dummy = FALSE)
head(d)
is.data.frame(d)
str(d)

# Hemos Convertido cada categoría (baja, media y alta) en su propia columna, con
# la finalidad de verificar que la RL Multinomial es un conjunto de RL binarias 
# que trabajan juntas. Usaremos un predictor (prog = programa de estudio)

# Visualizar la dummyfication
# devtools::install_github("strengejacke/strengejacke", force = TRUE)

d |>
  select(starts_with('ses')) |>
  mutate_all(factor) |>
  view_df(show.frq = TRUE, show.prc = TRUE)

gt_plt_summary(d, title = 'Factores Sobre el Nivel Socio Econòmico')

# La RL Multinomial es un conjunto de RL binarias que trabajan juntas (ejemplo:
# ses_low = 1 para bajo y 0 para otros(medio o alto))

d %>% 
  plot_frq(ses)

d %>% 
  group_by(ses) %>%
  plot_frq(schtyp) %>%
  plot_grid()

tbl_summary(d,
            by=ses,
            statistic = all_continuous() ~ '{mean} ({sd})')

tbl_summary(d,
            by=ses)

# Visualizar que tres RL binarias y una RL multinomial (con un predictor categòrico)
# producen resultados iguales

plot_ses = ggbarstats(
  d |> mutate(ses = factor(ses, levels = c('high', 'middle', 'low'))),
  x = ses, y = prog, label = 'both',
  results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_low = ggbarstats(d, ses_low, prog, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_middle = ggbarstats(d, ses_middle, prog, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_high = ggbarstats(d, ses_high, prog, label = 'both',
                      results.subtitle = FALSE) +
  theme(legend.position = 'top')

plot_ses + plot_low + plot_middle + plot_high

# Visualizar resultados de RL binaria y RL multinomial con un predictor 
# categórico

sl = glm(ses_low ~ prog, data = d, family = binomial())
low_plot = plot_model(sl, type = 'eff', terms = 'prog', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

sm = glm(ses_middle ~ prog, data = d, family = binomial())
middle_plot = plot_model(sm, type = 'eff', terms = 'prog', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

sh = glm(ses_high ~ prog, data = d, family = binomial())
high_plot = plot_model(sh, type = 'eff', terms = 'prog', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Modelo Multinomial
library(nnet)
m = multinom(ses ~ prog, data = d)

ses_plot = plot_model(m, type = 'eff', terms = 'prog', title = '') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))

# Gráficos idénticos
(low_plot + middle_plot + high_plot) / ( ses_plot)

# guardar el gráfico
ggsave('model_plots.png', device = png, dpi = 999, width = 9, height = 5)

# Modelo Multinomial con un predictor categòrico y uno cuantitativo

m1 = multinom(ses ~ prog + math, data = d)

# Como la variable dependiente es categórica ordinal, veamos si también se puede
# usar una regresión logística ordinal

# Modelo de regresion logistica ordinal
mord <- polr(ses ~ prog + math, data = d)

# H0: Las regresiones son paralelas (los Odds son proporcionales).
# Como el pvalor de "Omnibus" es 0.18 es mejor hacer una regresión ordinal
# pero haremos la Multinomial
brant(mord)
gofcat::brant.test(mord)
rms::plot.xmean.ordinaly(ses ~ prog + math, data=d)

# Supuestos

# 1. Independencia de las observaciones

# 2. Ausencia de multicolinealidad severa (VIF)
vif(m1)

# 3. Muestra suficientemente grande. Si converge no hay problemas
m2 = multinom(ses ~ prog + math, data = d, trace = TRUE)

# 4. Relaciòn lineal entre predictores continuos y el logit (Z). math_log resulta 
#    significativo → la relaciòn no es lineal (usar transformaciones o splines)
d$math_log = d$math * log(d$math)
m3 = multinom(ses ~ prog + math + math_log, data = d)
s = summary(m3)
z = s$coefficients / s$standard.errors
(pvalor = 2*(1-pnorm(abs(z))))

# 5. La razòn de probabilidades entre dos categorìas no cambia si se agrega o 
#    elimina otra categorìa (Prueba de: Hausman-McFadden o de Small-Hsiao)
library(mlogit)
d_mlogit = mlogit.data(d, choice = 'ses', shape = 'wide')
m4 = mlogit(ses ~ prog + math, data = d_mlogit)

# 6. No pueden haber predictores que separen completamente una categorìa del 
#    resto (no se podrìan estimar los coeficientes o las estimaciones tienden
#   al infinito)

# 7. Las variables categòricas deben tener frecuencias suficientes para cada
#    categorìa


ggeffect(m1)
ggeffect(m1) |> plot() %>% plot_grid

#Contraste de significatividad global.
# H0: Beta2=Beta3=...=Betak=0
modelo_ls0 <- multinom(ses ~ 1, data = d)

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

EMM = emmeans(m, ~ prog|ses) # con la RL multinomial
REMM = regrid(EMM, transform = 'logit')
contREMM = contrast(REMM, method = 'pairwise', type = 'response', df = Inf)

EMMsh = emmeans(sh, pairwise ~ prog, type = 'response') # RL binaria nivel alto = 1

EMM
EMMsh$emmeans

contREMM
EMMsh$contrasts

# No cometer las siguientes fallas, porque los resultados no son iguales
emmeans(m, ~ prog|ses)

contrast(REMM, method = 'pairwise', type = 'response')

# Se recomienda hacer varias RL binarias porque la interpretación es más fácil. 
# Pero si la variable dependiente tiene muchas categorías, hacer muchas RL 
# binarias puede se abrumador, en este caso es mejor usar la RL multinomial, 
# ya que tiene tres ventajas: En primer lugar da una clara descripción general 
# de los datos, en particular si se tiene múltiples predictores. En segundo 
# lugar, permite comparar las categorías directamente unas contra otras. 
# En tercer lugar tiene flexibilidad para realizar gráficos.

# Visualización de los datos de diferentes maneras. Con RLog Múltiple

emmip(m1, ~ ses | prog, CIs = TRUE, dodge = 0.5, type = 'response')

emmip(m1, ses ~ prog, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top') + # + para {patchwork}

emmip(m1, prog ~ ses, CIs = TRUE, dodge = 0.5,
      type = 'response') +
  theme(legend.position = 'top')

emtrends(m1, ~ math | ses, 'math', weights = 'prop', infer = TRUE)

EMM = emmeans(m1, ~ math | ses,
              at = list(math = c(40, 70)),
              weights = 'prop')

EMM

REMM = regrid(EMM, transform = 'logit')
contrast(REMM, method = 'revpairwise', type = 'response', df = Inf)

# Gráfico
emmip(m1, ses ~ math, CIs = TRUE, dodge = 2,
      at = list(math = c(40, 70)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, ses ~ math | prog, CIs = TRUE, dodge = 2,
      at = list(math = c(40, 70)), type = 'response') +
  theme(legend.position = 'top')

emmip(m1, ses ~ prog | math, CIs = TRUE, dodge = 2,
      at = list(math = c(40, 70)), type = 'response') +
  theme(legend.position = 'top')

# Predicciones
summary(d)
summary(m1)
nuevo1 = data.frame(prog = 'academic', math = 68)
b11 = -2.227982; b21 = 0.3432541; b31 = 0.9223466; b41 = 0.04961603
b12 = -4.202097; b22 = 0.9526703; b32 = 0.3181833; b42 = 0.07252600

e1 = exp(b11 + b21*1 + b31*0 + b41*68)
e2 = exp(b12 + b22*1 + b32*0 + b42*68)
(p1 = e1/(1+e1 +e2)) # probabilidad de ses middle
(p2 = e2/(1+e1 +e2)) # probabilidad de ses high
(p3 = 1/(1+e1 +e2))  # probabilidad de ses low
p1 +p2 +p3

# de otra manera
p=predict(m1,newdata=nuevo1,type="probs")
p

# Modelo multinomial con interacciones
head(d)
m5 = multinom(ses ~ prog * female, data = d)

emmip(m5, female ~ prog | ses, CIs = TRUE, dodge = 0.5,
      type = 'response')

emmip(m5, prog ~ female | ses, CIs = TRUE, dodge = 0.5,
      type = 'response')

EMM = emmeans(m5, ~ female | ses | prog)
REMM = regrid(EMM, transform = 'logit')
contrast(REMM, interaction = 'revpairwise', type = 'response', df = Inf)

# Modelo multinomial multivariado (con varios predictores)
library(glmulti)
glmulti(
  ses ~ schtyp + prog + write + awards + math +
    honors + read + socst + science + cid,
  data = d,
  level = 1,
  method = 'd',
  fitfunction = multinom)

multinom_glmulti = glmulti(
  ses ~ schtyp + prog + write + awards + math +
    honors + read + socst + science + cid,
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

# Historial de los datos (tomemos el  segundo modelo recomendado) y grafiquemos
m6 = multinom(ses ~ schtyp + prog + write + socst + honors + science,
             data = d)

ggeffect(m6) |> plot() |> plot_grid()

check_model(m6)
model_performance(m6)
tbl_regression(m6)
tbl_regression(m6, exponentiate = T)

# importancia de las variables.
rf = randomForest(ses ~ schtyp + prog + write + socst + honors + science,
                  data = d)

vip(rf)
