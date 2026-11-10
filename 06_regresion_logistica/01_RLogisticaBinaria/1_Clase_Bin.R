######################################################
#                                                    #
#  ANÁLISIS DE REGRESIÓN LOGÍSTICA BINARIA           #
#  Data: Empresas                                    #   
#                                                    #
#  Profesor: Clodomiro Fernando Miranda Villagómez   #                                      #
######################################################


#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o información cargada
rm(list = ls())

#---------------------------------------------------------
#  Librerias
library(tinytable)
library(MASS) 
library(ROCR)
library(pROC)
library(car)
library(ResourceSelection)
library(vcdExtra)
library(pscl)
library(tidyverse)
library(gtsummary)
library(ggstatsplot)
library(patchwork)
library(sjPlot)
library(jtools)
library(forestmodel)
library(emmeans)
library(performance)
library(DHARMa)
library(cardx)
library(ggeffects)
library(equatiomatic)
library(effectsize)
library(multcomp)

# El ARLBinaria se usa para modelar respuestas binarias, tales como el diagnóstico
# de una enfermedad (positivo o negativo), votar (a favor o en contra) o 
# la admisión a una universidad (admitido o no).

# Importar datos
tt(head(datos, 16))
dim(datos)
str(datos)

attach(datos)
is.data.frame(datos)

datos <- datos %>%
  mutate(across(c(resultado, gradnove, tipo, investigacion,
                  sectecng), as.factor)) %>%
  glimpse()

head(datos,8)

# Para ver la categoria de referencia
# con la función levels: la primera es la categoría de referencia
levels(c(datos$resultado, datos$gradnove, datos$tipo, 
         datos$investigacion, datos$sectecng))

# en forma detallada:
contrasts(datos$resultado)  
contrasts(datos$gradnove)

contrasts(datos$tipo)
contrasts(datos$investigacion)
contrasts(datos$sectecng)
?contrasts

# Para cambiar la categoria de referencia.
datos$resultado=relevel(datos$resultado,ref = "fracaso")
datos$gradnove = relevel(datos$gradnove,ref="prod_nuevos")
datos$sectecng = relevel(datos$sectecng,ref="alta")
datos$investigacion = relevel(datos$investigacion,ref="si")

# Para ver la categoria de referencia
contrasts(datos$resultado)  
contrasts(datos$gradnove)

contrasts(datos$tipo)
contrasts(datos$investigacion)
contrasts(datos$sectecng)

# Estadìstica descriptiva

view_df(datos, show.frq=T, show.prc = T, show.na = T)

datos %>% 
  plot_frq(sectecng)

datos %>% 
  group_by(resultado) %>%
  plot_frq(sectecng) %>%
  plot_grid()

datos %>% 
  group_by(resultado) %>% 
  plot_frq(
    publicid,
    type = 'histogram',
    show.mean = T,
    normal.curve = T
  ) %>% 
  plot_grid()

datos %>% 
  group_by(resultado) %>% 
  plot_frq(
    publicid,
    type = 'violin',
    show.mean = T,
    normal.curve = T
  ) %>% 
  plot_grid()

datos %>% 
  group_by(resultado) %>% 
  plot_frq(
    personal,
    type = 'bar',
    show.mean = T,
    normal.curve = T
  ) %>% 
  plot_grid()

tbl_summary(datos,
            by=resultado,
            statistic = all_continuous() ~ '{mean} ({sd})')

tbl_summary(datos,
            by=resultado)

# Tabla combinada con comparación por resultado

tabla_combinada = datos %>%
  select(resultado , gradnove, tipo, investigacion, sectecng) %>%
  tbl_summary(
    by = resultado, # variable de comparación
    include = c(gradnove, tipo, investigacion, sectecng),
    percent = 'row', # porcentaje por fila
    statistic = list(all_categorical() ~ '{n} ({p})'),
    digits = all_categorical() ~ c(0, 2),
  ) %>%
  add_overall(last = TRUE) %>% # agrega la columna 'total'
  add_p(
    test = everything() ~ 'chisq.test', # aplica chi-cuadrado a todas las variables
    pvalue_fun = ~ style_pvalue(.x, digits = 3) # formato del p-valor 
  ) %>%
  modify_header(label ~ '**Variable**') %>% # cambia el encabezado de la columna de variables
  modify_caption('Tabla 1: Comparación de Gradnove,\ntipo, investigación y sectecng  por resultado')
tabla_combinada


# CREAR BASE DE DATOS PARA ENTRENAMIENTO Y PARA PRUEBA.
# Con la data de entrenamiento se hace el modelo y
# con la data de prueba se mide la performance del modelo.
# Para no ser juez y parte.

library(caret)
tt(head(datos))

set.seed(8)
ind_train = createDataPartition(y = datos$resultado, p=0.7, list=FALSE)
base_entrenamiento = datos[ind_train,]
dim(base_entrenamiento)
base_prueba = datos[-ind_train,]
dim(base_prueba)
dim(datos)[1]
dim(base_entrenamiento)[1]/dim(datos)[1]
dim(base_prueba)[1]/dim(datos)[1]

prop.table(table(datos$resultado))
prop.table(table(base_entrenamiento$resultado))
prop.table(table(base_prueba$resultado))

library(MASS)
#  La seleccion de Variables se hace con la data de entrenamiento

tt(head(base_entrenamiento))
modelo <- glm(resultado ~ ., family = binomial(link = logit), 
            data = base_entrenamiento)

step <- stepAIC(modelo, direction = "forward", trace=FALSE)
step$anova # nos quedamos con todas las variables predictoras

# Estadistico òmnibus y Estadistico de Wald
#------------------------------------------

## Prueba Omnibus

# Primera manera
#Contraste de significatividad global.
# H0: Beta2=Beta3=...=Betak=0
# En Model Likelihood Ratio Test con un pvalor<0.0001 se rechaza la H0
# Se muestra la conclusiòn con todas las variables independientes.

library(rms)
mod_lrm = lrm(formula = resultado ~ ., data = base_entrenamiento, 
              x = TRUE, y = TRUE)
mod_lrm

# Segunda manera
#Contraste de significatividad global. Manualmente
# H0: Beta2=Beta3=...=Betak=0
deviance_model = modelo$deviance
modelo$deviance
deviance_base = modelo$null.deviance
modelo$null.deviance
ji = deviance_base - deviance_model
ji
ji_gl = modelo$df.null - modelo$df.residual
ji_gl
pvalor = pchisq(ji, df = ji_gl, lower.tail = FALSE)
pvalor

# Tercera manera
library(lmtest)
m0 = glm(resultado ~ 1, family = binomial(link = logit), data = base_entrenamiento)
lrtest0 = lrtest(m0, modelo)
lrtest0 # El modelo se ajusta bien a los datos (al menos un beta != cero)

#----------------------------------------------------------
# Devianza de cada variable

# publicid, investigación y personal tienen devianza significativa, eso indica 
# que cada una de ellas reduce significativamente la devianza del modelo al 
# incluirse (cada una de ellas mejora el ajuste del modelo)
anova(modelo, test="Chisq")

#------------------------------------------------------------------#
#                      Eleccion del modelo final                   #
#------------------------------------------------------------------#

modelo1<-glm(resultado ~ ., family = binomial(link = logit),data=base_entrenamiento)

# para exportar la base_entrenamiento a excel y usarla en Markdown
# library(writexl)
# write_xlsx(base_entrenamiento, "base_entrenamiento.xlsx")
#getwd()
# ver el R Markdown

summary(modelo1)

#CONTRASTE PARA LOS COEFICIENTES INDIVIDUALES
# Estadistico W de Wald

library(survey)
regTermTest(modelo1,"publicid",method="Wald")
regTermTest(modelo1,"gradnove",method="Wald")
#regTermTest(modelo1,"tipo",method="Wald")
regTermTest(modelo1,"investigacion",method="Wald")
regTermTest(modelo1,"sectecng",method="Wald")
regTermTest(modelo1,"personal",method="Wald")

#CALIBRACION DEL MODELO
#PRUEBA DE HOSMER LEMESHOW

hl <- hoslem.test(modelo1$y, fitted(modelo1),g=10)
hl
cbind(hl$observed,hl$expected)

# de otra forma
HL=HLtest(modelo1,g=10)
HL
HL$table

#EVALUANDO EL AJUSTE DEL MODELO
# Pseudos R cuadrado. Muchos autores coinciden que un pseudo R^2 es óptimo si 
# está entre 0.4 y 0.6
library(pscl)
pR2(modelo1)

library(rcompanion)
nagelkerke(modelo1)

library(blorr)
blr_rsq_count(modelo1)
blr_rsq_adj_count(modelo1)
blr_rsq_nagelkerke(modelo1)
blr_rsq_mckelvey_zavoina(modelo1)
blr_rsq_mcfadden(modelo1)
blr_rsq_mcfadden_adj(modelo1)
blr_rsq_cox_snell(modelo1)
blr_rsq_effron(modelo1)

# Calibración del modelo
dim(base_entrenamiento)
calibracion = rms::calibrate(mod_lrm, B = 110)

# se puede observar una buena calibracion del modelo, es decir el modelo predice
# correctamente probabilidades cercanas a lo real es decir está bien calibrado,
# ello se ve cuando la linea negra(bias-corrected) se acerca bastante a la 
# diagonal(ideal) las lineas sombreadas muestran intervalos de confianza, si el
# intervalo es estrecho el modelo es estable, es decir no cambian mucho los 
# resultados si cambiamos el boostrap.

# El MAE (Mean Absolute Error) no evalúa la discriminación (qué tan bien separa
# clases) como lo hacen las métricas AUC, sensibilidad, MCC, etc. En el ejemplo 
# El MAE es una medida de calibración (que tan bien coinciden probabilidades 
# predichas con la realidad). El MAE = 0.028 indica la distancia, en promedio, 
# de la probabilidad predicha del valor real (0 o 1) el modelo se equivoca en 
# 0.028 unidades de probabilidad
plot(calibracion)

# chequeo visual de los supuestos del modelo de RL binomial
#windows()
check_model(modelo1)
#dev.off()

# check_normality(modelo1) # no aplica a modelos glm
# plot(check_normality(modelo1))

check_autocorrelation(modelo1)

check_dag(modelo1)
plot(check_dag(modelo1))

check_residuals(modelo1)
plot(check_residuals(modelo1))

check_distribution(modelo1)
plot(check_distribution(modelo1))

check_predictions(modelo1)

check_collinearity(modelo1)
plot(check_collinearity(modelo1))

# check_heteroscedasticity(modelo1) # no aplica a modelos glm
# plot(check_heteroscedasticity(modelo1))

check_outliers(modelo1)
plot(check_outliers(modelo1))

# La overdispersion ocurre cuando la variabilidad observada en los datos es 
# mayor de lo que el modelo predice.
check_overdispersion(modelo1)

#----------------------------------------------------------
# Importancia de las variables
library(caret)
import=varImp(modelo1)

orden = dplyr::arrange(import, desc(import))
orden

# importancia de las variables. Criterio de Anova
# publicid es màs importante e investigación en segundo lugar
car::Anova(modelo1) 

library(vip)
vip(modelo1)
library(tidyverse)
vip(modelo1, geom = 'point', aesthetics = list(color = "tomato", size = 4))
vip(modelo1, aesthetics = list(color = "green2", fill = 'yellow2', size = 1.4))

library(randomForest)
rf = randomForest(resultado ~ publicid + gradnove + investigacion + sectecng + 
                    personal, data = base_entrenamiento)

vip(rf, aesthetics = list(color = 'darkorchid3', fill = '#A2CD5A', size = 1.4))

summary(modelo1)$coef


# visualisar resultados (gráfico de efectos parciales)
ggeffect(modelo1) %>%
  plot() %>%
  plot_grid()

marginalModelPlots(modelo1, col.line = c('#CD8162', '#3A5FCD'))
mmps(modelo1, col.line = c('tomato', 'green2'))

# En este ejemplo ocurrió que con el método step forward quedaron todas 
# las variables, por lo tanto, modelo y modelo1 son iguales.
model_performance(modelo)
model_performance(modelo1)
compare_performance(modelo, modelo1, rank=T)
plot(compare_performance(modelo, modelo1))
anova(modelo, modelo1)
anova(modelo1, modelo, test = "Chisq")
AIC(modelo)
AIC(modelo1)

# COEFICIENTES

str(modelo1)

summary(modelo1)
modelo1$coefficients
coef(modelo1)
coef(modelo1)[1]
coef(modelo1)[2]

rl_table = tbl_regression(modelo1)
rl_table

# ODDS RATIO

#---------------------------------------------------------
# Cociente de ventajas (Odd Ratio)
exp(coef(modelo1))
cbind(Coeficientes=modelo1$coef,ExpB=exp(modelo1$coef))
#----------------------------------------------------------
# Cociente de ventajas e IC 95% 
exp(cbind(OR = coef(modelo1),confint.default(modelo1,level = 0.95)))
?confint.default
#----------------------------------------------------------

rl_table1 = tbl_regression(modelo1, exponentiate = T)
rl_table1
exp(2.62059)
exp(-7.20829)
interpret_oddsratio(c(13.7, 2.62, 0.10, 0.00074, 1.03, 2.20))
?interpret_oddsratio

forest_model(modelo1, exponentiate = TRUE)

summ(modelo1, confint = TRUE, exp = TRUE)

plot_model(modelo1)
plot_model(modelo1,show.values = T,width=0.1)+
  ylab('Éxito o Fracaso de Empresas-"odds ratio"')

tab_model(modelo1,
          show.reflvl = T,
          show.intercept = F,
          p.style = 'numeric_stars')

library(ggstats)
ggcoef_model(modelo1)
ggcoef_model(modelo1, exponentiate = TRUE)

#####################################
# MISCELÄNEA con el paquete emmeans #
#####################################

# Considerando el sector tecnològico, hallar la PROBABILIDAD de obtener éxito 
# por categoría y comparar los ODDS en el sector tecnológico, cambiando la 
# referencia.

# $emmeans
# La probabilidad exacta de tener éxito al lanzar un nuevo producto al mercado 
# es 0.435 cuando el sector tecnológico es alto. La H0 (null = prob = 0.5), esto
# indica que la prob. de tener éxito es igual a la prob. de fracaso 
# (están 50 a 50).
# Ver el gráfico de abajo, la prueba es NS.

# $contrasts (OR)
# tener éxito al lanzar un nuevo producto al mercado es 3.370 más ventajoso
# cuando el ST es alto y no bajo (237% más ventajoso). 
# La H0: (null: OR = 1) está el el intervalo, la prueba es NS
# tener éxito al lanzar un nuevo producto al mercado es 0.289 menos ventajoso
# cuando el ST es bajo y no medio (0.711 menos ventajoso)

emmeans(modelo1, pairwise ~ sectecng,
        type = 'response', infer = T)

# Invertir el orden de los odds ratio si se necesita (con el sector tecnológico)

emmeans(modelo1, ~ sectecng, type = 'response') %>%
  pairs(reverse = T, infer = T)

# Hallar ODDSi = pi / (1-pi), por , i = 1, 2, 3 (con el sector tecnológico).
cld(emmeans(modelo1, pairwise ~ sectecng,
            type = 'response', infer = T))

# Para una empresa del ST alto es menos probable el èxito que el fracaso.
# Para una empresa del ST bajo es menos probable el èxito que el fracaso.
# Para una empresa del ST medio es menos probable el èxito que el fracaso.

emmeans(modelo1, ~ sectecng, infer = T) %>%
  as_tibble() %>%
  dplyr::select(sectecng, emmean) %>%
  mutate(odds = exp(emmean)) %>%
  mutate_if(is.numeric, ~ round(., 2))

0.435 / (1-0.435)
0.186 / (1-0.186)
0.441 / (1-0.441)

# MISCELÁNEA CON LA LIBRERÍA: emmeans

(emm1 = emmeans(modelo1, specs = 'sectecng', type = 'response'))
(emm2 = emmeans(modelo1, specs = 'sectecng', counterfactuals = 'sectecng'))

# Creando una cuadrícula de referencia
rg = ref_grid(modelo1)

# Prueba de hipótesis ordinaria H0: Pi = 0.6 
test(emm1, null = 0.6)
summary(emm1, null = 0.6, infer = c(TRUE, TRUE))
test(emm1, null = 0.6, side = 'left')

# Prueba de hipótesis ordinaria H0: Pi-Pj = 0.22 
emm4 = emmeans(modelo1, specs = 'sectecng')
# Creando una cuadrícula de referencia
rg.emm4 = regrid(emm4)
(pares = pairs(rg.emm4))
(pares1 = pairs(rg.emm4, null = 0.22))
(pares2 = pairs(rg.emm4, null = 0.22, side = 'great'))

# Pruebas de equivalencia para cada una de las 3 probabilidades
# Probar que las 3 probabilidades son equivalentes a 0.6 con una tolerancia de 
# ± 0.10 (es decir H0: Pi<=0.5 o Pi>=0.7  vs  H1: 0.5<Pi<0.7)
# H1 indica que la probabilidad es equivalente a 0.6

test(emm1, null = 0.6, delta = 0.10, side = 'equivalence')

# Pruebas de equivalencia para las diferencias de probabilidades por pares.
# Supongamos que una diferencia inferior a 0.58 es irrelevante
# es decir H0: abs(Pi-Pj) >= 0.58   vs  H1: abs(Pi-Pj) < 0.58
# H1 indica equivalencia

test(pares, delta = 0.58, side = 'equivalence', adjust = 'holm')

# con pairs se hace la equivalencia con diferencias de diferencias
pairs(pares, delta = 0.58, side = 'equivalence', adjust = 'holm')

# medias geométricas
(emm3 = emmeans(rg, specs = 'sectecng', type = 'response')) 
pairs(emm3) # se obtienen ratios
pairs(emm3, reverse = TRUE)

# medias aritméticas
(emm5 = emmeans(regrid(rg, transform = 'response'), specs = 'sectecng'))
pairs(emm5) # se obtienen diferencias
pairs(emm5, reverse = TRUE)
# medias aritméticas
(emm6 = emmeans(modelo1, specs = 'sectecng', regrid = TRUE))
pairs(emm6) # se obtienen diferencias
pairs(emm6, reverse = TRUE)

# Simulando una muestra de coeficientes de regresión.
# Resultados similares a emm1
set.seed(3.1416)
rgb = regrid(rg, N.sim = 400, transform = 'pass')
emmeans(rgb, specs = 'sectecng', type = 'response')
emm1

# Pendientes
emmeans(modelo1, pairwise ~ sectecng,
        infer = T)

# considerando algunos valores de las variables numéricas (se gana en 
# interpretación)
summary(datos)
tt(head(datos))
modelo1.rg = ref_grid(modelo1, at = list(publicid = c(1.612, 1.844, 4.224),
                                         personal = c(4, 10, 40)))
cld(emmeans(modelo1.rg, ~ personal, type = 'response'))
cld(emmeans(modelo1.rg, ~ publicid, type = 'response'))
cld(emmeans(modelo1.rg, ~ personal | investigacion, type = 'response'))
 
# Una pequeña miscelánea de gráficos

emmip(modelo1, sectecng ~ gradnove, type = 'response')
emmip(modelo1, investigacion ~ sectecng, type = 'response')
emmip(modelo1, tipo ~ gradnove, type = 'response')
emmip(modelo1.rg, investigacion ~ personal, type = 'response')
emmip(modelo1.rg, investigacion ~ personal, type = 'response', style = 'factor')
emmip(modelo1.rg, personal ~ investigacion | publicid, type = 'response')

emm = emmeans(modelo1, ~ investigacion | sectecng)
cld(emm)
plot(emm)
plot(emm, comparisons = TRUE, colors = '#EE3B3B')

cld(emm, type = 'response')
plot(emm, type = 'response')
plot(emm, type = 'response', comparisons = TRUE, colors = '#53868B')

#------------------------------------------------------------------#
#              Probabilidades y grupo estimados                    #
#------------------------------------------------------------------#

proba.pred = predict(modelo1, type="response")

#Probabilidades estimadas de pertenecer a la clase 1 (no a la clase de referencia)

head(proba.pred, 10)
clase.pred <- ifelse(proba.pred >= 0.5, 'exito', 'fracaso')

head(clase.pred,10)

finaldata = cbind(base_entrenamiento, proba.pred, clase.pred)
tt(head(finaldata,10))

library(RColorBrewer)
RColorBrewer::display.brewer.all()
library(ggthemes)

a=ggplot(finaldata, aes(x = proba.pred, fill = resultado)) + 
  geom_density(alpha = 0.7)+
  scale_fill_discrete(h=c(10,360)*25,c=30,l=20,h.start=0,
                      direction=1,na.value='grey50') +
  theme_economist()
a

b=ggplot(finaldata, aes(x = proba.pred, fill = resultado)) + 
  geom_histogram(alpha = 0.7)+
  scale_fill_discrete(h=c(10,360)*25,c=30,l=20,h.start=0,
                      direction=1,na.value='grey50') +
  theme_fivethirtyeight()
b

a + b

write.csv(finaldata,"Compras-Logistica-Probabilidades.csv")

#------------------------------------------------------------------
# Prediccion para nuevos individuos
# A continuacion se da el resultado de que una empresa tendra exito (1) con
# la probabilidad de que esa empresa tenga exito

# Manualmente

summary(modelo1)
b1 = -9.97725; b2 = 2.62059; b3 = 0.96447; b4 = -2.28516; b5 = -7.20829;
b6 = -1.21479; b7 = 0.02682; b8 = 0.79050

nuevo1 <- data.frame(publicid = 1.2, gradnove = "prod_nuevos", tipo = 'cons_indus',
                   investigacion = "si", sectecng = "alta", personal = 8)

z = b1 + b2*1.2 + b3*0 + b4*1 + b5*0 +  b6*0 + b8*8
ez = exp(z)
(Pe = ez / (1+ez)) 

# Usando la función predict
p = predict(modelo1, newdata = nuevo1, type = "response")# La empresa no tendra  exito
p 
#Considerando type="link"
log_odds <- predict(modelo1, newdata = nuevo1, type = "link")
exp(log_odds)/(1+exp(log_odds))

# Usando la distribución Logística
plogis(log_odds)
plogis(z)
 
odds <-p/(1-p)
odds

# Cuando una empresa tiene (publicid=1.2,gradnove="prod_nuevos", 
# tipo = "cons_indus", investigacion='si',secteng="alta",personal=8),la 
# probabilidad de tenga exito es 0.06119558 menos probable que no tenerlo, 
# al lanzar un nuevo producto al mercado

## Otro caso
nuevo2 <- data.frame(publicid = 3, gradnove = "prod_nuevos", 
                     tipo = "cons_final", investigacion = "si", 
                     sectecng = "media", personal = 4)

z1 = b1 + b2*3 + b3*0 + b4*0 + b5*0 +  b7*1 + b8*4
ez = exp(z1)
(Pe = ez / (1+ez))

# Usando la función predict
p = predict(modelo1, newdata = nuevo2, type = "response")# La empresa no tendra exito
p
#Considerando type="link"
log_odds <- predict(modelo1,newdata=nuevo2,type="link")
exp(log_odds)/(1+exp(log_odds))

# Usando la distribución Logística
plogis(log_odds)
plogis(z1)

odds <-p/(1-p)
odds

# Cuando una empresa tiene (publicid=3,gradnove="prod_nuevos", 
# tipo = "cons_final", investigacion='si',secteng="media",personal=4),la 
# probabilidad de que tenga exito es 2.925152 mas probable a que no lo tenga, 
# al lanzar un nuevo producto al mercado

# Guardar cualquier archivo 
#------------------------------
#setwd("G:/MODELO LOGIT")

# Guardar modelo en UCV O DISCO para utilizarlo en otra oportunidad
#------------------------------------------------------------------
saveRDS(modelo1, "./modelo_logit.rds")

# Leer el modelo final
#---------------------
super_model <- readRDS("./modelo_logit.rds")
print(super_model)


# Hacer predicciones
nuevos <- read.delim("clipboard")

final_predictions <- predict(super_model, nuevos, type = "response")
final_predictions
data.frame(nuevos, final_predictions)

# MODELO Log Binomial
# Como veremos, el Risk Ratio (RR) es más fácil de interpretar que el Odds Ratio

# Veamos el efecto no lineal en publicid
m = glm(resultado ~ poly(publicid, 3), 
        family = binomial(link = 'log'),
        data = base_entrenamiento,
        start = c(log(0.5),0 ,0, 0))

library(ggeffects)
ggeffect(m)
ggeffect(m) |> plot() # para ver intervalos de publicidad con efecto significativo

# Veamos el efecto no lineal en personal
m1 = glm(resultado ~ poly(personal, 3), 
         family = binomial(link = 'log'),
         data = base_entrenamiento,
         start = c(log(0.5),0 ,0, 0))

library(ggeffects)
ggeffect(m1)
ggeffect(m1) |> plot() # para ver intervalos de edades con efecto significativo

# Categorizando publicid y personal. Se pierde información pero se gana en
# interpretación
d = base_entrenamiento |>
  mutate(pub_cats = case_when(
    publicid < 4.5 ~ '<4.5',
    between(publicid, 4.5, 6) ~ '4.5-6',
    publicid > 6 ~ '>6'
  )) |>
  mutate(per_cats = case_when(
    personal < 20 ~ '<20',
    between(personal, 20, 50) ~ '20-50',
    personal > 50 ~ '>50'
  ))

# Para calcular los Risk Ratio (RR) ponemos 'log' en el link en lugar de logit
# 
m2 = glm(resultado ~ pub_cats + gradnove + tipo + investigacion + sectecng + 
           per_cats,
         data = d, family = binomial(link = 'log'),
         start = c(log(0.5), 0, 0, 0, 0, 0, 0, 0, 0, 0))

library(performance)
windows()
check_model(m2)

emmeans(m2, pairwise ~ pub_cats,
        type = 'response', infer = T)

# $emmeans
# Las probabilidades de obtener éxito al lanzar un nuevo producto al mercado son
# 0.471 (para publicidad menor de 4.5), 0.749 (para publicidad >6) y 0.855 (para 
# publicidad entre 4.5 y 6). El pvalor < 0.0001 indican que la probabilidad 
# de obtener éxito no es 1 (100%) (en la categoría <4.5). En las otras categorías
# no se rechaza la H0 de que la probabilidad de obtener éxito es 1.

# $contrasts
# Las empresas con publicidad <4.5 tenían una probabilidad de obtener éxito 
# 0.628 veces menor que las de publicidad >6 (menor en 37.2%). El pvalor es no 
# significativo (H0: RR = 1). Lo mismo ocurre con las empresas con 
# publicidad >6 tenían una probabilidad de éxito menor en 12.4% respecto a las 
# empresas con publicidad (4.5-6), el pvalor es significativo. Las empresas 
# con publicidad <4.5 tenían una probabilidad de tener éxito menor en 45% que
# las empresas con publicidad (4.5-6), el pvalor es significativo (H0: RR = 1).

# La librería logbin soluciona el problema de hallar valores iniciales pero a 
# veces, como ahora, se presentan problemas de convergencia.
library(logbin)
m3 = logbin(resultado ~ pub_cats + gradnove + tipo +investigacion + sectecng + 
              per_cats, data = d)

emmeans::emmeans(m3, pairwise ~ pub_cats,
        type = 'response', infer = T)

# Visualisar resultados del modelo
library(sjPlot)
fancy_plot = ggeffect(m2) |> plot() |> plot_grid()
fancy_plot

fancy_plot = ggeffect(m3) |> plot() |> plot_grid()
fancy_plot

nice_table = tbl_regression(m2,
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

############################################
#  Indicadores para Evaluación de Modelos  #
############################################

# Se utiliza la data de prueba para validar el modelo
# Tabla de confusión con los datos de prueba

modelo2<-glm(resultado ~ publicid + gradnove + tipo + investigacion + sectecng + 
               personal, family = binomial(link = logit),data=base_prueba)

#------------------------------------------------------------------#
#              Probabilidades y grupo estimados                    #
#------------------------------------------------------------------#

proba.pred1 = predict(modelo2, type = "response")

#Probabilidades estimadas de pertenecer a la clase 1 (no referencia)

head(proba.pred1, 10)
clase.pred1 <- ifelse(proba.pred1 >= 0.5, 'exito', 'fracaso')
head(clase.pred1,10)

finaldata1 = cbind(base_prueba, proba.pred1, clase.pred1)
head(finaldata1,10)

c=ggplot(finaldata1, aes(x = proba.pred1, fill = resultado)) + 
  geom_density(alpha = 0.9)
c

d=ggplot(finaldata1, aes(x = proba.pred1, fill = resultado)) + 
  geom_histogram(alpha = 0.9)
d

c | d

# Comparando lo que se obtuvo con la data de entrenamiento
(a + b)/(c + d)


##################################
# 1. Tabla o Matriz de Confusion #
##################################

table(clase.pred1, base_prueba$resultado)

d2 = base_prueba %>% 
  mutate(pred_probs = predict(modelo2, base_prueba, 'response')) %>%
  mutate(pred_classes = ifelse(pred_probs > 0.5, 'exito', 'fracaso')) %>%
  mutate(resultado = factor(resultado, levels = c('exito', 'fracaso')),
         pred_classes = factor(pred_classes, levels = c('exito', 'fracaso')))

# Tabla o matriz de confusión (esta es la forma apropiada, convencional)

cm2 = xtabs(data = d2, ~ pred_classes + resultado)
confusionMatrix(cm2, mode = 'everything', positive = 'exito')

# Calcular el accuracy
#-----------------------
accuracy <- mean(base_prueba$resultado == clase.pred1)
accuracy# Es la proporción de todas las predicciones correctas

# Calcular el error de mala clasificacion
#----------------------------------------
error <- mean(base_prueba$resultado != clase.pred1)
error# Es la proporción de todas las predicciones incorrectas
is.factor(clase.pred1)
clase = as.factor(clase.pred1)
library(caret)
cm = confusionMatrix(clase, base_prueba$resultado, positive="exito")
print(cm)

# No Information Rate (NIR) (Tasa de falta de información): es el Accuracy que
# se obtendría si siempre predijera la clase mayoritaria. Con esta data la clase
# mayoritaria es éxito (2+25= 27). El NIR ayuda a verificar si el Accuracy del
# modelo (0.913) es mejor que el  del modelo base sin información (0.587). Para
# verificar si la diferencia (0.913-0.587=0.326) es significativa se halla el 
# p-value [Acc > NIR] con una prueba binomial
(NIR = 27/46)
0.913-0.587

# H0: La precisión del modelo <= NIR
binom.test(
  x = 42,     # Número de clasificaciones correctas
  n = 46,     # Número de ensayos
  p = 0.587,  # Probabilidad del suceso esperado (en las hipótesis)
  alternative = 'greater')

# Sensibilidad: # Probabilidad de que un una empresa que tuvo éxito, sea bien 
# identificada
# Especificidad: Probabilidad de que una empresa que tuvo fracaso sea 
# correctamente clasificada
# Valor predictivo positivo (precision): Probabilidad de que un resultado 
# positivo (éxito) del modelo sea verdadero con las siguientes empresas que
# lanzarán un nuevo producto.
# Valor predictivo negativo: Probabilidad de que un resultado negativo (fracaso)
# del modelo sea verdadero con las siguientes empresas que lanzarán un nuevo 
# producto.
# Prevalencia: Proporción de casos positivos reales (éxitos). Es la proporción
# real de empresas que tuvieron éxito. 

# Prueba de Mcnemar
# H0: Las frecuencias de éxitos no ha variado (No hay diferencias significativas
#     entre los valores predichos y los reales)
# H0: Las frecuencias de éxitos si ha variado

bp_r = base_prueba$resultado
df = tibble(clase.pred1, bp_r)
# Los siguientes resultados no coinciden con la librería caret
mcnemar.test(clase.pred1,bp_r,correct = F)
mcnemar.test(clase.pred1,bp_r,correct = T)

ggbarstats(
  data=df,
  x=clase.pred1,
  y=bp_r,
  paired=T,
  label = 'both'
)

# Los siguientes resultados (con la corrección) si coinciden con la librería caret
tabla = matrix(c(17,2,2,25), 2, 2)
mcnemar.test(tabla, correct = F)
mcnemar.test(tabla, correct = T)

# Otra presentacion 

a = d2$pred_classes
is.factor(a)
b = d2$resultado
is.factor(b)

bin_cm=ConfusionTableR::binary_class_cm(a,b)
bin_cm$record_level_cm
library(dplyr)
glimpse(bin_cm$record_level_cm)
ConfusionTableR::binary_visualiseR(train_labels = a,
                                   truth_labels = b,
                                   class_label1 = "Exito",
                                   class_label2 = "Fracaso",
                                   quadrant_col1 = "#28ACB4",
                                   quadrant_col2 = "#4397D2",
                                   custom_title = "Matriz de Confusion-Empresas",
                                   text_col = "black",
                                   cm_stat_size = 1.2,positive="exito")

par(mfrow=c(1,1))

# Coeficiente de correlación de Mathews (MCC): mide la correlación 
# dentre las clases verdaderas y las clases predichas
library(mltools)

mcc(preds = d2$pred_classes, actuals = base_prueba$resultado)

#####################################
# Curva ROC y Area bajo la curva    #
#####################################

library(ROSE)
#Capacidad discriminante del modelo
roc.curve(d2$resultado, d2$pred_probs, lty=2,lwd=1.8,col="blue" ,
          main="ROC curves")

#Capacidad discriminante del modelo

library(caTools)
colAUC(d2$pred_probs, d2$resultado,plotROC = TRUE)
abline(0, 1,col="tomato",lty=4)

library(pROC)
df = tibble(a, b)
roc(b ~ fitted.values(modelo2), data = df,
    plot = TRUE, legacy.axes = TRUE,
    print.auc = TRUE, ci = TRUE)


################################################################################
## Descubra el poder secreto de DOR, LR+ y LR- en una Matriz de Confusión.    ##
##                         Más allá que el  Accuracy                          ##
################################################################################

## Positive likelihood ratio (LR+): La razón de verosimilitud positiva Responde
# a la pregunta ¿cuánto se debe confiar en un resultado positivo de una prueba.
# Nos dice que tan probable es un resultado positivo en una persona enferma en
# comparación con una sana.
# LR+ = TPR/FPR = Sensitivity/(1 - Specifity) = (TP/(TP+FN)) / (1 - TN/(TN+FP))

# Nota: TPR, FPR, FNR y TNR se definirán más abajo.

# Rango de la LR+              INTERPRETACIÓN
#     < 1                      La prueba es peor que adivinar al azar. Evitarla.
#     = 1                      La prueba no ofrece información útil. Ignorarla.
#     Entre 1 y 5              La prueba es algo útil pero no definitiva.
#                              Usarla con precaución
#     Enter 5 y 10             La prueba es moderadamente fuerte
#     > 10                     La prueba es muy fuerte, significa que un 
#                              resultado positivo es altamente confiable.

## Negative likelihood ratio (LR-): La razón de verosimilitud negativa Responde
# a la pregunta ¿cuánto se debe confiar en un resultado negativo de una prueba.
# Nos dice cuánto más probable es un resultado negativo en una persona sana en
# comparación con una persona enferma. Qué tan bien un modelo predice los 
# negativos cuando son realmente negativos en comparación con cuando son 
# positivos.
# LR- = FNR/TNR = (1-Sensitivity)/Specifity = (1-TP/(TP+FN)) / TN/(TN+FP)

# Rango de la LR-                   INTERPRETACIÓN
#     < 0.1                         La prueba es excelente en detectar la 
#                                   condición de negativo.
#     De 0.1 a 0.5                  La prueba es decente pero no perfecta.
#     De 0.5 a 1                    La prueba tiene un valor limitado
#     > 1                           La prueba no es confiable para detectar 
#                                   negativos.

# Una buena prueba de diagnóstico o modelo de aprendizaje automático debe tener 
# LR+ > 1 y LR- < 1. Si LR+ < 1 y LR- > 1 la prueba o modelo es esencialmente 
# inútil.

## Diagnostic Odds Ratio (DOR): Al igual que el LR+ y el LR- la razón de odds 
# (chance) diagnóstica (DOR) sirve para evaluar que tan bien se desempeña una 
# prueba o un modelo. La ventaja del DOR es que no direcciona hacia positivos o 
# negativos. El DOR es una métrica que resume el desempeño general de una 
# prueba o modelo sin dividirlo en direcciones positivas o negativas (cuánto 
# mejor discrimina entre positivos y negativos).

# DOR = [TP / FN] / [FP / TN] o también
# DOR = LR+ / LR-

# DOR = (TP/FN)/(FP/TN) = (TP/FP)(FN/TN) = (TP*TN)/(FP*FN) 
# DOR = (Sensitivity*Specifity) / ((1-Sensitivity)*(1-Specifity))

# Rango de la DOR (DOR)     INTERPRETACIÓN
#     < 1                   La prueba es peor que adivinar al azar. Evitarla.
#     = 1                   La prueba no ofrece información útil. Ignorarla.
#     Entre 1 y 5           La prueba es algo útil pero no definitiva. Usarla 
#                           con precaución
#     Enter 5 y 10          La prueba es moderadamente fuerte
#     > 10                  La prueba es muy fuerte, significa que discrimina 
#                           muy bien entre positivos y negativos.

# El DOR es robusto ante ligeros desequilibrios en los datos. Puede volverse 
# inestable con tamaños de muestra pequeños. Un desbalance extremo también 
# (por ejemplo 99% sanos y 1% enfermos) vuele inestable al DOR.

# USO del R

library(tidyverse)
library(caret)
library(epiR)      # para la función epi.tests()
library(tinytable)
theme_set(theme_bw())

# obtener datos
d = carData::TitanicSurvival|> filter(!(is.na(age)))
tt(head(d))
dim(d)

# dividir los datos en entrenamiento (80%) y prueba (20%)
set.seed(42)
train_index = createDataPartition(d$survived, p = 0.80, list = FALSE)
titanic_train = d[train_index, ]
titanic_test = d[-train_index, ]

dim(titanic_train)[1]/dim(d)[1]
dim(titanic_test)[1]/dim(d)[1]

prop.table(table(d$survived))
prop.table(table(titanic_train$survived))
prop.table(table(titanic_test$survived))

# modelo de regresión logística dicotómica considerando interacciones
model = glm(survived ~ sex * age * passengerClass, titanic_train,
            family = binomial)

# realizar predicciones en el conjunto de datos de prueba
d2 = titanic_test |>
  mutate(pred_probs = predict(model, titanic_test, type = 'response')) |>
  mutate(predicted = ifelse(pred_probs > 0.5, 'yes', 'no')) |>
  mutate(survived = factor(survived, levels = c('yes', 'no')),
         predicted = factor(predicted, levels = c('yes', 'no')))

tt(head(d2))

# crer la matriz de confusión y obtener resultados
cm = table(d2$predicted, d2$survived)
confusionMatrix(cm, mode = 'everything')
epi.tests(cm)

# Cálculo de LR+, LR- y DOR

TP = 46; FN = 39; FP = 8; TN = 115 ; Sensitivity = 0.5412 ; Specifity = 0.935
VPnegativo = 0.7468; Precision = 0.8519

# Razón de verosimilitud positiva y negativa

# un resultado positivo (sobrevivir) es moderadamente confiable
(LRmas = Sensitivity/(1-Specifity) ) 
(LRmas = (TP/(TP+FN)) / (1 - TN/(TN+FP)))

# un resultado negativo (no sobrevivir) es moderadamente confiable
(LRmenos = (1-Sensitivity)/Specifity)
(LRmenos = (1-TP/(TP+FN)) / (TN/(TN+FP)))

summary(epi.tests(cm))

# Razón de chance diagnóstica (DOR)

# La prueba discrimina bien entre positivos y negativos (> 10)
(DOR = LRmas / LRmenos)
(DOR = (TP/FN)/(FP/TN)) 
(DOR = (TP/FP)/(FN/TN)) 
(DOR = (TP*TN)/(FP*FN))
(DOR = (Sensitivity*Specifity) / ((1-Sensitivity)*(1-Specifity)))

# False Positive Rate (FPR): Tasa de falsos positivos.
# Es el porcentaje de verdaderos negativos clasificados incorrectamente como
# positivos. Es el porcentaje de personas sanas etiquetadas erroneamente como
# enfermas. En Estadística se conoce com Error tipo I o probabilidad de una 
# falsa alarma.

# False Negative Rate (FNR): Tasa de falsos negativos.
# Es el porcentaje de verdaderos positivos clasificados incorrectamente como 
# negativos. Es el porcentaje de personas enfermas etiquetadas erroneamente como
# sanas. En Estadística se conoce como error tipo II o probabilidad de un error.
# La FNR es peor que la FPR: por ejemplo si a alguien se le diagnostica cáncer 
# pero después de 4 meses se le dice que no tiene cáncer, esta es una situación
# que no es tan grave si lo comparamos con que a ese paciente se le dice que 
# está sano y después de 4 meses se le dice que tiene cáncer, la enfermedad puede haber
# empeorado.

# False Discovery Rate (FDR): Tasa de falsos descubrimientos.
# Es porcentaje de predicciones positivas que en realidad son falsos positivos
# entre todas las predicciones positivas. Es el porcentaje de individuos sanos
# identificados incorrectamente como enfermas.

# False Omission Rate (FOR): Tasa de falsas omisiones.
# Mide la proporción de predicciones negativas que son incorrectas. Porcentaje
# de personas enfermas etiquetadas incorrectamente como sanas

# USO del R

# Cálculos

# El 6.5% de los que realmente no sobrevivieron han sido clasificados como 
# sobrevivientes
(FPR = 1-Specifity)

# El 45.88% de los que realmente sobrevivieron han sido clasificados como no 
# sobrevivientes
(FNR = 1-Sensitivity)

# el 14.81% de todas las predicciones positivas fueron incorrectas
# el 14.81% de todas las predicciones de sobrevivencia fueron incorrectas
(FDR = 1-Precision) 

# el 25.32% de todas las predicciones negativas fueron incorrectas
# el 25.32% de todas las predicciones de no sobrevivencia fueron incorrectas
(FOR = 1-VPnegativo) 

#####################################################################
##  Métricas del error para el aprendizaje automático en el modelo ##
#####################################################################

# En la siguiente situación el accuracy falla, porque los datos 
# están desbalanceados. El accuracy no siempre es útil.
# Cada prueba médica o modelo de clasificación tiene un objetivo: Cometer el
# menor la menor cantidad de errores posible, allí entra la tasa de error o 
# tasa de clasificación errónea que es la proporción de predicciones erróneas:
# error_rate = (FP+FN)/total

# Número total de casos
total_cases = 1000

# Porcentaje de casos positivos
positive_percentage = 0.05

# Cálculo del número de casos positivos y negativos
num_positive = total_cases*positive_percentage # 50 +
num_negative = total_cases-num_positive        # 950 -

# Crear etiquetas actuales: 1 para positivo, 0 para negativo
actual = c(rep(1, num_positive), rep(0, num_negative))

# Predictor que siempre predice 'negativo' (0)
predicted = rep(0, total_cases)

library(caret)

conf_matrix = confusionMatrix(factor(predicted), factor(actual))

print(conf_matrix)

# Las siguientes tasas muestran exactamente dónde tiene dificultades el modelo

# Misclassification Rate (Error_Rate): Tasa de clasificación errónea
# error_rate = (FP+FN)/total
# error_rate = 1-accuracy
TP = 950; TN = 0; FP = 50; FN = 0; total = 1000
(Sensitivity = TP/(TP+FN))
(Specifity = TN/(TN+FP))
(Precision = TP/(TP+FP)) # VP+: frecuencia a la que una predicción positiva es correcta
(VPnegativo = TN/(TN+FN))

error_rate = function(TP, TN, FP, FN) {
  (FP+FN)/(TP+TN+FN+FP)
}

error_rate(TP, TN, FP, FN)


