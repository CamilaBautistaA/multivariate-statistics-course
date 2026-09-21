###################################################
#                                                 #
#            MANOVA Y MANCOVA                     #
#            EJEMPLOS: Varios                     #
# Profesor: Clodomiro Fernando Miranda Villagomez #
#         cfmiranda@lamolina.edu.pe               #
#                                                 #
###################################################

######################
#   MANOVA EN DCA   ##
######################

rm(list = ls())   

library(foreign)
library(tinytable)
library(tidyverse)
library(mvnormtest)
library(heplots)
library(biotools)
library(covTestR)
library(DFA.CANCOR)
library(psych)
library(MVTests)
library(emmeans)

datos <-read.spss("1.Maiz.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
str(datos)

tt(datos)
head(datos, 20)

# Supuesto de normalidad multivariada
# La siguiente funcion (mshapiro.test) se aplica grupo a grupo, 
# por lo tanto, primero es necesario
# dividir la base de datos en los grupos, seis en el ejemplo.

trat1 = datos %>% filter(Tratamientos == "Hibrido1") %>%
  dplyr::select(Y1, Y2, Y3, Y4, Y5)
trat2 = datos %>% filter(Tratamientos == "Hibrido2") %>%
  dplyr::select(Y1, Y2, Y3, Y4, Y5)
trat3 = datos %>% filter(Tratamientos == "Hibrido3") %>%
  dplyr::select(Y1, Y2, Y3, Y4 ,Y5)
trat4 = datos %>% filter(Tratamientos == "Hibrido4") %>%
  dplyr::select(Y1, Y2, Y3, Y4, Y5)
trat5 = datos %>% filter(Tratamientos == "Hibrido5") %>%
  dplyr::select(Y1, Y2, Y3, Y4, Y5)
trat6 = datos %>% filter(Tratamientos == "Hibrido6") %>%
  dplyr::select(Y1, Y2, Y3, Y4, Y5)

#Ejecutamos el test de normalidad multivariada

mshapiro.test(t(trat1))
mshapiro.test(t(trat2))
mshapiro.test(t(trat3))
mshapiro.test(t(trat4))
mshapiro.test(t(trat5))
mshapiro.test(t(trat6))

#Con cada hibrido no hay normalidad entonces deberian plantearse 
#transformaciones o usar las pruebas
#robustas que se describen más abajo.

# Supuesto de homogeneidad de matrices variancia covariancia

res <- boxM(datos[, -1], datos[, "Tratamientos"])
res
summary(res)

heplots::boxM(cbind(Y1, Y2, Y3, Y4, Y5) ~ Tratamientos, data = datos)

boxM(datos[-1], datos[, 1])

maiz <- unique(datos$Tratamientos)
maiz1 <- lapply(maiz,
         function(x){as.matrix(datos[datos$Tratamientos == x, 2:6])}
)

names(maiz1) <- maiz
Ahmad2017(maiz1)

## Prueba Wrapper
homogeneityCovariances(datos, group = Tratamientos, covTest = BoxesM)

HOMOGENEITY(data = datos,groups = 'Tratamientos', 
            variables = c('Y1','Y2','Y3','Y4','Y5'))

# Supuesto de variables dependientes correlacionadas. 
# Prueba de esfericidad de Bartlett

options(scipen = 0)
cortest.bartlett(cor(datos[, -1]), n = nrow(datos[, -1]))

res1 = Bsper(datos[, -1])
summary(res1)

#Trabajando con el modelo de MANOVA en DCA
modelo = manova(cbind(Y1, Y2, Y3, Y4, Y5) ~ Tratamientos, data = datos)

#Determinación de la matriz residual y la matriz factorial del MANOVA.
str(modelo)
Matrices = summary(modelo)$SS
F = Matrices$Tratamientos
W = Matrices$Residuals

#Variabilidad explicada por el factor (Tratamientos). Matriz suma de 
#cuadrados y productos cruzados del factor (SCOCF)
F

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)
#del factor 
T = F + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte
#de la variabilidad total puede atribuirse al factor, mientras que un 
#valor proximo a 0 significa que el factor explica muy poco
#de esa variabilidad total.

eta2 = 1 - det(W)/det(T)
eta2
det(F)/det(T) # La razon de determinantes no es aditiva como las
#suma de cuadrados.

# Pruebas de hipotesis del modelo. Calculo de contrastes del modelo

#Contrastes del modelo en relacion a los supuestos. Todos los 
#estadisticos son bastante robustos ante violaciones de normalidad,
#y la prueba de Roy es muy sensible a violaciones de la hipotesis
#de la matriz de covariancias. Cuando las muestras son iguales por
#grupo, la prueba de Pillai es el estadistico más robusto ante
#violaciones de los supuestos.

k = 6 #numero de grupos
p = 5 #numero de variables
n = 6 #numero de observaciones por grupo
datosc1 = datos
head(datosc1)
datosc1$Tratamientos <- as.numeric(datosc1$Tratamientos)
str(datosc1)

VMPG = matrix(NA, k, p) #vector de medias por grupo
for(i in 1:k){
  VMPG[i,]=colMeans(datos[datosc1$Tratamientos == i, -1])
}
VMPG #cada fila es un vector de medias
   
#Computar B
(B=n*(k-1)*cov(VMPG))
# Tambien
n*(t(VMPG)-colMeans(VMPG))%*%t(t(VMPG)-colMeans(VMPG))
  
# Computar W
W = (n-1)*cov(datos[datosc1$Tratamientos == 1, -1])
for(i in 2:k){
    W = W + (n-1)*cov(datos[datosc1$Tratamientos == 1, -1])
}
W
W+B

## autovalores de W^{-1}B
(lambdas = eigen(solve(W)%*%B)$values)

##traza de pillai
sum(lambdas/(1+lambdas))
#tambien
sum(diag(solve(W+B)%*%B))

#El software R obtiene un valor calculado de pillai diferente al que se calculó 
#manualmente, esto se debe a que hay varias formas de calcular este valor
summary(modelo, test = "Pillai")
1-pf(1.6833, 25,150)

##Lambda de Wilks
det(W)/det(W + B)
#tambien
prod(1/(1 + lambdas))

#El software R obtiene un valor calculado diferente al que se calculó 
#manualmente, esto se debe a que hay varias formas de calcular este valor
summary(modelo, test = "Wilks")
1-pf(1.8417, 25,98.088)

## Lawley Hotelling
LH = sum(lambdas)
LH

summary(modelo, test = "Hotelling-Lawley")
1-pf(1.9047, 25,122)

## Raiz mayor de Roy
lambdas[1]/(1 + lambdas[1])
summary(modelo, test = "Roy")
1-pf(6.7982, 5, 30)

#A continuacion se aprecia que las 5 variables respuesta son significativas
#lo que indica que las 5 variables contribuyeron para que se rechace la H0 
#que dice que los 6 vectores de promedios son iguales y se concluye que al 
#menos uno de los vectores de promedios es diferente
summary.aov(modelo)

# El manova resultó Significativo a un nivel de significación de 0.05, entonces
# las comparaciones por pares se deben hacer

#Comparar dos hibridos

modelo1 = manova(cbind(Y1, Y2, Y3, Y4, Y5) ~ Tratamientos, data = datos,
               subset = Tratamientos %in% c("Hibrido3","Hibrido4"))

summary(modelo1,test="Pillai")
summary(modelo1,test="Wilks")
summary(modelo1,test="Hotelling-Lawley")
summary(modelo1,test="Roy")
summary.aov(modelo1)

#---------------------------------------------------#
#           Comparación General por pares           #
#---------------------------------------------------#

# H0:los 2 vectores de promedios son iguales
# H1:los 2 vectores de promedios difieren

Tratam <- c("Hibrido1", "Hibrido2", "Hibrido3", "Hibrido4", "Hibrido5", "Hibrido6")
comb<-t(combn(length(Tratam) , 2))

for(i in 1:nrow(comb)){
  modelo.comp = manova(cbind(Y1, Y2, Y3, Y4, Y5) ~ Tratamientos, data=datos,
                     subset = Tratamientos %in% Tratam[comb[i,]])
  print(paste("Trat: ", Tratam[comb[i,]][1], "y", Tratam[comb[i,]][2]))
  print(summary(modelo.comp, test = "Pillai"))
  cat("\n")
  
}
trat6

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(Y1, Y2, Y3, Y4, Y5) ~ Tratamientos, data = datos)
emm = emmeans(fit, ~ Tratamientos | rep.meas)
emm
names(datos)
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'none')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'bonferroni')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'holm')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'fdr')

# Asumiendo que el tratamiento Hibrido1 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas')

# Asumiendo que el tratamiento Hibrido4 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas', 
           ref = 'Hibrido4', adjust = 'bonferroni')

# mult.name toma por defecto rep.meas
mvcontrast(emm, method = 'consec')
mvcontrast(emm, method = 'consec', show.ests = TRUE)

#########################
##   MANOVA EN DBCA    ##
#########################

# Leer datos de panetones

datos <- read.delim("clipboard")
tt(datos)
str(datos)

datos$Preparacion <- factor(datos$Preparacion,levels=c(1,2,3,4),
  labels=c("prep1","prep2","prep3","prep4"))
datos$Publico <- factor(datos$Publico,levels=c(1,2,3,4,5),
  labels=c("publ1","publ2","publ3","publ4","publ5"))
str(datos)
tt(head(datos))

#-------------------------------------------
#Supuesto de normalidad multivariada       -
#-------------------------------------------

# H0: Si hay normalidad
# H1: No hay normalidad

#La siguiente funcion (mshapiro.test) se aplica grupo a grupo,
#por lo tanto, primero es necesario dividir la base de datos en
#los grupos, cuatro en el ejemplo.

# sólo con las variables respuesta

trat1 = datos %>% filter(Preparacion == "prep1") %>%
  dplyr::select(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor)
trat2 = datos %>% filter(Preparacion == "prep2") %>%
  dplyr::select(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor)
trat3 = datos %>% filter(Preparacion == "prep3") %>%
  dplyr::select(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor)
trat4 = datos %>% filter(Preparacion == "prep4") %>%
  dplyr::select(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor)

#Ejecutamos el test

library(mvnormtest)
mshapiro.test(t(trat1))
mshapiro.test(t(trat2))
mshapiro.test(t(trat3))
mshapiro.test(t(trat4))

#Con las formas de preparacion no hay normalidad entonces
#deberian plantearse transformaciones o usar las pruebas
#robustas que se describen mas abajo.

# Supuesto de homogeneidad de matrices variancia covariancia

res <- boxM(datos[, 2:5], datos[, "Preparacion"])
res
summary(res)

heplots::boxM(cbind(PuntajeMasa, PuntajeSabor, PuntajeRelleno,
        PuntajeOlor) ~ Preparacion, data = datos)

boxM(datos[2:5], datos[, 1])

paneton <- unique(datos$Preparacion)
paneton1 <- lapply(paneton,
    function(x){as.matrix(datos[datos$Preparacion == x, 2:5])}
)

names(paneton1) <- paneton
Ahmad2017(paneton1)

## Prueba Wrapper
homogeneityCovariances(datos[, -6], group = Preparacion, covTest = BoxesM)

library(DFA.CANCOR)
HOMOGENEITY(data = datos[, -6], groups = 'Preparacion', 
            variables = c('PuntajeMasa', 'PuntajeSabor',
            'PuntajeRelleno', 'PuntajeOlor'))

#----------------------------------------------------------------------------------------------------
# Supuesto de variables dependientes correlacionadas. 
# Prueba de esfericidad de Bartlett
#----------------------------------------------------------------------------------------------------

# H0:no estan correlacionas  det(Rp)=1
# H1:si estan correlacionas

options(scipen = 0)
cortest.bartlett(cor(datos[,c(-1, -6)]), n = nrow(datos[, c(-1, -6)]))

res1 = Bsper(datos[, c(-1, -6)])
summary(res1)

# pvalor = 4.875483e-10<0.05, RH0, las variables respuesta 
# si estan correlacionas


#--------------------------------#
# Trabajando con el modelo      ##
#--------------------------------#

modelo = manova(cbind(PuntajeMasa, PuntajeSabor, PuntajeRelleno,
                    PuntajeOlor) ~ Preparacion + Publico, data = datos)

#Determinacion de la matriz residual y la matriz factorial del MANOVA.
str(modelo)
Matrices = summary(modelo)$SS
Matrices
F = Matrices$Preparacion
B = Matrices$Publico
W = Matrices$Residuals

#Variabilidad explicada por el factor (Preparacion). 
#Matriz suma de cuadrados y productos cruzados del factor (SCOCF)
F

#Variabilidad explicada por los Publico. Matriz suma de cuadrados
#y productos cruzados de Publico (SCOCB)
B

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados
#total (SCOCT) del factor 
T = F + B + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte
#de la variabilidad total puede atribuirse al factor, mientras que
#un valor proximo a 0 significa que el factor explica muy poco
#de esa variabilidad total.

eta2 = 1 - det(B + W)/det(T)
eta2
# [1]  0.8208076
#La bondad de ajuste, quiere decir que el 82%
#de la variabilidad total se puede atribuir al factor
#mientras mas cercano a 1 es mucho mejor


#--------------------------------------
# Pruebas de hipótesis del modelo.   -- 
# Cálculo de contrastes del modelo   --
#--------------------------------------

#Para: Preparacion
# H0: los vectores de promedios son iguales
# H1: al menos un vector es diferente

#Para: Publico
# H0: no existe diferencia entre Publico
# H1: si existe diferencia entre Publico, se justifica el bloqueo

summary(modelo, test = "Pillai")

#los vectores son iguales
#no existe diferencia entre Publico, no se justifica el bloqueo

summary(modelo, test = "Wilks")

#no se ha hecho bien en bloquear, ya que no hoy diferencias entre
#las categorias de Publico
#se quiere que haya homogeneidad dentro de Publico y diferencia
#entre Publico

summary(modelo, test = "Hotelling-Lawley")
summary(modelo,test = "Roy")

#Para cada variable respuesta
summary.aov(modelo)

# El manova resultó NS a un nivel de significación de 0.05, entonces las 
# comparaciones por pares no se deben hacer, pero se harán para ilustración.

#---------------------------------#
#      Comparacion por pares      #
#---------------------------------#

# H0:los 2 vectores son iguales
# H1:los 2 vectores difieren

modelo1 = manova(cbind(PuntajeMasa, PuntajeSabor, PuntajeRelleno,
        PuntajeOlor) ~ Preparacion + Publico, data = datos,
        subset = Preparacion %in% c("prep1", "prep2"))

summary(modelo1, test = "Pillai")
summary(modelo1, test = "Wilks")
summary(modelo1, test = "Hotelling-Lawley")
summary(modelo1, test = "Roy")
summary.aov(modelo1)

#---------------------------------------------------#
#           Comparación General por pares           #
#---------------------------------------------------#

# H0:los 2 vectores son iguales
# H1:los 2 vectores difieren

Prepar <- c("prep1", "prep2", "prep3", "prep4")
comb <- t(combn(length(Prepar), 2))

for(i in 1:nrow(comb)){
  modelo.comp = manova(cbind(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor) ~
                       Preparacion + Publico, data = datos,
                     subset = Preparacion %in% Prepar[comb[i,]])
  print(paste("Trat: ",Prepar[comb[i,]][1], "y", Prepar[comb[i,]][2]))
  print(summary(modelo.comp, test = "Pillai"))
  cat("\n")
  
}

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(PuntajeMasa, PuntajeSabor, PuntajeRelleno, PuntajeOlor) ~ 
           Preparacion + Publico, data = datos)
emm = emmeans(fit, ~ Preparacion | rep.meas)
emm
names(datos)

mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'none')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'bonferroni')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'holm')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'fdr')

# Asumiendo que el tratamiento prep1 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas')

# Asumiendo que el tratamiento prep2 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas', 
           ref = 'prep2', adjust = 'holm')

# mult.name toma por defecto rep.meas
mvcontrast(emm, method = 'consec', adjust = 'holm')
mvcontrast(emm, method = 'consec', show.ests = TRUE, adjust = 'holm')

##########################################
##    MANOVA CON FACTORIALES EN DCA     ##
##########################################

datos <-read.spss("1.Pimenton.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
str(datos)
tt(datos)

# Supuesto de normalidad multivariada

# Desagregamos la base por grupos por grupos: Temperatura

trat1 = datos %>% filter(Temperatura == "40C", Tiempo == "30 min") %>%
  dplyr::select(Caroteno, Vitamina1)
trat2 = datos %>% filter(Temperatura == "80C", Tiempo == "30 min") %>%
  dplyr::select(Caroteno, Vitamina1)
trat3 = datos %>% filter(Temperatura == "40C", Tiempo == "60 min") %>%
  dplyr::select(Caroteno, Vitamina1)
trat4 = datos %>% filter(Temperatura == "80C", Tiempo == "60 min") %>%
  dplyr::select(Caroteno, Vitamina1)

#Ejecutamos el test

mshapiro.test(t(trat1))
mshapiro.test(t(trat2))
mshapiro.test(t(trat3))
mshapiro.test(t(trat4))

#No hay normalidad con los trat1 y trat2 entonces deberían
#plantearse transformaciones o usar las pruebas mas robustas 
#descritas mas arriba.

# Supuesto de homogeneidad de matrices variancia covariancia

heplots::boxM(cbind(Caroteno, Vitamina1) ~ Temperatura * Tiempo, data = datos)

# Supuesto de variables dependientes correlacionadas. Prueba de esfericidad de Bartlett

datos1=datos[-1:-2]
head(datos1)

options(scipen = 0)
cortest.bartlett(cor(datos1), n = nrow(datos1))

#Determinacion de la matriz residual y las matrices factoriales 
#(Temperatura, Tiempo y Temperatura*Tiempo) del MANOVA.

modelo = manova(cbind(Caroteno, Vitamina1) ~ Temperatura * Tiempo, data = datos)
str(modelo)
Matrices = summary(modelo)$SS
Matrices
FTemperatura = Matrices$Temperatura
FTiempo = Matrices$Tiempo
FTemperaturaTiempo = Matrices$`Temperatura:Tiempo`
W = Matrices$Residuals

#Variabilidad explicada por el factor (Temperatura). Matriz suma de cuadrados
#y productos cruzados del factor (SCOCFTemperatura)
FTemperatura

#Variabilidad explicada por el factor Tiempo. Matriz suma de cuadrados
#y productos cruzados del factor Tiempo (SCOCTiempo)
FTiempo

#Variabilidad explicada por la interaccion Temperatura*Tiempo. Matriz suma de
#cuadrados y productos cruzados de la interaccion Temperatura*Tiempo 
#(SCOCTemperatura*Tiempo)

FTemperaturaTiempo

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)
#del factor 
T = FTemperatura + FTiempo + FTemperaturaTiempo + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte de la variabilidad
#total puede atribuirse al factorial, mientras que un valor proximo a 0 significa
#que el factorial explica muy poco de esa variabilidad total.

eta2 = 1 - det(W)/det(T)
eta2

# Pruebas de hipotesis del modelo

summary(modelo, test = "Pillai")
summary(modelo, test = "Wilks")
summary(modelo, test = "Hotelling-Lawley")
summary(modelo, test = "Roy")
summary.aov(modelo)

# Como la interacciòn resultò significativa hacemos el anàlisis de efectos simples

# Modelo MANOVA
modelo_manova <- manova(cbind(Caroteno, Vitamina1) ~ Temperatura * Tiempo, data = datos)

# Resumen multivariado
summary(modelo_manova, test = "Pillai")

# Efectos simples de A dentro de cada nivel de B

# Para cada nivel de B, analizamos el efecto de A
levels(datos$Tiempo) %>% lapply(function(nivel_Tiempo) {
  cat("\n--- Análisis para Tiempo =", nivel_Tiempo, "---\n")
  subdatos <- filter(datos, Tiempo == nivel_Tiempo)
  summary(manova(cbind(Caroteno, Vitamina1) ~ Temperatura, data = subdatos), test = "Pillai")
})

# Efectos simples de B dentro de cada nivel de A

levels(datos$Temperatura) %>% lapply(function(nivel_Temperatura) {
  cat("\n--- Análisis para Temperatura =", nivel_Temperatura, "---\n")
  subdatos <- filter(datos, Temperatura == nivel_Temperatura)
  summary(manova(cbind(Caroteno, Vitamina1) ~ Tiempo, data = subdatos), test = "Pillai")
})

# Visualizaciòn opcional

# Promedios para graficar
promedios <- datos %>%
  group_by(Temperatura, Tiempo) %>%
  summarise(Caroteno = mean(Caroteno), Vitamina1 = mean(Vitamina1), .groups = "drop")

ggplot(promedios, aes(x = Tiempo, y = Caroteno, group = Temperatura, color = Temperatura)) +
  geom_line() + geom_point(size = 2) +
  labs(title = "Interacción Temperatura*Tiempo sobre Caroteno")

ggplot(promedios, aes(x = Tiempo, y = Vitamina1, group = Temperatura, color = Temperatura)) +
  geom_line() + geom_point(size = 2) +
  labs(title = "Interacción Temperatura*<tiempo sobre Vitamina1")

# Tabla de efectos simples

# Inicializamos listas para guardar resultados
efecto_Temperatura_en_Tiempo <- list()
efecto_Tiempo_en_Temperatura <- list()

# Efecto de Temperatura dentro de cada nivel de Tiempo
for (nivel_Tiempo in levels(datos$Tiempo)) {
  subdatos <- filter(datos, Tiempo == nivel_Tiempo)
  manova_mod <- manova(cbind(Caroteno, Vitamina1) ~ Temperatura, data = subdatos)
  pval <- summary(manova_mod, test = "Pillai")$stats["Temperatura", "Pr(>F)"]
  efecto_Temperatura_en_Tiempo[[nivel_Tiempo]] <- pval
}

# Efecto de Tiempo dentro de cada nivel de Temperatura
for (nivel_Temperatura in levels(datos$Temperatura)) {
  subdatos <- filter(datos, Temperatura == nivel_Temperatura)
  manova_mod <- manova(cbind(Caroteno, Vitamina1) ~ Tiempo, data = subdatos)
  pval <- summary(manova_mod, test = "Pillai")$stats["Tiempo", "Pr(>F)"]
  efecto_Tiempo_en_Temperatura[[nivel_Temperatura]] <- pval
}

# Convertir a data frames
tabla_efecto_Temperatura_en_Tiempo <- enframe(efecto_Temperatura_en_Tiempo, name = "Nivel de Tiempo", value = "p-valor efecto de Temperatura")
tabla_efecto_Tiempo_en_Temperatura <- enframe(efecto_Tiempo_en_Temperatura, name = "Nivel de Temperatura", value = "p-valor efecto de Tiempo")

# Convertir listas a vectores numéricos antes de enmarcar
tabla_efecto_Temperatura_en_Tiempo <- tibble(
  `Nivel de Tiempo` = names(efecto_Temperatura_en_Tiempo),
  `p-valor efecto de Temperatura` = unlist(efecto_Temperatura_en_Tiempo)
)

tabla_efecto_Tiempo_en_Temperatura <- tibble(
  `Nivel de Temperatura` = names(efecto_Tiempo_en_Temperatura),
  `p-valor efecto de Tiempo` = unlist(efecto_Tiempo_en_Temperatura)
)

# Mostrar
print(tabla_efecto_Temperatura_en_Tiempo)
print(tabla_efecto_Tiempo_en_Temperatura)

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(Caroteno, Vitamina1) ~ Temperatura*Tiempo, data = datos)

# Efecto simple de Temperatura dentro de cada nivel de Tiempo
emmTemperatura = emmeans(fit, ~ Temperatura | Tiempo | rep.meas)
mvcontrast(emmTemperatura, method = 'pairwise', mult.name = 'rep.meas', 
           by = 'Tiempo', adjust = 'holm')
mvcontrast(emmTemperatura, method = 'pairwise', mult.name = 'rep.meas', 
           by = 'Tiempo', adjust = 'none')

# Efecto simple de Tiempo dentro de cada nivel de Temperatura
emmTiempo = emmeans(fit, ~ Tiempo | Temperatura | rep.meas)
mvcontrast(emmTiempo, method = 'pairwise', mult.name = 'rep.meas', 
           by = 'Temperatura', adjust = 'bonferroni')
mvcontrast(emmTemperatura, method = 'pairwise', mult.name = 'rep.meas', 
           by = 'Temperatura', adjust = 'none')

# Si la interacción hubiese resultado NS se haría el análisis de efectos
# principales
# Comparar las cuatro combinaciones de tratamiento

emmTempTiem = emmeans(fit, ~ Temperatura*Tiempo | rep.meas)
mvcontrast(emmTempTiem, method = 'pairwise', mult.name = 'rep.meas',
           adjust = 'bonferroni')

##############################################
##      MANOVA CON FACTORIALES EN DBCA      ##
##############################################

# Ingreso de datos (REDUFACIL)

datos <- read.delim("clipboard",T)
tt(datos)
str(datos)

datos <- datos %>%
  mutate(across(c(Tipo_Dieta, Tipo_Ejercicio, Sexo), as.factor)) %>%
  glimpse()

# Supuesto de normalidad multivariada

# Desagregamos la base por tratamientos

trat1 = datos %>% filter(Tipo_Dieta == "Dieta1",Tipo_Ejercicio == "Ejercicio1") %>%
  dplyr::select(Perdida_Peso, Grasa_Corporal, Glucosa)
trat2 = datos %>% filter(Tipo_Dieta == "Dieta1",Tipo_Ejercicio == "Ejercicio2") %>%
  dplyr::select(Perdida_Peso, Grasa_Corporal, Glucosa)
trat3 = datos %>% filter(Tipo_Dieta == "Dieta2", Tipo_Ejercicio == "Ejercicio1") %>%
  dplyr::select(Perdida_Peso, Grasa_Corporal, Glucosa)
trat4 = datos %>% filter(Tipo_Dieta == "Dieta2", Tipo_Ejercicio == "Ejercicio2") %>%
  dplyr::select(Perdida_Peso, Grasa_Corporal, Glucosa)

#Ejecutamos el test

library(mvnormtest)
mshapiro.test(t(trat1))
mshapiro.test(t(trat2))
mshapiro.test(t(trat3))
mshapiro.test(t(trat4))

#No hay normalidad en 2 tratamientos entonces deberían
#plantearse transformaciones o usar las pruebas más robustas 
#descritas mas arriba.

# Supuesto de homogeneidad de matrices variancia covariancia

heplots::boxM(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ 
                Tipo_Dieta * Tipo_Ejercicio, data = datos)

# Supuesto de variables dependientes correlacionadas. Prueba de esfericidad de Bartlett

datos1=datos[3:5]
head(datos1)

options(scipen=0)
cortest.bartlett(cor(datos1), n = nrow(datos1))

#Determinacion de la matriz residual y las matrices factoriales 
#(Tipo_Dieta, Tipo_Ejercicio, Tipo_Dieta*Tipo_Ejercicio y Sexo) del MANOVA.

modelo = manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa)~
                Tipo_Dieta * Tipo_Ejercicio + Sexo, data = datos)
str(modelo)
coef(modelo)
Matrices = summary(modelo)$SS
Matrices
FDieta = Matrices$Tipo_Dieta
FEjercicio = Matrices$Tipo_Ejercicio
FDietaEjercicio = Matrices$`Tipo_Dieta:Tipo_Ejercicio`
FBloque = Matrices$Sexo
W = Matrices$Residuals

#Variabilidad explicada por el factor (Dieta). Matriz suma de cuadrados y 
#productos cruzados del factor (SCOCFDieta)
FDieta

#Variabilidad explicada por el factor Ejercicio. Matriz suma de cuadrados y 
#productos cruzados del factor Velocidad (SCOCEjercicio)
FEjercicio

#Variabilidad explicada por la interaccion Dieta*Ejercicio. Matriz suma de
#cuadrados y productos cruzados de la interaccion Dieta*Ejercicio 
#(SCOCLevadura*Velocidad)

FDietaEjercicio

#Variabilidad explicada por el Bloque. Matriz suma de cuadrados y productos cruzados
#del Bloque (SCOCBloque)
FBloque

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)
#del factor 
T = FDieta + FEjercicio + FDietaEjercicio + FBloque + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte de la variabilidad
#total puede atribuirse al factorial, mientras que un valor proximo a 0 significa 
#que el factor explica muy poco de esa variabilidad total.

eta2 = 1 - det(FBloque+W)/det(T)
eta2

# Pruebas de hipotesis del modelo

summary(modelo, test = "Pillai")
summary(modelo, test = "Wilks")
summary(modelo, test = "Hotelling-Lawley")
summary(modelo, test = "Roy")
summary.aov(modelo)

# Como la interacciòn resultò significativa hacemos el anàlisis de efectos
# simples

# Modelo MANOVA
modelo_manova <- manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ Tipo_Dieta * Tipo_Ejercicio + Sexo, data = datos)

# Resumen multivariado
summary(modelo_manova, test = "Pillai")

# Efectos simples de Tipo_Dieta dentro de cada nivel de Tipo_Ejercicio

library(dplyr)

# Para cada nivel de B, analizamos el efecto de A
levels(datos$Tipo_Ejercicio) %>% lapply(function(nivel_Tipo_Ejercicio) {
  cat("\n--- Análisis para Tipo_Ejercicio =", nivel_Tipo_Ejercicio, "---\n")
  subdatos <- filter(datos, Tipo_Ejercicio == nivel_Tipo_Ejercicio)
  summary(manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ Tipo_Dieta, data = subdatos), test = "Pillai")
})

# Efectos simples de Tipo_Ejercicio dentro de cada nivel de Tipo_Dieta

levels(datos$Tipo_Dieta) %>% lapply(function(nivel_Tipo_Dieta) {
  cat("\n--- Análisis para Tipo_Dieta =", nivel_Tipo_Dieta, "---\n")
  subdatos <- filter(datos, Tipo_Dieta == nivel_Tipo_Dieta)
  summary(manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ Tipo_Ejercicio, data = subdatos), test = "Pillai")
})

# Visualizaciòn opcional

# Promedios para graficar
promedios <- datos %>%
  group_by(Tipo_Dieta, Tipo_Ejercicio) %>%
  summarise(Perdida_Peso = mean(Perdida_Peso), Grasa_Corporal = mean(Grasa_Corporal), Glucosa = mean(Glucosa), .groups = "drop")

ggplot(promedios, aes(x = Tipo_Ejercicio, y = Perdida_Peso, group = Tipo_Dieta, color = Tipo_Dieta)) +
  geom_line() + geom_point(size = 2) +
  labs(title = "Interacción Tipo_Dieta*Tipo_ERjercicio sobre Perdida_Peso")

ggplot(promedios, aes(x = Tipo_Ejercicio, y = Grasa_Corporal, group = Tipo_Dieta, color = Tipo_Dieta)) +
  geom_line() + geom_point(size = 2) +
  labs(title = "Interacción Temperatura*<tiempo sobre Grasa corporal")

ggplot(promedios, aes(x = Tipo_Ejercicio, y = Glucosa, group = Tipo_Dieta, color = Tipo_Dieta)) +
  geom_line() + geom_point(size = 2) +
  labs(title = "Interacción Temperatura*<tiempo sobre Glucosa")

# Tabla de efectos simples

# Inicializamos listas para guardar resultados
efecto_Dieta_en_Ejercicio <- list()
efecto_Ejercicio_en_Dieta <- list()

# Efecto de la Dieta dentro de cada nivel de Ejercicio
for (nivel_Ejercicio in levels(datos$Tipo_Ejercicio)) {
  subdatos <- filter(datos, Tipo_Ejercicio == nivel_Ejercicio)
  manova_mod <- manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ Tipo_Dieta, data = subdatos)
  pval <- summary(manova_mod, test = "Pillai")$stats["Tipo_Dieta", "Pr(>F)"]
  efecto_Dieta_en_Ejercicio[[nivel_Ejercicio]] <- pval
}

# Efecto de Ejercicio dentro de cada nivel de Dieta
for (nivel_Dieta in levels(datos$Tipo_Dieta)) {
  subdatos <- filter(datos, Tipo_Dieta == nivel_Dieta)
  manova_mod <- manova(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ Tipo_Ejercicio, data = subdatos)
  pval <- summary(manova_mod, test = "Pillai")$stats["Tipo_Ejercicio", "Pr(>F)"]
  efecto_Ejercicio_en_Dieta[[nivel_Dieta]] <- pval
}

# Convertir a data frames
tabla_efecto_Dieta_en_Ejercicio <- enframe(efecto_Dieta_en_Ejercicio, name = "Tipo de Ejercicio", value = "p-valor efecto de dieta")
tabla_efecto_Ejercicio_en_Dieta <- enframe(efecto_Ejercicio_en_Dieta, name = "Tipo de Dieta", value = "p-valor efecto de Ejercicio")

# Convertir listas a vectores numéricos antes de enmarcar
tabla_efecto_Dieta_en_Ejercicio <- tibble(
  `Tipo de Ejercicio` = names(efecto_Dieta_en_Ejercicio),
  `p-valor efecto de Dieta` = unlist(efecto_Dieta_en_Ejercicio)
)

tabla_efecto_Ejercicio_en_Dieta <- tibble(
  `Tipo de Dieta` = names(efecto_Ejercicio_en_Dieta),
  `p-valor efecto de Ejercicio` = unlist(efecto_Ejercicio_en_Dieta)
)

# Mostrar
print(tabla_efecto_Dieta_en_Ejercicio)
print(tabla_efecto_Ejercicio_en_Dieta)

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(Perdida_Peso, Grasa_Corporal, Glucosa) ~ 
           Tipo_Dieta*Tipo_Ejercicio + Sexo, data = datos)

# Efecto simple de Tipo_Dieta dentro de cada nivel de Tipo_Ejercicio
emmDieta = emmeans(fit, ~ Tipo_Dieta | Tipo_Ejercicio | rep.meas)
mvcontrast(emmDieta, 'pairwise', mult.name = 'rep.meas', by = 'Tipo_Ejercicio',
           adjust = 'holm')
mvcontrast(emmDieta, 'pairwise', mult.name = 'rep.meas', by = 'Tipo_Ejercicio',
           adjust = 'none')

# Efecto simple de Tipo_Ejercicio dentro de cada nivel de Tipo_Dieta
emmEjercicio = emmeans(fit, ~ Tipo_Ejercicio | Tipo_Dieta | rep.meas)
mvcontrast(emmEjercicio, 'pairwise', mult.name = 'rep.meas', by = 'Tipo_Dieta',
           adjust = 'bonferroni')
mvcontrast(emmEjercicio, 'pairwise', mult.name = 'rep.meas', by = 'Tipo_Dieta',
           adjust = 'none')

# Si la interacción hubiese resultado NS se haría el análisis de efectos
# principales
# Comparar las cuatro combinaciones de tratamiento

emmDietEjer = emmeans(fit, ~ Tipo_Dieta*Tipo_Ejercicio | rep.meas)
mvcontrast(emmDietEjer, 'pairwise', mult.name = 'rep.meas',
           adjust = 'holm')

# Algunos resultados
library(car)
linearHypothesis(modelo, "Tipo_DietaDieta2:Tipo_EjercicioEjercicio2")
coef(modelo)
lh.out <- linearHypothesis(modelo, hypothesis.matrix =
                             c("Tipo_DietaDieta2 = 0", 
                               "Tipo_EjercicioEjercicio2 = 0"))
lh.out

#############################
##     MANCOVA EN DCA      ##
#############################

datos <-read.spss("1.Algebra.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
tt(datos)
datos

# Supuesto de normalidad multivariada
#La siguiente funcion (mshapiro.test) se aplica grupo a grupo, por lo tanto,
#primero es necesario dividir la base de datos en los grupos, dos en el ejemplo.

#desagregamos la base por grupos.

trat1 = datos %>% filter(Grupo == "Grupo1") %>%
  dplyr::select(Despues_Teoria ,Despues_Practica)
trat2 = datos %>% filter(Grupo == "Grupo2") %>%
  dplyr::select(Despues_Teoria, Despues_Practica)


#Ejecutamos el test

mshapiro.test(t(trat1))
mshapiro.test(t(trat2))

#En cada Grupo hay normalidad.


# Supuesto de homogeneidad de matrices variancia covariancia

res <- boxM(datos[, 4:5], datos[, "Grupo"])
res
summary(res)

heplots::boxM(cbind(Despues_Teoria, Despues_Practica) ~
                Grupo, data = datosc)

boxM(datos[, 4:5], datos[, 1])

algebra <- unique(datos$Grupo)
algebra1 <- lapply(algebra,
                function(x){as.matrix(datos[datos$Grupo == x, 4:5])}
)

names(algebra1) <- algebra
Ahmad2017(algebra1)

## Prueba Wrapper
homogeneityCovariances(datos, group = Grupo, covTest = BoxesM)

HOMOGENEITY(data = datos[c(1, 4, 5)],groups = 'Grupo', 
            variables = c('Despues_Teoria', 'Despues_Practica'))

# Supuesto de variables dependientes correlacionadas. 
# Prueba de esfericidad de Bartlett

datos1=datos[, -1:-3]
datos1
options(scipen=0)
cortest.bartlett(cor(datos1), n = nrow(datos1))

#Trabajando con el modelo de MANCOVA en DCA.
modelo = manova(cbind(Despues_Teoria, Despues_Practica) ~ 
                  Grupo + Antes_Teoria + Antes_Practica, data = datos)

#Determinacion de la matriz residual y la matriz factorial del MANCOVA.
str(modelo)
Matrices = summary(modelo)$SS
Matrices
F = Matrices$Grupo
X1 = Matrices$Antes_Teoria
X2 = Matrices$Antes_Practica
W = Matrices$Residuals

#Variabilidad explicada por el factor (grupo). Matriz suma de cuadrados y 
#productos cruzados del factor (SCOCF)
F

#Variabilidad de la covariable 1. Matriz suma de cuadrados y productos cruzados
#de la covariable 1 (SCOCX1)
X1

#Variabilidad de la covariable 2. Matriz suma de cuadrados y productos cruzados
#de la covariable 2 (SCOCX2)
X2

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)
 
T = F + X1 + X2 + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte de la variabilidad
#total puede atribuirse al factor, mientras que un valor proximo a 0 significa que 
#el factor explica muy poco de esa variabilidad total.

eta2 = 1 - det(X1 + X2 + W)/det(T)
eta2

# Pruebas de hipotesis del modelo

summary(modelo, test = "Pillai")
summary(modelo, test = "Wilks")
summary(modelo, test = "Hotelling-Lawley")
summary(modelo, test = "Roy")
#summary.aov(modelo)

#De otra manera
library(jmv)
modelo2 = mancova(data = datos,deps = vars(Despues_Teoria, Despues_Practica),
                factors = Grupo, covs = c(Antes_Teoria, Antes_Practica), boxM = T,
                shapiro = T, qqPlot = T)
modelo2

# Se tiene que hacer un MANOVA en DCA porque las Xs son no significativas
modelo1 = manova(cbind(Despues_Teoria, Despues_Practica) ~ 
                   Grupo, data = datos)
summary(modelo1, test = "Pillai")
summary(modelo1, test = "Wilks")
summary(modelo1, test = "Hotelling-Lawley")
summary(modelo1, test = "Roy")

#Comparar dos Grupos

modelo1 = manova(cbind(Despues_Teoria, Despues_Practica) ~ 
                   Grupo, data = datos,
               subset = Grupo %in% c("Grupo1", "Grupo2"))

summary(modelo1,test="Pillai")
summary(modelo1,test="Wilks")
summary(modelo1,test="Hotelling-Lawley")
summary(modelo1,test="Roy")
summary.aov(modelo1)

#---------------------------------------------------#
#           Comparación General por pares           #
#---------------------------------------------------#

# H0:los 2 vectores son iguales
# H1:los 2 vectores difieren

Grup<-c("Grupo1", "Grupo2")
comb<-t(combn(length(Grup), 2))

for(i in 1:nrow(comb)){
  modelo.comp = manova(cbind(Despues_Teoria, Despues_Practica) ~ 
                         Grupo, data = datos,
                     subset=Grupo %in% Grup[comb[i,]])
  print(paste("Trat: ",Grup[comb[i,]][1], "y",Grup[comb[i,]][2]))
  print(summary(modelo.comp,test = "Pillai"))
  cat("\n")
  
}

# NOTA: En este caso las dos covariables resultaron NS es por esa razón que se
# hizo un Manova en DCA. Ahora, sólo con fines de explicación, se asumirá que al
# menos una covariable resultó significativa y se seguirá haciendo el Mancova en
# DCA. Por lo tanto, para ver si hay diferencias entre grupos se usará la 
# librería emmeans.

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(Despues_Teoria, Despues_Practica) ~ Grupo + Antes_Teoria + 
           Antes_Practica, data = datos)
emm = emmeans(fit, ~ Grupo | rep.meas)
emm
names(datos)
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'none')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'bonferroni')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'holm')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'fdr')

# Asumiendo que el tratamiento Grupo1 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas')

# Asumiendo que el tratamiento Grupo2 es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas', 
           ref = 'Grupo2', adjust = 'bonferroni')

# mult.name toma por defecto rep.meas
mvcontrast(emm, method = 'consec')
mvcontrast(emm, method = 'consec', show.ests = TRUE)

#################################
##      MANCOVA EN DBCA        ##
#################################

datos <-read.spss("1.Marranas.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)

tt(datos)

# Supuesto de normalidad multivariada
#La siguiente funcion (mshapiro.test) se aplica grupo a grupo, por lo tanto
#primero es necesario dividir la base de datos en los grupos, tres en el ejemplo.

#desagregamos la base de datos por Racion.

trat1 = datos %>% filter(Racion == "18%") %>%
  dplyr::select(PESCAM, NUMLECH)
trat2 = datos %>% filter(Racion == "20%") %>%
  dplyr::select(PESCAM, NUMLECH)
trat3 = datos %>% filter(Racion == "22%") %>%
  dplyr::select(PESCAM, NUMLECH)

#Ejecutamos el test

mshapiro.test(t(trat1))
mshapiro.test(t(trat2))
mshapiro.test(t(trat3))

#Hay normalidad solo el el primer grupo.

# Supuesto de homogeneidad de matrices variancia covariancia

res <- boxM(datos[, 3:4], datos[, "Racion"])
res
summary(res)

heplots::boxM(cbind(PESCAM, NUMLECH) ~ Racion, data = datos)

boxM(datos[, 3:4], datos[, 5])

marranas <- unique(datos$Racion)
marranas1 <- lapply(marranas,
            function(x){as.matrix(datos[datos$Racion == x, 3:4])}
)

names(marranas1) <- marranas
Ahmad2017(marranas1)

## Prueba Wrapper
data=datos[, 3:5]
homogeneityCovariances(data, group = Racion, covTest = BoxesM)

HOMOGENEITY(data = datos[c(3,4,5)],groups = 'Racion', 
            variables = c('PESCAM', 'NUMLECH'))

# Supuesto de variables dependientes correlacionadas.
# Prueba de esfericidad de Bartlett

datos1=datos[, 3:4]
datos1
options(scipen=0)
cortest.bartlett(cor(datos1), n = nrow(datos1))
cor.test(datos$PESCAM, datos$NUMLECH)

#Trabajando con el modelo de MANCOVA.
modelo = manova(cbind(PESCAM, NUMLECH) ~ 
                  Racion + Repeticion + PESSERV + GRASDOR, data = datos)
modeloDCA = manova(cbind(PESCAM, NUMLECH) ~ 
                  Racion + PESSERV + GRASDOR, data = datos)
MatricesDCA = summary(modeloDCA)$SS
Wdca = MatricesDCA$Residuals

#Determinacion de la matriz residual y las matrices factoriales del MANCOVA.
str(modelo)
Matrices = summary(modelo)$SS
Matrices
F = Matrices$Racion
Bloque = Matrices$Repeticion
X1 = Matrices$PESSERV
X2 = Matrices$GRASDOR
W=Matrices$Residuals

#Variabilidad explicada por el factor (Racion). Matriz suma de cuadrados y productos cruzados
#del factor (SCOCF)
F

#Variabilidad explicada por el Bloque. Matriz suma de cuadrados y productos cruzados
#del Bloque (SCOCBloque)
Bloque

#Variabilidad de la covariable 1. Matriz suma de cuadrados y productos cruzados
#de la covariable 1 (SCOCX1)
X1

#Variabilidad de la covariable 2. Matriz suma de cuadrados y productos cruzados
#de la covariable 2 (SCOCX2)
X2

#Variabilidad residual. Matriz suma de cuadrados y productos cruzados
#del residual (SCOCR)
W

#Variabilidad Total. Matriz suma de cuadrados y productos cruzados total (SCOCT)

T = F + Bloque + X1 + X2 + W
T

#Bondad de ajuste. Un valor proximo a 1 indica que la mayor parte de la variabilidad
#total puede#atribuirse a la Racion, mientras que un valor proximo a 0 significa 
#que el factor explica muy poco de esa variabilidad total.

eta2 = 1 - det(Bloque + X1 + X2 + W)/det(T)
eta2

det(Wdca)/det(W) # El mancova en DCA es más eficiente

# Pruebas de hipotesis del modelo

summary(modelo, test = "Pillai")
summary(modelo, test = "Wilks")
summary(modelo, test = "Hotelling-Lawley")
summary(modelo, test = "Roy")
summary.aov(modelo)

# NOTA: En este caso al menos una covariable resultó S es por esa razón que se
# se compararán los niveles de Racion.

# Comparar dos raciones

modelo1 = manova(cbind(PESCAM, NUMLECH) ~ 
                   Racion + Repeticion + PESSERV + GRASDOR, data = datos,
               subset = Racion %in% c("20%", "22%"))

summary(modelo1, test = "Pillai")
summary(modelo1, test = "Wilks")
summary(modelo1, test = "Hotelling-Lawley")
summary(modelo1, test = "Roy")
summary.aov(modelo1)

#---------------------------------------------------#
#           Comparación General por pares           #
#---------------------------------------------------#

# H0:los 2 vectores son iguales
# H1:los 2 vectores difieren

Raci <- c("18%", "20%", "22%")
comb<-t(combn(length(Raci), 2))

for(i in 1:nrow(comb)){
  modelo.comp = manova(cbind(PESCAM, NUMLECH) ~ 
                         Racion + Repeticion + PESSERV + GRASDOR,
                       data=datos,
                     subset=Racion %in% Raci[comb[i,]])
  print(paste("Trat: ",Raci[comb[i,]][1], "y",Raci[comb[i,]][2]))
  print(summary(modelo.comp, test = "Pillai"))
  cat("\n")
  
}

#-------------------------------------------------------------------#
# Comparación General por pares usando la librería emmeans          #
# donde la función mvcontrast se basa en la distribución Hotelling  # 
#-------------------------------------------------------------------#

fit = lm(cbind(PESCAM, NUMLECH) ~ Racion + Repeticion + PESSERV + GRASDOR,
         data = datos)
emm = emmeans(fit, ~ Racion | rep.meas)
emm

mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'none')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'bonferroni')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'holm')
mvcontrast(emm, method = 'pairwise', mult.name = 'rep.meas', adjust = 'fdr')

# Asumiendo que el tratamiento 18% es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas')

# Asumiendo que el tratamiento 20% es el control
mvcontrast(emm, method = 'trt.vs.ctrl1', mult.name = 'rep.meas', ref = '20%', 
           adjust = 'bonferroni')

# mult.name toma por defecto rep.meas
mvcontrast(emm, method = 'consec')
mvcontrast(emm, method = 'consec', show.ests = TRUE)

