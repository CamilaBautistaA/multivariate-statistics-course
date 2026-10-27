
#EJEMPLO DE LAS EMPRESAS

library(psych)
library(MASS)
library(klaR)
library(gains)
library(caret)
library(Boruta) 
library(gmodels)
library(GGally)
library(reshape2)
library(heplots)
library(MVN)
library(dplyr)


# Lectura de datos con 5 variables
library(foreign)
datosd <-read.spss("Empresas.sav",
                   use.value.labels=TRUE, 
                   to.data.frame=TRUE)
str(datosd)
head(datosd)

# Para ver la categoria de referencia
contrasts(datosd$resultad)  

str(datosd)
attach(datosd)

#########################################################
# 2. Analisis Descriptivo de las variables predictoras  #
#########################################################

library(psych)
summary(datosd)
summarytools::descr(datosd)

tapply(inc_prom,resultad,mean)
tapply(inc_publ,resultad,mean)
tapply(patrocin,resultad,mean)
tapply(promocio,resultad,mean)
tapply(publicid,resultad,mean)

###########################
# Seleccion de Variables  #
###########################

# Criterio de lambda de Wilks para seleccion de variables
library(klaR)
greedy.wilks(resultad ~ .,data=datosd)

#TRABAJANDO SEGUN EL CRITERIO DE LAMDA DE WILKS

datosd=datosd[,-2:-3]
head(datosd)

# Supuesto de homogeneidad de matrices variancia covariancia

library(biotools)
boxM(datosd[,-1],datosd[,1])

#Supuesto del poder discriminante de dada variable

mod.dca = lm(patrocin ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # An?lisis de varianza

mod.dca = lm(promocio ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # An?lisis de varianza

mod.dca = lm(publicid ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # An?lisis de varianza

# Supuesto de normalidad multivariada de cada grupo

#Desagregamos la base por grupos, transponemos y efectuamos la prueba.
datosd
grupo1=datosd[datosd$resultad=="Bajos",2:4]
grupo2=datosd[datosd$resultad=="Medios",2:4]
grupo3=datosd[datosd$resultad=="Altos",2:4]

grupo1=t(grupo1)
grupo2=t(grupo2)
grupo3=t(grupo3)

library(mvnormtest)
mshapiro.test(grupo1)
mshapiro.test(grupo2)
mshapiro.test(grupo3)

#########################################################
#  Estimacion de la Funcion Discriminante Regularizada  #
#########################################################
?rda
head(datosd)
library(klaR)
modelo<-rda(resultad ~.,datosd,gamma=0.15,lambda=1)
str(modelo)
print(modelo)

#Significacion conjunta de las funciones discriminantes
#Las funcines discriminantes son estadisticamente significativas
head(datosd)
fit.manova=manova(data=datosd,cbind(patrocin,promocio,publicid)~resultad)
summary((fit.manova),test="Wilks")
#GL Residuales=46

#Autovalores
summary((fit.manova),test="Wilks")$Eigenvalues
lambda1=4.005055
lambda2=0.2005152
eta2.1=lambda1/(1+lambda1)
eta2.2=lambda2/(1+lambda2)
eta2.1
eta2.2

#El 80% de la variacion total es explicada por la funcion discriminante 1.
#El 16.7% de la variacion total es explicada por la funcion discriminante 2.


#Correlacion canonica
eta1=eta2.1^0.5
eta2=eta2.2^0.5
eta1
eta2
#Eta1=0.895.Como se aproxima a 1 entonces la FD1 tiene gran capacidad explicativa.
#Eta2=0.409.Como no se aproxima a 1 entonces la FD2 no tiene gran capacidad explicativa.

#GL Residuales=46

#########################################################
#  Clasificacion y Probabilidades de las observaciones  #
#########################################################

#--------------------------
# Clase predicha
head(datosd)
clase.pred <-predict(modelo,datosd[,-1])$class
head(clase.pred)


#---------------------------
# Probabilidad predicha
proba.pred  <-predict(modelo,datosd[,-1])$posterior
head(proba.pred)
str(proba.pred)
predict(modelo,datosd)


#-------------------------------------------------------------
# Almacenamiento de datos con clase y probabilidad predecida
datosf=cbind(datosd,clase.pred,proba.pred)
head(datosf)
str(datosf)
write.csv(datosf,"resultados-adl-scores.csv")

library(ggplot2)
proba.pred1 <- proba.pred[,1]
ggplot(datosf, aes(x = proba.pred1, fill = resultad)) + 
  geom_histogram(alpha = 0.25)
ggplot(datosf, aes(x = proba.pred1, fill = resultad)) + 
  geom_density(alpha = 0.25)


############################################
#  Indicadores para Evaluacion de Modelos  #
############################################

#############################
# 1. Tabla de clasificacion #
#############################

library(gmodels)
# Tabla de clasificacion
CrossTable(x = resultad, y = clase.pred,
           prop.t=FALSE, prop.c=FALSE, prop.chisq = FALSE)

addmargins(table(Real=resultad,Clase_Predicha=clase.pred))
prop.table(table(Real=resultad,Clase_Predicha=clase.pred),1)

# Calcular el accuracy
accuracy <- mean(resultad==clase.pred) ; accuracy
(18+16+14)/49
# Calcular el error de mala clasificacion
error <- mean(resultad!=clase.pred) ; error


#-----------------------------------
# Usando la libreria caret
library(caret)
confusionMatrix(clase.pred,resultad,positive="Si")


#####################################
# 2. Curva ROC y Area bajo la curva #
#####################################

head(proba.pred)
proba.pred1 <- proba.pred[,1] # Probabilidad de prediccion Bajos
proba.pred2 <- proba.pred[,2] # Probabilidad de prediccion Medios
proba.pred3 <- proba.pred[,3] # Probabilidad de prediccion Altos
head(proba.pred1)
head(proba.pred2)
head(proba.pred3)

library(caTools)

# Considerando probabilidades de prediccion bajos
colAUC(proba.pred1,resultad,plotROC = T)
abline(0, 1,col="red",lty=3)
colAUC(proba.pred1,resultad,plotROC = F)

# Considerando probabilidades de prediccion medios
colAUC(proba.pred2,resultad,plotROC = T)
abline(0, 1,col="red",lty=3)
colAUC(proba.pred2,resultad,plotROC = F)

# Considerando probabilidades de prediccion altos
colAUC(proba.pred3,resultad,plotROC = T)
abline(0, 1,col="red",lty=3)
colAUC(proba.pred3,resultad,plotROC = F)

#######################################
#  Prediccion para nuevos individuos  #
#######################################

# publicid promocio patrocin
#   243      153      21.3

nuevo1<-data.frame(publicid=243,promocio=153,patrocin=21.3)
predict(modelo,nuevo1)$class
predict(modelo,nuevo1)$posterior

nuevo2<-data.frame(publicid=236,promocio=152,patrocin=20)
predict(modelo,nuevo2)$class
predict(modelo,nuevo2)$posterior


