
#EJEMPLO DE VICUÑAS

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

# Lectura de datos con 6 variables
library(foreign)
datosd <-read.spss("Vicunas.sav",
                   use.value.labels=TRUE, 
                   to.data.frame=TRUE)
#datosd<-read.spss(file.choose(),header=T,to.data.frame = T)
str(datosd)
head(datosd)
summary(datosd)

# No considerar las columnas 1, 3 y 4.
datosd <- datosd[,-c(1,3,4)]
head(datosd,40)
attach(datosd)

# Para ver la categoria de referencia (positivo)
contrasts(datosd$Zona)  
# Para cambiar la categoria de referencia (positivo)
# datosd$Zona = relevel(datosd$Zona,ref="Puna Seca")

##############################################################################
# 2. Analisis exploratorio con algunas medidas de las variables predictoras  #
##############################################################################

library(psych)
summary(datosd)
summarytools::descr(datosd)

tapply(Hemoglobina,Zona,mean)
tapply(Proteina,Zona,mean)
tapply(Albumina,Zona,mean)
tapply(Creatinina,Zona,mean)
tapply(Colesterol,Zona,mean)

##################################
#  Analisis Exploratorio gráfico #
##################################

# Correlaciones
# -------------
ggscatmat(data = datosd,color="Zona",alpha = 0.4) +
  theme_bw() +
  labs(title = "Correlaci?n por pares") +
  theme( plot.title = element_text(size = 16, face = "bold"),
         axis.text = element_blank(),
         strip.text = element_text(colour = "black", size = 8, face = 2))

# Boxplot
# -------------
df.m <- melt(datosd, id.var = "Zona")
p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Zona))+ 
  facet_wrap( ~ variable, scales="free") 
p

###########################
# Seleccion de Variables  #
###########################

# Criterio de lambda de Wilks para seleccion de variables
library(klaR)
greedy.wilks(Zona ~ .,data=datosd)

# Seleccion de variables con la libreria Boruta
library(Boruta) 
set.seed(123)
boruta <- Boruta(Zona~., data = datosd, doTrace = 2)
print(boruta)

plot(boruta,cex.axis=0.5)

#TRABAJANDO SEGUN EL CRITERIO DE LAMDA DE WILKS
head(datosd)
datosd=datosd[,c(1,4,5)]
head(datosd)

# Supuesto de homogeneidad de matrices variancia covariancia

library(biotools)
boxM(datosd[,-1],datosd[,1])
#Resulta que las matrices variancia covariancia no son homogeneas.
#Se debe usa el analisis discriminante cuadratico.

#Supuesto del poder discriminante de cada variable.

mod.dca = lm(Albumina ~ Zona,datosd) # modelo DCA
summary(aov(mod.dca)) # An?lisis de varianza

mod.dca = lm(Creatinina ~ Zona,datosd) # modelo DCA
summary(aov(mod.dca)) # An?lisis de varianza

# Supuesto de normalidad multivariada de cada grupo

#Desagregamos la base por grupos, transponemos y efectuamos la prueba.
datosd
grupo1=datosd[datosd$Zona=="Puna Humeda",2:3]
grupo1
grupo2=datosd[datosd$Zona=="Puna Seca",2:3]
grupo2

grupo1=t(grupo1)
grupo1
grupo2=t(grupo2)

library(mvnormtest)
mshapiro.test(grupo1)
mshapiro.test(grupo2)

#######################################################
#  Estimacion de la Funcion Discriminante Cuadratica  #
#######################################################

head(datosd)
library(MASS)
table(datosd$Zona)
modelo<-qda(Zona ~.,datosd,prior=c(36,37)/73)
str(modelo)
print(modelo)

#Significacion conjunta de las funciones discriminantes
#Las funciones discriminantes son estadisticamente significativas
head(datosd)
fit.manova=manova(data=datosd,cbind(Albumina,Creatinina)~Zona)
summary((fit.manova),test="Wilks")
#GL Residuales=71

#Autovalores
summary((fit.manova),test="Wilks")$Eigenvalues
lambda1=1.950961
eta2.1=lambda1/(1+lambda1)
eta2.1

#El 66.1% de la variacion total es explicada por la funcion discriminante.

#Correlacion canonica
eta1=eta2.1^0.5
eta1

#Eta1=0.8130974.Como se aproxima a 1 entonces la FD tiene gran capacidad
#explicativa.

#GL Residuales=71

summary((fit.manova),test="Wilks")$SS
modelo$scaling

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
head(datosf,20)
datosf
str(datosf)
write.csv(datosf,"resultados-adl-scores.csv")

library(ggplot2)
proba.pred1 <- proba.pred[,2]
ggplot(datosf, aes(x = proba.pred1, fill = Zona)) + geom_histogram(alpha = 0.25)
ggplot(datosf, aes(x = proba.pred1, fill = Zona)) + geom_density(alpha = 0.25)


############################################
#  Indicadores para Evaluacion de Modelos  #
############################################

#############################
# 1. Tabla de clasificacion #
#############################

library(gmodels)
# Tabla de clasificaci?n
CrossTable(x = Zona, y = clase.pred,
           prop.t=FALSE, prop.c=FALSE, prop.chisq = FALSE)

addmargins(table(Real=Zona,Clase_Predicha=clase.pred))
prop.table(table(Real=Zona,Clase_Predicha=clase.pred),1)

# Calcular el accuracy
accuracy <- mean(Zona==clase.pred) ; accuracy
(34+35)/73
# Calcular el error de mala clasificacion
error <- mean(Zona!=clase.pred) ; error

#-----------------------------------
# Usando la libreria caret
library(caret)
confusionMatrix(clase.pred,Zona,positive="Puna Seca")


#####################################
# 2. Curva ROC y Area bajo la curva #
#####################################

# Primera forma

library(ROSE)
roc.curve(datosd$Zona, proba.pred[,2],lty=2,lwd=1.8,col="blue" ,
          main="ROC curves")

library(caTools)
colAUC(proba.pred,Zona,plotROC = TRUE)
abline(0, 1,col="green2",lty=3)

#######################################
#  Prediccion para nuevos individuos  #
#######################################

nuevo1<-data.frame(Albumina=4.4,Creatinina=1.2)
predict(modelo,nuevo1)$class
predict(modelo,nuevo1)$posterior

nuevo2<-data.frame(Albumina=6.7,Creatinina=1.9)
predict(modelo,nuevo2)$class
predict(modelo,nuevo2)$posterior


