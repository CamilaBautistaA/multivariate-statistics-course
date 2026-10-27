##################################################
#                                                #
# ANALISIS DISCRIMINANTE LINEAL MULTIPLE         #
#      EJEMPLO EMPRESAS                          #
#     Clodomiro Fernando Miranda Villagomez      #   
#   cfmiranda@lamolina.edu.pe                    #
#                                                #
##################################################
#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o informacion cargada
rm(list = ls())

##############
#  Paquetes  #
##############

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

######################
#  Lectura de Datos  #
######################

# Lectura de datos con 5 variables independientes
library(foreign)
datosd <-read.spss("2.Discriminante_Multiple_Empresas.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
str(datosd)
head(datosd)
summary(datosd)
attach(datosd)

#########################################################
# 2. Analisis Descriptivo de las variables predictoras  #
#########################################################

library(psych)
summary(datosd)

tapply(inc_prom,resultad,mean)
tapply(inc_publ,resultad,mean)
tapply(patrocin,resultad,mean)
tapply(promocio,resultad,mean)
tapply(publicid,resultad,mean)

## Caras de Chernov (promedio y desviacion estandar)
head(datosd)

library(tidyverse)

prev1 <- datosd %>%
  group_by(resultad) %>%
  summarise(incproprom = mean(inc_prom), publprom = mean(inc_publ),  
            patrocinprom = mean(patrocin), promocioprom = mean(promocio),
            publicidprom = mean(publicid))
prev1

prev2 <- datosd %>%
  group_by(resultad) %>%
  summarise(incprodesv = sd(inc_prom), publdesv = sd(inc_publ),  
            patrocindesv = sd(patrocin), promociodesv = sd(promocio),
            publiciddesv = sd(publicid))
prev2

library(aplpack)
windows()
faces(prev1[, 2:6], labels = prev1$resultad)  #Media
faces(prev2[, 2:6], labels = prev2$resultad)  #Desv

###########################
#  Analisis Exploratorio  #
###########################

# Correlaciones
# -------------
ggscatmat(data = datosd,color="resultad",alpha = 0.4) +
  theme_bw() +
  labs(title = "Correlacion por pares") +
  theme( plot.title = element_text(size = 16, face = "bold"),
         axis.text = element_blank(),
         strip.text = element_text(colour = "black", size = 8, face = 2))

library(dplyr)
require(GGally)
ggpairs(datosd,mapping = ggplot2::aes(color = resultad),alpha = 0.4, 
        lower = list(continuous = "smooth"), 
        diag = list(continuous = "bar"), 
        axisLabels = "none")

# Boxplot
# -------------
df.m <- melt(datosd, id.var = "resultad")
p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=resultad))+ 
  facet_wrap( ~ variable, scales="free") 
p

###########################
# Seleccion de Variables  #
###########################

# Criterio de lambda de Wilks para seleccion de variables
library(klaR)
greedy.wilks(resultad ~ .,data=datosd)

# Seleccion de variables con la libreria Boruta
library(Boruta) 
set.seed(123)
boruta <- Boruta(resultad~., data = datosd, doTrace = 2)
print(boruta)

plot(boruta,cex.axis=0.5)

#TRABAJANDO SEGUN EL CRITERIO DE LAMDA DE WILKS
head(datosd)
datosd=datosd[,-2:-3]
head(datosd)

# Supuesto de homogeneidad de matrices variancia covariancia

bajos<-subset(datosd,datosd$resultad=="Bajos")
medios<-subset(datosd,datosd$resultad=="Medios")
altos<-subset(datosd,datosd$resultad=="Altos")
bajos=bajos[,-1]
medios=medios[,-1]
altos=altos[,-1]

var(bajos) # Matriz variancia covariancia muestral de bajos
var(medios) # Matriz variancia covariancia muestral de medios
var(altos) # Matriz variancia covariancia muestral de altos

library(biotools)
boxM(datosd[,-1],datosd[,1])

#Supuesto del poder discriminante de cada variable

mod.dca = lm(patrocin ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza

mod.dca = lm(promocio ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza

mod.dca = lm(publicid ~ resultad,datosd) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza

# Supuesto de normalidad multivariada de cada grupo

#Desagregamos la base por grupos, transponemos y efectuamos la prueba.
head(datosd)
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
#  Estimacion de las Funciones Discriminantes Lineales  #
#########################################################
?lda
head(datosd)
table(datosd[,1])
dim(datosd)
library(biotools)
modelo<-lda(resultad ~.,datosd,prior=c(18,17,14)/49)
str(modelo)
print(modelo)

# Para ver la categoria de referencia (positivo)
contrasts(datosd$resultad)  
# Para cambiar la categoria de referencia (positivo)
# datosd$resultad = relevel(datosd$resultad,ref="Medios")

# Pruebas individuales
# --------------------
library(dplyr)
wilk.ind <- function(modelo, data){
  vars <- all.vars(modelo)
  dif <- list()
  for(i in 2:length(vars)){
    dif[[i-1]] <- greedy.wilks(as.formula(paste0(vars[1], "~", vars[i])),
                               data=data,niveau = 1)$results
  }
  wilk <- bind_rows(dif, .id = "column_label")[,1:3]
  wilk
}

wilk.ind(resultad ~.,datosd)

# Las 3 variables en forma individual contribuyen a formar las funciones
# discriminantes del modelo
greedy.wilks(resultad ~., data= datosd, niveau = 0.1)

# importancia de las variables

library(randomForest)
rf = randomForest(resultad ~., data= datosd)

library(vip)
vip(rf)


#Significacion conjunta de las funciones discriminantes
#Las funciones discriminantes son estadisticamente significativas 
#(discriminan bien)
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
#Eta1=0.895.Como se aproxima a 1 entonces la FD1 tiene gran 
#capacidad explicativa.
#Eta2=0.409.Como no se aproxima a 1 entonces la FD2 no tiene gran
#capacidad explicativa.

#GL Residuales=46
#Calculo de los coeficientes estandarizados y de la MATRIZ de ESTRUCTURA

summary((fit.manova),test="Wilks")$SS
modelo$scaling
std.b.patrocinLD1=(sqrt(summary(fit.manova)$SS$Residuals[1,1]/46))*modelo$scaling[1,1]
std.b.promocioLD1=(sqrt(summary(fit.manova)$SS$Residuals[2,2]/46))*modelo$scaling[2,1]
std.b.publicidLD1=(sqrt(summary(fit.manova)$SS$Residuals[3,3]/46))*modelo$scaling[3,1]
std.b.publicidLD1

std.b.patrocinLD2=(sqrt(summary(fit.manova)$SS$Residuals[1,1]/46))*modelo$scaling[1,2]
std.b.promocioLD2=(sqrt(summary(fit.manova)$SS$Residuals[2,2]/46))*modelo$scaling[2,2]
std.b.publicidLD2=(sqrt(summary(fit.manova)$SS$Residuals[3,3]/46))*modelo$scaling[3,2]
std.b.publicidLD2

SCPC.residual=summary((fit.manova),test="Wilks")$SS$Residuals
SCPC.residual
SCPC.residual.varianzas=SCPC.residual/46
SCPC.residual.correlaciones=cov2cor(SCPC.residual.varianzas)
SCPC.residual.correlaciones

#Matriz de Estructura
l.patrocinLD1=SCPC.residual.correlaciones[1,1]*std.b.patrocinLD1+SCPC.residual.correlaciones[1,2]*std.b.promocioLD1+
  SCPC.residual.correlaciones[1,3]*std.b.publicidLD1
l.promocioLD1=SCPC.residual.correlaciones[2,1]*std.b.patrocinLD1+SCPC.residual.correlaciones[2,2]*std.b.promocioLD1+
  SCPC.residual.correlaciones[2,3]*std.b.publicidLD1
l.publicidLD1=SCPC.residual.correlaciones[3,1]*std.b.patrocinLD1+SCPC.residual.correlaciones[3,2]*std.b.promocioLD1+
  SCPC.residual.correlaciones[3,3]*std.b.publicidLD1

l.patrocinLD2=SCPC.residual.correlaciones[1,1]*std.b.patrocinLD2+SCPC.residual.correlaciones[1,2]*std.b.promocioLD2+
  SCPC.residual.correlaciones[1,3]*std.b.publicidLD2
l.promocioLD2=SCPC.residual.correlaciones[2,1]*std.b.patrocinLD2+SCPC.residual.correlaciones[2,2]*std.b.promocioLD2+
  SCPC.residual.correlaciones[2,3]*std.b.publicidLD2
l.publicidLD2=SCPC.residual.correlaciones[3,1]*std.b.patrocinLD2+SCPC.residual.correlaciones[3,2]*std.b.promocioLD2+
  SCPC.residual.correlaciones[3,3]*std.b.publicidLD2

#MATRIZ DE ESTRUCTURA
l.patrocinLD1
l.promocioLD1
l.publicidLD1

l.patrocinLD2
l.promocioLD2
l.publicidLD2

variables=c("patrocinio","promocion","publicidad")
LD1=c(l.patrocinLD1,l.promocioLD1,l.publicidLD1)
LD2=c(l.patrocinLD2,l.promocioLD2,l.publicidLD2)
data.frame(variables,LD1,LD2)

# MATRIZ DE ESTRUCTURA. De otra forma

head(datosd)
scores <- predict(modelo)$x # Scores discriminantes
M_Estructura <- cor(datosd[, -1], scores)
M_Estructura

#########################################################
#  Clasificacion y Probabilidades de las observaciones  #
#########################################################

#--------------------------
# Clase predicha
head(datosd)
clase.pred <-predict(modelo,datosd[,-1])$class
head(clase.pred, 20)


#---------------------------
# Probabilidad predicha
proba.pred  <-round(predict(modelo,datosd[,-1])$posterior,4)
head(proba.pred)
str(proba.pred)


#-------------------------------------------------------------
# Almacenamiento de datos con clase y probabilidad predecida
datosf=cbind(datosd,clase.pred,proba.pred)
head(datosf)
str(datosf)
write.csv(datosf,"resultados-adl-scores.csv")

plot(modelo,dimen = 1)
plot(modelo,dimen = 2)

z=predict(modelo,datosd[-1])$x
z=predict(modelo)$x
ldahist(z,datosd[,1])

library(ggplot2)
proba.pred2 <- proba.pred[,1]

ggplot(datosf, aes(x = proba.pred2, fill = clase.pred)) + 
  geom_histogram()

ggplot(datosf, aes(x = proba.pred2, fill = clase.pred)) + 
  geom_density()



# Grafico de particion todas las variables
#-----------------------------------------

# Si se consideraran las variables por pares, los rojos son mal clasificados
partimat(resultad~.,data=datosd,method="lda",nplots.vert=2,
         image.colors = c("darkgoldenrod1", "skyblue2","green2"))

############################################
#  Indicadores para Evaluacion de Modelos  #
############################################

#############################
# 1. Tabla de clasificacion #
#############################
dim(datosd)
library(gmodels)
# Tabla de clasificacion
CrossTable(x = clase.pred, y = resultad,
           prop.t=FALSE, prop.c=FALSE, prop.chisq = FALSE)

addmargins(table(Clase_Predicha=clase.pred,Real=resultad))
prop.table(table(Clase_Predicha=clase.pred,Real=resultad),1)

# Calcular el accuracy o tasa de buena clasificacion
accuracy <- mean(resultad==clase.pred) ; accuracy
(18+15+13)/49
#La probabilidad de que una empresa sea
#bien clasificada es 0.939
# Calcular la tasa de error o tasa de mala clasificacion
error <- mean(resultad!=clase.pred) ; error
3/49
#-----------------------------------
# Usando la libreria caret
library(caret)
cm=confusionMatrix(clase.pred,resultad)
print(cm)

# Implementing function to collapse data
library(ConfusionTableR)
mc_df <- ConfusionTableR::multi_class_cm(clase.pred, resultad)
# Access the reduced data for storage in databases
print(mc_df$record_level_cm)
glimpse(mc_df$record_level_cm)

mc_df$confusion_matrix
mc_df$cm_tbl

# Coeficiente Kappa de Cohen
# ?que magnitud de kappa refleja un acuerdo adecuado? Las directrices
# serian utiles, pero otros factores distintos del acuerdo pueden influir
# en su magnitud, lo que hace problematica la interpretacion de una
# determinada magnitud

head(datosf)
datos_kappa=datosf[,c(1,5)]
head(datos_kappa)
library(psych)
cohen.kappa(datos_kappa, w=NULL,n.obs=NULL,alpha=.05,levels=NULL)

# Con BAJOS vs OTROS
# Calculo de la sensibilidad=VP/(FN+VP)
18/(0+18)
# Calculo de la especificidad=VN/(VN+FP)
28/31
# Caculo de la precision=valor predictivo(+)=VP/(FP+VP)
18/21
# Calculo del valor predictivo(-)=VN/(VN+FN)

# Prevalencia=(VP+FN)/Total


# Con MEDIOS vs OTROS
# Calculo de la sensibilidad=VP/(FN+VP)

# Calculo de la especificidad=VN/(VN+FP)

# Calculo de la precision=valor predictivo(+)=VP/(FP+VP)

# Calculo del valor predictivo(-)=VN/(VN+FN)

# Prevalencia=(VP+FN)/Total


# Con ALTOS vs OTROS
# Calculo de la sensibilidad=VP/(FN+VP)

# Calculo de la especificidad=VN/(VN+FP)

# Calculo de la precision=valor predictivo(+)=VP/(FP+VP)

# Calculo del valor predictivo(-)=VN/(VN+FN)

# Prevalencia=(VP+FN)/Total

# Coeficiente de correlación de Mattews multiclase (MCC multiclase): mide la correlación 
# entre las clases verdaderas y las clases predichas

library(yardstick)
head(datosf)

mcc(datosf,
    truth = resultad,
    estimate = clase.pred)

# F1-score macro: Es especialmente útil cuando las clases están desbalanceadas
# porque otorga el mismo peso a cada clase, independientemente de cuántas 
# observaciones tenga.
# Es el promedio de k F1-score binarios
# Mide que tan buena es la clasificación (debe ser >0.6)

library(yardstick)

f_meas_vec(
  truth = datosd$resultad,
  estimate = clase.pred,
  estimator = "macro"
)

f_meas(
  datosf,
  truth = resultad,
  estimate = clase.pred,
  estimator = "macro"
)

#####################################
# 2. Curva ROC y Area bajo la curva #
#####################################

# Área bajo la curva multiclase: mide la capacidad del modelo de distinguir
# correctamente entre las tres categorías.
# Puede ser el promedio de los AUC de cada categoría vs otros, o el promedio de
# los AUC de cada categoría vs cada una de las otras categorías.

library(pROC)
head(datosd)
head(proba.pred)
multiclass.roc(datosd$resultad, proba.pred)

# Primera forma

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

# Segunda forma
head(datosd)
library(Epi)
ROC(form=resultad~patrocin+promocio+publicid,plot="ROC")

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

(n=table(datosd[,1]))
(VMEG=apply(datosd[,-1],2,function(x) tapply(x,datosd[,1],mean))) #medias de grupos
(Sp=Reduce("+",by(datosd[,2:4],datosd[,1],function(x) (nrow(x)-1)*cov(x)))/(sum(n)-length(n)))

mahalanobis(VMEG,c(21.3,153,243),Sp)#clasifica como Bajos (menor distancia=6.293933)
mahalanobis(VMEG,c(20,152,236),Sp)#clasifica como Altoss (menor distancia=2.081771)

#Para ver la clasificacion de todas las empresas
(tmp=apply(datosd[,-1],1,function(x) mahalanobis(VMEG,as.numeric(x),Sp)))

prd.empresas=apply(tmp,2,which.min)
adl.empresas=predict(modelo)$class
#La siguiente tabla muestra un gran de acuerdo
table(prd.empresas,adl.empresas)

#La siguiente tabla muestra errores de clasificacion. Esta tabla coincide con
#las que se vieron anteriormente
table("Empresa predicha"=adl.empresas,"Empresa actual"=datosd[,1])













