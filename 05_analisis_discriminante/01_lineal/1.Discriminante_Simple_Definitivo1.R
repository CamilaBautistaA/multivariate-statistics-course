##################################################
#                                                #
# ANÁLISIS DISCRIMINANTE LINEAL SIMPLE           #
#          EJEMPLO SUSCRIPCIÓN                   #
# Clodomiro Fernando Miranda Villagómez          #   
#       cfmiranda@lamolina.edu.pe                #
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
library(klaR)
library(gains)
library(caret)
library(Boruta) 
library(rsample)
library(janitor)
library(gmodels)
library(GGally)
library(reshape2)
library(heplots)
library(MVN)
library(dplyr)
library(MASS)
library(tidyverse)
library(gtsummary)
library(DataExplorer)
library(patchwork)

# El ADLS se usa para modelar respuestas binarias como el diagnóstico de una
# enfermedad (positivo o negativo), votar (a favor o en contra) o la admisión
# a una universidad (admitido o no).

#######################
# 1. Lectura de datos #
#######################

datosd <- read.delim("clipboard",T)
head(datosd)
dim(datosd)
str(datosd)

datosd$Suscripcion=factor(datosd$Suscripcion)
str(datosd)

# Para ver la categoria de referencia (positivo)
contrasts(datosd$Suscripcion)  
# Para cambiar la categoria de referencia (no). Positivo (si)
datosd$Suscripcion = relevel(datosd$Suscripcion,ref="si")

#------------------------------------------------------------------#
#           Analisis descriptivo y exploratorio                    #
#------------------------------------------------------------------#

tbl_summary(datosd,
            by=Suscripcion,
            statistic = all_continuous() ~ '{mean} ({sd})') %>%
  add_ci()

tbl_summary(datosd,
            by=Suscripcion)

plot_histogram(datosd)

plot_density(datosd)

plot_qq(datosd)

plot_qq(datosd,by='Suscripcion')
plot_qq(datosd,by='Suscripcion',ncol = 2,ggtheme = theme_dark())

plot_boxplot(datosd,by='Suscripcion')

plot_scatterplot(datosd,by='Edad')

plot_prcomp(datosd) # Componentes principales


# Analisis descriptivo
# --------------------
describeBy(datosd[,-6],datosd$Suscripcion)

# Correlaciones
# -------------
ggscatmat(data = datosd,color="Suscripcion",alpha = 0.4) +
  theme_bw() +
  labs(title = "Correlacion por pares") +
  theme( plot.title = element_text(size = 16, face = "bold"),
         axis.text = element_blank(),
         strip.text = element_text(colour = "black", size = 8, face = 2))

ggpairs(datosd,mapping = ggplot2::aes(color = Suscripcion),alpha = 0.4, 
        lower = list(continuous = "smooth"), 
        diag = list(continuous = "bar"), 
        axisLabels = "none")

# Boxplot
# -------------
df.m <- melt(datosd, id.var = "Suscripcion")
p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Suscripcion))+ 
  facet_wrap( ~ variable, scales="free") 
p

###########################
# Seleccion de Variables  #
###########################

# Criterio de lambda de Wilks para seleccion de variables
library(klaR)
greedy.wilks(Suscripcion ~ .,data=datosd)

# Seleccion de variables con la libreria Boruta
library(Boruta) 
set.seed(123)
boruta <- Boruta(Suscripcion~., data = datosd, doTrace = 2)
print(boruta)

plot(boruta,cex.axis=0.5)

#Trabajando con las mejores variables segun lambda de Wilks

head(datosd)
datosd=datosd[,-c(4,5)]
head(datosd)

# Pruebas individuales
# --------------------

# Dividir la base de datos en datos: para el entrenamiento (E) y 
# para la prueba (P). Con la base de entrenamiento se crea el modelo y con la
# base de prueba para validar el modelo (y no ser juez y parte)

set.seed(10)
suscrip_split = initial_split(datosd, prop = 0.70, strata = Suscripcion)

suscrip_train = training(suscrip_split)
suscrip_test = testing(suscrip_split)

# se mantienen los porcentajes
tabyl(datosd$Suscripcion)
tabyl(suscrip_train$Suscripcion)
tabyl(suscrip_test$Suscripcion)

dim(suscrip_test)
dim(suscrip_train)
dim(datosd)[1]
dim(suscrip_test)[1]/dim(datosd)[1]
dim(suscrip_train)[1]/dim(datosd)[1]

prop.table(table(datosd$Suscripcion))
prop.table(table(suscrip_test$Suscripcion))
prop.table(table(suscrip_train$Suscripcion))

# Trabajando con la data de entrenamiento para crear el modelo

######################
##       MODELO     ##
######################
table(suscrip_train$Suscripcion)
modelo<-lda(Suscripcion ~ Educacion+Edad+Tvdiario ,prior=c(230/309,79/309),data=suscrip_train)  
print(modelo)


#------------------------------------------------------------------#
#                         Supuestos basicos                        #
#------------------------------------------------------------------#

# Prueba de homogeneidad de matrices de variancias-covariancias
#--------------------------------------------------------------
G1 <- subset(suscrip_train,Suscripcion=="no")
head(G1)
dim(G1)
G2 <- subset(suscrip_train,Suscripcion=="si")
dim(G2)
var(G1[,-4]) # Matriz variancia covariancia muestral de no
var(G2[,-4]) # Matriz variancia covariancia muestral de si
library(biotools)
head(suscrip_train)
boxM(suscrip_train[, -4], suscrip_train[,4])

boxM(suscrip_train[,c(1,2,3)],suscrip_train[,4])

# Prueba de Normalidad multivariada
#------------------------------------

library(mvnormtest)
mshapiro.test(t(G1[, -4]))
mshapiro.test(t(G2[, -4]))

library(MVN)
royston_test <- mvn(data = G1[,-4], mvn_test = "royston")
royston_test
royston_test <- mvn(data = G2[,-4], mvn_test = "royston")
royston_test

#Deteccion de outliers multivariados
head(datosd)
# Distancia de Mahalanobis
outliers_quan <- MVN::mvn(
  data = as.data.frame(datosd[,-4]),
  multivariate_outlier_method = "quan",
  alpha = 0.05
)

outliers_quan
summary(outliers_quan, 'outliers')

# Distancia ajustada de Mahalanobis
outliers_adj <- MVN::mvn(
  data = as.data.frame(datosd[,-4]),
  multivariate_outlier_method = "adj",
  alpha = 0.05
)

summary(outliers_adj, 'outliers')

# Distancia de Mahalanobis

centro <- colMeans(datosd[,-4])
S <- cov(datosd[,-4])

dist_mahal <- mahalanobis(datosd[,-4], center = centro, cov = S)
(corte_mahal <- qchisq(0.975, df = ncol(datosd[,-4])))

tabla_mahal <- data.frame(
  ID = 1:nrow(datosd[,-4]),
  Distancia_Mahalanobis = dist_mahal,
  Atipico = dist_mahal > corte_mahal
)

ggplot(tabla_mahal, aes(x = ID, y = Distancia_Mahalanobis)) +
  geom_point(size = 3) +
  geom_hline(yintercept = corte_mahal, linetype = "dashed") +
  labs(
    title = "Distancias de Mahalanobis de los datos",
    subtitle = "La linea horizontal representa el punto de corte exploratorio",
    x = "Vinos (ID)",
    y = expression(D^2)
  ) +
  theme_minimal()
#Supuesto del poder discriminante de cada variable

mod.dca = lm(suscrip_train$Educacion ~ suscrip_train$Suscripcion) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza
head(suscrip_train)
mod.dca = lm(suscrip_train$Edad ~ suscrip_train$Suscripcion) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza

mod.dca = lm(suscrip_train$Tvdiario ~ suscrip_train$Suscripcion) # modelo DCA
summary(aov(mod.dca)) # Analisis de varianza

#------------------------------------------------------------------#
#       Estimacion de la Funcion Discriminante Lineal              #
#------------------------------------------------------------------#

# Prueba global
# --------------

#Significacion de la funcion discriminante
#La funcion discriminante es estadisticamente significativa (clasifica bien)
dim(suscrip_train)
modelo.adj <- manova(data=suscrip_train, cbind(Educacion,Edad,Tvdiario) ~ Suscripcion)
summary(modelo.adj,test="Wilks")

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

wilk.ind(Suscripcion ~ Educacion+Edad+Tvdiario, data= suscrip_train)

# Las 3 variables en forma individual contribuyen a formar la funcion
# discriminante del modelo
greedy.wilks(Suscripcion ~ Educacion+Edad+Tvdiario, data= suscrip_train, niveau = 0.1)
?greedy.wilks

# importancia de las variables

library(randomForest)
rf = randomForest(Suscripcion ~ Educacion+Edad+Tvdiario, data= suscrip_train)

library(vip)
vip(rf)

#Significacion de la funcion discriminante
#La funcion discriminante es estadisticamente significativa (clasifica bien)
modelo.adj1 <- manova(data=suscrip_train, cbind(Educacion,Edad,Tvdiario) ~ Suscripcion)
summary(modelo.adj1,test="Wilks")

#Autovalor
summary((modelo.adj1),test="Wilks")$Eigenvalues
lambda=0.5861708
eta2=lambda/(1+lambda)
eta2
#El 36.96% de la variabilidad total es explicada por la 
#funcion discriminante 1.

#Correlacion canonica
eta = sqrt(eta2)
eta
#Eta=0.6079.Como se aproxima a 1 entonces la FD 
# tiene gran capacidad explicativa.

#GL Residuales=307
#Calculo de los coeficientes estandarizados y de la MATRIZ de ESTRUCTURA

summary((modelo.adj1),test="Wilks")$SS
std.b.Educacion=(sqrt(summary(modelo.adj1)$SS$Residuals[1,1]/307))*modelo$scaling[1,1]
std.b.Edad=(sqrt(summary(modelo.adj1)$SS$Residuals[2,2]/307))*modelo$scaling[2,1]
std.b.Tvdiario=(sqrt(summary(modelo.adj1)$SS$Residuals[3,3]/307))*modelo$scaling[3,1]
std.b.Educacion
modelo$scaling

SCPC.residual=summary((modelo.adj1),test="Wilks")$SS$Residuals
SCPC.residual
SCPC.residual.varianzas=SCPC.residual/307
SCPC.residual.correlaciones=cov2cor(SCPC.residual.varianzas)
SCPC.residual.correlaciones

#Matriz de Estructura
l.Educacion=SCPC.residual.correlaciones[1,1]*std.b.Educacion+
  SCPC.residual.correlaciones[1,2]*std.b.Edad+
  SCPC.residual.correlaciones[1,3]*std.b.Tvdiario
l.Edad=SCPC.residual.correlaciones[2,1]*std.b.Educacion+
  SCPC.residual.correlaciones[2,2]*std.b.Edad+
  SCPC.residual.correlaciones[2,3]*std.b.Tvdiario
l.Tvdiario=SCPC.residual.correlaciones[1,3]*std.b.Educacion+
  SCPC.residual.correlaciones[2,3]*std.b.Edad+
  SCPC.residual.correlaciones[3,3]*std.b.Tvdiario

#MATRIZ DE ESTRUCTURA (la edad es la mas correlacionada
# a la funcion discriminante: -0.6098455)
rbind(l.Educacion,l.Edad,l.Tvdiario)

# MATRIZ DE ESTRUCTURA. De otra forma

head(suscrip_train)
scores <- predict(modelo)$x # Scores discriminantes
M_Estructura <- cor(suscrip_train[, -4], scores)
M_Estructura

#------------------------------------------------------------------#
#                    Prediccion de nuevos casos                    #
#------------------------------------------------------------------#
head(suscrip_train)
# Educacion   Edad    Tvdiario
#     12       25         4

nuevo<-data.frame(Educacion=12,Edad=25,Tvdiario=4)
predict(modelo,nuevo)$class
predict(modelo,nuevo)$posterior

datosn<-read.delim("clipboard")
clase=predict(modelo,datosn)$class
proba=predict(modelo,datosn)$posterior
datonsnp=cbind(datosn,clase,proba)
datonsnp

############################################################
#  Clase y Probabilidad predicha de lOs datos reales       #
############################################################

# Centroides
# ----------
D1 <- sum(modelo$means[1,]*modelo$scaling)
D2 <- sum(modelo$means[2,]*modelo$scaling)

# Punto de corte discriminante
# ----------------------------
C=(D1+D2)/2
C

# Prediccion para la observacion 10
# --------------------------------
head(suscrip_train,11)
(D <- sum(modelo$scaling*suscrip_train[10,-4]))
D-C    #Valor positivo. Se le clasifica como "no"

clase.pred <-predict(modelo,suscrip_train[10,-4])$class 
clase.pred
clase.pred <-predict(modelo,suscrip_train[10,-4])$x 
str(clase.pred)
clase.pred

# Valores estimados con la ecuacion discriminante
#--------------------------
fd.pred <- predict(modelo,suscrip_train[,-4])$x
head(fd.pred, 40)
tail(fd.pred, 40)
# Clase predicha
# ---------------
clase.pred <-predict(modelo,suscrip_train[,-4])$class 
clase.pred[50:100]
suscrip_train$Suscripcion[50:100]

# Probabilidad predicha
#---------------------------
proba.pred  <-round(predict(modelo,suscrip_train[,-4])$posterior,5)
head(proba.pred,30)

# Almacenamiento de datos con clase y probabilidad predicha
#-------------------------------------------------------------
datosf=cbind(suscrip_train,fd.pred,clase.pred,proba.pred)
head(datosf,30)
write.csv(datosf,"suscripcion_Expo1.csv")

# Grafico de la densidad e histograma
#---------------------------------
ggplot(datosf, aes(x = proba.pred[, 1], fill = Suscripcion)) + 
  geom_density(alpha = 0.6) 

plot(modelo,dimen = 1)

z=predict(modelo,suscrip_train[-4])$x
z=predict(modelo)$x
ldahist(z,suscrip_train[,4])

library(ggplot2)
ggplot(datosf, aes(x = proba.pred[,1], fill = Suscripcion)) + 
  geom_histogram(alpha = 0.6)

library(RColorBrewer)
RColorBrewer::display.brewer.all()

library(ggplot2)
a=ggplot(datosf) +
  aes(x = proba.pred[,1], fill = Suscripcion) +
  geom_histogram(bins = 30L) +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal()
a

b=ggplot(datosf, aes(x = proba.pred[,1], fill = Suscripcion)) + 
  geom_density(alpha = 0.6)+
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal()
b
a + b

# Grafico de particion todas las variables
#-----------------------------------------

# Si se consideraran las variables por pares, los rojos son mal clasificados
partimat(Suscripcion~.,data=suscrip_train,method="lda",nplots.vert=2,
         image.colors = c("darkgoldenrod1", "skyblue2"))

#----------------------------------------------------#
#         Matriz o Tabla de Confusión                #
#           Para validar el modelo                   #
#----------------------------------------------------#

# Se llama matriz de confusión porque el modelo puede hacer predicciones 
# correctas (TP y TN) pero el modelo se puede confundir (FP y FN).
# Error tipo I: cuando el modelo predice incorrectamente la clase positiva (FP),
# el resultado real es negativo.
# Error tipo II: cuando el modelo predice incorrectamente la clase negativa (FN),
# el resultado real es positivo.
head(suscrip_test)
d2 = suscrip_test %>% 
  mutate(pred_probs = predict(modelo, suscrip_test[,-4])$posterior[,1]) %>%
  mutate(pred_classes = ifelse(pred_probs > 0.5, 'si', 'no')) %>%
  mutate(Suscripcion = factor(Suscripcion, levels = c('si', 'no')),
         pred_classes = factor(pred_classes, levels = c('si', 'no')))

# Tabla o matriz de confusión (esta es la forma apropiada, convencional)

cm2 = xtabs(data = d2, ~ pred_classes + Suscripcion)
confusionMatrix(cm2, mode = 'everything', positive = 'si')


# Trabajando con la data de prueba

clase.pred1 <-predict(modelo,suscrip_test[,-4])$class
proba.pred1  <-round(predict(modelo,suscrip_test[,-4])$posterior,5)

# Tabla de Confusión (de otra manera)

clase.pred2 <-predict(modelo,suscrip_test[,-4])$class 
data1 = suscrip_test %>% mutate(clase.pred2)
head(data1)
cm1 = xtabs(data = data1, ~ clase.pred2 + Suscripcion)
confusionMatrix(cm1, mode = 'everything', positive = 'si')

TP = 91; FN = 8; FP = 12; TN = 22; Total = 133
TP + FN + FP + TN
table(datosd$Suscripcion)
113/442
# Accuracy ----------------------
# Mide la proporción de todas las predicciones correctas, independientemente de
# la clase (Ver VPpos y VPneg). Con datos desbalanceados (desequilibrados) es 
# mejor usar el Balanced Accuracy (BA = (sensibilidad+ especificidad)/2)
accuracy <- mean(suscrip_test$Suscripcion==clase.pred1)
accuracy
Accuracy <- (TP+TN)/Total
Accuracy
# Calcular el error de mala clasificacion (Tasa de error aparente)
#----------------------------------------
error <- mean(suscrip_test$Suscripcion!=clase.pred1)
error
error <- (FP+FN)/Total
error

# Calculo de la sensibilidad o Recall = TP/(FN+TP)
# Probabilidad de que un si, sea bien identificado
# Tasa de verdaderos positivos identificados correctamente por el modelo
# El modelo identificó el 91.92% de los positivos reales
(sensibilidad = TP/(FN+TP))

# Calculo de la especificidad=TN/(TN+FP)
# Tasa de verdaderos negativos identificados correctamente por el modelo
# Probabilidad de que un no, sea bien identificado
# Probabilidad de rechazo correcto
# El modelo identificó el 64.71% de los negativos reales
(especificidad = TN/(TN+FP))

# Balanced Accuracy
# Es una medida más robusta que el Accuracy (Ambas clases (si y no) reciben
# igual consideración). En nuestro ejm. el BA = 0.7831 es sólo un 6.7% inferior que el 
# Accuracy. Esto sugiere que las clases están relativamente balanceadas.
(BA = (sensibilidad + especificidad)/2)
.8496-.7831

# Índice J de Youden
# Maximiza la efectividad general del modelo
# Mide cuánto mas informado está el modelo respecto a resultados aleatorios
# Maximiza los verdaderos positivos y minimiza los falsos positivos
# J = 0 indica que la prueba o modelo no funciona mejor que el azar (no tiene información útil)
# J = 1 indica que la prueba o modelo distingue perfectamente entre positivos y negativos
# J nos dice que tan bien el modelo mejora la toma de decisiones sobre el azar
# J = 0.5663 es aceptable pero no fuerte
(J = sensibilidad + especificidad - 1)
library(cutpointr)
youden(tp = TP, fp = FP, tn = TN, fn = FN)
library(epiR)
summary(epi.tests(cm1))

# Con la función cutpointr se puede equilibrar mejor la sensibilidad, 
# la especificidad o la J de Youden. Esa parte no se verá en este curso.

# Calculo de la Precision = valor predictivo(+) = TP/(FP+TP)
# Mide que tan precisas son las predicciones positivas
# Probabilidad de que un resultado positivo del modelo sea verdadero con los
# los siguientes clientes de la revista
# Porcentaje de resultados positivos que son correctos.
# Un VPpos = 1 indica que todas las predicciones positivas fueron correctas.
# no hay FP. Un VPpos = 0 indica que las predicciones positivas fueron incorrectas.
# En nuestro caso: VPpos = 0.8835 quiere decir que el 88.35% de los clientes
# marcados por el modelo como suscritos realmente se suscribieron. 
# se suscribieron 
(VPpos = TP/(FP+TP))

# Calculo del valor predictivo(-) = TN/(TN+FN)
# Porcentaje de resultados negativos que son correctos.
# Un VPneg = 0.7333 quiere decir que el 73.33% de clientes marcados por el modelo
# como no suscritos realmente no se suscribieron,
(VPneg = TN/(TN+FN))

# Prevalencia: Proporción de casos positivos reales. Es la proporción real de 
# clientes que se suscribieron
(Prevalencia = (FN+TP)/Total)

# Detection Rate (DR): mide la proporción de positivos que el modelo identificó
# correctamente
(DR = TP/Total)

# Detection Prevalence (DP) o Apparent Prevalence: La prevalencia de detección 
# es la proporción de casos positivos predichos por el modelo. El modelo predice
# que el 77.44% se suscribirán.
(DP = (TP+FP)/Total)

# Prueba de Mc Nemar. 
# H0: No existe variacion entre los valores reales y predichos
# de suscripcion.

tabla=table(Clase_Predicha=clase.pred1,Real=suscrip_test$Suscripcion)
tabla
mcnemar.test(tabla)

# Coeficiente Kappa de Cohen
# ¿que magnitud de kappa refleja un acuerdo adecuado? Las directrices
# serian utiles, pero otros factores distintos del acuerdo pueden
# influir en su magnitud, lo que hace problematica la interpretacion
# de una determinada magnitud
### Prueba de Kolmogorov Smirnov
# Esta prueba no es recomendable cuando hay empates con una muestra. 
# Aqui se hace la prueba con dos muestras


datosf1=cbind(suscrip_test,clase.pred1,proba.pred1)

head(datosf1)
datos_kappa=datosf1[,c(4,5)]
head(datos_kappa,15)
library(psych)
cohen.kappa(datos_kappa, w=NULL,n.obs=NULL,alpha=.05,levels=NULL)

## H0: Las distribuciones de las suscripciones observadas y predichas 
## no difieren

head(datosf1)
str(datosf1)
library(tidyverse)
datos_ks=datosf1[,c(4,5)]
head(datos_ks)
str(datos_ks)
datos_ks1=mutate(datos_ks,Suscripcion=as.numeric(Suscripcion),
                 clase.pred1=as.numeric(clase.pred1))
head(datos_ks1)
datos_ks2=mutate(datos_ks1,Suscripcion1=Suscripcion-1,
                 clase.pred2=clase.pred1-1)
head(datos_ks2)
ks.test(x=datos_ks1$Suscripcion, y=datos_ks1$clase.pred1)
ks.test(x=datos_ks2$Suscripcion1, y=datos_ks2$clase.pred2)


#-----------------------------------
# Usando la libreria caret
library(caret)
cm = confusionMatrix(clase.pred1,suscrip_test$Suscripcion,positive="si")
(TP+FN)/133  # Prevalencia
(TP+FP)/133 # Detection Prevalence

# No Information Rate (NIR) (Tasa de falta de información): es el Accuracy que se obtendría
# si siempre predijera la clase mayoritaria. Con esta data la clase mayoritaria es 'si'
# (91+8= 99). El NIR ayuda a verificar si el Accuracy del modelo (0.8496) es mejor que el del
# modelo base sin información (0.7444). Para verificar si la diferencia (0.8496-0.7444=-0.1052) 
# es significativa se halla el p-value [Acc > NIR] con una prueba binomial
0.8496-0.7444
table(suscrip_test$Suscripcion)
(NIR = 99/133)

# H0: El accuracy del modelo <= NIR
# H1: El accuracy del modelo > NIR
binom.test(
  x = 113,     # Número de predicciones positivas
  n = 133,     # Número de ensayos
  p = 0.7443609,  # Probabilidad del suceso esperado (en las hipótesis)
  alternative = 'greater')

# Matriz de confusión. Otra presentación

# Un alto Recall (o sensibilidad) puede reducir la Precision (Valor predictivo +)
# Por eso a menudo se usa el F1 Score, porque balancea ambos (es la media armònica)
(F1 = 2 * sensibilidad * VPpos / (sensibilidad + VPpos))

confusionMatrix(clase.pred1,suscrip_test$Suscripcion,positive="si",
                mode = 'everything')

# Matriz de confusión. Otra presentación

bin_cm=ConfusionTableR::binary_class_cm(clase.pred1,suscrip_test$Suscripcion)
bin_cm$record_level_cm
library(dplyr)
glimpse(bin_cm$record_level_cm)
ConfusionTableR::binary_visualiseR(train_labels = clase.pred1,
                                   truth_labels = suscrip_test$Suscripcion,
                                   class_label1 = "Si",
                                   class_label2 = "No",
                                   quadrant_col1 = "violet",
                                   quadrant_col2 = "tomato",
                                   custom_title = "Matriz de Confusion-Suscripcion",
                                   text_col = "black",
                                   cm_stat_size = 1.2,positive="si")


library(caret)
confusionMatrix(clase.pred1,suscrip_test$Suscripcion,positive="si")

#####################################
#  Curva ROC y Area bajo la curva   #
#####################################

library(ROSE)
roc.curve(suscrip_test$Suscripcion, proba.pred1[,1],lty=2,lwd=1.8,col="blue" ,
          main="ROC curves")


library(caTools)
head(datosf1)
colAUC(proba.pred1[,1],datosf1$Suscripcion,plotROC = TRUE)
abline(0, 1,col="gold",lty=7)

library(Epi)
ROC(proba.pred1[,2],datosf1$Suscripcion,lwd = 3)
head(proba.pred1[,2])

ROC(data = datosf1, form = Suscripcion ~ si)

2*0.8611111-1 #Coeficiente de Gini
329/442
#####################################################################
#   OBSERVACIÓN: Validación de un modelo con datos desbalanceados  ##
#####################################################################

## ¡Deje de utilizar el Accuracy! Cuatro métricas esenciales para datos 
## desequilibrados

# El Accuracy, la Precission (VPpos) y el Recall (sensibilidad) fallan en medir 
# el éxito de los modelos de aprendizaje automático  en escenarios de alto riesgo
# como la predicción de tornados, inundaciones, detección de fraudes, diagnóstico
# de enfermedades raras etc. en los cuales los datos no están balanceados. En esas
# situaciones se puede usar las métricas:

# 1. F1 Score = Puntaje F1. Es el promedio armónico del VPpos y la sensibilidad

# Equilibra la Precission y el Recall, por lo tanto es mejor que ambos. F1 Score 
# y el Accuracy miden la performance del modelo, pero F1 es más eficaz con datos 
# no balanceados

# 2. Critical Success Index (CSI) = Índice Crítico de Éxito 
# (o Threat Score (TS) = Puntaje de Amenaza)

# Se usa para predecir eventos raros. Un CSI = 0 indica que el modelo no hizo 
# predicciones correctas del evento, un CSI = 1 quiere decir que el modelo 
# hizo predicciones perfectas. Si el CSI es de 0.8 o más indica un gran 
# rendimiento del modelo, si el CSI es a lo más 0.2 el modelo tiene un 
# rendimiento malo. Para valores de CSI mayores de 0.2 pero menores de 0.8 se
# considera que el modelo tiene un rendimiento aceptable.
# El CSI depende del evento raro que se quiere predecir, por ejemplo, si un 
# evento raro E1 es más "fácil" de predecir que otro E2 se puede esperar un CSI
# entre 0.6-0.8, mientras que para E2 entre 0.3-0.5, para ser aceptables.
# El CSI ignora los TN (ver la fórmula más abajo). El CSI funciona mejor que el F1
# con datos extremadamente desbalanceados. Con eventos muy raros el CSI puede hacer
# predicciones correctas (TP) lo que hace que el CSI sea más grande de lo necesario,
# lo que dificulta la comparación del CSI de este estudio con otros eventos raros
# como huracanes vs inundaciones. Para solucionar este inconveniente es mejor usar
# el ETS.

# 3. Equitable Threat Score (ETS) = Puntaje de Amenaza Equitativa

# El ETS tampoco considera los TN (ver la fórmula más abajo).
# El ETS elimina los efectos aleatorios (Hc: impactos debidos al azar o 
# TP_ebc (True Positives expected by chance)) tanto en
# el numerador como en el denominador del CSI, ver la fórmula más abajo. Esto hace que 
# el ETS sea más preciso y justo para comparar diferentes eventos raros.

# 4. Matthews Corr. Coeff. (MCC) = Coeficiente de Correlación de Matthews
# El MCC, a diferencia del CSI y del ETS, se calcula con los 4 valores de la
# matriz de confusión (incluyendo los TN). El MCC es más confiable que el
# Accuracy especialmente con datos desbalanceados. Da una mejor visión del
# modelo tanto con datos balanceados o no balanceados.
# La MCC mide la correlación entre las clasificaciones observadas y predichas.
# Se interpreta como una correlación

# HIPOCALCEMIA EN BOVINOS

library(readxl)
library(tinytable)
library(janitor)
library(caret)
library(MASS)
library(tidyverse)

datos <- read_excel("Dataset_Hipocalcemia_4.xlsx")
tt(head(datos))
tabyl(datos, Hipocalcemia)
datos$Hipocalcemia = as.factor(datos$Hipocalcemia)
contrasts(datos$Hipocalcemia)
datos$Hipocalcemia = relevel(datos$Hipocalcemia,ref="si")


# Entrenamiento y Prueba (70/30)

set.seed(84)
idx_train <- createDataPartition(datos$Hipocalcemia, p = 0.70, list = FALSE)

train_data <- datos[ idx_train, ]
test_data  <- datos[-idx_train, ]

# Estimación del Modelo LDA (Con ajuste a priori)

prop_train <- prop.table(table(train_data$Hipocalcemia))

table(train_data$Hipocalcemia)
modelo <- lda(Hipocalcemia ~ ., data = train_data, prior = c(70, 630)/700)
print(modelo)

# Clasificación y probabilidades

d2 = test_data %>% 
  mutate(pred_probs = predict(modelo, test_data[,-7])$posterior[,1]) %>%
  mutate(pred_classes = ifelse(pred_probs > 0.5, 'si', 'no')) %>%
  mutate(Hipocalcemia = factor(Hipocalcemia, levels = c('si', 'no')),
         pred_classes = factor(pred_classes, levels = c('si', 'no')))

tt(head(d2))

# Tabla o matriz de confusión (esta es la forma apropiada, convencional)

cm2 = xtabs(data = d2, ~ pred_classes + Hipocalcemia)
confusionMatrix(cm2, mode = 'everything', positive = 'si')

# Evaluación del Modelo (Métricas Robustas para Desbalance)

# Cálculo de F1, PRECISSION y SENSITIVITY

TP = 9; FN = 21; FP = 4; TN = 266
Recall = 0.30 # o sensibilidad

(Precision = TP / (TP + FP)) # o valor predictivo +

(F1 = 2 * Precision * Recall / (Precision + Recall))

(Sensitivity = TPR = TP / (TP + FN)) # o Recall

# Critical Success Index (CSI) = Índice Crítico de Éxito
# Cálculo de CSI o TS. El CSI ignora los TN

csi =function(TP, FP, FN) {
  TP / (TP + FP + FN)
}

# csi = 0.2647, esto indica que el modelo no hizo predicciones correctas del evento

csi(TP, FP, FN) 
# CSI (%) & ETS (%) & CSI - ETS (%)

# Equitable Threat Score (ETS) = Puntaje de Amenaza Equitativa
# Cálculo del ETS y CSI
csi_ets = function(TP, FP, FN, total) {
  
# cálculo de TP_ebc (True Positives expected by chance) o Hc (impactos debidos al azar)
TP_ebc = ((TP + FP)  * (TP + FN)) / total
  
# cálculo de ETS
ETS = (TP -TP_ebc) / (TP + FP + FN - TP_ebc)
  
# cálculo del CSI
CSI = (TP) / (TP + FP + FN)
  
# resultados
results = list(
    'CSI (%)' = CSI * 100,
    'ETS (%)' = ETS * 100,
    'CSI - ETS (%)' = (CSI - ETS) * 100) |>
    as_tibble() |>
    mutate_all( ~round(., 1))
  return(results)
}

csi_ets(TP, FP, FN, total = 300)

# Coeficiente de correlación de Matthews (MCC)

tt(head(test_data))
MCC = function(TP, TN, FP, FN) {
  numerador = (TP * TN) - (FP * FN)
  denominador = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  numerador / denominador
}

# Hay una correlación negativa entre las clasificaciones observadas y predichas
MCC(TP, TN, FP, FN)

# Coeficiente de correlación de Mathews (MCC): mide la correlación 
# entre las clases verdaderas y las clases predichas
library(mltools)

clase.pred1 <-predict(modelo, test_data[,-7])$class

mcc(preds = clase.pred1, actuals = test_data$Hipocalcemia)

#Capacidad discriminante del modelo

tt(head(d2))
library(ROSE)

roc.curve(d2$Hipocalcemia, d2$pred_probs, lty=4, lwd=2, col="blue2" ,
          main="ROC curves")

