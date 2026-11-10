################################################
#                                              #
#  ANALISIS DE REGRESION LOGISTICA Multinomial #
#      Clodomiro Fernando Miranda Villagomez   #   
#    cfmiranda@lamolina.edu.pe                 #
#                                              #
################################################


#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o informacion cargada
rm(list = ls())


#---------------------------------------------------------
#  Librerias

library(nnet) 
library(MASS) 
library(brant)
library(tidyverse)
library(gtsummary)
library(patchwork)
library(sjPlot)
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

###############################
# CASO 1: Peso de Ninos       #
###############################

# Se quiere estudiar el peso de recien nacidos (bajo, normal y alto)
# en funcion de 6 variables independientes.


#-------------------------------------------------
# Lectura de datos
library(foreign)
datos <-read.spss("PesoNinos1.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
head(datos)
datosd <- datos[,-1]
head(datosd)
dim(datosd)
str(datosd)

#View(datosd)

attach(datosd)

round(tapply(Edad_Madre,Peso_RN_1,mean),3)

boxplot(Edad_Madre ~ Peso_RN_1, 
        main= "BoxPlot de Edad de la Madre vs Peso",
        xlab = "Peso", ylab="Edad de la Madre",
        col = c("red","blue","green"))

view_df(datosd, show.frq=T, show.prc = T, show.na = T)

datosd %>% 
  plot_frq(Peso_RN_1)

datosd %>% 
  group_by(Peso_RN_1) %>%
  plot_frq(Tipo_de_parto) %>%
  plot_grid()

tbl_summary(datosd,
            by=Peso_RN_1,
            statistic = all_continuous() ~ '{mean} ({sd})')

tbl_summary(datosd,
            by=Peso_RN_1)

# Como la variable dependiente es categórica ordinal, veamos si corresponde
# hacer una regresión logística ordinal o una regresión logística multinomial

# Modelo de regresion logistica ordinal
model <- polr(Peso_RN_1~.,data=datosd)

# H0: Las regresiones son paralelas.
# Como el pvalor de "Omnibus" es cero se debe hacer una regresión multinomial
brant(model)

#-----------------------------------------------------#
# 2. Estimacion del Modelo Logistico Multinomial      #
#-----------------------------------------------------#

modelo_ls <- multinom(Peso_RN_1~.,data=datosd)

# Para ver la categoria de referencia
contrasts(datosd$Peso_RN_1)
contrasts(datosd$Tipo_de_parto)
contrasts(datosd$Sexo_RN)

# Para cambiar la categoria de referencia
datosd$Peso_RN_1 = relevel(datosd$Peso_RN_1,ref="Peso Normal")
datosd$Sexo_RN = relevel(datosd$Sexo_RN,ref="Masculino")
datosd$Tipo_de_parto = relevel(datosd$Tipo_de_parto,ref="Parto  Normal")

contrasts(datosd$Peso_RN_1)  
contrasts(datosd$Sexo_RN)
contrasts(datosd$Tipo_de_parto)

#  Seleccion de Variables  

step <- stepAIC(modelo_ls,direction="forward", trace=FALSE)
step$anova # Seguimos con las seis variables independientes
?stepAIC

# Dividir la base de datos en datos: para el entrenamiento (E) y 
# para la prueba (P)

head(datosd)

set.seed(10)
peso_split = initial_split(datosd, prop = 0.70, strata = Peso_RN_1)

peso_train = training(peso_split)
peso_test = testing(peso_split)

# se mantienen los porcentajes
tabyl(datosd$Peso_RN_1)
tabyl(peso_train$Peso_RN_1)
tabyl(peso_test$Peso_RN_1)

dim(peso_test)
dim(peso_train)
dim(datosd)[1]
dim(peso_test)[1]/dim(datosd)[1]
dim(peso_train)[1]/dim(datosd)[1]

prop.table(table(datosd$Peso_RN_1))
prop.table(table(peso_test$Peso_RN_1))
prop.table(table(peso_train$Peso_RN_1))

# Trabajando con la data de entrenamiento

modelo_E <- multinom(Peso_RN_1 ~ ., data = peso_train)

# visualisar resultados (gráfico de efectos parciales)
ggeffect(modelo_E) %>%
  plot() %>%
  plot_grid()

model_performance(modelo_E)

rl_table = tbl_regression(modelo_E)
rl_table

rl_table1 = tbl_regression(modelo_E, intercept = T)
rl_table1

rl_table2 = tbl_regression(model, exponentiate = T)
rl_table2
interpret_oddsratio(c(1.19, 0.89, 0.53, 4.10, 1.85, 1.81))
?interpret_oddsratio

plot_model(modelo_E)
plot_model(modelo_E, show.intercept = F, show.values = T, width=0.1)+
  ylab('ODDS RATIO')
exp(86) # alpha

tab_model(modelo_E,
          show.reflvl = T,
          show.intercept = F,
          p.style = 'numeric_stars')

# chequeo visual de todos los supuestos del modelo
check_model(modelo_E)

#Calculamos el estadistico ji cuadrado
#SIGNIFICATIVIDAD GLOBAL DEL MODELO
#Estimamos el modelo solo con la constante
modelo_ls0 <- multinom(Peso_RN_1~1,data=peso_train)

ji=modelo_ls0$deviance-modelo_E$deviance
gl.ji=modelo_E$edf-modelo_ls0$edf
pvalor=1-pchisq(ji,gl.ji)
print(cbind(ji,gl.ji,pvalor))

str(modelo_E)
summary(modelo_E)
summary(modelo_E)$coefficients
coef(modelo_E)

#---------------------------------------------------------
#Intervalos para B

confint(modelo_E, level=0.95)

# Cociente de ventajas: Exp(B) (Odd Ratio)
exp(coef(modelo_E))
exp(confint(modelo_E, level=0.95))

# Estadistico W de Wald

regTermTest(modelo_E,"Edad_Madre",method="Wald")
regTermTest(modelo_E,"Num_Partos",method="Wald")
regTermTest(modelo_E,"Tipo_de_Parto",method="Wald")
regTermTest(modelo_E,"Sexo_RN",method="Wald")
regTermTest(modelo_E,"Talla_RN",method="Wald")
regTermTest(modelo_E,"Edad_Gestacional",method="Wald")

# Pseudo R cuadrado
pR2(modelo_E)

#----------------------------------------------------------
# Importancia de las variables
import=varImp(modelo_E)
orden=dplyr::arrange(import, desc(import))
orden

# importancia de las variables. Criterio de Anova
car::Anova(modelo_E) # Talla_RN es màs importante y Num_Partos la menos importante

rf = randomForest(Peso_RN_1 ~ ., data = peso_train)

vip(rf)

#Significatividad de los coeficientes
stargazer(modelo_E,type="text")

summary(modelo_E)
coef(modelo_E)[1,]

BajoPeso=coef(modelo_E)[1,]
BajoPeso
SobrePeso=coef(modelo_E)[2,]
SobrePeso

# PREDICCIONES
peso_train[1,]
peso_train[1,-1]
head(peso_train)
str(peso_train)
x1=peso_train[1,-1]
x1

predict(modelo_E,x1)

#   Calculo para valores predichos

Xp=c(1,20,2,0,1,50,30)
A=coef(modelo_E)[1,]; B=coef(modelo_E)[2,]
#A tiene los coeficientes de bajo peso y B los de sobre peso

# Probabilidad de bajo peso

Xp=c(1,20,2,0,1,50,30); P1=exp(sum(A*Xp))/(1+exp(sum(A*Xp))+exp(sum(B*Xp))); P1
log_oddsA=log(exp(sum(A*Xp)))
plogis(log_oddsA)

# Probabilidad de sobre peso
Xp=c(1,20,2,0,1,50,30); P2=exp(sum(B*Xp))/(1+exp(sum(A*Xp))+exp(sum(B*Xp))); P2

# Probabilidad de peso normal
P3=1/(1+exp(sum(A*Xp))+exp(sum(B*Xp))); P3

P1+P2+P3

#Se le clasificara como de bajo peso con probabilidad 0.9950425


# EVALUACION DEL MODELO DE ENTRENAMIENTO CON LA DATA DE PRUEBA

head(peso_test,11)

# Tabla de confusion con los datos de prueba

Prediccion <-predict(modelo_E, peso_test)
head(Prediccion)
#Pred=rep("Datos4_P$Peso_RN_1", length(Prediccion))
Tabla_Conf=table(Prediccion,peso_test$Peso_RN_1)
Tabla_Conf
# 

dim(peso_test)
tbc=(58+2+4)/70
tbc
tmc=1-tbc
tmc

clase.pred <- predict(modelo_E,peso_test)
head(clase.pred)
a=cbind(p_test = peso_test$Peso_RN_1, clas = clase.pred)
head(a,18)
head(peso_test)
Peso_RN_1.Ajustado=modelo_E$fitted.values
head(Peso_RN_1.Ajustado)
Residuales=modelo_E$residuals
head(Residuales)
head(cbind(Peso_RN_1.Ajustado,Residuales))

## De otro modo

cm <- caret::confusionMatrix(clase.pred, peso_test$Peso_RN_1)
print(cm)

## De otra manera

# Implementing function to collapse data
library(ConfusionTableR)
mc_df <- ConfusionTableR::multi_class_cm(clase.pred, peso_test$Peso_RN_1)
# Access the reduced data for storage in databases
print(mc_df$record_level_cm)
library(tibble)
glimpse(mc_df$record_level_cm)

mc_df$confusion_matrix
mc_df$cm_tbl

