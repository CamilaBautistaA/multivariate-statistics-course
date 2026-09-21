###################################################
#                                                 #
#    EJEMPLOS: Varios                             #
# Profesor: Clodomiro Fernando Miranda Villagomez #
# cfmiranda@lamolina.edu.pe                       #
#                                                 #
###################################################

# Ejemplo de Personalidad

library(tinytable)
library(mvnormtest)
library(psych)
library(car)

bd <- read.delim("clipboard")
tt(head(bd))
dim(bd)

# Supuesto de normalidad multivariada
pers2 <- bd[4:8]
tt(head(pers2))
mshapiro.test(t(pers2))

# Supuesto de variables dependientes correlacionadas. 
# Prueba de esfericidad de Bartlett

options(scipen=0)
cortest.bartlett(cor(pers2), n=nrow(pers2))

# Para hacer la prueba en conjunto (Significacion del modelo)
# H0: Ninguna de las variables independientes influye
# sobre las dependientes
# Ha: Al menos una de las variables independientes influye
# sobre las dependientes
mod1 = lm(cbind(Apertura, Responsabilidad, Extraversion, Amabilidad, Neuroticismo)~
          edad + depresion_puntaje + consumo_alcohol, data = bd)

alias(mod1)
# Los siguientes estadisticos son bastante robustos ante violaciones
# de normalidad

ev_mod1 <- linearHypothesis(mod1, hypothesis.matrix =
                            c("edad = 0", "depresion_puntaje = 0", 
                              "consumo_alcohol = 0"))
ev_mod1

# Pillai
1-pf(19.97782, 15,1944)
#Wilks
1-pf(23.11089, 15,1783.723)

# los cuatro criterios rechazan H0.

# Pruebas de hipotesis para verificar cuál o cuáles de las tres variables
# independientes influyen sobre las dependientes

summary(manova(mod1),test="Pillai")
summary(manova(mod1),test="Wilks")
summary(manova(mod1),test="Hotelling-Lawley")
summary(manova(mod1),test="Roy")

car::Anova(mod1)
car::Anova(mod1, type = 'III')
car::Anova(mod1, test = 'Wilks')
car::Anova(mod1, type = 'III', test = 'Wilks')
Manova(mod1)
Manova(mod1,type="III")

#A continuacion procedemos a utilizar las 20 primeras filas 
#de la data a manera de data de prueba para predecir
#valores de las variables dependientes y comparar con los resultados reales.
head(bd)
test.data = bd[c(1:20), c(1:3)]
predictions <- predict(mod1, test.data)
reales <- bd[c(1:20), c(4:8)]
comparacion <- cbind(predictions, reales)
colnames(comparacion) <- c("Apertura_predicha", 
  "Responsabilidad_predicha", "Extraversion_predicha", "Amabilidad_predicha", 
  "Neuroticismo_predicho", "Apertura_real","Responsabilidad_real", 
  "Extraversion_real", "Amabilidad_real", "Neuroticismo_real")
comparacion

# Calculando el error (metrica = mape)
mape = function(y_pred, y_real){
  n = length(y_pred)
  mape_res = abs(y_pred - y_real)/(y_real*n)
  return(sum(mape_res))
}

mape(comparacion$Apertura_predicha, comparacion$Apertura_real)
mape(comparacion$Responsabilidad_predicha, comparacion$Responsabilidad_real)
mape(comparacion$Extraversion_predicha, comparacion$Extraversion_real)
mape(comparacion$Amabilidad_predicha, comparacion$Amabilidad_real)
mape(comparacion$Neuroticismo_predicho,comparacion$Neuroticismo_real)

# Podemos observar que la estimación de la Responsabilidad no es tan
# buena como la estimación para la Apertura. Así mismo, podemos ver el R2 de los 
# dos modelos

res = summary(mod1)

res$`Response Apertura`$r.squared
res$`Response Responsabilidad`$r.squared
res$`Response Extraversion`$r.squared
res$`Response Amabilidad`$r.squared
res$`Response Neuroticismo`$r.squared

# Ejemplo de Colesterol


# Ingreso de datos
# Lectura de datos SPSS

library(foreign)
datos <-read.spss("Colesterol.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
head(datos)

# Supuesto de normalidad multivariada

datos1 = datos[c(-4:-7)]
head(datos1)

library(mvnormtest)
mshapiro.test(t(datos1))

# Supuesto de variables dependientes correlacionadas.
# Prueba de esfericidad de Bartlett

library(psych)
library(rela)
options(scipen = 0)
cortest.bartlett(cor(datos1), n=dim(datos1))
cortest.bartlett(cor(datos1), n=nrow(datos1))

# Para hacer la prueba en conjunto (Significacion del modelo)
# H0: Ninguna de las variables independientes influye
# sobre las dependientes
# Ha: Al menos una de las variables independientes influye 
# sobre las dependientes

datosc = datos
head(datosc)

#Los siguientes estadisticos son bastante robustos ante violaciones
#de normalidad

modelo1=lm(cbind(Y1,Y2,Y3) ~ X1 + X2 + X3 + X4, data = datosc)

library(car)
alias(modelo1)

library(car)
lh.out <- linearHypothesis(modelo1, 
                           hypothesis.matrix = c("X1 = 0", "X2 = 0",
                                                 "X3 = 0","X4 = 0"))
lh.out
?linearHypothesis

# Los cuatro criterios rechazan H0.

# Pruebas de hipotesis para verificar cual o cuales de las cuatro 
# variables independientes son las que influyen sobre las dependientes

datosc = datos
head(datosc)

summary(manova(modelo1), test="Pillai")
summary(manova(modelo1), test="Wilks")
summary(manova(modelo1), test="Hotelling-Lawley")
summary(manova(modelo1), test="Roy")

library(car)
Manova(modelo1, type="II")  # Sin intercepto
Manova(modelo1, type="III") # Con intercepto

#A continuacion procedemos a utilizar las 4 primeras filas 
#de la data a manera de data de prueba para predecir
#los valores de Y1, Y2 y Y3 y comparar con los resultados reales.
head(datos)
test.data = datos[c(1:4), c(4:7)]
predictions <- predict(modelo1, test.data)
resultados_reales <- datos[c(1:4), -c(4:7)]
comparacion <- cbind(predictions, resultados_reales)
colnames(comparacion) <- c("Y1_predicho", "Y2_predicho",
  "Y3_predicho", "Y1_real", "Y2_real", "Y3_real")
comparacion

str(modelo1)
modelo1$coefficients

# Calculando el error (metrica=mape)
mape = function(y_pred, y_real){
  n = length(y_pred)
  mape_res = abs(y_pred - y_real)/(y_real*n)
  return(sum(mape_res))
}

mape(comparacion$Y1_predicho, comparacion$Y1_real)
mape(comparacion$Y2_predicho, comparacion$Y2_real)
mape(comparacion$Y3_predicho, comparacion$Y3_real)
# Podemos observar que la estimación de la Y3 es la mejor.
# Así mismo, podemos ver el R2 de los 
# tres modelos

res = summary(modelo1)

res$`Response Y1`$r.squared
res$`Response Y2`$r.squared
res$`Response Y3`$r.squared


