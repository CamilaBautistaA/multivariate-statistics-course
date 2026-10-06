#######################################################
#                                                     #
#                ANALISIS FACTORIAL                   #
#       EJEMPLO: ALPACAS CON SEIS VARIABLES           #
#   Mg.Sc. Clodomiro Fernando Miranda Villagomez      #
#           cfmiranda@lamolina.edu.pe                 #
#                                                     #                                          #
#######################################################

#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o informacion cargada
rm(list = ls())


# Paquetes
library(psych)
library(ade4)
library(psych)
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(mvnormtest)
library(tidyverse)
library(GGally)
library(PerformanceAnalytics)
library(performance)

######################
#  Lectura de Datos  #
######################

# Lectura de datos con 6 variables
library(foreign)
datos <-read.spss("Alpacas.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
head(datos)
str(datos)

# Coeficientes de Correlacion 
data2=datos
head(data2)

# Abajo se ve que las variables (Peso, Mecha, Diametro, DesvE y Picazon) por
# un lado y CV por otro tienen alta correlacion esto es un indicio de que podemos
# hallar 2 factores
(correl=round(cor(data2),3))
round(solve(correl),3)#tiene inversa (es una matriz no singular). 
# El AF se puede trabajar con varios metodos para estimar parametros.
# En algunos de ellos no se requiere que la matriz correlaciones 
# sea no singular

# Prueba estadistica
corr.test(data2)

#################  
##             ##
##  SUPUESTOS  ##
##             ##
#################


# Supuesto de variables correlacionadas. 
# Prueba de Esfericidad de Bartlett

# H0: |Rp|=1  , Las variables no estan correlacionadas

options(scipen=0)
esfer=cortest.bartlett(cor(datos),n=dim(datos))
str(esfer)
esfer$p.value[1]

# Supuesto de normalidad multivariada: Prueba de Shapiro Wilk
dim(datos)
mshapiro.test(t(datos))

# Indicador Kaiser-Meyer-Olkin KMO y MSAj 

kmo = psych::KMO(datos)
str(kmo)
kmo$MSA
kmo$MSAi
kmo

KMO(datos)

# Prueba de esfericidad y kmo de otra manera
dim(datos)
check_factorstructure(datos)
check_factorstructure(correl, n = 94)
check_sphericity_bartlett(datos)
check_kmo(datos)
check_kmo(correl, n = 94)
check_sphericity_bartlett(correl, n = 94)

?check_kmo # para ver criterios respecto a los valores de kmo.

##########################################################
# Analisis Factorial con 6 variables con el Metodo de    # 
# Componentes Principales                                #
##########################################################

# Analisis Factorial sin rotacion con funcion principal 
# (tambien existe la funcion fa y factanal)
.72^2+(-.11)^2
facto.sin.rota <- principal(r=datos,nfactors=2,rotate="none")
?principal
# Con el pvalor=0.96 se concluye que 2 factore son suficientes
print(facto.sin.rota)
str(facto.sin.rota)
a=(facto.sin.rota$loadings[3,])^2 
# comunalidad de la variable Peso, considerando los 6 factores
sum(a)
sum(facto.sin.rota$values) # suma de los autovalores
#Proportion Var
facto.sin.rota$values[1]/6
facto.sin.rota$values[2]/6

#Proportion Explained
0.56/1
0.17/1
0.12/1

# Autovalores
facto.sin.rota$values
sum(facto.sin.rota$values)
#Se deben retener los 2 primeros factores porque sus autovalores 
#son mayores que 1.

# Grafica de Valores propios

# Primera forma
plot(facto.sin.rota$values,type="b",pch=20,col="steelblue",lwd=2)
abline(h=1,lty=3,col="tomato",lwd=3)

# Comunalidades y unicidad considerando los 6 factores
variable=c('peso','mecha','diámetro','desv','cv','picazón')
comunal=facto.sin.rota$communality
especif=facto.sin.rota$uniquenesses
d=tibble(variable,comunal,especif)
d %>% arrange(desc(comunal))

# Cargas Factoriales, Correlaciones Factor-variable
facto.sin.rota$loadings# Se omiten los valores menores de 0.10
.717^2+(-.105)^2+(-0.472)^2 #Comunalidad de "Peso" con 4 factores
head(facto.sin.rota$scores) # puntuaciones o scores de los 6 factores
scor=facto.sin.rota$scores
round(cor(scor),3)

# Analisis Factorial con rotacion utilizando la funcion principal
# Primero veamos el siguiente ejemplo

# Ejemplo: Una chica de 12 años hizo 5 calificaciones sobre sus conocidos.
Y=data.frame('gente'=c('fsm1','hermana','fsm2','padre','profesor','msm','fsm3'),
             'amabilidad'=c(1,8,9,9,1,9,9),
             'inteligencia'=c(5,9,8,9,9,7,7),
             'feliz'=c(5,7,9,9,1,7,9),
             'simpatia'=c(1,9,9,9,1,9,9),
             'justo'=c(1,8,8,9,9,9,7))
Y

spec=eigen(cor(Y[,2:6]))
str(spec)
round(spec$values,3)#autovalores
# Proportion var
round(spec$values/sum(spec$values),3)
# Cumulative var
round(cumsum(spec$values)/sum(spec$values),3)
# First 2 eigenvectors associated with 2 largest eigenvalues
round(spec$vectors[,1:2],3)

loadings=spec$vectors[,1:2]%*%diag(sqrt(spec$values[1:2]))
colnames(loadings)=paste0('Factor',1:2)
rownames(loadings)=colnames(Y)[-1]
round(loadings,3)

# A continuacion se aprecia que los factores rotados se asocian mejor con las variables
par(mfcol=c(1,2))
plot(loadings,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
text(loadings,row.names(loadings),pos=c(1,2,2,2,1))
abline(v=0,h=0,lty=2)
rot.loadings=loadings%*%matrix(c(cos(30.8),sin(30.8),-sin(30.8),cos(30.8)),nrow=2)
plot(rot.loadings,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),xlab='Factor1 Rot',ylab='Factor2 Rot')
text(rot.loadings,row.names(loadings),pos=c(1,2,2,2,1))
abline(v=0,h=0,lty=2)
par(mfcol=c(1,1))

# Volviendo al problema. Con el criterio del promedio se encontró que
# reteniendo 2 factores, ellos explican el 72.5% de la variabilidad
# total y viendo las correlaciones(loadings) entre factores y variables
# se concluye que no es necesario rotar, pero haremos una rotacion
# sólo para mostrar el procedimiento que se tendría que hacer de
# ser necesaria la rotación.

facto.con.rota <- principal(r=datos,nfactors=2,rotate="varimax")
?principal
dim(datos)
str(facto.con.rota)
print(facto.con.rota)
str(facto.sin.rota)
print(facto.sin.rota)
facto.con.rota$values #Autovalores sin rotar
#Notar que los autovalores rotados si coinciden con los anteriores.
#Si se rota con diferentes #s de factores cambian los Vaccounted
facto.con.rota$ Vaccounted 
facto.sin.rota$ Vaccounted
sum(facto.con.rota$ Vaccounted[1,])
sum(facto.sin.rota$ Vaccounted[1,])
facto.sin.rota$values
facto.con.rota$communality
facto.sin.rota$communality
facto.con.rota$loadings
facto.sin.rota$loadings
head(facto.con.rota$scores)
head(facto.sin.rota$scores)
round(cor(facto.con.rota$scores),4)

# Grafica de variables sobre el primer plano de componentes

# Dividiendo la ventana de graficos
#x <- c(1, 2, 3, 4) # cuatro gráficos diferentes
x=c(1,2,3,3)
m <- matrix(x, ncol = 2)
m
layout(m)
nf <- layout(m)
layout.show(nf)

load <- facto.con.rota$loadings[,1:2]
plot(load, pch=16, xlim=c(-1,1), ylim=c(-1,1),col="chartreuse3",
     xlab="Factor 1",ylab="Factor 2") 
abline(h=0,lty=3,col="brown1",lwd=2)
abline(v=0,lty=3,col="burlywood1",lwd=2)
text(load,pos=1,labels=names(datos),cex=1.1)#agrega los nombres a las variables

# Grafica de circulo de correlaciones

s.corcircle(load,grid=FALSE)
fa.diagram(facto.con.rota)
par(mfcol=c(1,1))

# Almacena los datos y los resultados de los scores en un
# archivo CSV
salida.facto <- cbind(datos,facto.sin.rota$scores)
head(salida.facto)
str(salida.facto)
write.csv(salida.facto,"Alpacas-factorial-con-acp.csv")


# Grafica de la Matriz de Correlacion entre Factores

cor.plot(cor(salida.facto[7:8]),
         main="Correlaciones entre Factores sin Rotar", 
         diag=TRUE,
         show.legend = TRUE)  

?fa

#Otra forma de retener el numero de factores
#Analisis Paralelo

library(paran)
paran(datos,iterations=5000,graph=TRUE,color=4)
3.350073-2.994749
#El numero de componentes que deberia retenerse es 1.
?fa.parallel
#Test Estadistico para retener m CP (H0:lamda(m+1)=lamda(m+2)=...=lamda(p)=0)
#Esta prueba tiene la limitacion de aconsejar la retencion de demasiadas Componentes Principales.
library(nFactors)
nBartlett(datos,N=94,alpha=0.01,cor=TRUE,details=TRUE)
#con este test se recomienda retener 5 factores.
head(datos)


