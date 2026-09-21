
#####################################
##      INFERENCIA MULTIVARIADA    ##
#####################################

rm(list = ls())

library(nortest)
library(gtExtras)
library(sjPlot)
# devtools::install_github("strengejacke/strengejacke")
library(summarytools)
library(GGally)
library(BSDA)
library(inferr)
library(mvnormtest)
library(MVTests)
library(MVN)
library(tidyverse)
library(ggstatsplot)
library(effectsize)
library(car)
library(reshape2) 
library(ggthemes)
library(ggpubr)
library(ICSNP)
library(heplots)
library(biotools)
library(ergm)
library(rrcov)
library(Hotelling)
library(echarts4r)
library(gganimate)
library(hrbrthemes)
library(png)
library(gifski)
library(tinytable)
library(rstatix)

############################################################
##                     PRUEBA Z^2                         ##
##   Inferencia Respecto a un Vector de Medias Simple.    ##         # 
##       Con Matriz Varianza Covarianza Conocida          ##             #
############################################################


## EJEMPLO 1

### CEBICHE DE PESCADO

# Importar datos
head(datos)
tt(head(datos, 20))
datos=datos[-1]
tt(tail(datos))
glimpse(datos)
datos$distrito = as.factor(datos$distrito)


# Análisis exploratorio

gt_plt_summary(datos, title = 'Ceviche de Pescado')
view_df(datos, show.frq=T, show.prc = T, show.na = T)
summary(datos)
summarytools::descr(datos)

table(datos$distrito)
datos %>% freq_table(distrito)

# Correlaciones

datos1 = datos[-5]
tt(head(datos1))

# Cálculo de una correlación
var(datos1)
cor(datos1$energia,datos1$agua)
-4.366129/(7.218653*4.339015)^.5

datos1 %>% cor_test()
datos1 %>% cor_test(method = 'spearman')

datos |>
  group_by(distrito) |>
  cor_test(
    vars = energia,
    vars2 = c(agua, carbohidrato, proteina), method = 'spearman') |>
  add_significance()

cor_mat(datos1, method = 'pearson') |>
  cor_mark_significant()

datos %>% ggpairs(.,
                   mapping = ggplot2::aes(colour = distrito),
                   lower = list(continuous = wrap('smooth', alpha = 0.3, size = 0.1)),
                   upper = list(continuous = wrap('cor', size = 4.5))) +
  scale_color_manual(values = c('Rimac' = 'yellow', 'Ate' = 'red' ,'Surco' = 'darkgreen')) +
  scale_fill_manual(values = c('Rimac' = 'yellow', 'Ate' = 'red' ,'Surco' = 'darkgreen'))

# CASO UNIVARIADO: Prueba de hipotesis para la media de agua con sigma=2.38
# conocida
# H0: mu = 76.2        H1: mu != 76.2
# H0: El contenido de agua proviene de una distribucion normal
d1 = datos$agua
shapiro.test(d1)
lillie.test(d1)

z.test(datos$agua,alternative="two.sided",mu=76.2,sigma.x = 2.38,conf.level=0.95)

# Prueba de hipotesis respecto a una varianza (H0:sigma2 = 2.38^2,   
# H1:sigma2 != 2.38^2)
ifr_os_var_test(datos, agua, sd=2.38, alternative = 'both')
2*pchisq(37.535,49)
qchisq(0.025,49)
qchisq(0.975,49)

## CASO MULTIVARIADO: Prueba de hipotesis para un vector de medias con la
## matriz varianza-covarianza conocida (CEBICHE)

# Supuesto: H0: La muestra proviene de una distribucion normal tetravariada
tt(head(datos1))

mshapiro.test(t(datos1))

# Parámetros
n <- 50       # Tamaño de la muestra
p <- 4         # Número de variables
mu <- c(90,75.3,13,7.3)  # Vector de medias
sigma <- matrix(c(9,-4.5,3.5,3.5,
                  -4.5,4,-2.5,-2.5,
                  3.5,-2.5,2,1.5,
                  3.5,-2.5,1.5,2),nrow=p)  # Matriz de varianza-covarianza

# Calcular la media muestral
media_muestral <- colMeans(datos1)

# Calcular la estadística Z^2
Z2 <- t(media_muestral - mu) %*% solve(sigma / n) %*% (media_muestral - mu)

# Comparar con la distribución chi-cuadrado con p grados de libertad
p_valor <- 1 - pchisq(Z2, df = p)

# Resultados
cat("Estadística Z^2:", Z2, "\n")
cat("p-valor:", p_valor, "\n")

# Decisión
if (p_valor < 0.05) {
  cat("Rechazar H0: el vector de medias es diferente del vector de medias hipotético.\n")
} else {
  cat("No se puede rechazar H0: no hay suficiente evidencia para decir que el vector de medias es diferente del vector de medias hipotético.\n")
}

## Prueba: H0:la matriz de covarianza es igual a la Matriz dada (sigma)
result <- Bcov(datos1, Sigma=sigma)
summary(result)

#############################################################
##                    PRUEBA T^2 de HOTELLING              ##
##    Inferencia Respecto a un Vector de Medias Simple.    ##
##       Con Matriz Varianza Covarianza Desconocida        ##
#############################################################

## EJEMPLO 2

# pH, alcohol, densidad y cloruros en vinos

#INFERENCIAS RESPECTO A UN VECTOR DE MEDIAS SIMPLE

#Extension multivariada de la prueba t de Student con un vector de medias

# Importar datos
dim(datos)
tt(head(datos,10))
summarytools::descr(datos)
sd(datos$densidad)

# CASO UNIVARIADO: Prueba de hipotesis para la media de alcohol en 100ml
# de vino con sigma2 desconocido
# H0: mu = 10.5        H1: mu != 10.5

# Supuesto
# H0: alcohol tiene distribucion normal

shapiro.test(datos$alcohol)
lillie.test(datos$alcohol)

t.test(datos$alcohol, alternative = "two.sided" , mu = 10.5, conf.level = 0.95)
t.test(alcohol ~ 1, data = datos, alternative = "two.sided", mu = 10.5,
       conf.level = 0.95)

datos %>% t_test(alcohol ~ 1, alternative = 'two.sided', mu = 10.5)

# Prueba t univariada H0: mu=10.5
gghistostats(data = datos$alcohol %>% as_tibble(),
             x=value,
             test.value = 10.5,  # 0 es por default
             type = 'p',      # 'np' para Wilcoxon, 'r' prueba robusta, 'bs' bayesiana
             normal.curve = T,
             normal.curve.args = list(linewidth=2,col='purple'))

# el efecto del alcohol es muy pequeño para rechazar H0
interpret_hedges_g(-0.15)
?interpret_hedges_g

# CASO MULTIVARIADO
#Deteccion de outliers multivariados

par(mfrow = c(1, 2))
# Distancia de Mahalanobis
outliers <- mvn(data = datos, multivariateOutlierMethod = "quan")
# Distancia ajustada de Mahalanobis
outliers.adj <- mvn(data = datos, multivariateOutlierMethod = "adj")

par(mfrow = c(1, 1))

# Supuesto de normalidad multivariada de datos (la matriz de datos): 
# Prueba de Shapiro Wilk

# Para muestras grandes, T^2 de Hotelling puede ser robusto frente a desviaciones
# leves de normalidad multivariada debido al TLC.

mshapiro.test(t(datos))

# Test MVN de Royston.  no se aconseja usar este test si los datos
# cuentan con mas de 5000 o menos de 3 observaciones, ya que depende
# del test de Shapiro Wilk.
royston <- mvn(data = datos, mvnTest = "royston", multivariatePlot = "qq")
str(royston)
royston
royston$multivariateNormality$`p value`

#Test MVN de Mardia
mardia <- mvn(data = datos, mvnTest = "mardia")
mardia$multivariateNormality

#Test MVN de Henze-Zirkler
hz <- mvn(data = datos, mvnTest = "hz")
hz$multivariateNormality
hz$multivariateNormality$`p value`

#Test MVN de Doornik-Hansen
dh <- mvn(data = datos, mvnTest = "dh")
dh$multivariateNormality

#Test MVN de Szekely-Rizzo
energy <- mvn(data = datos, mvnTest = "energy")
energy$multivariateNormality

# T2 de Hotelling
colMeans(datos)
is.data.frame(datos)
mu=c(0,0,0,0)
result=OneSampleHT2(datos,mu0=mu,alpha = 0.05)
summary(result)
mu=c(3.32, 10.37, 0.997, 0.089)
result=OneSampleHT2(datos,mu0=mu,alpha = 0.05)
summary(result)

#Usando la funcion lm

y <- as.matrix(datos)
head(y)

y <- as.matrix(datos) - matrix(c(0,0,0,0),nrow(datos),ncol(datos),byrow=T)
anova(lm(y ~ 1))
y <- as.matrix(datos) - matrix(c(3.32, 10.37, 0.997, 0.089),nrow(datos),ncol(datos),byrow=T)
anova(lm(y ~ 1))

# de otra manera
mu1=c(0,0,0,0)
mu2=c(3.32, 10.37, 0.997, 0.089)

HotellingsT2(datos, mu = mu1)
HotellingsT2(datos, mu = mu2)

#Prueba robusta (aplicar cuando no se cumplen los supuestos y muestra grande)
#Las constantes alpha, d y q se deben elegir convenientemente.

RHT2(data=datos,mu0=mu1,alpha=0.72,d=1396.59,q=1132.99)
RHT2(data=datos,mu0=mu2,alpha=0.72,d=1396.59,q=1132.99)

#Extendiendo los intervalos de confianza a las regiones

#Regiones de confianza para el vector de medias normal

n <- nrow(datos)
p <- ncol(datos)
xbar <- colMeans(datos)
xbar
S <- cov(datos)
S

tconst <- sqrt((p/n)*((n-1)/(n-p)) * qf(0.99,p,n-p))
id <- c(1,2)
id
plot(ellipse(center=xbar[id], shape=S[id,id], radius=tconst, draw=F),
     xlab="pH",ylab="alcohol",col="tomato",lwd = 10)

id <- c(1,3)
id
plot(ellipse(center=xbar[id], shape=S[id,id], radius=tconst, draw=F),
     xlab="pH",ylab="densidad",col="green3",lwd = 10)

id <- c(2,4)
id
plot(ellipse(center=xbar[id], shape=S[id,id], radius=tconst, draw=F),
     xlab="alcohol",ylab="cloruro",col="darkblue",lwd = 10)


#Intervalos confidenciales simultaneos basados en T^2

#Funcion R para Intervalos confidenciales simultaneos basados en T^2.

T.ci <- function(mu, Sigma, n, avec=rep(1,length(mu)), level=0.95){
  p <- length(mu)
  if(nrow(Sigma)!=p) stop("Need length(mu) == nrow(Sigma).")
  if(ncol(Sigma)!=p) stop("Need length(mu) == ncol(Sigma).")
  if(length(avec)!=p) stop("Need length(mu) == length(avec).")
  if(level <=0 | level >= 1) stop("Need 0 < level < 1.")
  cval <- qf(level, p, n-p) * p * (n-1) / (n-p)
  zhat <- crossprod(avec, mu)
  zvar <- crossprod(avec, Sigma %*% avec) / n
  const <- sqrt(cval * zvar)
  c(lower = zhat - const, upper = zhat + const)
}

n <- nrow(datos)
p <- ncol(datos)
xbar <- colMeans(datos)
S <- cov(datos)

T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,1,0,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,0,1,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,0,0,1))

#Comparando I de C simultaneos T^2 con los I de C clasicos t.

TCI <- tCI <- NULL
for(k in 1:4){
  avec <- rep(0, 4)
  avec[k] <- 1
  TCI <- c(TCI, T.ci(xbar, S, n, avec))
  tCI <- c(tCI,
           xbar[k] - sqrt(S[k,k]/n) * qt(0.975, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(0.975, df=n-1))
}

rtab <- rbind(TCI, tCI)
round(rtab, 6)

#I de C simultaneos mas precisos

#ICs simultaneos a traves del metodo de Bonferroni con R

TCI <- tCI <- bon <- NULL
alpha <- 1 - 0.05/(2*4)
for(k in 1:4){
  avec <- rep(0, 4)
  avec[k] <- 1
  TCI <- c(TCI, T.ci(xbar, S, n, avec))
  tCI <- c(tCI,
           xbar[k] - sqrt(S[k,k]/n) * qt(0.975, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(0.975, df=n-1))
  bon <- c(bon,
           xbar[k] - sqrt(S[k,k]/n) * qt(alpha, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(alpha, df=n-1))
}

rtab <- rbind(TCI, tCI, bon)
round(rtab, 2)

#Formar intervalos de confianza para muestras grandes 
# (T^2  aprox.   Chi(p))

chi <- NULL
for(k in 1:4){
  chi <- c(chi,
           xbar[k] - sqrt(S[k,k]/n) * sqrt(qchisq(0.95, df=p)),
           xbar[k] + sqrt(S[k,k]/n) * sqrt(qchisq(0.95, df=p)))
}

rtab <- rbind(TCI, tCI, bon,chi)
rtab
round(rtab, 2)

#Regiones de prediccion para futuras observaciones

#Formando regiones de prediccion con R
n <- nrow(datos)
p <- ncol(datos)
xbar <- colMeans(datos)
S <- cov(datos)

pconst <- sqrt((p/n)*((n^2-1)/(n-p)) * qf(0.99,p,n-p))
id <- c(1,2)
plot(ellipse(center=xbar[id], shape=S[id,id], radius=pconst, draw=F),
     xlab="pH",ylab="alcohol",col="pink3",lwd = 10)

id <- c(1,3)
plot(ellipse(center=xbar[id], shape=S[id,id], radius=pconst, draw=F),
     xlab="pH",ylab="densidad",col="red",lwd = 10)

id <- c(2,4)
plot(ellipse(center=xbar[id], shape=S[id,id], radius=pconst, draw=F),
     xlab="alcohol", ylab="cloruro",col="green3",lwd = 10)

######################################################
##              PRUEBA T^2 de HOTELLING             ##
##     Inferencias sobre dos vectores de medias     ##
##          Con dos muestras dependientes           ##
######################################################

# Extensión multivariable de la prueba t de dos muestras
# Al igual que en el caso univariado, la prueba T^2 para muestras
# dependientes multivariables realiza la prueba de una muestra en
# una puntuacion de DIFERENCIA.

## EJEMPLO 3

# 25 individuos califican dos modelos de calzado

# Importar datos
tt(head(datos, 12))
d = datos[-1]
tt(head(d, 12))
dim(datos)

# CASO UNIVARIADO: Comparando solo la variable confort
# H0:mu1 = mu2         H1:mu1 != mu2
# H0:mu1-mu2 = 0       H1:mu1-mu2 != 0

ggplot(datos, aes(x = MODELO, y = confort))+
  geom_boxplot(fill=c('#00FFFF','#00FF00'))+
  labs(x='Modelos de Calzado',y='Calificación',
    title='Calificación del Confort de dos Modelos de Calzado')+
  theme_economist()

modelo1 = datos$confort[1:25]
modelo2 = datos$confort[26:50]

t.test(modelo1, modelo2, paired = TRUE)

t_test(data = datos, confort ~ MODELO, alternative = "two.sided", paired = TRUE,
       conf.level = 0.95)

t_test(data = datos, confort ~ MODELO, alternative = "two.sided", paired = TRUE,
       conf.level = 0.95, detailed = TRUE)

dif = modelo1 - modelo2

shapiro_test(dif)
shapiro.test(dif)

id = 1:25
df = tibble(id, modelo1, modelo2)
df_largo = melt(df, id = 'id')

comparaciones=list(c('modelo1','modelo2'))

ggplot(data=df_largo,aes(x=variable,y=value))+
  geom_boxplot(fill=c('#FFD700','#00BFFF'))+
  labs(x='Confort',y='Calificación',
       title='Dos Modelos de Calzado')+
  theme_dark()+
  stat_compare_means(comparisons=comparaciones,method='t.test',paired=T)

ggwithinstats(
  data = df_largo,
  x=variable,
  y=value,
  type = 'p'
)

interpret_hedges_g(0.98) # El efecto del confort es grande para rechazar H0

# Lo anterior es equivalente a lo que sigue
df1 = tibble(dif)
gghistostats(data = df1,
             x=dif,
             test.value = 0,  # 0 es por default
             type = 'p',      # 'np' para Wilcoxon
             normal.curve = T,
             normal.curve.args = list(linewidth=2,col='tomato'))

# CASO MULTIVARIADO
m1 = datos[1:25, 3:7]
m2 = datos[26:50, 3:7]
X=m1-m2
X
str(X)

# Distancia de Mahalanobis
outliers <- mvn(data = X, multivariateOutlierMethod = "quan")
# Distancia ajustada de Mahalanobis
outliers.adj <- mvn(data = X, multivariateOutlierMethod = "adj")

# Supuesto de normalidad multivariada de la matriz de diferencias: 
# Prueba de Shapiro Wilk y otras

mshapiro.test(t(X))

royston <- mvn(data = X, mvnTest = "royston", multivariatePlot = "qq")
royston$multivariateNormality

#Test MVN de Mardia
mardia <- mvn(data = X, mvnTest = "mardia")
mardia$multivariateNormality

#Test MVN de Henze-Zirkler
hz <- mvn(data = X, mvnTest = "hz")
hz$multivariateNormality

#Test MVN de Doornik-Hansen
dh <- mvn(data = X, mvnTest = "dh")
dh$multivariateNormality

#Test MVN de Szekely-Rizzo
energy <- mvn(data = X, mvnTest = "energy")
energy$multivariateNormality

## PRUEBA T2 de Hotelling

# Con la libreria MVTests
result <- Mpaired(T1=m1,T2=m2)
summary(result)

#Usando la libreria ICSNP 

xbar <- colMeans(X)
xbar
HotellingsT2(X)
HotellingsT2(X,mu=c(-0.8,2.6,-0.1,0.82,4.2))

#Usando la funcion lm

y <- as.matrix(X)
anova(lm(y ~ 1),test="Hotelling-Lawley")
y <- as.matrix(X) - matrix(c(0,0,0,0,0),nrow(X),ncol(X),byrow=T)
anova(lm(y ~ 1),test="Hotelling-Lawley")
y <- as.matrix(X) - matrix(c(-0.8,2.6,-0.1,0.82,4.2),nrow(X),ncol(X),byrow=T)
anova(lm(y ~ 1),test="Hotelling-Lawley")

######################################################
##              PRUEBA T^2 de HOTELLING             ##
##     Inferencias sobre dos vectores de medias     ##
##          Con dos muestras independientes         ##
######################################################

## EJEMPLO 4

# CON MUESTRAS INDEPENDIENTES. 
# Dos especies de escarabajo

# importar los datos
tt(head(datos, 18))
str(datos)
datos$Especie = as.factor(datos$Especie)
torax1=datos %>% filter(Especie=='Especie1') %>% pull(Torax)
torax2=datos %>% filter(Especie=='Especie2') %>% pull(Torax)

summary(datos)

# CASO UNIVARIADO:  = la variable Torax (H0: mu1 = mu2      H1:mu1 != mu2)

# Prueba de homogeneidad de varianzas (H0:sigma1^2=sigma2^2)
var.test(torax1,torax2)
var(torax1)/var(torax2)
datos %>% levene_test(Torax ~ Especie)

shapiro.test(torax1)
shapiro.test(torax2)

datos3 <- datos[,-1]
head(datos3)
str(datos3)
summarytools::descr(datos3)
(Ala=25.39/280.09)
(Torax=16.25/185.61)

# El siguiente gráfico sólo considera varianzas heterogéneas
comparaciones=list(c('Especie1','Especie2'))

ggplot(data=datos,aes(x=Especie,y=Torax))+
  geom_boxplot(fill=c('#E6E6FA','#FFE4B5'))+
  labs(x='Especies',y='Tórax',
       title='Longitud del Tórax de Escarabajos')+
  theme_economist()+
  stat_compare_means(comparisons=comparaciones,
                     method='t.test')

## La prueba t de student para comparar promedios se usa
## cuando ambas muestras independientes son 'normales', de
## esa manera el promedio representa a cada grupo. También
## las varianzas tienen que ser homogéneas.
## La prueba t de Welch también compara promedios y asume 
## que los grupos son normales pero las varianzas no son
## homogéneas.

t.test(torax1,torax2,alternative="two.sided",
       paired=F,var.equal=T,conf.level = 0.95)

datos %>% t_test(Torax ~ Especie, alternative="two.sided",
       paired=F,var.equal=T,conf.level = 0.95)

t.test(torax1,torax2,alternative="two.sided",
       paired=F,var.equal=F,conf.level = 0.95)

head(datos)

ggbetweenstats(
  data = datos,
  x=Especie,
  y=Torax,
  type = 'p', # probar con var.equal = F
  var.equal = T
)

interpret_hedges_g(1.05) #El efecto de la especie es grande para rechazar H0

## Estadistica Descriptiva.

# Full summary statistics
head(datos)
datos %>% get_summary_stats(Torax)

datos %>% get_summary_stats(Ala)

# Summary statistics of grouped data
# Show only common summary
datos2 %>%
  group_by(Especie) %>%
  get_summary_stats(Torax, type = "common")

datos %>%
  group_by(Especie) %>%
  get_summary_stats(Ala, type = "common")

summarytools::descr(datos)

# Grouped statistics
library(summarytools)
with(datos, stby(Torax, Especie, descr))

with(datos, stby(Ala, Especie, descr))


# Supuesto de normalidad multivariada de la matriz X e Y: 
# Prueba de Shapiro Wilk

X=filter(datos,
         Especie == "Especie1")[, 2:3]
X
mshapiro.test(t(X))

Y=filter(datos,
         Especie == "Especie2")[, 2:3]
Y
mshapiro.test(t(Y))

# Prueba de Homogeneidad de matrices variancia covariancia
# H0: Sigma1 = Sigma2
results <- BoxM(data=datos[,2:3],datos$Especie)
summary(results)

X=(datos %>% filter(Especie == "Especie1"))[, 2:3]
Y=(datos %>% filter(Especie == "Especie2"))[, 2:3]

var(X)# matriz varianza covarianza de la especie 1
var(Y)# matriz varianza covarianza de la especie 2

#De acuerdo a la prueba M de Box no se puede asumir que 
#las matrices de varianzas y covarianzas son iguales, por lo 
#anterior se buscan otras funciones de R ya que las librerias MVTests,
#ergm y rrcov nos permiten realizar la prueba T2 de Hotelling
#indicando este aspecto

datos
G<-c(rep(1,10),rep(2,13))
results1 <- TwoSamplesHT2(data=datos[,2:3],group=G,alpha=0.05,
                          Homogenity = F)
summary(results1)

X
Y
summary(X)
summary(Y)
approx.hotelling.diff.test(X, Y, mu0 = 0,
                           assume.indep = T, var.equal = F)

approx.hotelling.diff.test(X, Y, mu0 = c(16.5,-16.8),
                           assume.indep = T, var.equal = F)

#Prueba ROBUSTA de Hotelling (solo hace pruebas con mu=c(0,0) )
T2.test(X, Y, mu=0,method = "c")

## COMO EJEMPLO: ASUMIENDO QUE LAS MATRICES VARIANZA COVARIANZA SON IGUALES
## se puede trabajar de la siguiente manera (En nuetro caso
## no se cumple esta asunción) 
head(datos)
HotellingsT2(filter(datos,
                    Especie == "Especie1")[, 2:3],
             filter(datos,
                    Especie == "Especie2")[, 2:3])

HotellingsT2(X, Y)
HotellingsT2(X, Y,mu=c(16.5,-16.8))

fit=hotelling.test(X, Y, shrinkage = FALSE, perm = FALSE,
                   B = 10000, progBar = (perm && TRUE))
fit

######################
#                    #
#     MISCELANEA     #
#                    #
######################

#Lectura de ambas especies en una data
datos
head(datos)

## Analisis Exploratorio Grafico de las dos especies de escarabajo

library(paletteer)
paletteer_d("ggthemes::Color_Blind")
mi_color <- paletteer_d("ggthemes::Color_Blind")[c(4,5,6)] 

## Grafico PIE
table(datos$Especie)

Mi_df <- data.frame(n = c(10, 13), 
                    x = c("Especie1", "Especie2")) %>%
  mutate(percent = round(n/sum(n), 2) )
p=Mi_df %>%
  e_charts(x)  %>% 
  e_pie(n) %>%
  e_tooltip()

p %>%
  e_color(mi_color)

## Grafico de barras
p1=datos |>
  e_charts(x=Especie)|>
  e_bar(serie=Torax) |>
  e_tooltip(trigger = 'item')

p1 %>%
  e_color(mi_color)


# grafica de dispersion basica#

# Un scatterplot basico de Torax y Alas con el color que depende de
# la Especie
p1=ggplot(datos, aes(x=Torax, y=Ala, color=Especie)) + 
  geom_point(size=3) +  theme_ipsum()+
  xlab("Longitud Torax en micras")+
  ylab( "Longitud Alas en 0.01 mm " )

p1

# B
ggplot(datos, aes(x=Torax, y=Ala, shape=Especie, 
                   size=Especie, color=Especie)) + 
  geom_point() +  theme_ipsum()+
  xlab("Longitud Torax en micras")+
  ylab( "Longitud Alas en 0.01 mm " )+
  labs(colour="ESPECIE")

# BOXPLOT

boxplot(datos$Torax~datos$Especie,
        main="Escarabajos",xlab="Especie",ylab = "Longitud Torax",
        col=c("deepskyblue4","darkolivegreen2"))

df.m <- melt(datos, id.var = "Especie")
p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Especie))+ 
  facet_wrap( ~ variable, scales="free") 
p

is.data.frame(datos)
datos |>
  e_charts() |>
  e_boxplot(Torax, outliers = TRUE) |>
  e_boxplot(Ala, outliers = TRUE)


##HISTOGRAMAS

hist(datos$Torax,col="steelblue",main = "Longitud de Torax",
     xlab = "Long. Torax en micras",ylab = "Frecuencia")
shapiro.test(datos$Torax)

datos |>
  e_charts() |>
  e_histogram(Torax, name = "Histograma") |>
  e_density(Torax, areaStyle = list(opacity = .4), smooth = TRUE, 
            name = "Densidad", y_index = 1) |>
  e_tooltip(trigger = "axis")

# Histogramas espejo#
p4 <- datos %>%
  ggplot( aes(x=Torax, fill=Especie)) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="ESPECIE")+
  xlab("Long. Torax en micras")+
  ylab( "Frecuencia de Escarabajos" )

p4

library(RColorBrewer)
RColorBrewer::display.brewer.all()

p5 <- datos %>%
  ggplot( aes(x=Ala, fill=Especie)) +
  geom_histogram(alpha=0.7, position = 'identity') +
  scale_fill_brewer(palette="Spectral", 
                    name="Especie") +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="ESPECIE")+
  xlab("Long. Ala en micras")+
  ylab( "Frecuencia de Escarabajos" )

p5

# Correlaciones
head(datos)
library(GGally)
datos %>% ggpairs(.,
                   mapping = ggplot2::aes(colour = Especie),
                   lower = list(continuous = wrap('smooth', alpha = 0.3, size = 0.1)),
                   upper = list(continuous = wrap('cor', size = 4.5))) +
  scale_color_manual(values = c('Especie1' = 'darkgreen', 'Especie2' = 'yellow')) +
  scale_fill_manual(values = c('Especie1' = 'darkgreen', 'Especie2' = 'yellow'))

### Uso de la libreria ExPanDaR

### ESCARABAJOS
datos
library(ExPanDaR)
ExPanD(datos)


#### USO DEL ggplot builder

#install.packages("esquisse")

attach(datos)
str(datos)

esquisse::esquisser()



















