
###################################################
# FAMD: ANÁLISIS FACTORIAL CON DATOS MIXTOS       #
#                                                 #
# Profesor: Clodomiro Fernando Miranda Villagómez #
###################################################

# Carga de librerías ------------------------------------------------------

library(FactoMineR)
library(factoextra) # Visualización de resultados
library(tidyverse) # Manipulación, transformación y visualización de datos
library(sjPlot)#fácil visualización
library(performance)
library(tidymodels) # Herramientas para estimación modelos
library(jsonlite) # Leer código json
library(visdat) # Vistazo a los datos
library(naniar) # Análisis de datos perdidos
library(GGally) # Complemento de visualización
library(tidytext)
library(effectsize)


# Lectura de datos --------------------------------------------------------

library(foreign)
Y <-read.spss("AFMixto_Diabetes.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
head(Y)
Y=Y[,-1]

Y$diabetes = fct_recode(Y$diabetes,
                    Si_Diab = 'Si',
                    No_Diab = 'No')

head(Y)
dim(Y)
summary(Y)
# Análisis Exploratorio de los Datos

view_df(Y,show.frq=T,show.prc=T,show.na=T)

Y %>% 
  plot_frq(niveles_peso)

p=Y %>% 
  group_by(diabetes) %>%
  plot_frq(niveles_peso) %>%
  plot_grid()
p

save_plot(
  filename = 'diabetes_vs_peso.jpg',
  fig=p,width=30
)

plot_grpfrq(
  var.cnt = Y$niveles_peso,
  var.grp = Y$antec.fa
)

plot_xtab(
  x = Y$diabetes,
  grp = Y$antec.fa,
  margin = 'row',
  bar.pos = 'stack',
  show.summary = T,
  coord.flip = T
)

tab_xtab(
  var.row = Y$diabetes,
  var.col = Y$antec.fa,
  show.row.prc = T
)

Y %>% 
  group_by(antec.fa) %>% 
  plot_frq(
    colester,
    type = 'histogram',
    show.mean = T,
    normal.curve = T
  ) %>% 
  plot_grid()

head(Y)
m=lm(colester~edad+diabetes,data=Y)
summary(m)
contrasts(Y$diabetes)
plot_model(m)
plot_model(m,show.values = T,width=0.1)+
  ylab('íncremento de colesterol comparado con la edad')

tab_model(m,
          show.reflvl = T,
          show.intercept = F,
          p.style = 'numeric_stars')

# El efecto de la edad y diabetes sobre el colesterol es muy débil

interpret_r2(0.048,rules='cohen1988')
?interpret_r2

windows()
check_model(m, verbose = TRUE)

dev.off()
check_normality(m)
plot(check_normality(m))
check_heteroscedasticity(m)
plot(check_heteroscedasticity(m))
check_collinearity(m)
plot(check_collinearity(m))
check_outliers(m)
plot(check_outliers(m))

# FAMD --------------------------------------------------------------------

dim(Y)
summary(Y)

# Se muestran 5 gráficos en los ejes 1 y 2 (buscar con la flecha inversa)
famd_pre_result <- Y %>% 
  FAMD(ncp = 7, graph = T)
famd_pre_result
famd_pre_result$eig
famd_pre_result$var

# Se muestran 5 gráficos en los ejes 2 y 3 (buscar con la flecha inversa)
famd_pre_result1 <- Y %>% 
  FAMD(ncp = 7, graph = T, axes = c(2,3))
famd_pre_result1

str(famd_pre_result)

# famd_pre_result %>% fviz_famd()
famd_pre_result %>% fviz_famd_var()
?fviz_famd
# Analisis de valores propios ---------------------------------------------

(eigenvalues <- famd_pre_result %>% get_eigenvalue())
famd_pre_result %>% fviz_screeplot()
famd_pre_result %>% fviz_screeplot(addlabels = TRUE, ylim = c(0, 30),
                                   barfill="white", barcolor ="darkblue",
                                   linecolor ="red")
a = famd_pre_result$eig
a
str(a)
b = a[1:7,1]
b

# Suma de autovalores
sum(b)

# Para los porcentajes sin acumular y acumulados se divide entre 7
(b/7)*100
(b[[1]]/7)*100
(b[[2]]/7)*100

c = a[1:7,2]
c

d = (b/sum(b))*100
d
cumsum(d)

e = (b/7)*100
e
cumsum(e)

# Analisis de variables

# General

var_analysis <- famd_pre_result %>% get_famd_var()
var_analysis$coord
var_analysis$cos2
var_analysis$contrib

famd_pre_result %>% fviz_contrib(choice = "var", axes=1,)
famd_pre_result %>% fviz_contrib(choice = "var", axes=1:2)
famd_pre_result %>% fviz_contrib(choice = "var", axes=1:5)

# Cuantitativas

var_analysis_quanti <- famd_pre_result %>% get_famd_var(element = "quanti.var")
var_analysis_quanti$coord
var_analysis_quanti$cos2
aaa = var_analysis_quanti$contrib
aaa
str(aaa)
(s = sum(aaa[1:2, 1])) # total contrib de la dim.1 con variables cuantitativas

head(Y)
summary(Y)
famd_pre_result %>% 
  fviz_famd_var(choice = "quanti.var", col.var = "contrib",
                gradient.cols = c("red","yellow","green"))
?fviz_famd_var

famd_pre_result %>% 
  fviz_famd_var(choice = "quanti.var", col.var = "contrib",axes=c(2,3),
                gradient.cols = c("red","yellow","green"))

famd_pre_result %>% 
  fviz_famd_var(choice = "quanti.var", col.var = "cos2",
                gradient.cols = c("red","yellow","green"))

# Cualitativas

var_analysis_quali <- famd_pre_result %>% get_famd_var(element = "quali.var")
var_analysis_quali
var_analysis_quali$coord
# v.test: pruebas Z. Si el valor es mayor en valor absoluto a 1.96 entonces
# la categoria aporta a la dimension (es significativa)

var_analysis_quali$v.test
aa = var_analysis_quali$contrib
aa
str(aa)
(s1 = sum(aa[1:9, 1])) # La suma de contrib (con variables cualitativas) respecto a la Dim.1 no es 100
s + s1 # contrib total a la dim.1 de variables cuantitativas y cualitativas

(bb = var_analysis_quali$cos2)
str(bb)
sum(bb[9, 1:7]) # La suma de cos2 de obesidad da 1

famd_pre_result %>% 
  fviz_famd_var("quali.var", col.var = "contrib",
                gradient.cols = c("red","yellow","green"))

var_analysis_quali$cos2

famd_pre_result %>% 
  fviz_famd_var("quali.var", col.var = "cos2",
                gradient.cols = c("red","yellow","green"))


# Analisis de individuos --------------------------------------------------
dim(Y)
ind_analysis <- famd_pre_result %>% get_famd_ind()
ind_analysis

str(ind_analysis)

ind_analysis$coord[1:7,] %>% View()

ind_analysis$contrib[1:10,] %>% View()
sum(ind_analysis$contrib[1:1750, 4])#La suma de las contrib para la Dim.4 es 100

ind_analysis$cos2[1:8,] %>% View()
sum(ind_analysis$cos2[3, 1:7]) # La suma de cos2 del paciente 3 es 1

famd_pre_result %>% 
  fviz_famd_ind(col.ind = "contrib", 
                gradient.cols = c("red", "yellow", "green"))

famd_pre_result %>% 
  fviz_famd_ind(col.ind = "cos2", 
                gradient.cols = c("red", "yellow", "green"))

#A continuacion entrar a Factor Analysis on Mixed Data del Factoshiny, 

library(Factoshiny)
result=Factoshiny(Y)
Factoshiny(result)

#==================================================================#
# VALIDACIÓN DE CONSTRUCTO (DATOS ORDINALES)                       #
#==================================================================#
Constructo<-read.delim("clipboard")
head(Constructo)

test=Constructo
head(test)
test = test[, -1]

#------------------------------------------------------------------#
#  Análisis Factorial Policórico                                   #
#------------------------------------------------------------------#
library(psych)

#Probando supuestos

#Prueba de Normalidad multivariante (No tiene sentido con datos ordinales)
# ----------------------------
#library(mvnormtest)
#mshapiro.test(t(test))

# ============================
#Esfericidad de Bartlet
# ----------------------------
r.poly=polychoric(test)
R=r.poly$rho    # La matriz de correlaciones policórica
cor(test)      # correlaciones de Pearson (solo para comparar)
n = nrow(test)  # Tamaño de la muestra
cortest.bartlett(R,n)


#KMO(Kaiser-Meyer-Olkin)
# ----------------------------
KMO(R)
test1=test[-13]
head(test1)
r.poly1=polychoric(test1)
R1=r.poly1$rho    # La matriz de correlaciones policórica
KMO(R1)
library(performance)
check_factorstructure(R, n)

#Identificando el N° de factores
# ----------------------------
fap <- fa.parallel(R,n.obs=n,fa="fa")

# ----------------------------

#### METODOS DE ESTIMACION

# ml=Máxima verosimilitud
# minres=Mínimo residuals
# pa=Ejes principal Components
# gls=Mínimo cuadrados generalizados
# wls=Mínimo cuadrados ponderados (Apropiado para variables dicotomicas)
# uls=Mínimo cuadrado no ponderado (Apropiado para variables ordinales)
# --------------------------------------------------------------
?fa
factorial1=fa(test, nfactors=3,n.obs=n,rotate='none',fm="uls",cor="poly")
factorial1
e1=fa.diagram(factorial1, e.size=.05,rsize=3.5,digits=2,col="blue")

factorial2=fa(test, nfactors=3,n.obs=n,rotate="varimax",fm="uls",cor="poly")
factorial2

#Gráfica de segmentación de factores rotados
# ----------------------------
e=fa.diagram(factorial2, e.size=.05,rsize=3.5,digits=2,col="blue")

#Cargas Factoriales
# ----------------------------
dat = data.frame(ifelse(abs(factorial2$loadings)==apply(abs(factorial2$loadings),
                                                       1,max), apply(abs(factorial2$loadings),1,max),NA))
library(tidyverse)
dat1 = data.frame(Item=rownames(dat),dat)
fact = gather(dat1,Factor,cargas,-Item)%>% drop_na()
Cargas=split(fact, fact$Factor)
Cargas

#==================================================================#
# VALIDACIÓN DE CONSTRUCTO (DATOS DICOTOMICOS)                       #
#==================================================================#
Constructo = read.delim("clipboard")
test = Constructo
head(test)

#------------------------------------------------------------------#
#  Análisis Factorial Tetracorico                                   #
#------------------------------------------------------------------#
library(psych)
summary(test)
#Probando supuestos

#Prueba de Normalidad multivariante (No tiene sentido con datos dicotomicos)
# ----------------------------
# library(mvnormtest)
# mshapiro.test(t(test))

# ============================
#Esfericidad de Bartlet
# ----------------------------
r.tetra=tetrachoric(test)
R=r.tetra$rho    # La matriz de correlaciones tetracoricas
cor(test)
n = nrow(test)  # Tamaño de la muestra
cortest.bartlett(R,n)


#KMO(Kaiser-Meyer-Olkin)
# ----------------------------
KMO(R)
library(performance)
check_factorstructure(R, n)

#Identificando el N° de factores
# ----------------------------
fap <- fa.parallel(R,n.obs=n,fa="fa")

#Corriendo el modelo
# ----------------------------

# Ya que con las variables son dicotomicas se tiene que usar correlaciones tetracoricas

## Usar el Metodo de Minimos Cuadrados ponderados

factorial2=fa(test, nfactors=2,n.obs=n,rotate="none",fm="wls",cor="tet")
factorial2
e=fa.diagram(factorial2, e.size=.05,rsize=3.5,digits=2,col="blue")


factorial3=fa(test, nfactors=2,n.obs=n,rotate="varimax",fm="wls",cor="tet")
factorial3

#Gráfica de segmentación de factores
# ----------------------------
e1=fa.diagram(factorial3, e.size=.05,rsize=3.5,digits=2,col="blue")

#Cargas Factoriales
# ----------------------------
dat = data.frame(ifelse(abs(factorial3$loadings)==apply(abs(factorial3$loadings),
                                                        1,max), apply(abs(factorial3$loadings),1,max),NA))
library(tidyverse)
dat1 = data.frame(Item=rownames(dat),dat)
fact = gather(dat1,Factor,cargas,-Item)%>% drop_na()
Cargas=split(fact, fact$Factor) # se obtiene el valor absoluto de las cargas
Cargas



