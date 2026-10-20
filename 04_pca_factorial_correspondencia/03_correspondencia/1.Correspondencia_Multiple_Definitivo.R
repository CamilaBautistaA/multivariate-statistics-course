
# ==================================== #
#     Análisis de Correspondencia      #
#              Múltiple                #
# Caso: Admitidos a Ciencias en una    #
# Universidad de Colombia en el 2013-I #
# ==================================== #

library(FactoClass)
library(tidyverse)
library(GGally)
library(gtsummary)
library(ggrepel)
library(sjPlot)
library(ggpubr)
browseVignettes("sjPlot")

# NOTA: En la literatura clásica (Benzécri, Greenacre y Husson) el ACM se 
# considera principalmente un método exploratorio y descriptivo, por lo que el
# énfasis no está en una prueba de hipótesis global, sino en cuantificar y 
# representar la estructura de asociación existente entre las variables
# categóricas. El supuesto es la asociación entre el conjunto de variables 
# categóricas y para verificar esto se puede hacer:
# 1. La prueba chi cuadrado para cada par de variables.
# 2. Revisar medidas de asociación (como la V de Cramér).
# 3. Evaluar la inercia y los autovalores obtenidos en el ACM.
# 4. Interpretar los mapas factoriales para identificar estructuras de
#    asociación.


data(admi)
Z1=admi[,1:15]
head(Z1)
dim(Z1)
Y<-admi[,c(8,11,9,10)] 
head(Y)
dim(Y)
summary(Y) 
view_df(Y, show.frq=T, show.prc = T, show.na = T)

Y %>% 
  plot_frq(edad)

Y %>% 
  group_by(gene) %>%
  plot_frq(estr) %>%
  plot_grid()

plot_xtab(
  x = Y$edad,
  grp = Y$estr,
  margin = 'row',
  bar.pos = 'stack',
  show.summary = T,
  coord.flip = T
)

Y %>%  
  ggpairs(cardinality_threshold = 20)

# Tabla combinada con comparación por estrato

tabla_combinada = Y %>%
  select(gene, edad, estr, orig) %>%
  tbl_summary(
    by = estr, # variable de comparación
    include = c(gene, edad, orig),
    percent = 'row', # porcentaje por fila
    statistic = list(all_categorical() ~ '{n} ({p})'),
    digits = all_categorical() ~ c(0, 2),
  ) %>%
  add_overall(last = TRUE) %>% # agrega la columna 'total'
  add_p(
    test = everything() ~ 'chisq.test', # aplica chi-cuadrado a todas las variables
    pvalue_fun = ~ style_pvalue(.x, digits = 3) # formato del p-valor 
  ) %>%
  modify_header(label ~ '**Variable**') %>% # cambia el encabezado de la columna de variables
  modify_caption('Tabla 1: Comparación de Genero,\nEdad y Origen  por Estrato')
tabla_combinada

table(Y$gene,Y$edad)
table(Y$gene)

# PIE de las variables activas

par ( mfrow =c(2 ,2) , mai =c (0.1 ,0 ,0.2 ,0))
for (i in 1:4 ){
  cat <- attributes (Y[,i])$ levels
  per <- tabulate (Y[,i])/ nrow (Y)* 100
  paste (cat , round (per ,1) , sep =" ")->eti
  pie ( summary (Y[,i]),eti ,
        col = gray ( seq (1.0 ,0.6 , length = length ( cat ))) ,
        main = colnames (Y)[i])
}

par(mfrow = c(1, 1))

library(RColorBrewer)
RColorBrewer::display.brewer.all()

par ( mfrow =c(2 ,2) , mai =c (0.1 ,0 ,0.2 ,0))
for (i in 1:4 ){
  cat <- attributes (Y[,i])$ levels
  per <- tabulate (Y[,i])/ nrow (Y)* 100
  paste (cat , round (per ,1) , sep =" ")->eti
  pie ( summary (Y[,i]),eti ,
        col = hcl.colors ( 4,"YlGnBu") ,
        main = colnames (Y)[i])
}

par(mfrow = c(1, 1))

table(Y$edad)
(a = 100*c(118, 171, 56, 100)/445)
sum(a)

df <- data.frame(
  Edad = c("a16m", "a17", "a18", 'a19M'),
  Frec = c(26.52, 38.43, 12.58, 22.47)
) 

labs <- paste0(df$Edad, " (", df$Frec, "%)")

ggpie(data = df, x = "Frec",
      label = labs,
      fill = "Edad", color = "white",
      palette = 'npg')

ggpie(data = df, x = "Frec",
      label = labs,
      lab.pos = "in", lab.font = "black",
      fill = "Edad", color = "white",
      palette = c("#00AFBB", "yellow", "#FC4E07", 'green2')
)

ggdonutchart(data = df, x = "Frec",
      label = labs,
      fill = "Edad", color = "white",
      palette = 'npg')

ggdonutchart(data = df, x = "Frec",
      label = labs,
      lab.pos = "in", lab.font = "white",
      fill = "Edad", color = "white",
      palette = 'lancet'
)

# Para la tabla , registros multiplos de 25 
Z<-acm.disjonctif(Y); data.frame(Y,Z)[seq(0,nrow(Y),25),] 
xtable(data.frame(Y,Z)[seq(0,nrow(Y),25),],digits=rep(0,17))
data.frame(Y,Z)[seq(0,nrow(Y),25),]
head(Z,14)
dim(Z)


# Tabla de Burt del ejemplo Admitidos a Ciencias
# la tabla de Burt cruza las 12 categorias de las 4 variables de
# los admitidos a Ciencias. Algunas lecturas son:
# entre los admitidos hay 128 mujeres y 317 hombres;
# 46 mujeres tienen 16 años o menos y 81 de los admitidos son de 
# estrato alto.

B<-acm.burt(Y,Y); xtable(B,digits=rep(0,13))
B

# Sea p=Num de categorias consideradas, s=Num de variables categoricas y 
# n=el total de individuos.
# Las masas marginales fila fij =1/n ; i=1,...,445 y las  marginales
# columna f.j =nj/(n*s); j=1,...,12 donde nj es el numero de individuos
# que asumen la categoria j. En el ejemplo fij =1/445=0.00225  ;
# i=1, 2, ..., 445. f.2 =nj/(n*s)=317/(445*4)=0.1781 y 
# f.11=38/(445*4)=0.02135

summary(Y)
1/445
317/(445*4)
38/(445*4)
445*4-(128+317+118+171+56+100+179+185+81+311+38+96)

# un subespacio de dimension p-s, es el numero de valores propios
# mayores que cero.En el ejemplo hay 12-4 = 8 valores propios mayores que cero,
# que se muestran a continuacion

acm<-dudi.acm(Y,scannf=FALSE,nf=8) 
barplot(acm$eig,cex.axis=0.6,col="darkolivegreen1") 
dev.print(device = xfig,file="ACMadmiValP.fig") 
eiglst<-data.frame(vp=acm$eig,porce=acm$eig*100/sum(acm$eig), 
                   acupor=cumsum(acm$eig)*100/sum(acm$eig)) 
xtable(eiglst,digits=c(1,3,1,1)) #tabla en formato LaTeX
eiglst
eiglst$vp
sum(eiglst$vp)
0.2293303/sum(eiglst$vp)
sum(eiglst$porce)
# vp quiere decir valores propios
# La decisión del número de ejes a retener se basa sobre todo en
# la forma del histograma puesto que el ACM tiene ejes parásitos, 
# es decir, que aparecen y no contienen informacion. 

# El porcentaje de inercia no es un criterio apropiado en el caso del ACM.
# Del histograma de la figura anterior se puede concluir que tres ejes
# son suficientes, a pesar de que retienen menos del 50% de la inercia (46.6%).

#  Histograma de valores propios del ACM de admitidos. 
# Notese que, los tres primeros ejes se destacan sobre los demás
# y que retienen el 46.6% de la inercia

# Codigo para obtener el primer plano factorial

# con gg=TRUE se usa ggplot2 y ggrepel para la grafica 

selin<-seq(25,445,25) # Individuos cuya numeracion es multiplo de 25
plot(acm,Tcol=FALSE,roweti=as.character(selin),cframe=1,cex.row=0.6,
     cex.global=0.8,gg=TRUE)

# Criterios para interpretar el siguiente gráfico:
# 1. La cercania de los atributos indica mayor asociacion. Ejemplo: estr.bajo y
# orig.Cund.De otro lado, orig.Otro no se asocia con otros atributos.
# 2. Cuanto mas lejano del origen mayor asociacion. Ejemplo estr.medio y 
# edad.a18 estan mas asociadas que edad.a17 y gene.M
# 3. Cuando las categorias estan opuestas por el origen se dice que están 
# asociadas inversamente. Ejemplo: gene.F y edad.a19M o edad.a16m y orig.Bogo
plot(acm,Tcol=TRUE,roweti=as.character(selin),cframe=1,cex.row=0.6,
     cex.global=0.8,gg=TRUE)

# Las contribuciones de cada individuo: absolutas, relativas e inercia
# están en porcentaje (%)

# Gráfico

head(Y,10)
summary(Y)
317/(445*4)
library(ca)
cm=mjca(Y[,1:4],lambda='indicator',nd=NA)
?mjca
str(cm)
cm
cm$Burt

plot(cm)
head(cm$ rowmass,8) # masa de fila
1/(445)
sum(cm$ rowmass)
cm$ colmass # masa de columna
sum(cm$ colmass)
cm$ colinertia # inercia de columna
sum(cm$ colinertia)

# DE OTRA MANERA

library(FactoMineR)
?MCA
# estrato y origen estan asociados a la dim 1 y edad con genero a la dim 2
# si la ubicacion de una variable estuviera cerca al origen no 
# mostraria relacion con las dimensiones
fit2=MCA(Y)
str(fit2)
# v.test: verifica si la contribución absoluta de la categoría a la dimension
# es significativa (en otras palabras si aporta a la dimension). Si 
# abs(v.test)>1.96 la prueba es significativa (es una prueba de Z)

summary(fit2,nb.dec=1,ncp=8)
qnorm(0.975)

library(factoextra)
fviz_mca_var(fit2,col.var = "darkorchid")
?fviz_mca_var
fviz_contrib(fit2,choice="var",axes=1,fill = "mediumspringgreen",
             color = "mediumorchid4")
?fviz_contrib
fviz_contrib(fit2,choice="var",axes=2,fill = "lemonchiffon",
             color = "lightcoral")
fviz_contrib(fit2,choice="var",axes=1:2,fill = "lightskyblue1",
             color = "indianred1")
fviz_contrib(fit2,choice="var",axes=5,fill = "lightskyblue1",
             color = "indianred1")
fviz_screeplot(fit2, addlabels = TRUE, ylim = c(0, 18),hjust=0.5)
fviz_eig(fit2, addlabels=TRUE, hjust = 0.5,
         barfill="white", barcolor ="darkblue",
         linecolor ="red") + ylim(0,18) + theme_minimal()

library(Factoshiny)
result=Factoshiny(fit2)
Factoshiny(result)

library(Factoshiny)
result=Factoshiny(Y)
Factoshiny(result)

# Si con la data se quiere hacer la correspondencia simple entre
# Edad y Estrato hay que hacer una tabla.
# Edad en filas

head(Y,20)
tabla=table(Y$edad,Y$estr)
tabla
is.matrix(tabla)

############################################
# Prueba de Independencia Chi-Cuadrado     #
############################################

# H0: La edad y el estrato son independientes
prueba <- chisq.test(tabla)
prueba
pchisq(18.991,6,lower.tail=FALSE)#pvalor

qchisq(0.05,6,lower.tail=FALSE)
library(fastGraph)
shadeDist( 12.59159, "dchisq", 6, 
           lower.tail=FALSE,col = c("limegreen","red4") )

# De otra manera

library(ggstatsplot)
library(effectsize)
library(rstatix)

ggbarstats(
  data=Y,
  x=estr,
  y=edad,
  label = 'both'
)

# Se observa que el número de postulantes de estrato alto disminuye 
# conforme aumenta la edad (lo contrario ocurre en el estrato bajo).
# El pvalor de 0.00418 indica que las variables son dependientes 
# (relacionadas) pero no dice que tan fuerte es esa relación, para 
# eso tenemos el tamaño del efecto que en este caso es el valor v de 
# Cramer (0.12) que dice que esa relación es pequeña y su IC va de
# 0 a 0.18. También se aprecia la v de Cramer bayesiana (0.13) que 
# también afirma que la relación es pequeña, pero el tamaño del 
# efecto no prueba la hipótesis, aparte del pvalor también se muestra
# el factor de sesgo o factor de Bayes (ln(BF01) = 0.02) que indica 
# que la prueba muestra evidencia anecdótica a favor de la H1,

interpret_cramers_v(0.12)
interpret_cramers_v(0.14)
?interpret_cramers_v

interpret_bf(exp(0.02))
?interpret_bf

interpret_p(0.004179)
?interpret_p

# Si se quiere probar H0:La proporción de postulantes de los tres
# estratos dentro de cada grupo de edad es igual, se tienen los 4
# pvalores: 0.05, 0.00458, 0.0035 y 1.92e-08.

# Perfiles y masas
library(gmodels)
CrossTable(tabla,
           prop.t=T,
           prop.r=TRUE,
           prop.c=TRUE,
           prop.chisq=FALSE)

library(FactoMineR) 

res.ca = CA(tabla,ncp=4,graph=TRUE)

print(res.ca)

summary(res.ca)
#etc

#Si se quiere considerar Estrato en filas.

tabla=table(Y$estr,Y$edad)
is.matrix(tabla)
tabla


############################################
# Prueba de Independencia Chi-Cuadrado     #
############################################

prueba <- chisq.test(tabla)
prueba
# etc

# MISCELÁNEA

library(janitor)
tabyl(Y$edad) # show_na = FALSE

Y1 = Y %>%
  tabyl(edad, estr)

# funciones agradabls
Y %>%
  tabyl(edad, estr) %>%
  #adorn_totals()
  adorn_totals(c('col', 'row'))

# para obtener los porcentajes
Y %>%
  tabyl(edad, estr) %>%
  adorn_totals('row') %>%              # también se puede usar 'col'
  adorn_percentages('row') %>%         # también se puede usar 'col'
  adorn_pct_formatting() %>%
  adorn_ns() %>%                       # position = 'front
  adorn_title() %>%
  adorn_title('combined')

# porcentajes independientes de las filas y columnas
Y %>%
  tabyl(edad, estr) %>%
  adorn_totals(c('col', 'row')) %>%
  adorn_percentages('all') %>%
  adorn_pct_formatting() %>%
  adorn_ns()
