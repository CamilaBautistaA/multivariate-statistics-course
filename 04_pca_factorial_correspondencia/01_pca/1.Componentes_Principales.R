#####################################################
#                                                   #
#        ANALISIS DE COMPONENTES PRINCIPALES        # 
#              EJEMPLO Postulantes                  #
#   Mg.Sc. Clodomiro Fernando Miranda Villagomez    #
#           cfmiranda@lamolina.edu.pe               #
#                                                   #                     
#####################################################


#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o informacion cargada
rm(list = ls())


###############
#  Paquetes   #
###############

library(car)
library(tinytable)
library(GGally)
library(MVN)
library(epiDisplay)
library(scales)
library(rgl)
library(psych)
library(hrbrthemes)
library(gganimate)
library(png)
library(gifski)
library (dplyr)
library(viridis)
library(tidyverse)
library(forcats)
library(BSDA)
library(dlookr)
library(ggpubr)
library(summarytools)
library(pastecs)
library(corrplot)
library(mvnormtest)
library(PerformanceAnalytics)
library(ggcorrplot)
library(ade4)
library(factoextra)
library(patchwork)
library (adegenet)
library(FactoMineR)

####################
# Lectura de datos #
####################

library(foreign)
datos <-read.spss("1.Postulantes.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
tt(head(datos))
str(datos)
dim(datos)
datos1=datos[,-1]
tt(head(datos1))
attach(datos1)

#Verificando si hay datos perdidos
per.miss.col=100*colSums(is.na(datos1))/dim(datos1)[1]
per.miss.col

#Aplicando el criterio de la puntuacion Z (para 
#Observaciones outliers de cada variable)
is_outlier2 <- function(x,k = 2) {
  return(abs(scale(x)) > k)
}

datosP2=datos1[,-10]
head(datosP2)
datosP2[is_outlier2(datosP2$RV,3),]
datosP2[is_outlier2(datosP2$RM,3),]
datosP2[is_outlier2(datosP2$MAT,3),]
datosP2[is_outlier2(datosP2$PSI,3),]
datosP2[is_outlier2(datosP2$FIS,3),]
datosP2[is_outlier2(datosP2$LOG,3),]
datosP2[is_outlier2(datosP2$BIO,3),]
datosP2[is_outlier2(datosP2$HIS,3),]
datosP2[is_outlier2(datosP2$QUI,3),]

# No considerar la primera ni la columna once: Id y condicion.

datosacp <- datos[,c(-1,-11)]
tt(head(datosacp))
str(datosacp)

########################
# Analisis Descriptivo #
########################

summary(datosacp)

summarytools::descr(datosacp)

X=datosacp$RM
Xbar=mean(datosacp$RM)
Sd=sd(X)
mean(((X-Xbar)/Sd)^3)# Asimetria Muestral
mean(((X-Xbar)/Sd)^4)-3# Curtosis muestral

# Grouped statistics

head(datos1)
with(datos1, stby(RM, CON, descr))

with(datos1, stby(MAT, CON, descr))

describe(datosacp)

round(stat.desc(datosacp),2)
round(stat.desc(datosacp,basic=FALSE),2)

###########################
# Analisis de Correlacion #
###########################

#-----------------------------
# Coeficientes de Correlacion 
round(cor(datosacp),3)
diag(cor(datosacp))
sum(diag(cor(datosacp)))

# Prueba estadistica

corr.test(datosacp)

#----------------------------------------------
# Graficos de Correlacion

i=cor(datosacp,method="pearson")
corrplot(i,sig.level=0.05,type="lower")
?corrplot
corrplot(i,method = "number",order = "original")
corrplot(i,method = "ellipse",order = "original")
corrplot(i,method = "ellipse",order = "original",addCoef.col = "magenta")
corrplot(i,method = "ellipse",order = "original",addCoef.col = "magenta",type = "upper")
corrplot(i,method = "ellipse",order = "original",addCoef.col = "magenta",type = "lower")
corrplot(i,method = "square",order = "original",tl.pos = "d",addCoef.col = "magenta")
corrplot(i,method = "square",order = "original",tl.pos = "d")
corrplot(i,method = "pie",order = "AOE",tl.pos = "d",addCoef.col = "green2")
corrplot(i,method = "circle",order = "FPC",tl.pos = "d",addCoef.col = "green2")
corrplot(i,method = "circle",order = "FPC",tl.pos = "d")
corrplot(i,method = "color",order = "original",tl.pos = "d",addCoef.col = "green2")
corrplot(i,method = "color",order = "original",tl.pos = "d")
corrplot(i,method = "ellipse",order = "AOE",type="upper",tl.pos = "d")
corrplot(i,add=TRUE,type="lower",method = "number",order = "AOE",
         diag = FALSE,tl.pos = "n",cl.pos = "n")

res1=cor.mtest(datosacp,conf.level=0.05)
res1
a=datosacp$RV
b=datosacp$PSI
c=cbind(a,b)
head(c,10)

mshapiro.test(t(c))

corrplot(i, p.mat = res1$p, sig.level = 0.05)
corrplot(i, p.mat = res1$p, order = "hclust", insig = "pch", addrect = 3)
corrplot(i, p.mat = res1$p, insig = "p-value",sig.level = -1)#agrega los pvalores

# Primera forma
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))
corrplot(cor(datosacp),
         title = "Matriz de correlacion", mar=c(0,0,1,0),
         method = "color", outline = T, addgrid.col = "darkgray",
         order = "hclust", addrect = 3, col=col1(100),
         tl.col='black', tl.cex=.75)

# Segunda forma
pairs(datosacp,col="green2")

# Tercera forma

chart.Correlation(datosacp, histogram=TRUE, pch=20)

# Cuarta forma - Mapas de Calor
cor.plot(cor(datosacp),
         main="Mapa de Calor", 
         diag=TRUE,
         show.legend = TRUE)  

data2=datosacp

head(data2)
corr=cor(data2)

## Correlogramas

ggcorrplot(corr) +
  ggtitle("Correlograma de Postulantes") +
  theme_minimal()

ggcorrplot(corr, method = 'circle') +
  ggtitle("Correlograma de Postulantes") +
  theme_minimal()

ggcorrplot(corr, method = 'circle', type = 'lower') +
  ggtitle("Correlograma de Postulantes") +
  theme_minimal()

ggcorrplot(corr, method = 'circle', type = 'lower', lab = TRUE) +
  ggtitle("Correlograma de Postulantes") +
  theme_minimal() +
  theme(legend.position="none")

## Incluyendo pvalores

# Calculando los pvalores
p.mat <- cor_pmat(datosacp)
p.mat

# Agregando la No significacion (X) de las correlaciones
# --------------------------------
corr = cor(datosacp)
ggcorrplot(corr, hc.order = TRUE,
           type = "lower", p.mat = p.mat)

##################################################################
# 1. Analisis de Componentes Principales usando la libreria ade4 #
#                   Matriz de Correlaciones                      #
##################################################################
head(datosacp)

acp <- dudi.pca(datosacp,
                scannf=FALSE, scale=TRUE,
                nf=ncol(datosacp))# Con scale se tipifican las variables
summary(acp)
3.6617/9
acp[["eig"]]
str(acp)
print(acp)
3.6617/9
(3.6617+1.5135)/9
# Valores propios (autovalores)
acp$eig
sum(acp$eig)

inertia.dudi(acp)

# Vectores propios
acp$c1

# Correlaciones entre las variables originales y las componentes principales
acp$co
?cor
# Grafica de Valores propios - ScreePlot

# Primera forma
plot(acp$eig,type="b",pch=20,col="blue",lwd=2)
abline(h=1,lty=3,col="red",lwd=3.5)

a=fviz_eig(acp,choice='eigenvalue',geom="line",linecolor = '#3A5FCD',xlab = 'Componentes Principales')+
  geom_hline(yintercept = 1,color='#EE6363')+
  theme_grey()
a
?fviz_eig

# Segunda forma

eig.val <- get_eigenvalue(acp)
eig.val

barplot(eig.val[, 2], names.arg=1:nrow(eig.val), 
        main = "Autovalores",
        xlab = "Componentes Principales",
        ylab = "Porcentaje de variancias",
        col ="steelblue")
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type="b", pch=19, col = "red")

# Tercera forma

fviz_screeplot(acp)
fviz_screeplot(acp, ncp=6)
fviz_eig(acp, addlabels=TRUE, hjust = 0.5)

fviz_eig(acp, addlabels=TRUE, hjust = 0.5,
              barfill="white", barcolor ="darkblue",
              linecolor ="red") + ylim(0,50) + theme_minimal()

b=fviz_screeplot(acp, ncp=9, addlabels=TRUE,hjust = 0.5,linecolor = "#FC4E07",
               barfill = "#00AFBB",xlab = "Componentes Principales")
b

# Grafica de Variables sobre el circulo de correlaciones

# Primera forma 
s.corcircle(acp$co,grid=FALSE,xax = 1, yax = 2)
?s.corcircle
# Segunda forma
fviz_pca_var(acp,col.var = '#EE8262',axes = c(1, 3))
?fviz_pca_var
c=fviz_pca_var(acp, col.var="#FF3030")+theme_minimal()
c
dim(datosacp)

# Scores o Puntuaciones de cada individuo
acp$li[1:10,]

options(scipen=999)
round(cov(acp$li),4)
acp$eig
round(cor(acp$li),4)

describe(acp$li)

# Grafica de individuos sobre el primer plano de componentes

# Primera forma
s.label(acp$li,xax=1,yax=2,clabel=0.7,grid=FALSE,boxes=FALSE)

# Segunda forma 
fviz_pca_ind(acp,col.ind = "steelblue")

# Grafica de individuos sobre los componentes 2 y 3
s.label(acp$li,xax=2,yax=3,clabel=0.7,grid=FALSE,boxes=FALSE)


# Grafica de individuos sobre el primer plano con biplot

# Primera forma
s.label(acp$li,clabel=0.7,grid=FALSE,boxes=FALSE)
s.corcircle(acp$co,grid=FALSE,add.plot = TRUE,clabel=0.7)

# Segunda forma
d=fviz_pca_biplot(acp, repel = F,
                col.var = "#EE3A8C",
                col.ind = "green" )
d

(a + b)/(c + d)
a | b / c | d

# Grabar los datos y los resultados de los scores en un archivo CSV
salidaacp=cbind(datosacp,acp$li[,c(1,2,3)])
head(salidaacp)
str(salidaacp)
write.csv(salidaacp,"P.csv")

fviz_eig(acp, ncp = 9, addlabels=TRUE, hjust = 0.5,barfill = "violet",
         barcolor = "blue")
fviz_pca_ind(acp, repel = F,col.ind = 'steelblue')# Evitar superposicion de texto

#otra opcion de grafico

#Pareciera que hay dos grupos de postulantes
colorplot (acp$li, acp$li,transp = F, cex = 3,
           xlab="PC1", ylab = "PC2")
title ("Analisis PCA de Admision")
abline (v=0, h=0, col="black", lty=2)

# Formando grupos de postulantes

head(datosacp)
datosacp1=scale(datosacp)
round(head(datosacp1),2)
round(colMeans(datosacp1),8)
cov(datosacp1)

#Darle 3 componentes principales retenidos y 2 grupos
grp <- find.clusters(datosacp1, max.n.clust=8)

head(grp)

#DAPC: Analisis discriminante de Componentes Principales
#Darle 3 componentes y 1 funcion discriminante
dapc.WIDIV <- dapc(datosacp1, grp$grp)
?dapc
scatter (dapc.WIDIV, posi.da = "bottomright", bg = "white", pch = 17:22,
         cstar = 0)

# Cursos que discriminan mejor. Barras mas altas indican mejor discriminacion
contrib <- loadingplot(dapc.WIDIV$var.contr, axis=1,
                       thres=.07, lab.jitter=1)
is.data.frame(contrib)
contrib1=tibble(nombre=contrib$var.names,valor=contrib$var.values)

orden=contrib1 %>%
  arrange(desc(valor))

orden

######################################################
#OTRA ALTERNATIVA PARA HACER COMPONENTES PRINCIPALES #
######################################################

acp = PCA(datosacp,scale.unit = TRUE,ncp=9,graph = TRUE)
summary(acp)
str(acp)
sum(acp$ind$contrib[,1])#suma de ctr de los individuos sobre la componente 1
sum(acp$ind$cos2[1,])#suma de los cos2 del individuo 1 sobre todas las componentes
sum(acp$var$contrib[,3])#suma de los ctr de los cursos sobre la componente 3
sum(acp$var$cos2[3,])#suma de los cos2 del curso MAT sobre todas las componentes

head(datosacp)

respca = PCA(datosacp, scale.unit=TRUE, ncp=9, graph=TRUE)
summary(respca)
plot(respca, label = "none")
?plot
fviz_pca_var(respca, col.var = "steelblue")
respca$eig
respca$var

fviz_pca_ind(respca, label="none", habillage=datos1$CON)

fviz_pca_ind(respca, label="none", habillage=datos1$CON,
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_biplot(respca, label = "var", habillage=datos1$CON,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_minimal())

fviz_pca_biplot(respca, 
                # Individuals
                geom.ind = "point",
                fill.ind = datos1$CON, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "RdBu",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "v1", color = "Contrib",
                                    alpha = "Contrib"))

fviz_pca_ind(respca, col.ind = "#00AFBB", repel = TRUE)

#Con la libreria ExPosition
library(ExPosition)
ePCA <- epPCA(datosacp)
summary(datosacp)
#Editar las etiquetas
library(explor)
explor(respca)

#Realizando el PCA en Factoshiny
library(Factoshiny)
result=Factoshiny(datosacp)
#res.shiny=PCAshiny(respca)

#Determinando conglomerados(clusters) jerarquicos con componentes principales
res.hcpc <- HCPC(respca,nb.clust = 3)
?HCPC

str(res.hcpc)
res.hcpc$data.clust
Grupos=res.hcpc$data.clust$clust
table(Grupos)

##################################
# 2. DESCRIPCION DE LOS CLUSTERS #
##################################

datosf=cbind(datos1,Grupos)
head(datosf)
str(datosf)

# Diagrama de Cajas de variable RM segun Cluster
ggboxplot(datosf, x = 'Grupos', y = 'RM',
          fill = 'Grupos',
          palette = 'aaas',
          xlab = 'Grupos', ylab = 'Razonamiento Matemático') +
  stat_kruskal_test(        # También stat_anova_test, stat_welch_anova_test
    label = 'as_detailed_italic',
    label.y = 35,
    size = 7) +
  geom_pwc(
    # También 't_test', 'sign_test', 'dunn_test', 'emmeans_test', 'tukey_hsd', 'games_howell_test' 
    method = 'wilcox_test',  
    label = 'Wilcoxon, p = {p.adj.format}{p.adj.signif}',
    label.size = 6,
    p.adjust.method = 'fdr',
    y.position = c(22, 26, 30),
    bracket.nudge.y = 0.03,
    tip.length = 0.03,
    step.increase = 0.1
  )
boxplot(datosf$RM ~ datosf$Grupos, 
        main= "BoxPlot de RM  vs CLUSTER",
        xlab = "Cluster", 
        names=c("Cluster 1", "Cluster 2", "Cluster 3"),
        col = c("red","blue","peru"))


# Diagrama de Cajas de la variable MAT segun Cluster
ggboxplot(datosf, x = 'Grupos', y = 'MAT',
          fill = 'Grupos',
          palette = 'npg',
          xlab = 'Grupos', ylab = 'Matemáticas') +
  stat_kruskal_test(
    label = 'as_detailed_italic',
    size = 7,
    label.y = 35) +
  geom_pwc(
    method = 'wilcox_test',
    label = 'Wilcoxon, p = {p.adj.format}{p.adj.signif}',
    label.size = 6,
    p.adjust.method = 'bonferroni',
    y.position = c(22, 26, 30),
    bracket.nudge.y = 0.03,
    tip.length = 0.03,
    step.increase = 0.1
  )

boxplot(datosf$MAT ~ datosf$Grupos, 
        main= "BoxPlot de MAT  vs CLUSTER",
        xlab = "Cluster", 
        names=c("Cluster 1", "Cluster 2", "Cluster 3"),
        col = c("red","blue","peru"))

# Perfil en base a las medias de los resultados
attach(datosf)
rv <- tapply(RV,Grupos,mean) ; rv
rm <- tapply(RM,Grupos,mean)  ; rm
mat <- tapply(MAT,Grupos,mean) ; mat
psi <- tapply(PSI,Grupos,mean) ; psi
fis <- tapply(FIS,Grupos,mean) ; fis
log <- tapply(LOG,Grupos,mean)   ; log
bio <- tapply(BIO,Grupos,mean) ; bio
his <- tapply(HIS,Grupos,mean)  ; his
qui <- tapply(QUI,Grupos,mean) ; qui

medias <- rbind(rv,rm,mat,psi,fis,log,bio,his,qui)   ; medias
general <- c(mean(RV),mean(RM),mean(MAT),
             mean(PSI),mean(FIS),mean(LOG),
             mean(BIO),mean(HIS),mean(QUI))   ; general
medias <- cbind(medias,general)
str(medias)
medias

matplot(medias,
        main = "Grafico de promedios de Variables segun Cluster",
        xlab = "Variables",
        ylab = "Promedios",
        type="l",
        xaxt="n",         # Permite eliminar los nombres del eje X
        ylim=c(-2,20), 
        col=c("blue","red","green2","black"))
axis(1,at=1:9,labels=c("RV","RM","MAT","PSI","FIS","LOG","BIO","HIS","QUI"))

legend("topright", c("Cluster 1", "Cluster 2", "Cluster 3","General"), 
       pch=c(5,5,5,5), ncol=4, cex=0.8, 
       col=c("blue","red","green2","black"), bty="n")

#Para obtener conglomerados con los componentes principles tambien se puede
#usar la libreria Factoshiny. Entrar a "Principal Component Analysis", hacer 
#check en "perform clustering after leaving PCA app?", escoger el numero de 
#cluster en "Number of dimensions kept for clustering" y despues 
#"quit the app", finalmente explorar los cluster.
#los cluster que se haya elegido.
library(Factoshiny)
result1=Factoshiny(datosacp)

## MISCELANEA

#Metodo Paralelo para la retencion de componentes principales.
#Cuando hay subjetividad en el grafico de sedimentacion (Scree Plot) respecto al numero de CP a retener
#se puede recurrir al Metodo Paralelo

library(paran)
paran(datosacp,iterations=5000,graph=TRUE,color=2)
#se confirma que se debe retener 3 CP.

#Test Estadistico para retener m CP (H0:lamda(m+1)=lamda(m+2)=...=lamda(p)=0)
#Esta prueba tiene la limitacion de aconsejar la retencion de demasiadas Componentes Principales.
library(nFactors)
nBartlett(cor(datosacp),N=541,alpha=0.01,cor=TRUE,details=TRUE)
#con este test se recomienda retener 7 u 8 CP.

