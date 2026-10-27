#######################################################
#                                                     #
#         ANALISIS DE CORRESPONDENCIAS SIMPLE         #
#   EJEMPLO: Horas Chateando                          #
#   Profesor: Clodomiro Fernando Miranda Villagomez   #
#         cfmiranda@lamolina.edu.pe                   #
#                                                     #
#######################################################

#---------------------------------------------------------
# Para limpiar el workspace, por si hubiera algun dataset 
# o informacion cargada
rm(list = ls())

#################
#  1. Paquetes  #
#################

library(MASS)
library(ca)
library(anacor)
library(FactoMineR)
library(gplots)
library(vcd)
library(graphics)
library(factoextra)
library(vegan)


########################
# 2. Ingreso de Datos  #
########################
 
# Antes de realizar el Analisis de Correspondencias se ingresan 
# los datos, se calculan las frecuencias y se hallan los perfiles
# de perfiles.
#

datos.acs <- matrix(c(13,10,11,16,5,
                     20,27,27,19,2,
                     9,27,71,16,32,
                     8,11,41,24,11),nrow=4,byrow=T)
datos.acs
is.matrix(datos.acs)

# Asignacion de nombres a las filas y columnas de la tabla

dimnames(datos.acs)<-list(Horas=c("H1", "H2", "H3", "H4")
                          ,Calificacion=c("A","B","C","D","E"))

datos.acs
addmargins(datos.acs)

############################################
# 4. Prueba de Independencia Chi-Cuadrado  #
############################################

# H0: Las dos variables son independientes

prueba <- chisq.test(datos.acs)
prueba

# Chi-squared distribution with 12 d.f.

qchisq(0.05,12,lower.tail=FALSE)
library(fastGraph)
shadeDist( 21.02607, "dchisq", 12, lower.tail=FALSE )

prueba$observed
prueba$expected
55*50/400
# Tabla con el paquete gmodels y función CrossTable()
# Perfiles fila, perfiles columna, masas de fila y masas de columna
library(gmodels)
CrossTable(datos.acs,
           prop.t=T,
           prop.r=TRUE,
           prop.c=TRUE,
           prop.chisq=F)

############################################
#  5. Analisis de Correspondencias Simple  #
#     con el paquete FactoMiner            #
############################################

library(FactoMineR)

?CA
res.ca = CA(datos.acs,ncp=5,graph=TRUE)

print(res.ca)

summary(res.ca)
str(res.ca)
res.ca$eig
res.ca$eig[,1]
sum(res.ca$eig[,1])

#----------------------------------
# Scree Plot de los Autovalores

library(factoextra)
eig.val <- get_eigenvalue(res.ca)
eig.val



fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 90))

# Graficos- Biplot

# Primera forma - usando plot.CA de FactoMineR
plot.CA(res.ca) # Mapa Simetrico
plot.CA(res.ca, axes = c(1,2), col.row = "blue", col.col = "red")
plot.CA(res.ca,mass=c(T,T))

# Segunda forma - usando fviz_ca_biplot de factoextra
library(factoextra)
fviz_ca_biplot(res.ca, repel = T)
fviz_ca_biplot(res.ca, repel = T) + theme_minimal()
fviz_ca_biplot(res.ca, repel = T) + theme_light()
fviz_ca_biplot(res.ca, repel = T) + theme_void()
fviz_ca_biplot(res.ca, repel = T) + theme_gray()

library(ca)
res.ca1 = ca(datos.acs)

plot(res.ca1) #simetrica
plot(res.ca1, map="rowprincipal") #asimetrica
plot(res.ca1, map="symbiplot")

#--------------------------------------------
# Interpretacion de los Indicadores del ACS

summary(res.ca1,nb.dec = 3, ncp = 2) 

# Interpretacion de la Contribucion Absoluta (ctr)
# Por ejemplo: para la fila H1 y la dimension 1 se tiene 
# una ctr = 21.205
# El 21.205% de la inercia de la dimension 1 es explicada 
# por la fila H1
# 
# Interpretacion de la Contribucion Relativa (cos2)
# Por ejemplo: para la fila H1 y la dimension 1 se tiene 
# una cos2 = 0.776
# El 77.6% de la inercia de la fila H1 es explicada 
# por la dimension 1

library(Factoshiny)
result=Factoshiny(res.ca)
Factoshiny(result)

#####################################################################
# 6. Analisis de Correspondencias Simples con puntos suplementarios #
#####################################################################

#             E1	E2	E3	E4	E5	E6	E7	E8	E9	Ideal
# Precios	    16	17	18	19	16	45	15	19	18	45
# Variedad	   8	15	18	17	27	20	 2	14	53	53
# Rapidez	    20	20	23	21	29	20	18	19	25	29
# Informacion	11	13	12	17	20	16	15	10	44	44
# Trato	      28	25	25	22	30	26	24	22	26	30
# Condiciones	21	21	20	24	27	22	18	21	24	27
# Acceso	    21	21	21	23	26	15	16	18	21	26


datos_s.acs <- matrix(c(16,17,18,19,16,45,15,19,18,45,
                       8,15,18,17,27,20, 2,14,53,53,
                      20,20,23,21,29,20,18,19,25,29,
                      11,13,12,17,20,16,15,10,44,44,
                      28,25,25,22,30,26,24,22,26,30,
                      21,21,20,24,27,22,18,21,24,27,
                      21,21,21,23,26,15,16,18,21,26),
                      nrow=7,byrow=T)

datos_s.acs

# Asignacion de nombres a las filas y columnas de la tabla

dimnames(datos_s.acs)<-list(Atributos=c("Precios", "Variedad", "Rapidez", 
                                  "Informacion","Trato","Condiciones","Acceso")
                          ,Empresa=c("Empresa 1","Empresa 2","Empresa 3",
                                     "Empresa 4","Empresa 5","Empresa 6",
                                     "Empresa 7","Empresa 8","Empresa 9",
                                     "Ideal"))

datos_s.acs
addmargins(datos_s.acs)

#---------------------------------------
# Prueba de Independencia Chi-Cuadrado  

chisq.test(datos_s.acs[,-10])
qchisq(0.05,48,lower.tail=FALSE)
library(fastGraph)
shadeDist( 65.17077, "dchisq", 48, lower.tail=FALSE,col = c("black", "steelblue") )
?shadeDist
#---------------------------------
# ACS con el paquete FactoMiner   

library(FactoMineR) 

library(factoextra)
res.ca2 = CA(datos_s.acs,
            ncp=2,
            graph=T,
            col.sup = 10)

#----------------------------------
# Scree Plot de los Autovalores
get_eigenvalue(res.ca2)


fviz_screeplot(res.ca2, addlabels = TRUE, ylim = c(0, 80))

#-------------------------------------------------
# Biplot filas, columnas y columna suplementaria
fviz_ca_biplot(res.ca2, repel = T) + theme_light()


summary(res.ca2,nb.dec = 2, ncp = 2) 
res.ca2$eig
sum(res.ca2$eig[,1])
#str(res.ca2)

library(Factoshiny)
result=Factoshiny(res.ca2)
Factoshiny(result)


