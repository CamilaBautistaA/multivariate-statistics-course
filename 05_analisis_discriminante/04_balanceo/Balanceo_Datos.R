
### BALANCEO DE CLASES DESBALANCEADAS

## Se considera que las clases estan desbalanceadas cuando la clase
## minoritaria tiene un porcentaje de 10% o menos.

## NOTA IMPORTANTE: El balanceo se hace con fines predictivos. Nunca 
## se debe usar con fines descriptivos porque los datos han cambiado.

## El balanceo se hace con la data de entrenamiento.

library(caret)

library(foreign)
datos <-read.spss("Balanceo_Datos.sav",
                  use.value.labels=TRUE, 
                  to.data.frame=TRUE)
head(datos)
str(datos)
dim(datos)
predictoras=datos[,-c(1,6)]
dependiente=datos[,6]
head(predictoras)
head(dependiente)
table(dependiente)
prop.table(table(dependiente))

## Metodos para tratar datos desbalanceados: Undersampling, Oversampling,
## Datos Sinteticos (Smote y Rose y variaciones del Smote y Rose). 
## Tambien se pueden usar algoritmos robustos que trabajan con datos
## desbalanceados como: Random Forest, Maquinas de Soporte Vectorial, 
## algunos algoritmos de Redes Neuronales, Modelos con Enlaces 
## Asimetricos.

## UNDERSAMPLING: Se obtienen clases de igual tamano de acuerdo a la clase
## minoritaria

?downSample

MuestraDescendente=downSample(predictoras, dependiente)
MuestraDescendente 
dependiente1=MuestraDescendente[,6]
table(dependiente1)
prop.table(table(dependiente1))

## OVERSAMPLING: Se obtienen clases de igual tamano de acuerdo a la clase
## mayoritaria

MuestraAscendente=upSample(predictoras, dependiente)
MuestraAscendente
dependiente2=MuestraAscendente[,6]
table(dependiente2)
prop.table(table(dependiente2))









