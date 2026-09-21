

# Comparacion de algunas Pruebas de Normalidad Multivariada (Shapiro, Royston y 
# Henze-Zirkler)

library(nortest)

CPN3<-function(data,n,r){
  library(mvnormtest)
  library(MVN)
  pval.shapiro<-NULL
  pval.royston<-NULL
  pval.henze<-NULL
  for(i in 1:r){
    datos<-data[sample(nrow(data),size = n,replace=T),]
    pval.shapiro[i]<-mshapiro.test(t(datos))$p.value
    pval.royston[i]<-mvn(datos,mvnTest = 'royston')$multivariateNormality$`p value`
    pval.henze[i]<-mvn(datos,mvnTest = "hz")$multivariateNormality$`p value`
  }
  p.shapiro<-mean(pval.shapiro<0.05)*100
  p.royston<-mean(pval.royston<0.05)*100
  p.henze<-mean(pval.henze<0.05)*100
  return(list(p.shapiro=p.shapiro,p.royston=p.royston,p.henze=p.henze))
}

# Con la data del ejemplo 2 de inferencia multivariada funciona con muestras
# de tamano 14 para adelante. Con valores menores resultan matrices singulares.
# Con el ejemplo 1 si funciona para cualquier tamano de muestra.

datos=read.delim("clipboard",T)
is.data.frame(datos)
datos
CPN3(datos,32,2000)
CPN3(datos,200,2000)
