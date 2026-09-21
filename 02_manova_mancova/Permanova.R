
# Permutational Multivariate Analysis of Variance (Análisis de Variana 
# Multivariado por Permutaciones).

# El Permanova es un método no paramétrico alternativo al Manova y se basa en.
# distancias.

#  Permanova en DCA

library(vegan)

# Crear datos simulados: abundancias de especies en diferentes sitios y hábitats
set.seed(123)
species_data <- data.frame(
  Species1 = rpois(15, lambda = c(rep(5, 5), rep(10, 5), rep(15, 5))),
  Species2 = rpois(15, lambda = c(rep(3, 5), rep(8, 5), rep(12, 5))),
  Species3 = rpois(15, lambda = c(rep(2, 5), rep(6, 5), rep(9, 5)))
)
habitat <- factor(rep(c("Forest", "Grassland", "Wetland"), each = 5))

# Calcular la matriz de distancias (Bray-Curtis)
dist_matrix <- vegdist(species_data, method = "bray")

# Aplicar PERMANOVA con adonis2
permanova_result <- adonis2(dist_matrix ~ habitat, data = data.frame(habitat = habitat))

# Mostrar resultados
print(permanova_result)
# conclusión: Hay diferencias entre habitats.

# Visualización opcional con NMDS
ordination <- metaMDS(species_data, distance = "bray", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(habitat))
legend("topright", legend = levels(habitat), col = 1:3, pch = 19)

# Permanova en DBCA

library(vegan)

# Datos simulados
set.seed(123)
# Cuatro tratamientos, tres bloques, y dos variables de respuesta
data <- data.frame(
  Block = factor(rep(1:3, each = 4)),  # Bloques
  Treatment = factor(rep(1:4, times = 3)),  # Tratamientos
  Variable1 = c(rnorm(4, mean = 5, sd = 1), 
                rnorm(4, mean = 7, sd = 1), 
                rnorm(4, mean = 6, sd = 1)),  # Respuesta 1
  Variable2 = c(rnorm(4, mean = 3, sd = 0.5), 
                rnorm(4, mean = 4, sd = 0.5), 
                rnorm(4, mean = 3.5, sd = 0.5))  # Respuesta 2
)

# Matriz de datos multivariados
response_matrix <- data[, c("Variable1", "Variable2")]

# Calcular la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Treatment + Block, by = 'terms', 
                            data = data)

# Mostrar resultados
print(permanova_result)

# Visualización opcional
# Ordenación NMDS para visualizar las distancias
ordination <- metaMDS(response_matrix, distance = "euclidean", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(data$Treatment))
legend("topright", legend = levels(data$Treatment), col = 1:4, pch = 19)

# Permanova: Alternativa no paramétrica de un factorial en DCA.

# Instalar y cargar el paquete vegan
library(vegan)

# Crear datos simulados
set.seed(123)
data <- data.frame(
  FactorA = factor(rep(c("Baja", "Media", "Alta"), each = 6)),  # Factor A con 3 niveles
  FactorB = factor(rep(c("Presente", "Ausente"), times = 9)),   # Factor B con 2 niveles
  Variable1 = c(rnorm(6, mean = 5, sd = 1), 
                rnorm(6, mean = 6, sd = 1), 
                rnorm(6, mean = 7, sd = 1)),  # Respuesta 1
  Variable2 = c(rnorm(6, mean = 3, sd = 0.5), 
                rnorm(6, mean = 4, sd = 0.5), 
                rnorm(6, mean = 5, sd = 0.5))  # Respuesta 2
)

# Matriz de datos multivariados
response_matrix <- data[, c("Variable1", "Variable2")]

# Calcular la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar PERMANOVA para el diseño factorial
permanova_result <- adonis2(dist_matrix ~ FactorA * FactorB, by = 'terms', data = data)

# Mostrar resultados
print(permanova_result)

# Visualización opcional

# Ordenación NMDS para visualizar las distancias
ordination <- metaMDS(response_matrix, distance = "euclidean", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(data$FactorA))
text(ordination, display = "sites", labels = as.numeric(data$FactorB), cex = 0.8, pos = 3)
legend("topright", legend = levels(data$FactorA), col = 1:3, pch = 19, title = "FactorA")

# Permanova: Alternativa no paramétrica de un factorial en DBCA

library(vegan)

# Crear datos simulados
set.seed(123)
data <- data.frame(
  Block = factor(rep(1:3, each = 6)),  # Tres bloques
  FactorA = factor(rep(c("Baja", "Media", "Alta"), each = 2, times = 3)),  # Factor A (3 niveles)
  FactorB = factor(rep(c("Presente", "Ausente"), times = 9)),  # Factor B (2 niveles)
  Variable1 = c(rnorm(6, mean = 5, sd = 1), 
                rnorm(6, mean = 6, sd = 1), 
                rnorm(6, mean = 7, sd = 1)),  # Respuesta 1
  Variable2 = c(rnorm(6, mean = 3, sd = 0.5), 
                rnorm(6, mean = 4, sd = 0.5), 
                rnorm(6, mean = 5, sd = 0.5))  # Respuesta 2
)

# Matriz de datos multivariados
response_matrix <- data[, c("Variable1", "Variable2")]

# Calcular la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar PERMANOVA con bloques y factores
permanova_result <- adonis2(dist_matrix ~ Block + FactorA * FactorB, 
                            by = 'terms', data = data)

# Mostrar resultados
print(permanova_result)

# Visualización opcional

# Ordenación NMDS para visualizar las distancias
ordination <- metaMDS(response_matrix, distance = "euclidean", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(data$Block))
text(ordination, display = "sites", labels = interaction(data$FactorA, data$FactorB), cex = 0.8, pos = 3)
legend("topright", legend = levels(data$Block), col = 1:3, pch = 19, title = "Bloques")

# Permanova: Alternativa no paramétrica del Mancova en DCA

library(vegan)

# Crear datos simulados
set.seed(123)
data <- data.frame(
  Treatment = factor(rep(1:3, each = 10)),  # Tres tratamientos
  Covariate = runif(30, 50, 100),  # Covariable continua
  Variable1 = c(rnorm(10, mean = 5, sd = 1), 
                rnorm(10, mean = 7, sd = 1), 
                rnorm(10, mean = 6, sd = 1)),  # Respuesta 1
  Variable2 = c(rnorm(10, mean = 3, sd = 0.5), 
                rnorm(10, mean = 4, sd = 0.5), 
                rnorm(10, mean = 3.5, sd = 0.5))  # Respuesta 2
)

# Matriz de datos multivariados
response_matrix <- data[, c("Variable1", "Variable2")]

# Calcular la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar PERMANOVA con covariable
permanova_result <- adonis2(dist_matrix ~ Covariate + Treatment,
                            by = 'terms', data = data)

# Mostrar resultados
print(permanova_result)

# Visualización opcional

# Ordenación NMDS para visualizar las distancias
ordination <- metaMDS(response_matrix, distance = "euclidean", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(data$Treatment))
legend("topright", legend = levels(data$Treatment), col = 1:3, pch = 19)

# Permanova: Alternativa no paramétrica del Mancova en DBCA

library(vegan)

# Crear datos simulados
set.seed(123)
data <- data.frame(
  Block = factor(rep(1:3, each = 9)),  # Tres bloques
  Treatment = factor(rep(1:3, times = 9)),  # Tres tratamientos
  Covariate = runif(27, 50, 100),  # Covariable continua
  Variable1 = c(rnorm(9, mean = 5, sd = 1), 
                rnorm(9, mean = 6, sd = 1), 
                rnorm(9, mean = 7, sd = 1)),  # Respuesta 1
  Variable2 = c(rnorm(9, mean = 3, sd = 0.5), 
                rnorm(9, mean = 4, sd = 0.5), 
                rnorm(9, mean = 5, sd = 0.5))  # Respuesta 2
)

# Matriz de datos multivariados
response_matrix <- data[, c("Variable1", "Variable2")]

# Calcular la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar PERMANOVA con bloques y covariable
permanova_result <- adonis2(dist_matrix ~ Block + Covariate + Treatment,
                            by = 'terms',data = data)

# Mostrar resultados
print(permanova_result)

# Visualización

# Ordenación NMDS
ordination <- metaMDS(response_matrix, distance = "euclidean", k = 2, trymax = 50)
plot(ordination, type = "n")
points(ordination, display = "sites", pch = 19, col = as.numeric(data$Block))
text(ordination, display = "sites", labels = as.numeric(data$Treatment), cex = 0.8, pos = 3)
legend("topright", legend = levels(data$Block), col = 1:3, pch = 19)

## Alternativa no paramétrica a la Regresión Multivariada.
## Regresión Multivariada basada en Distancias (o Distance-Based
## Redundancy Analysis, db-RDA).

# Objetivo: Evaluar si las variables independientes explican las variaciones
# en las respuestas multivariadas.

library(vegan)

# Crear datos simulados
set.seed(123)
data <- data.frame(
  Temp = runif(30, 15, 25),         # Variable independiente 1 (temperatura)
  Precip = runif(30, 100, 200),     # Variable independiente 2 (precipitación)
  Resp1 = rnorm(30, mean = 50, sd = 10),  # Variable dependiente 1
  Resp2 = rnorm(30, mean = 30, sd = 5)    # Variable dependiente 2
)

# Matriz de variables dependientes
response_matrix <- data[, c("Resp1", "Resp2")]

# Crear la matriz de distancias (Euclidiana)
dist_matrix <- vegdist(response_matrix, method = "euclidean")

# Aplicar db-RDA
db_rda_result <- capscale(dist_matrix ~ Temp + Precip, data = data)

# Resumen de resultados
summary(db_rda_result)

# Prueba de significancia para los predictores
anova(db_rda_result, permutations = 999)
anova(db_rda_result, permutations = 999, by = 'terms')

# Graficar la ordenación de db-RDA
plot(db_rda_result, main = "db-RDA", scaling = 2)

# Incluir interacciones o covariables

db_rda_result <- capscale(dist_matrix ~ Temp * Precip, data = data)
summary(db_rda_result)
anova(db_rda_result, permutations = 999)
anova(db_rda_result, permutations = 999, by = 'terms')

# Probar diferentes métricas de distancia

dist_matrix <- vegdist(response_matrix, method = "bray")

# Aplicar db-RDA
db_rda_result <- capscale(dist_matrix ~ Temp * Precip, data = data)

# Resumen de resultados
summary(db_rda_result)

# Prueba de significancia para los predictores
anova(db_rda_result, permutations = 999)
anova(db_rda_result, permutations = 999, by = 'terms')

# Graficar la ordenación de db-RDA
plot(db_rda_result, main = "db-RDA", scaling = 2)




