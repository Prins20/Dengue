
## A.*Pipeline de clasificación de archivos.gz por códigos de dengue en R*
### Establecer directorio de trabajo
setwd("/home/sec-dengue")

### Mostrar contenidos iniciales
cat("Archivos en la carpeta 'prueba':\n")
print(list.files("prueba"))

### Dentro de la carpeta pruebas, codigos.txts dengue
contenido1 <- readLines("prueba/dengue1.txt")
contenido2 <- readLines("prueba/dengue2.txt")
contenido3 <- readLines("prueba/dengue3.txt")
length(contenido1)
length(contenido2)
length(contenido3)

### Mostrar códigos cargados
cat("===== dengue1.txt =====\n")
print(contenido1)
cat("===== dengue2.txt =====\n")
print(contenido2)
cat("===== dengue3.txt =====\n")
print(contenido3)

### Crear carpetas 
dir.create("prueba/dengue1", showWarnings = FALSE)
dir.create("prueba/dengue2", showWarnings = FALSE)
dir.create("prueba/dengue3", showWarnings = FALSE)

### Buscar archivos .gz en carpeta prueba/
archivos <- list.files("prueba", pattern = "\\.gz$", full.names = TRUE)

### Contadores
cont1 <- 0
cont2 <- 0
cont3 <- 0

### copiar y separar cada archivo
for (archivo in archivos) {
  nombre <- basename(archivo)
  codigo <- substr(nombre, 1, 10)  # primeros 10 caracteres del nombre

  if (codigo %in% contenido1) {
    file.copy(archivo, file.path("prueba/dengue1", nombre), overwrite = TRUE)
    cont1 <- cont1 + 1
  } else if (codigo %in% contenido2) {
    file.copy(archivo, file.path("prueba/dengue2", nombre), overwrite = TRUE)
    cont2 <- cont2 + 1
  } else if (codigo %in% contenido3) {
    file.copy(archivo, file.path("prueba/dengue3", nombre), overwrite = TRUE)
    cont3 <- cont3 + 1
  }
}

cat("Archivos copiados a 'dengue1':", cont1, "\n")
cat("Archivos copiados a 'dengue2':", cont2, "\n")
cat("Archivos copiados a 'dengue3':", cont3, "\n")
