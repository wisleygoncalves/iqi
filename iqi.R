################################## IQI ######################################

# Carregando Pacotes

if(!require(pacman))
  install.packages("pacman")
library(pacman)

pacman::p_load("ggplot2", 'dplyr',  "data.table", "raster", "gridExtra")

################################# Funções ######################################

# Banco de Dados
db = function(files) {
  
  db_drone = list()
  
  for(file in files) {
    raster_drone = paste(file, '.tif', sep = '')
    data_drone = raster(raster_drone)
    x = values(data_drone)
    x[is.na(x)] = 0
    x = x[which(x != 0)]
    x = ((x / 255) - (min(x) / 255)) / ((max(x) / 255) - (min(x) / 255))
    db_drone[[file]] = x
    glimpse(x)
  }
  
  return(db_drone)
}

# Indice de Qualidade de Imagem
IQI = function(x, y) {
  
  if (length(x) > length(y)) {
    x = x[1:length(y)]
  } else if (length(y) > length(x)) {
    y = y[1:length(x)]
  }
  
  cc = cov(x,y) / (sd(x) * sd(y)) # Coeficiente de Correlação
  vb = (2 * mean(x) * mean(y)) / ((mean(x) ^ 2) + (mean(y) ^ 2)) # variação do brilho
  sdc = (2 * sd(x) * sd(y)) / (var(x) + var(y)) # similaridade do contraste
  IQI = cc * vb * sdc # IQI
  Med1 = mean(x)
  Med2 = mean(y)
  Var1 = var(x)
  Var2 = var(y)
  
  data = data.frame(Variavel = c('CC', 'VB','SDC', 'IQI', 'Med1', 'Med2', 'Var1',
                                 'Var2'), ALT = c(cc, vb, sdc, IQI, Med1, Med2,
                                                  Var1, Var2))
  
  print(data)
  
  return(data)
}

# Histograma
gf = function(data, a) {
  
  hist(data, breaks = 50, col = "lightgray", border = 'black', freq = FALSE,
       main = "", xlab = "", ylab = "", axes = FALSE)
  
  axis(1, at = pretty(data, n = 10), cex.axis = 1.2, font.axis = 2)
  axis(2, las = 1, cex.axis = 1.2, font.axis = 2)
  
  title(xlab = a, ylab = "Frequência", line = 3, font.lab = 2, cex.lab = 1.5)
}

######################### Teste de Comparação ###############################

# Estresse Hídrico

DB = db(c('exr5', 'exr10', 'exr15', 'vdvi5', 'vdvi10', 'vdvi15'))

## Histograma

titles = c('ExR - 5m', 'ExR - 10m', 'ExR - 15m', 'VDVI - 5m', 'VDVI - 10m',
           'VDVI - 15m')

windows()
layout(matrix(c(1:6), nrow = 2, byrow = TRUE))
for (i in 1:length(titles)) {
  gf(DB[[i]], titles[i])
}
dev.off()

## Análise de IQI
IQI_5 = IQI(DB[[1]], DB[[4]]) # Altura de 5 m

IQI_10 = IQI(DB[[2]], DB[[5]]) # Altura de 10 m

IQI_15 = IQI(DB[[3]], DB[[6]]) # Altura de 15 m

# Tabela - IQI
iqi = data.frame(Variavel = IQI_5$Variavel, ALT5 = IQI_5$ALT, ALT10 = IQI_10$ALT,
                 ALT15 = IQI_15$ALT)
write.csv(iqi, file = 'iqi_estresse_hidrico.csv') # Salvando Arquivo

rm(DB, IQI_5, IQI_10, IQI_15, iqi)

# Teor de Clorofila

DB = db(c('ngrdi5', 'ngrdi10', 'ngrdi15', 'tgi5', 'tgi10', 'tgi15'))

## Histograma

titles = c('NGRDI - 5m', 'NGRDI - 10m', 'NGRDI - 15m', 'TGI - 5m', 'TGI - 10m',
           'TGI - 15m')

windows()
layout(matrix(c(1:6), nrow = 2, byrow = TRUE))
for (i in 1:length(titles)) {
  gf(DB[[i]], titles[i])
}
dev.off()

## Análise de IQI
IQI_5 = IQI(DB[[1]], DB[[4]]) # Altura de 5 m

IQI_10 = IQI(DB[[2]], DB[[5]]) # Altura de 10 m

IQI_15 = IQI(DB[[3]], DB[[6]]) # Altura de 15 m

# Tabela - IQI
iqi = data.frame(Variavel = IQI_5$Variavel, ALT5 = IQI_5$ALT, ALT10 = IQI_10$ALT,
                 ALT15 = IQI_15$ALT)
write.csv(iqi, file = 'iqi_teor_de_clorofila.csv') # Salvando Arquivo

rm(DB, IQI_5, IQI_10, IQI_15, iqi)

# Cobertura do Solo

DB = db(c('exg5', 'exg10', 'exg15', 'gli5', 'gli10', 'gli15'))

## Histograma

titles = c('ExG - 5m', 'ExG - 10m', 'ExG - 15m', 'GLI - 5m', 'GLI - 10m',
           'GLI - 15m')

windows()
layout(matrix(c(1:6), nrow = 2, byrow = TRUE))
for (i in 1:length(titles)) {
  gf(DB[[i]], titles[i])
}
dev.off()

## Análise de IQI
IQI_5 = IQI(DB[[1]], DB[[4]]) # Altura de 5 m

IQI_10 = IQI(DB[[2]], DB[[5]]) # Altura de 10 m

IQI_15 = IQI(DB[[3]], DB[[6]]) # Altura de 15 m

# Tabela - IQI
iqi = data.frame(Variavel = IQI_5$Variavel, ALT5 = IQI_5$ALT, ALT10 = IQI_10$ALT,
                 ALT15 = IQI_15$ALT)
write.csv(iqi, file = 'iqi_cobertura_solo.csv') # Salvando Arquivo

rm(DB, IQI_5, IQI_10, IQI_15, iqi)

# Biomassa Verde

DB = db(c('vari5', 'vari10', 'vari15', 'bgi5', 'bgi10', 'bgi15'))

## Histograma

titles = c('VARI - 5m', 'VARI - 10m', 'VARI - 15m', 'BGI - 5m', 'BGI - 10m',
           'BGI - 15m')

windows()
layout(matrix(c(1:6), nrow = 2, byrow = TRUE))
for (i in 1:length(titles)) {
  gf(DB[[i]], titles[i])
}
dev.off()

## Análise de IQI
IQI_5 = IQI(DB[[1]], DB[[4]]) # Altura de 5 m

IQI_10 = IQI(DB[[2]], DB[[5]]) # Altura de 10 m

IQI_15 = IQI(DB[[3]], DB[[6]]) # Altura de 15 m

# Tabela - IQI
iqi = data.frame(Variavel = IQI_5$Variavel, ALT5 = IQI_5$ALT, ALT10 = IQI_10$ALT,
                 ALT15 = IQI_15$ALT)
write.csv(iqi, file = 'iqi_biomassa_verde.csv') # Salvando Arquivo

rm(DB, IQI_5, IQI_10, IQI_15, iqi)
