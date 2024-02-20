###############################################################################
#                                                                             #
#----------------------- A B I E S   R E L I G I O S A -----------------------#
#--------------- F O R E S T   G O V E R N A N C E   M O D E L ---------------#
#                                                                             #
#                                                                             #
# Description: 
#
# Authors: , 2023
#                                                                             #
###############################################################################



###############################################################################
#                        F O R E S T   D Y N A M I C S                        #
#                                                                             #
# Description: Coupled, dynamic, model that models forest population dynamics #
#              in terms of growth, regeneration and mortality. This is        #
#              considered an hybrid model because it can simulate processes   #
#              as well as states                                              #
#                                                                             #
###############################################################################

###############################################################################
############# P R E P A R A T I O N   A N D   P R E P R O C E S S #############
###############################################################################

#---------------- P A C K A G E S   A N D   L I B R A R I E S ----------------#

library("dplyr")

#----------------- S E T   W O R K I N G   D I R E C T O R Y -----------------#

dir <- ('C:/Users/luisv/Documents/Academico/01_Maestria/00_Proyecto de Investigación/01_Proyecto/04_Modelo/')

#------------------- L O A D   A N D   C L E A N   D A T A -------------------#

rodales_area <- read.csv(file=paste0(dir,"Abies_R/_data/rodales.csv"), header=T)
trees <- read.csv(file=paste0(dir,"Abies_R/_data/numero_arboles_diametro.csv"), header=T)
colnames(trees) <- c("RODAL","DBH_CLASS","NUMBER")

#----------------- C O U N T   N U M B E R   O F   T R E E S -----------------#
# trees <- trees %>%
#   group_by(RODAL, DBH_CLASS) %>%
#   mutate(NUMBER = n())

#------------------ A G G R E G A T E   D A T A   F R A M E ------------------#
trees <- aggregate(x = trees, by = list(trees$RODAL, trees$DBH_CLASS, trees$NUMBER), FUN = SUM)
trees <- trees[,-c(1,2,3)]
trees <- trees[order(trees$RODAL, trees$DBH_CLASS),]

###############################################################################
############################## F U N C T I O N S ##############################
###############################################################################

################## S A M P L E   A   S I N G L E   S T A N D ##################

subpop <- function(d = rodales_area) {
  res <-  d[sample(1:nrow(d), 1), ]
  numberrodal <- res[1]
  arearodal <- res[2]
  return(list(numberrodal=numberrodal,arearodal=arearodal))
  }

############ O B T A I N S   T H E   N U M B E R   O F   T R E E S ############
######## A N D   T H E I R   D B H   F O R   S A M P L E D   S T A N D ########

GetTrees <- function(T, Plots) {
  res <- NULL
  for (i in Plots) {
    res <- rbind(res, subset(T, RODAL %in% i)) #Aquí se cambio "Rodal" por "RODAL" porque así está en trees
  }
  return(res)
}

# #################################################################### !!!!!!!!!!!!!!! VERIFICAR QUE ESTO SÍ VA AQUÍ !!!!!!!!!!!!!!!!
# ### Samples a single stand from the data 
# Sampled <- subpop ()
# ### Gets the stand number
# Plots <- Sampled$numberrodal

######################## G E N E R A T E S   S T A N D ########################
################### P O P U L A T I O N   S T R U C T U R E ###################

process<- function (trees,rodales_area){
  ### Get vectors for number of trees (stand) and DBH class (class)
  Sampled <- subpop ()
  Plots <- Sampled$numberrodal
  Tree.List <- GetTrees (trees,Plots)
  stand <- Tree.List$NUMBER 
  class <- Tree.List$DBH_CLASS
  ### 
  basal_big_class <- 0.424455347
  BAB <- rep(0,100)
  ### Partition the basal area of big trees >31 cm and add number of trees that the surplus of basal area represents
  ### Redistributes the basal area of DBH classes higher than 75 
  TBA <- 3.142*(Tree.List[Tree.List[,2]>75,2]/200)^2
  BAB <- round(TBA/basal_big_class,digits=0)
  y <- sum(BAB)
  ### Places the number of threes along a 15 length vector ???
  rodal <- c(0,stand[1:14])
  rodal[15] <- rodal[15]+y
  rodal[1] <- 30
  res <- rodal
  return (res)
}   

############## G E T S   H E I G H T   F R O M   D I A M E T E R ##############

### From Guzman Santiago et al. 2019. 
height_fir <- function (diametro){ 
 bo <-119.8798
 b1 <-154.2599
 alt <- 1.3+ (bo*(diametro/(b1+diametro))) 
 alt
}

################ G E T S   V O L U M E   F R O M   H E I G H T ################
######################### A N D   B A S A L   A R E A ######################### !!!!!!!!!!!!! FUENTE

volumen <- function (area_basal, alturas){ ##volume equation###suponiendo geometr?a aproximada c?nico 
  v <- 0.42*area_basal*alturas
  v
}

# ################# I N T E R P O L A T E S   E X P E C T E D ###################   !!!!!!!!!!!!!!!!!!!!!!!!!!
# ######### P R O D U C T I V I T Y   F R O M   K N O W N   V A L U E S #########

# ### Create data frame with known values
# knownproductivity <- data.frame(x <- c(30,150,300),
#                                 y <- c(0.1,0.5,1))
# ### Interpolate
# expected_productivity <- approx(knownproductivity$x, knownproductivity$y)
# ### Plot
# plot(expected_productivity)
# ### Exctract productivity
# prod<- expected_productivity$x ### DÓNDE SE USA ESTA????
# ### Extract effect
# efecto <- expected_productivity$y

############### T R A N S I T I O N   P R O B A B I L I T I E S ###############

### Empiric transition probability = 0.004 
exe <- function(Y,trees,rodales_area,prob_trans) {
  
  #------- G E N E R A T E S   P O P U L A T I O N   S T R U C T U R E -------#
  rodal<-process(trees,rodales_area) 
  n <- length(rodal)
  N1s <- rep(1, n)
  N0s <- rep(0, n)
  ### Sets upper and lower limits for DBH classes
  dbhl <- c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71)
  dbhu <- c(6,11,16,21,26,31,36,41,46,51,56,61,66,71,76)
  ### Gets quadratic DBH in order to assign a range value to every DBH class
  dbhq <- sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl) * 3))
  ### Gets quadratic basal area
  baq  <- (dbhq/2)^2 * pi/1e4
  ### Gets height mean per DBH class [meters]
  alturas <- height_fir(dbhq)
  ### Gets 5 centimeter DBH classes
  iw <- dbhu-dbhl 
  
  #------------- C R E A T E S   R E S U L T S   M A T R I C E S -------------#
  Rodal <- matrix(0, Y, n, byrow = TRUE)
  Recruits <- numeric(Y)
  BA <- numeric(Y)
  Senescence <- matrix(0, Y, n, byrow = TRUE)
  Crecimiento <- matrix(0, Y, n, byrow = TRUE)
  AnnualGrowth <-  matrix(0, Y, n, byrow = TRUE)
  mortalidades <- matrix(0, Y, n, byrow = TRUE)
  Productividad <- numeric(Y)
  
  #---------------- I N T E R P O L A T E S   E X P E C T E D ----------------#
  #--------------- T R A N S I T I O N   P R O B A B I L I T Y ---------------#
  ### Transition probabilities from Fir demography paper (Forest Ecology and Management) !!!!!!!!!!!!!! FUENTE
  knownpoints_t <- data.frame(x <- c(1,3,5,10,20,25,30,35,40,45,50,60,100),
                            y <- c(0.04,0.05,0.05,0.03,0.01,0.01,0.1,0.1,0.1,0.11,0.11,0.11,0.11))
  ### Survival probabilities 
  knownpoints <- data.frame(x <- c(1,3,5,10,20,25,30,35,40,45,50,60,100),
                            y <- c(0.6162,0.873,0.873,0.874,0.874,0.740,0.740,0.756,0.756,0.84,0.84,0.85,0.85))
  ### Survival vs diameter interpolation
  interpolation <- approx(knownpoints$x, knownpoints$y)
  diameter <- interpolation$x
  survival_diameter <- interpolation$y
  #plot(interpolation, xlab = "Diameter [cm]", ylab ="Survival probability")
  ### Annual diameter interpolation
  interpolation_t <- approx(knownpoints_t$x, knownpoints_t$y)
  diameter2 <- interpolation_t$x
  annualdiameter <- interpolation_t$y ##proveniente de las interopolaciones...constante
  #plot(interpolation_t, xlab = "Diameter [cm]", ylab = "Transition probability")
  
  #RegenLagTime <- 30   # assume it takes a 1-yr old natural seedling 30 years to grow to
  # a sapling with a DBH 2.0 cm (first diameter class) data from Jos? L. Gallardo-Salazar1
  #, Dante A. Rodr?guez-Trejo1*, Salvador Castro-Zavala. 2019
  ###Mortality Temi?o-Vi?ota et al. 2016. Modelling initial mortality of Abies religiosa
  ###produce altas mortalidades
  #a<- -(1.7033-0.0185*dbhq)
  #pmort <- 1/(1+exp(a))
  ####MORTALIDAD
  ###probabilidad de mortalidad de Abies (Vieilledent et al. 2009). 
 # knownmortality <- data.frame(x <- c(2,5,15,45,75,100),
                          #  y <- c(0.05,0.043,0.0197,0.0029,0.0009,0.0006))
  #yos_mort <- approx( knownmortality$x,  knownmortality$y)
  #diameter_mort <- yos_mort$x
  #mortality_prob <- yos_mort$y 
  #c1 <- 0.536#1-20 
  #c2 <- 0.157#20-40
  #c3 <- 0.371#>40 
  #probabilidad_mortalidad_1_20 <- runif(24000, c1-0.02,c1+0.02)
  #probabilidad_mortalidad_20_40 <- runif(24000, c2-0.02,c2+0.02)
  #probabilidad_mortalidad_40 <- runif(24000, c3-0.02,c3+0.02)
 
  #---- C A L C U L A T E S   T R A N S I T I O N   P R O B A B I L I T Y ----#
  for (y in 1:Y) {
    #714
    ### Randomly generates between 1000 and 1600 seedlings
    semillas <- floor(runif(2000,1000,1600)) #(#hacer ejercicio calibracion, 1500 a 2000 pl?ntulas por hectarea; Guzm?n-Aguilar et al. 2020)
    ### Empiric transition probability, from ************             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUENTE
    prob_trans <- 0.0044
    ### Gets number of seedlings that transition according to random sample
    plantulas <- round((sample(semillas,1,replace=TRUE))* prob_trans)
    ### Survival probabilities per DBH class
    surv <- rep(0,n-1)
    diame <- rep(0,n-1)
    
    for (i in 1:length(surv)){
      ### Gets survival diameters 
      vector <- which(diameter >= dbhq[i] & diameter < dbhq[i+1])
      randomvector <- sample(vector,1)
      tmp1 <- diameter[randomvector]
      tmp <- survival_diameter[randomvector] 
      surv[i] <- tmp
      diame[i] <- tmp1
    }
    
    
    vector15 <- which(diameter >= floor(dbhq[n-1]) & diameter < floor(dbhq[n]))
    surv15 <- sample(survival_diameter[vector15],1) 
    diam15 <- sample(diameter[vector15],1) 
    surv[15]<- surv15
    
    ### Obtains transition probabilities per DBH class
    anndiam <- rep(0,n-1)
    diame_tmp <- rep(0,n-1)
    for (i in 1:length(anndiam)){
      vector <- which(diameter2 >= dbhq[i] & diameter2 < dbhq[i+1])
      randomvector <- sample(vector,1)
      tmp1 <- diameter2[randomvector]
      tmp <-   annualdiameter[randomvector] 
      anndiam[i] <- tmp
      diame_tmp[i] <- tmp1
    }
    
    vector15 <- which(diameter2>=floor(dbhq[n-1]) & diameter2<floor(dbhq[n]))
    ann15 <- sample(annualdiameter[vector15],1) 
    diame_tmp15 <- sample(diameter2[vector15],1) 
    anndiam[15]<- ann15
    
    sobrevivencia <- surv
    adi <- anndiam
    graduating <- adi
    tmp_3 <-  sobrevivencia+adi
 
    ### Get number of rees that pass to the next DBH class
    survadi <- rep(0,length(tmp_3))
    for (x in 1:length(tmp_3)){
    if (tmp_3[x]>1){      
      survadi[x] <- 1
    }else{
     survadi[x] <-tmp_3[x]
    }
  }
    prob_mort <- 1-survadi
    
    rodal_crecimiento <- rep(0,length(rodal))
    for (x in 1:length(rodal_crecimiento)){
      if (rodal[x]==0){      
        rodal_crecimiento[x] <- 0
      }else{
        rodal_crecimiento[x] <- rbinom(1, rodal[x], graduating[x])
      }
    }
    
    ### Updates the structure according to graduation 
    rodal <- rodal - rodal_crecimiento # Resta de los arboles que se "graduan"
    rodal <- rodal + c(0, rodal_crecimiento[1:n-1]) # Suma de los arboles "graduados"
    LastDCT <- rodal_crecimiento[15]  # prevents loosing trees in last diameter class
    rodal [15] <- rodal[15] + LastDCT
   
    BAIncremento <- sum(rodal_crecimiento * baq)    # patch BA increment due to growth

    
    # ##Reclutamiento
    #paper"Plantaciones forestales vs. regeneración natural in situ: 
    #El caso de los pinos y la rehabilitación en el Parque Nacional Cofre de Perote"
    ##Se estimaron de 1200 a 1600 plántulas por ha
    
    
    ##Mortalidad (size-dependent mortality rates)
    ##prob_mort <- rep(0,n-1)
    ##diame_mort <- rep(0,n-1)
    ##for (i in 1:length(prob_mort)){
    ## vector <- which(diameter_mort >= dbhq[i] &  diameter_mort<dbhq[i+1])
    ##tmp <- sample(mortality_prob[vector],1) 
    ##tmp1 <- sample( diameter_mort[vector],1) 
      #randomanndiam <- sample(runif(300, tmp-0.02,tmp+0.02),1)##aleatoriedad den el modelo
    ##prob_mort[i] <- tmp
    ##diame_mort[i] <- tmp1
    ##}
    
    ##vector_mort20 <- which(diameter_mort>=floor(dbhq[n-1]) & diameter_mort<floor(dbhq[n]))
    ##mort20 <- sample(mortality_prob[vector_mort20],1) 
    ##diam_mort20 <- sample(diameter_mort[vector_mort20],1) 
    ##prob_mort[n]<- mort20
    
    ### Calculates mortality
    rodal_mortalidad <- rep(0,length(rodal))
    for (x in 1:length(rodal_mortalidad)){
      if (rodal[x]==0){      
        rodal_mortalidad[x] <- 0
      }else{
        rodal_mortalidad[x] <- rbinom(1, rodal[x], prob_mort[x]) # Arboles que mueren
      }
    }
    
    ### Updates stand according to mortality results
    rodal <- rodal-rodal_mortalidad #solo por senescencia
    
    ### Adds recruitments
    #RegenCohorts <- rpois(RegenLagTime, sample(plantulas, 1, replace = TRUE))
    Recruits[y] <- plantulas
    #RegenCohorts <- c(NewRegen, RegenCohorts[1:RegenLagTime - 1])##pensar en esto cuando la regenarcion artificial
    rodal[1] <- rodal[1] + Recruits[y]
    
    abasal_clase <- rodal*baq #m2/ha
    productividad <- sum(volumen(abasal_clase[3:15],alturas[3:15])) ##m3/ha
      
    ### Saves results in results matrices
    Rodal[y, ] <- rodal # stand after regeneration, captures regeneration pulses
    Senescence[y, ] <- rodal_mortalidad
    Crecimiento[y, ] <-   rodal_crecimiento
    BA[y] <- BAIncremento
    AnnualGrowth[y, ] <-  adi
    mortalidades[y, ] <- prob_mort
    Productividad[y] <- productividad
    
  }
  
  ### Results
  res <- list(Rodal=Rodal, Productividad=Productividad, AnnualGrowth=AnnualGrowth,mortalidades=mortalidades,  BA=BA, Senescence=Senescence,Crecimiento=Crecimiento, Recruits=Recruits)
  return(res)
}

###############################################################################
################## T R A N S I T I O N   E X P E R I M E N T S ################
###############################################################################

####primer experimento con proba_trans=0.06
set.seed(12)
### 
exp <- exe(1000,trees,rodales_area,0.06) 

###############################################################################
############################ S I M U L A T I O N S ############################
###############################################################################

#---------- N U M B E R   O F  S T A N D S   T O   E V A L U A T E -----------#
n.iter <- 1000
#----------- N U M B E R   O F   Y E A R S   T O   E V A L U A T E -----------#
Y <- 100 
#------------------- P R O D U C T I V I T Y   M A T R I X -------------------#
Productividad_s <- matrix(0, n.iter, Y, byrow = T) ## Each line is a plot

#Productividad_s_efecto <- numeric(n.iter)
#Trans_lineal_prod <- numeric(n.iter)

#----------------------- S I M U L A T I O N   L O O P -----------------------#
for (i in 1:n.iter) {
  Forest_model <- exe(Y,trees,rodales_area,0.005)

  Productividad_s[i, ] <- Forest_model$Productividad
 # Productividad_s_efecto[i] <- Productividad_s[i, 25]
  print(i)
}

##Escenarios
#1: 0.005
#2: 0.9
#3:0.006
#4: 0.01
#5: 0.0001
#6: 0.5
#7:0.3
#8:0
#9:0.001
#10: 0.05

prod_1<- colMeans(Productividad_s)
prod_2<- colMeans(Productividad_s)
prod_3<- colMeans(Productividad_s)
prod_4<- colMeans(Productividad_s)
prod_5<- colMeans(Productividad_s)
prod_6<- colMeans(Productividad_s)
prod_7<- colMeans(Productividad_s)
prod_8<- colMeans(Productividad_s)
prod_9 <- colMeans(Productividad_s)
prod_10 <- colMeans(Productividad_s)

experimentos <- as.data.frame(c(prod_1, prod_2, prod_3, prod_4, prod_5, prod_6, prod_7,prod_8,prod_9, prod_10))
colnames(experimentos)<-"productividad"
exp <-rep(c("0.005","0.9","0.006","0.01","0.0001","0.5","0.3", "0", "0.001","0.05"), each=100)
experimentos$tipo <- as.factor(exp)
simu <-rep(seq(1,100),10)
experimentos$tiempo <- simu
write.csv(experimentos, file = paste0(dir,"Abies_R/_results/experimentos_gobernanza.csv"))

#--------------------- P L O T   P R O D U C T I V I T Y ---------------------#

#Borrar siguiente linea
experimentos <- read.csv("C:/Users/luisv/Documents/Academico/01_Maestria/00_Proyecto de Investigación/01_Proyecto/07_Resultados/Modelo Forestal/productividad.csv")

library(ggplot2)
# Basic line plot with points
ggplot(data = experimentos, aes(x = tiempo, y = productividad)) +
  geom_line(aes(colour = tipo)) +
  geom_point(aes(colour = tipo)) +
  scale_colour_gradient2(low = "magenta", high = "red", mid = "turquoise", midpoint = 0.5)



##obtener vectores y transformarlos funci{on de valor lineal
vectoresprod_cond_normales <- (subset(Productividad_s,Productividad_s[,1]<=300))/max(prod_1)
write.csv(vectoresprod_cond_normales, "vectoresprod_cond_normales.csv")



#### D E R I V A T E ####


derivada1 <- (-vectores_cond_normales[,1]*log(vectores_cond_normales[,1]))*0.1 ##derivadas del valor inicial
derivadas_all <- matrix(0,length(vectores_cond_normales[,1]), 25, byrow = TRUE)


derivadas_all <- derivadas_all / 600
derivadas_all <- derivadas_all*(log(derivadas_all))
derivadas_all <- (-(derivadas_all)*(log(derivadas_all)))*0.1


for (i in length(vectores_cond_normales[,1]-1)){
#=(-B2*LN(B2))*SIGNO(B2-B1)*0.1 
derivadas_all[,i+1] <- (-vectores_cond_normales[,i+1]*log(vectores_cond_normales[,i+1])) 
*sign(vectores_cond_normales[,i+1]-vectores_cond_normales[,i]) *0.1
  
  
}

Productivity <- mean(Productividad_s[1:n.iter, 25])

###tmp1 <- sample( diameter_mort[vector],1) 

efecto_prod_transformado <- efecto[xx] ##transformado a una funcion de valor
prod_xit <- (-xit * xit_log)

###############################################################################
#                                   K S I M                                   #
#                                                                             #
# Description: Coupled, dynamic, model that models forest population dynamics #
#              in terms of growth, regeneration and mortality. This is        #
#              considered an hybrid model because it can simulate processes   #
#              as well as states                                              #
#                                                                             #
###############################################################################

############## E C U A C I O N   D I F E R E N C I A L   K S I M ##############

#  dxit = xit * ln(xit) (suma_alfa + suma_beta) * dt

# DONDE:
# xi = cantidad de una variable
# xi0 = cantidad inicial
# xit = valor de la variable en el tiempo anterior t
# suma_alfa = suma de productos entre los valores xi en la matriz alfa de la variable xj
# suma_beta = suma de productos entre los valores xi en la matriz beta de la variable xj
# dt = intervalo de tiempo

########################## K S I M   F U N C T I O N ##########################

func_dxit <- function(xit, ksim_evaluar, alfaiter, betaiter, dksim_evaluar, dt) {
  
  xit_log <- log(xit)
  prod_xit <- (-xit * xit_log)
  suma_alfa <- crossprod(alfaiter, ksim_evaluar)
  suma_beta <- crossprod(betaiter, dksim_evaluar)
  suma_alfabeta <- (suma_alfa + suma_beta)
  dxit <- (prod_xit * suma_alfabeta) * dt
  
  return(dxit)
}

################### P R E P A R A C I O N   D E   D A T O S ###################

# Variables
variables <- c("Regeneracion natural",
               "Unidades de recursos",
               "Volumen de corta",
               "Conflictos con mercados",
               "Acuerdos de venta",
               "Ingresos",
               "Trabajo comunitario",
               "Ganancias",
               "Disminucion fuerza de trabajo",
               "Participacion",
               "Toma de decisiones",
               "Asistencia tecnica y uso de la informacion",
               "Diversificación del ingreso")

n.vars <- length(variables)

# Tiempo
time <- 0:100
n.time <- length(time)

# Valores iniciales
ini <- c(0.30, # Regeneracion natural
         0.50, # Unidades de recursos
         0.20, # Volumen de corta
         0.25, # Conflictos con mercados
         0.45, # Acuerdos de venta
         0.05, # Ingresos
         0.60, # Trabajo comunitario
         0.10, # Ganancias
         0.35, # Disminucion fuerza de trabajo
         0.75, # Participacion
         0.70, # Toma de decisiones
         0.80, # Asistencia tecnica y uso de la informacion
         0.05) # Diversificación del ingreso

# Intervalo de cambio en el tiempo
dt <- 0.05 

########## C A R G A R   M A R I C E S   D E   I N T E R A C C I O N ##########

#-------------------- A L F A   ( A C U M U L A T I V O ) --------------------#
input_alfa <- read.csv(paste0(dir,"/KSIM_R/_datos/input_alfa.csv"))
input_alfa <- as.matrix(input_alfa)

#-------------------------- B E T A   ( F L U J O ) --------------------------#
input_beta <- read.csv(paste0(dir,"/KSIM_R/_datos/input_beta.csv"))
input_beta <- as.matrix(input_beta)

################# C R E A R   M A T R I C E S   D E S D E   R #################

# #-------------------- A L F A   ( A C U M U L A T I V O ) --------------------#
# input_alfa <- matrix(data = 0, nrow = n.vars, ncol = n.vars)
#   colnames(input_alfa) <- variables
#   rownames(input_alfa) <- variables
# 
# input_alfa <- edit(input_alfa) # Editar manualmente
# 
# #-------------------------- B E T A   ( F L U J O ) --------------------------#
# input_beta <- matrix(data = 0, nrow = n.vars, ncol = n.vars)
#   colnames(input_beta) <- variables
#   rownames(input_beta) <- variables
# 
# input_alfa <- edit(input_alfa) # Editar manualmente

############# M A R T I C E S   D E   C O M P O R T A M I E N T O #############

#---------------------------------- K S I M ----------------------------------#
ksim_evaluar <- matrix(c(
  ini,
  rep(rep(0, n.vars), n.time-1)
  
), nrow = n.time, ncol = n.vars, byrow = TRUE)

  colnames(ksim_evaluar) <- c(variables)
  r.names <- time
  rownames(ksim_evaluar) <- c(r.names)

#---------------------- K S I M   D I F E R E N C I A L ----------------------#
dksim_evaluar <- matrix(c(
  rep(rep(0, n.vars), n.time)
  
), nrow = n.time, ncol = n.vars, byrow = TRUE)

  colnames(dksim_evaluar) <- c(variables)
  r.names <- time
  rownames(dksim_evaluar) <- c(r.names)

################## S I M U L A R   I N T E R A C C I O N E S ##################

for (i in 1:100){ ###i son filas
 
  j<-1

     for (j in 1:13){
      xit <- ksim_evaluar[i,j] # Estado de la variable a evaluar en iteracion i
  
      alfa <- as.vector(input_alfa[j,]) #Estos tienen que cambiar por fila
      beta <- as.vector(input_beta[j,]) #Estos tienen que cambiar por fila
      

       alfaiter<- alfa
       betaiter <- beta
 
 
       if (i==1) { # The condition must return TRUE or FALSE
        ksimini <-  as.vector(ksim_evaluar[i,])
        dksimini<- rep(0,13)
        dksim_evaluar[i+1,j] <- func_dxit(xit, ksimini, alfaiter, betaiter,dksimini, dt)
          
        # if (ksim_evaluar[i,j]+dksim_evaluar[i+1,j] >= 1){
        #     ksim_evaluar[i+1,j] <- 1
        #   } else if (ksim_evaluar[i,j]+dksim_evaluar[i+1,j] <= 0){
        #     ksim_evaluar[i+1,j] <- 0
        #   } else {
        #   ksim_evaluar[i+1,j]  <-  ksim_evaluar[i,j]+dksim_evaluar[i+1,j]
        #   }
        
        } else {
        dksim_evaluar[i+1,j] <- func_dxit(xit, ksim_evaluar[i,], alfaiter, betaiter, dksim_evaluar[i,], dt)
        ksim_evaluar[i+1,j]  <-  ksim_evaluar[i,j]+dksim_evaluar[i+1,j]
        }
     }
}

##################### P L O T   I N T E R A C C I O N E S #####################
  
plot(ksim_evaluar)


  
  
  
  
  
  
  
  
################## S I M U L A R   I N T E R A C C I O N E S ##################
############################### O R I G I N A L ###############################

for (i in 1:25){ ###i son filas
  
  j<-1
  
  for (j in 1:12){
    xit <- ksim_evaluar[i,j] # Estado de la variable a evaluar en iteracion i
    
    alfa <- as.vector(input_alfa[j,]) #Estos tienen que cambiar por fila
    beta <- as.vector(input_beta[j,]) #Estos tienen que cambiar pr fila
    
    
    alfaiter <- alfa
    betaiter <- beta
    
    
    if (i==1) { # The condition must return TRUE or FALSE
      ksimini <-  as.vector(ksim_evaluar[i,])
      dksimini<- rep(0,13)
      dksim_evaluar[i+1,j] <- func_dxit(xit, ksimini, alfaiter, betaiter,dksimini, dt)
      ksim_evaluar[i+1,j]  <-  ksim_evaluar[i,j]+dksim_evaluar[i+1,j]
    } else {
      dksim_evaluar[i+1,j] <- func_dxit(xit, ksim_evaluar[i,], alfaiter, betaiter, dksim_evaluar[i,], dt)
      ksim_evaluar[i+1,j]  <-  ksim_evaluar[i,j]+dksim_evaluar[i+1,j]
    }
  }
}














































library(ggplot2)
rodala?o1 <- exp$Rodal[1,]
rodala?o100 <- exp$Rodal[1000,]
clasdiam<- seq(1:15)
numeroarboles <-  data.frame(clasdiam, rodala?o1, rodala?o100)

plot(clasdiam,exp$AnnualGrowth[1,], col="red", xlab="Clase diametral", ylab= "Probabildiad de transici?n" )
lines(clasdiam,exp$AnnualGrowth[1,], col = "gray")

plot(clasdiam,exp$mortalidades[1,], col="red", xlab="Clase diametral", ylab= "Probabildiad de mortalidad" )
lines(clasdiam,exp$mortalidades[1,], col = "gray")


library(reshape2)
numeroarboles2<- melt(numeroarboles, id.vars=c("clasdiam"))
ggplot(data=numeroarboles2, aes(x=clasdiam , y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  
  scale_fill_brewer(palette="Paired")+
  theme_minimal()


hist(rodala?o1)










#bsenescencia_1_20 <- rep(0, 10)
# for (a in 1:length(probsenescencia_1_20)){
#   rodaltmp1<- rodal[1:10]
#   if ( rodaltmp1[a]==0){
#     probsenescencia_1_20[a]<- 0
#   } else{
#     muertes_senescencia1 <- rbinom(1,  rodaltmp1[a], sample(probabilidad_mortalidad_1_20,1))
#     probsenescencia_1_20[a]<- muertes_senescencia1
#   }
# }
# 
# 
# probsenescencia_20_40 <- rep(0, 9)
# rodaltmp2<- rodal[11:19]
# for (b in 1:length(probsenescencia_20_40)){
#   if (rodaltmp2[b]==0){
#     probsenescencia_20_40[b]<- 0
#   }else{
#     muertes_senescencia2 <- rbinom(1,rodaltmp2[b], sample(probabilidad_mortalidad_20_40,1))
#     probsenescencia_20_40[b]<- muertes_senescencia2
#   }
# }
# 
# if (rodal[20]==0){      
#   probsenescencia_40 <- 0
# }else{
#   probsenescencia_40 <- rbinom(1, rodal[20], sample(probabilidad_mortalidad_40,1))
# }
# 
# senescencia <- c(probsenescencia_1_20,probsenescencia_20_40,probsenescencia_40 )
