#######################################################
#########  Juan Sebastian Nocua         ###############
#########  js.nocua@uniandes.edu.co     ###############
######### Last modified: 03/06/2021     ###############
#######################################################

##############################################################
######################### Functions ##########################
##############################################################

library(EnvStats)

heapsim <- function(heap, timestep, season, lengthseason, epocaInicial = FALSE) {
  
  ##################################### Input HeapSim #####################################
  #--settings of the model
  # non-decomposable fraction of the manure on the heap
  inertDM <- 0.3      # (-)
  # inertC  <- 0.3
  # inertN  <- 0.3
  # inertP  <- 0.3
  # inertK  <- 0.3
  # relative decomposition rate of the manure on the heap
  rdrDM <- 0.15        # month-1
  # rdrC  <- 0.3
  # rdrN  <- 0.3
  # rdrP  <- 0.3
  # rdrK  <- 0.3
  # non-decomposable fraction of the organic residues
  inertDM.OR <- 0.3   # (-)
  inertC.OR  <- 0.3   # (-)
  inertN.OR  <- 0.3   # (-)
  inertP.OR  <- 0.3   # (-)
  inertK.OR  <- 0.3   # (-)
  # C-content of the animal manure is estimated from data (see dataset_maize.xls)
  manureCContent <- 0.27    # kg C kg DM-1
  
  porcentajeCacaoOrganic <- 0.4
  
  # Literature information
  
  if (epocaInicial == TRUE){
    totalFaecalDM <-  1200*(1-runif(1,0.7,0.8))  # kg DM
    totalFaecalC  <-  totalFaecalDM*18.83*0.0291    #totalFaecalDM * manureCContent  # kg C
    totalFaecalN  <- totalFaecalDM*0.0291    # kg N
    totalFaecalP  <- totalFaecalDM*0.0308*0.44     # kg P
    totalFaecalK  <-  totalFaecalDM*0.0281*0.83    # kg K
    
    inputOrganicDM <- 1200*(1-porcentajeCacaoOrganic)*(1-0.1506)
    inputCocoaDM <- 1200*porcentajeCacaoOrganic*(1-0.1004)
    
    inputTotalDM <- totalFaecalDM + inputOrganicDM + inputCocoaDM                           # kg DM
    inputTotalC  <- totalFaecalC + inputOrganicDM*0.0199*1.38 +inputCocoaDM*(0.484)                               # kg C
    inputTotalN  <- totalFaecalN + inputOrganicDM*0.0199 +inputCocoaDM*(0.0136)  # agregar materia orgánica                 # kg N
    inputTotalP  <- totalFaecalP + inputOrganicDM*0.0002 +inputCocoaDM*(0.00017)                  # kg P
    inputTotalK  <- totalFaecalK + inputOrganicDM*0.02768 +inputCocoaDM*(0.0133)
    
    #NEW
    
    #percentages
    oxygenPercentageFaecal <- 0.0308*(1-0.44)
    hydrogenPercentagefaecal <- oxygenPercentageFaecal*0.052/0.688
    
    oxygenPercentageCocoa <- rtri(n = 1,min = 0.56,mode = 0.688,max =0.737)
    hydrogenPercentageCocoa <- rtri(n = 1,min = 0.051,mode = 0.052,max =0.060)
    
    oxygenPercentageOrganic <- 0.4735
    hydrogenPercentageOrganic <- 0.0544
    
    #quantites
    totalFaecalO <- totalFaecalDM*oxygenPercentageFaecal
    inputTotalO <- totalFaecalO+inputOrganicDM*oxygenPercentageOrganic+inputCocoaDM*oxygenPercentageCocoa
    
    totalFaecalH <- totalFaecalDM*hydrogenPercentagefaecal
    inputTotalH <- totalFaecalO+inputOrganicDM*hydrogenPercentageOrganic+inputCocoaDM*hydrogenPercentageCocoa
    
    #Heap
    heap$amountDM[timestep] <- inputTotalDM  # kg DM
    heap$amountC[timestep] <- inputTotalC # kg C
    heap$amountN[timestep] <- inputTotalN  # kg N
    heap$amountP[timestep] <- inputTotalP  # kg O
    heap$amountK[timestep] <- inputTotalK  # kg K
    heap$amountCocoa[timestep] <- inputCocoaDM  # kg Cocoa
    heap$amountO[timestep] <- inputTotalO  # kg O
    heap$amountH[timestep] <- inputTotalH  # kg H
    
    heap$temperature[timestep] <- 75  # °C
    heap$reactionRate[timestep] <- NA 
    heap$mass[timestep] <- 2400  # kg
    
    heap$amountOrganic[timestep] <- inputOrganicDM  # kg
    
    return(heap)
  }
  else{
    inputCocoaDM <- 0
    totalFaecalDM <-  (1-0.75)  # kg DM
    totalFaecalC  <- (1-0.75)*18.83*0.0291    #totalFaecalDM * manureCContent # kg C
    totalFaecalN  <- (1-0.75)*0.0291    # kg N
    totalFaecalP  <- (1-0.75)*0.0308*0.44 # kg P
    totalFaecalK  <-  (1-0.75)*0.0281*0.83
    
    totalFaecalO <- 0.1
    totalFaecalH <- 0.01
    
    inputTotalDM <- totalFaecalDM  # kg DM
    inputTotalC  <- totalFaecalC   # kg C
    inputTotalN  <- totalFaecalN   # kg N
    inputTotalP  <- totalFaecalP   # kg P
    inputTotalK  <- totalFaecalK
    inputTotalO <- totalFaecalO
    inputTotalH <- totalFaecalH
  }
  
  ################################################# HeapSim ##########################################################################
  # Currently the decomposition of DM is calculated directly, the decomposition rates of the other nutrients is based on the DM 
  # decomposition and the ratio of the other nutrients to the DM content of the manure.
  
  heapDecompositionDM <- rdrDM * (1 - inertDM) * heap$amountDM[timestep]        # kg DM month-1
  heapDecompositionC  <- heapDecompositionDM / (totalFaecalDM / totalFaecalC)   # kg C month-1
  heapDecompositionN  <- heapDecompositionDM / (totalFaecalDM / totalFaecalN)   # kg N month-1
  heapDecompositionP  <- heapDecompositionDM / (totalFaecalDM / totalFaecalP)   # kg P month-1
  heapDecompositionK  <- heapDecompositionDM / (totalFaecalDM / totalFaecalK)   # kg K month-1
  heapDecompositionO  <- heapDecompositionDM / (totalFaecalDM / totalFaecalO)   # kg O month-1
  heapDecompositionH  <- heapDecompositionDM / (totalFaecalDM / totalFaecalH)   # kg H month-1
  heapDecompositionOrganic  <- heapDecompositionDM*(1-porcentajeCacaoOrganic)*(1-0.1506)  # kg Organic month-1
  
  heapNetRateOfChangeDM <- -heapDecompositionDM + inputTotalDM
  heapNetRateOfChangeC  <- -heapDecompositionC + inputTotalC
  heapNetRateOfChangeN  <- -heapDecompositionN + inputTotalN
  heapNetRateOfChangeP  <- -heapDecompositionP + inputTotalP
  heapNetRateOfChangeK  <- -heapDecompositionK + inputTotalK
  heapNetRateOfChangeCocoa <- -heapDecompositionDM*porcentajeCacaoOrganic + inputCocoaDM
  heapNetRateOfChangeO  <- -heapDecompositionO + inputTotalO
  heapNetRateOfChangeH  <- -heapDecompositionH + inputTotalH
  heapNetRateOfChangeOrganic  <- -heapDecompositionOrganic
  
  if (is.nan(heapNetRateOfChangeC)) {heapNetRateOfChangeC <- 0} 
  if (is.nan(heapNetRateOfChangeN)) {heapNetRateOfChangeN <- 0} 
  if (is.nan(heapNetRateOfChangeP)) {heapNetRateOfChangeP <- 0} 
  if (is.nan(heapNetRateOfChangeK)) {heapNetRateOfChangeK <- 0}
  if (is.nan(heapNetRateOfChangeCocoa)) {heapNetRateOfChangeCocoa <- 0}
  if (is.nan(heapNetRateOfChangeO)) {heapNetRateOfChangeO <- 0}
  if (is.nan(heapNetRateOfChangeH)) {heapNetRateOfChangeH <- 0}
  
  heap$amountDM[timestep + 1] <- max(heap$amountDM[timestep] + heapNetRateOfChangeDM,0)  # kg DM
  heap$amountC[timestep + 1] <- max(heap$amountC[timestep] + heapNetRateOfChangeC,0)     # kg C
  heap$amountN[timestep + 1] <- max(heap$amountN[timestep] + heapNetRateOfChangeN,0)     # kg N
  heap$amountP[timestep + 1] <- max(heap$amountP[timestep] + heapNetRateOfChangeP,0)     # kg P
  heap$amountK[timestep + 1] <- max(heap$amountK[timestep] + heapNetRateOfChangeK,0)     # kg K
  heap$amountCocoa[timestep + 1] <- max(heap$amountCocoa[timestep] + heapNetRateOfChangeCocoa,0)  # kg DM
  heap$amountO[timestep + 1] <- max(heap$amountO[timestep] + heapNetRateOfChangeO,0)     # kg O
  heap$amountH[timestep + 1] <- max(heap$amountH[timestep] + heapNetRateOfChangeH,0)     # kg H
  heap$amountOrganic[timestep + 1] <- max(heap$amountOrganic[timestep] + heapNetRateOfChangeH,0)     # kg Organic
  
  ################ PETRIC METHODOLOGY ###############
  
  #Yield
  R <- 100*(2.66*heap$amountC[timestep]+7.94*heap$amountH[timestep]-heap$amountO[timestep])/398.9
  heap$reactionRate[timestep] <- R
  
  Q <- (127*R+400)*4.186
  
  #specific heats
  faecalAshPercentage <- 0.2755
  organicAshPercentage <- 0.0961
  cocoaAshPercentage <- rtri(n =1, min = 0.0105, mode = 0.041,max = 0.143)
  
  faecalHumidity <- runif(1,0.7,0.8)
  organicHumidity <- 0.1506
  cocoaHumidity <- 0.1004
  
  faecalCp <- 1.48-0.64*1200*faecalAshPercentage+4.18*1200*faecalHumidity
  organicCp <- 1.48-0.64*1200*(1-porcentajeCacaoOrganic)*organicAshPercentage+4.18*1200*(1-porcentajeCacaoOrganic)*organicHumidity
  cocoaCp <- 1.48-0.64*1200*porcentajeCacaoOrganic*cocoaAshPercentage+4.18*1200*porcentajeCacaoOrganic*cocoaHumidity
  
  deltaT <- Q/(faecalCp*heap$amountDM[timestep]+organicCp*heap$amountOrganic[timestep]+cocoaCp*heap$amountCocoa[timestep])
  
  #Temperature update
  
  heap$temperature[timestep+1] <- heap$temperature[timestep]-deltaT
  
  #Mass update
  
  Kt <- heap$temperature[timestep]*(80-heap$temperature[timestep])/1600
  deltaMass <- -Kt*heap$mass[timestep]
  heap$mass[timestep+1] <- heap$mass[timestep]+deltaMass
  
  return(heap)
}



library(plotly)

elementGraph <- function(element, heapPrueba){
  
  m = list(size = 10,
           color = 'rgba(255, 182, 193, .9)',
           line = list(color = 'rgba(152, 0, 0, .8)',
                       width = 2))
  
  if(element == "DM"){
    fig <- plot_ly(data = as.data.frame(heapPrueba), x =~amountCocoa[2:7], y = ~amountDM[2:7],
                   marker = m)
    tit <- "Dry matter Quantity vs Cocoa Quantity"
    yax <- "Dry matter Quantity (kg)"
  }
  if(element == "C"){
    fig <- plot_ly(data = as.data.frame(heapPrueba), x =~amountCocoa[2:7], y = ~amountC[2:7],
                   marker = m)
    tit <- "Carbon Quantity vs Cocoa Quantity"
    yax <- "Carbon Quantity (kg)"
  }
  if(element == "N"){
    fig <- plot_ly(data = as.data.frame(heapPrueba), x =~amountCocoa[2:7], y = ~amountN[2:7],
                   marker = m)
    tit <- "Nitrogen Quantity vs Cocoa Quantity"
    yax <- "Nitrogen Quantity (kg)"
  }
  if(element == "P"){
    fig <- plot_ly(data = as.data.frame(heapPrueba), x =~amountCocoa[2:7], y = ~amountP[2:7],
                   marker = m)
    tit <- "Phosphorus Quantity vs Cocoa Quantity"
    yax <- "Phosphorus Quantity (kg)"
  }
  if(element == "K"){
    fig <- plot_ly(data = as.data.frame(heapPrueba), x =~amountCocoa[2:7], y = ~amountK[2:7],
                   marker = m)
    tit <- "Potassium Quantity vs Cocoa Quantity"
    yax <- "Potassium Quantity (kg)"
  }
  
  
  fig <- fig %>% layout(title = 'Styled Scatter',
                        yaxis = list(zeroline = FALSE),
                        xaxis = list(zeroline = FALSE))
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  x <- list(
    title = "Cocoa Quantity (kg)",
    titlefont = f)
  
  y <- list(
    title = yax,
    titlefont = f)
  
  fig <- fig %>% layout(title = tit, xaxis = x, yaxis = y)
  fig
}

##############################################################
########################## Tests ###########################
##############################################################

#Heap with parameters is initialized based on number of months of analysis
numberOfMonths = 6
Heap <- list(amountDM = vector('numeric', length = numberOfMonths),
             amountC  = vector('numeric', length = numberOfMonths),
             amountN  = vector('numeric', length = numberOfMonths),
             amountP  = vector('numeric', length = numberOfMonths),
             amountK  = vector('numeric', length = numberOfMonths),
             amountCocoa  = vector('numeric', length = numberOfMonths),
             amountO  = vector('numeric', length = numberOfMonths),
             amountH  = vector('numeric', length = numberOfMonths),
             temperature  = vector('numeric', length = numberOfMonths),
             reactionRate  = vector('numeric', length = numberOfMonths),
             mass  = vector('numeric', length = numberOfMonths),
             amountOrganic  = vector('numeric', length = numberOfMonths),
             manureApplied = 0)

#Initilization of chemical components 
heapPrueba = heapsim(Heap,1,1,numberOfMonths, TRUE)

#Iterations
for(i in 1:(numberOfMonths-1)){
  heapPrueba = heapsim(heapPrueba,i,1,12, FALSE)
}

########################################################
#################### Graphs ############################
########################################################

#Chemical element graph
elementGraph(element = "K", heapPrueba)

#Mass as function of time (Modified Heapsim)
finalMass <- c()
for(i in 1:numberOfMonths){
  finalMass <- c(finalMass,heapPrueba$amountDM[i]+heapPrueba$amountC[i]+heapPrueba$amountK[i]+heapPrueba$amountN[i]+heapPrueba$amountP[i]+heapPrueba$amountCocoa[i])
}

m = list(size = 10,
         color = 'rgba(255, 182, 193, .9)',
         line = list(color = 'rgba(152, 0, 0, .8)',
                     width = 2))

fig <- plot_ly(x =c(1:numberOfMonths), y = finalMass,
               marker = m)
tit <- "Mass of pile versus Time"
yax <- "Mass of compost (kg)"

fig <- fig %>% layout(title = 'Styled Scatter',
                      yaxis = list(zeroline = FALSE),
                      xaxis = list(zeroline = FALSE))
f <- list(
  family = "Courier New, monospace",
  size = 18,
  color = "#7f7f7f"
)
x <- list(
  title = "Time (months)",
  titlefont = f)

y <- list(
  title = yax,
  titlefont = f)

fig <- fig %>% layout(title = tit, xaxis = x, yaxis = y)
fig


#Mass as function of time (Enhanced with Petric Heapsim)
m = list(size = 10,
         color = 'rgba(255, 182, 193, .9)',
         line = list(color = 'rgba(152, 0, 0, .8)',
                     width = 2))

fig <- plot_ly(data = as.data.frame(heapPrueba), x =c(1:numberOfMonths), y = ~mass,
               marker = m)
tit <- "Mass of pile versus Time"
yax <- "Mass of compost (kg)"

fig <- fig %>% layout(title = 'Styled Scatter',
                      yaxis = list(zeroline = FALSE),
                      xaxis = list(zeroline = FALSE))
f <- list(
  family = "Courier New, monospace",
  size = 18,
  color = "#7f7f7f"
)
x <- list(
  title = "Time (months)",
  titlefont = f)

y <- list(
  title = yax,
  titlefont = f)

fig <- fig %>% layout(title = tit, xaxis = x, yaxis = y)
fig



