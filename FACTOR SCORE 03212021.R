#*** Loading Required/Useful Libraries ##
setwd("L:/Twin Analysis/10092020")
require(OpenMx)   #Loads OpenMx
require(psych)   #Loads Psych package

library(devtools)
library(umx)
library (OpenMx)
library (haven)
library (dplyr)
library(car)
library(psych) 
source("GenEpiHelperFunctions.R")
source("polychoricMeansMatrix3.R")

library(haven)
twins <- read_sav("twindata12032020.sav")
View(twins)
twins$RATE_V1C <- ordered(twins$RATE_V1C, levels = 0:1,labels = c("0", "1")) # conversion
twins$RATE_V2C <- ordered(twins$RATE_V2C, levels = 0:1,labels = c("0", "1")) # conversion
twins$RATE_V3C <- ordered(twins$RATE_V3C, levels = 0:1,labels = c("0", "1")) # conversion


twins$WORRY_V1C <- ordered(twins$WORRY_V1C, levels = 0:1,labels = c("0", "1")) # conversion
twins$WORRY_V2C <- ordered(twins$WORRY_V2C, levels = 0:1,labels = c("0", "1")) # conversion
twins$WORRY_V3C <- ordered(twins$WORRY_V3C, levels = 0:1,labels = c("0", "1")) # conversion


twins$CONCERN_V1C <- ordered(twins$CONCERN_V1C, levels = 0:1,labels = c("0", "1")) # conversion
twins$CONCERN_V2C <- ordered(twins$CONCERN_V2C, levels = 0:1,labels = c("0", "1")) # conversion
twins$CONCERN_V3C <- ordered(twins$CONCERN_V3C, levels = 0:1,labels = c("0", "1")) # conversion


selVars2		<-	c("vetsaid","RATE_V1C", "WORRY_V1C", "CONCERN_V1C",
               "RATE_V2C", "WORRY_V2C", "CONCERN_V2C","RATE_V3C", "WORRY_V3C", "CONCERN_V3C")

testdata <- as.data.frame(twins[selVars2])
testdata2 <- na.omit(testdata)
require(OpenMx)

manif<-c("RATE_V1C","WORRY_V1C", "CONCERN_V1C", "RATE_V2C","WORRY_V2C", "CONCERN_V2C", "RATE_V3C","WORRY_V3C", "CONCERN_V3C")
manifestsv1 <- c("RATE_V1C","WORRY_V1C", "CONCERN_V1C")
manifestsv2 <- c("RATE_V2C","WORRY_V2C", "CONCERN_V2C")
manifestsv3 <- c("RATE_V3C","WORRY_V3C", "CONCERN_V3C")

latents <- c("v1","v2","v3")

umx::umx_set_optimizer(opt="SLSQP")
mxOption(NULL,"mvnRelEps",0.0055)
mxOption(NULL, 'Number of Threads', parallel::detectCores())


factorModel <- mxModel("OneFactor", 
                       type="RAM",
                       manifestVars=c("RATE_V1C","WORRY_V1C", "CONCERN_V1C", "RATE_V2C","WORRY_V2C", "CONCERN_V2C", "RATE_V3C","WORRY_V3C", "CONCERN_V3C"), latentVars=c("v1","v2","v3"),
                       # residual variances
                       
                       mxPath(from=c("RATE_V1C","WORRY_V1C", "CONCERN_V1C", "RATE_V2C","WORRY_V2C", "CONCERN_V2C", "RATE_V3C","WORRY_V3C", "CONCERN_V3C"), arrows=2, # 2-headed arrow for residual variance
                              
                              free=TRUE, # All variances to be estimated
                              
                              values=1 # Start value is 1
                              
                       ),
                       # factor variance
                       mxPath(from=c("v1","v2","v3"), arrows=2, free=FALSE, values = 1),
                       # factor covariance
                       mxPath(from="v1", to="v2",  arrows=1, free=TRUE, values=0.5),
                       mxPath(from="v2", to="v3",  arrows=1, free=TRUE, values=0.5),
                       #mxPath(from="v1", to="v3",  arrows=2, values=0.5),
                       # Expected means 
                       mxPath(from = 'one', to = manif),
                       # factors manifest vars
                       mxPath(from="v1", to=c("RATE_V1C","CONCERN_V1C", "WORRY_V1C"), arrows=1, free=TRUE, values=.6),
                       mxPath(from="v2", to=c("RATE_V2C","CONCERN_V2C", "WORRY_V2C"), arrows=1, free=TRUE, values=.6),
                       mxPath(from="v3", to=c("RATE_V3C","CONCERN_V3C", "WORRY_V3C"), arrows=1, free=TRUE, values=.6),
                       # error terms
                       mxPath(from=manif, arrows=2, values=1, free=TRUE),
                       mxThreshold(vars=manif, free=T, nThresh=1, values = .70),
                       mxData(observed=testdata2[manif], type="raw"))
model2 <- mxOption(factorModel, 'mvnRelEps', mxOption(factorModel, 'mvnRelEps')/5)
summary(factorRun <- mxTryHardOrdinal(model2))
options(max.print = 500)
uls 				<- mxAutoStart(factorModel)
summary( uls_fit	<- mxTryHardOrdinal(uls))    
factorscore3<-mxFactorScores(factorModel, type=c('WeightedML'),
                             minManifests=as.integer(3))


factorscore_v3 <- as.data.frame(factorscore3)

total <- cbind(testdata2,factorscore_v3)

write_xlsx(total,"factorscore03212021.xlsx")
