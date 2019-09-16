######################################################
### Code to Analyze SW-CRTs Using Methods from     ###
### in Kennedy-Shaffer et al. 2019                 ###
### Update: 9/15/2019                              ###
### Contact: lee_kennedyshaffer@g.harvard.edu      ###
### See Latest Update at:                          ###
### https://github.com/leekshaffer/SW-CRT-analysis ###
######################################################


######################################################
###  check/install/load needed packages            ###
######################################################

if (!require(MASS)) {
  install.packages("MASS")
  library(MASS)
}
if (!require(Synth)) {
  install.packages("Synth")
  library(Synth)
}
if (!require(lme4)) {
  install.packages("lme4")
  library(lme4)
}


######################################################
###  SC-SWT Functions                              ###
######################################################

### Helper Functions ###

expit <- function(x) exp(x)/(1+exp(x))

## Functions Defining Contrasts of Interests and Inverting those Contrasts ##
RDFunc <- function(a,b) a-b
RDInvApply <- function(out,effect) out-effect
logRRFunc <- function(a,b) log(a/b)
logRRInvApply <- function(out,effect) out*exp(-1*effect)
logORFunc <- function(a,b) log(a*(1-b)/((1-a)*b))
logORInvApply <- function(out,effect) expit(-1*effect+log(out/(1-out)))

## Create a vector with each cluster's first intervention period, by observation ##
StartTimeVec <- function(Periods, Clusters, Trts) {
  StartTimes <- tapply(Periods[Trts==1], Clusters[Trts==1], FUN=min)
  NumPds <- tapply(Periods, Clusters, FUN=length)
  return(rep(as.vector(StartTimes), times=NumPds))
}


### Methods for Estimation from SW-CRT ###

## NPWP Method from Thompson et al. 2018##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with one element:
### Est.NPWP: the estimated treatment effect on the specified contrast scale

NPWP.Effect.Est <- function(Periods, Outcomes, StartTimes, ContrastFunc) {
  PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0 & PercInt < 1])[[1]])
  ContrastPds <- ContrastPds[order(ContrastPds)]
  Contrasts <- matrix(rep(NA,2*length(ContrastPds)),nrow=2, ncol=length(ContrastPds))
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    Ints <- Outcomes[Periods == Pd & StartTimes <= Pd]
    Conts <- Outcomes[Periods == Pd & StartTimes > Pd]
    NumInt <- length(Ints)
    NumCont <- length(Conts)
    IntMean <- mean(Ints)
    ContMean <- mean(Conts)
    IntVar <- ifelse(NumInt==1,0,var(Ints))
    ContVar <- ifelse(NumCont==1,0,var(Conts))
    Contrast <- ContrastFunc(IntMean, ContMean)
    Weight <- ifelse(IntVar==0 & ContVar==0, (10^(-5)*(1/NumCont + 1/NumInt))^(-1), 
                     ((((NumCont-1)*ContVar + (NumInt-1)*IntVar)/(NumCont+NumInt-2))*(1/NumCont + 1/NumInt))^(-1))
    Contrasts[,i] <- c(Contrast,Weight)
  }
  TxEst <- sum(Contrasts[1,]*Contrasts[2,])/sum(Contrasts[2,])
  return(list(Est.NPWP=TxEst))
}


## Crossover Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
## CtrlType: one  of "Ctrl","Both","CtWt","BoWt", indicating whether only clusters
###          on control in both periods should be used as comparisons ("Ctrl"/"CtWt")
###          or cluster on intervention in both periods should be used as well ("Both"/"BoWt")
###          and whether the weights should be equal across periods ("Ctrl"/"Both")
###          or proportional to the harmonic mean of the number of intervention and comparison clusters  ("CtWt"/"BoWt").
###          Multiple values may be specified as a vector to get multiple effect estimates.
###          In the article, "Ctrl" is CO-1, "CtWt" is CO-2, and "Both" is CO-3.
# Outputs:
## A list with one element for each value in CtrlType:
### Est.CO.[CtrlType]: the estimated treatment effect for the specified CtrlType

CO.Effect.Est <- function(Periods, Outcomes, StartTimes, ContrastFunc, 
                          CtrlType=c("Ctrl","Both","CtWt","BoWt")) {
  Pds <- sort(unique(Periods))
  PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0])[[1]])
  ContrastPds <- sort(ContrastPds)
  if (ContrastPds[1] == Pds[1]) {ContrastPds <- ContrastPds[-1]}
  PdEsts <- matrix(rep(NA,length(CtrlType)*length(ContrastPds)),
                   nrow=length(ContrastPds),ncol=length(CtrlType))
  Weights <- matrix(rep(0,length(CtrlType)*length(ContrastPds)),
                    nrow=length(ContrastPds),ncol=length(CtrlType))
  for (i in 1:length(ContrastPds)) {
    Pdi <- ContrastPds[i]
    Type <- ifelse(StartTimes[Periods == Pdi] <= Pdi,1,0) + ifelse(StartTimes[Periods == Pdi] <= Pdi-1,1,0)
    n1 <- sum(ifelse(Type==1,1,0), na.rm=TRUE)
    n0 <- sum(ifelse(Type==0,1,0), na.rm=TRUE)
    n2 <- sum(ifelse(Type==2,1,0), na.rm=TRUE)
    Contrasts <- ContrastFunc(Outcomes[Periods == Pdi], Outcomes[Periods == Pdi-1])
    for (j in 1:length(CtrlType)) {
      if (CtrlType[j]=="Ctrl") {
        CtrlCont <- mean(Contrasts[Type==0], na.rm=TRUE)
        Weights[i,j] <- 1
      } else if (CtrlType[j]=="CtWt") {
        CtrlCont <- mean(Contrasts[Type==0], na.rm=TRUE)
        Weights[i,j] <- (1/n1 + 1/n0)^(-1)
      } else if (CtrlType[j]=="Both") {
        CtrlCont <- mean(Contrasts[Type %in% c(0,2)])
        Weights[i,j] <- 1
      } else if (CtrlType[j]=="BoWt") {
        CtrlCont <- mean(Contrasts[Type %in% c(0,2)])
        Weights[i,j] <- (1/n1 + 1/(n0+n2))^(-1)
      } else {
        stop("Please enter a vector with possible elements 'Both', 'Ctrl', 'BoWt', and 'CtWt' for variable CtrlType")
      }
        IntCont <- mean(Contrasts[Type==1], na.rm=TRUE)
        if (!is.nan(CtrlCont)) {PdEsts[i,j] <- IntCont - CtrlCont}
    }
  }
  outlist <- list(ContrastMat=Contrasts,PdEstsMat=PdEsts)
  for (j in 1:length(CtrlType)) {
    nam <- paste0("Est.CO.",CtrlType[j])
    TxEff <- sum(PdEsts[,j]*Weights[,j], na.rm=TRUE)/sum(ifelse(is.na(PdEsts[,j]),0,Weights[,j]), na.rm=TRUE)
    outlistint <- list(TxEff)
    names(outlistint) <- nam
    outlist <- append(outlist, outlistint)
  }
  return(outlist)
}


## Mixed Effects Model Method from Hussey & Hughes 2007 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: For binary outcomes, a vector of the number of individuals measured in each
###        cluster-period, in the same order as Periods; or, a single number if the number
###        of individuals measured is the same in each cluster-period.
###        For non-binary outcomes, leave as NULL and specify each individual's
###        outcome, cluster, period, and start time in the respective input vector.
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity"
# Outputs:
## A list with two elements:
### Est.MEM: the estimated treatment effect
### PVal.MEM: the asymptotic p-value for the null hypothesis of no treatment effect

MEM.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                           family, link) {
  Trts <- ifelse(Periods >= StartTimes, 1, 0)
  
  if (!is.null(Indivs)) {
    if (length(Indivs)==1) {
      Indivs = rep(Indivs, length(Periods))
    }
    NewTrts <- rep(Trts, times=Indivs)
    NewClusts <- rep(Clusters, times=Indivs)
    NewPds <- rep(Periods, times=Indivs)
    NewOuts <- NULL
    for (i in 1:length(Outcomes)) {
      NewOuts <- append(NewOuts, rep(1,round(Outcomes[i]*Indivs[i], digits=0)))
      NewOuts <- append(NewOuts, rep(0,Indivs[i]-round(Outcomes[i]*Indivs[i], digits=0)))
    }
  } else {
    NewTrts <- Trts
    NewClusts <- Clusters
    NewPds <- Periods
    NewOuts <- Outcomes
  }
  PdFactor <- relevel(factor(NewPds), ref=1)
  Indexes <- NewClusts*max(NewPds)+NewPds
  fam <- family(link=link)
  capture.output(res <- tryCatch({glmer(NewOuts~NewTrts+PdFactor+(1|NewClusts), 
                                        family=fam)},
                                      error=function(err) {1}), file="/dev/null")
  if (is.numeric(res)) {
    Est <- NA
    PVal <- NA
  } else {
    Est <- as.numeric(summary(res)$coefficients['NewTrts','Estimate'])
    SE <- as.numeric(summary(res)$coefficients['NewTrts','Std. Error'])
    PVal <- as.numeric(summary(res)$coefficients['NewTrts','Pr(>|z|)'])
  }
  return(list(Est.MEM=Est, PVal.MEM=PVal))
}

## Mixed Effects Model with Cluster-Period Interaction Method from Hooper et al. 2016 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: For binary outcomes, a vector of the number of individuals measured in each
###        cluster-period, in the same order as Periods; or, a single number if the number
###        of individuals measured is the same in each cluster-period
###        For non-binary outcomes, leave as NULL and specify each individual's
###        outcome, cluster, period, and start time in the respective input vector.
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity"
# Outputs:
## A list with two elements:
### Est.CPI: the estimated treatment effect
### PVal.CPI: the asymptotic p-value for the null hypothesis of no treatment effect

CPI.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                           family, link) {
  Trts <- ifelse(Periods >= StartTimes, 1, 0)
  
  if (!is.null(Indivs)) {
    if (length(Indivs)==1) {
      Indivs = rep(Indivs, length(Periods))
      }
    NewTrts <- rep(Trts, times=Indivs)
    NewClusts <- rep(Clusters, times=Indivs)
    NewPds <- rep(Periods, times=Indivs)
    NewOuts <- NULL
    for (i in 1:length(Outcomes)) {
      NewOuts <- append(NewOuts, rep(1,round(Outcomes[i]*Indivs[i], digits=0)))
      NewOuts <- append(NewOuts, rep(0,Indivs[i]-round(Outcomes[i]*Indivs[i], digits=0)))
    }
  } else {
    NewTrts <- Trts
    NewClusts <- Clusters
    NewPds <- Periods
    NewOuts <- Outcomes
  }
  PdFactor <- relevel(factor(NewPds), ref=1)
  Indexes <- NewClusts*max(NewPds)+NewPds
  fam <- family(link=link)
  capture.output(res <- tryCatch({glmer(NewOuts~NewTrts+PdFactor+(1|NewClusts)+(1|Indexes), 
                                        family=fam)},
                                 error=function(err) {1}), file="/dev/null")
  if (is.numeric(res)) {
    Est <- NA
    PVal <- NA
  } else {
    Est <- as.numeric(summary(res)$coefficients['NewTrts','Estimate'])
    SE <- as.numeric(summary(res)$coefficients['NewTrts','Std. Error'])
    PVal <- as.numeric(summary(res)$coefficients['NewTrts','Pr(>|z|)'])
  }
  return(list(Est.CPI=Est, PVal.CPI=PVal))
}

## Synthetic Control Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with three elements:
### Est.SCSWT1: the estimated treatment effect for equal weighting (SC-1 in the article)
### Est.SCSWT2: the estimated treatment effect for weighting by inverse-MSPE within equivalence groups (SC-2 in the article)
### ContrastMat: the matrix of cluster-period-specific treatment effect estimates, as well as the weights used in SC-1 and SC-2.
####             Can be used for alternate weighting schemes

SCSWT.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc) {
  
  Index <- Clusters*max(Periods)+Periods
  DataMat <- matrix(data=Outcomes[order(Index)], nrow=length(unique(Periods)), 
                    ncol=length(unique(Clusters)), byrow=FALSE,
                    dimnames=list(sort(unique(Periods)),sort(unique(Clusters))))
  ## Note: DataMat here only works if all clusters have same periods ##
  
  PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0 & PercInt < 1])[[1]])
  ContrastPds <- ContrastPds[order(ContrastPds)]
  Contrasts <- data.frame(Period = NULL, Cluster = NULL, StartTime = NULL, Contrast = NULL, MSPE = NULL,
                          ClustMin = NULL)
  
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    IntClusts <- Clusters[Periods == Pd & StartTimes <= Pd]
    ContClusts <- Clusters[Periods == Pd & StartTimes > Pd]
    NumInts <- length(IntClusts)
    NumConts <- length(ContClusts)
    ClustMinVec <- rep(min(NumInts,NumConts),NumInts)
    
    if (Pd == min(Periods)) {
      SCOutcomes <- rep(mean(Outcomes[Periods == Pd & Clusters %in% ContClusts]),NumInts)
      MSPEs <- rep(10^(-8), NumInts)
      StartTimeForDF <- rep(1, NumInts)
    } else {
      PrevPds <- sort(unique(Periods[Periods < Pd]))
      DataMati <- DataMat[as.character(PrevPds),,drop=FALSE]
    
      MSPEs <- rep(10^(-8), NumInts) ## Truncates MSPE at minimum of 10^(-8) to prevent inverse of zero issue
      StartTimeForDF <- rep(NA, NumInts)
      if (length(ContClusts) == 1) {
        SCOutcomes <- rep(Outcomes[Periods == Pd & Clusters == ContClusts], NumInts)
        for (j in 1:NumInts) {
          StartTimej <- min(StartTimes[Clusters==IntClusts[j]])
          StartTimeForDF[j] <- StartTimej
          PrevPdsj <- PrevPds[PrevPds < StartTimej]
          MSPEs[j] <- max(MSPEs[j], mean((DataMati[as.character(PrevPdsj),as.character(ContClusts)]-DataMati[as.character(PrevPdsj),as.character(IntClusts[j])])^2))
        }
      } else {
        SCOutcomes <- rep(0, NumInts)
        for (j in 1:NumInts) {
          StartTimej <- min(StartTimes[Clusters==IntClusts[j]])
          StartTimeForDF[j] <- StartTimej
          PrevPdsj <- PrevPds[PrevPds < StartTimej]
          X0ij <- DataMati[as.character(PrevPdsj),as.character(ContClusts),drop=FALSE]
          X1ij <- DataMati[as.character(PrevPdsj),as.character(IntClusts[j]),drop=FALSE]
          capture.output(SynthOut <- tryCatch({synth(X1=X1ij, X0=X0ij, Z0=X0ij, Z1=X1ij)},
                                              error=function(err) {1}), file="/dev/null")
          if (is.numeric(SynthOut)) {
            SCOutcomes[j] <- mean(Outcomes[Periods == Pd & Clusters %in% ContClusts])
            MSPEs[j] <- max(MSPEs[j],mean((X1ij - 
                                apply(X=X0ij,MARGIN=1,FUN=mean))^2))
          } else {
            MSPEs[j] <- max(MSPEs[j],as.numeric(SynthOut$loss.v))
            SCOutcomes[j] <- as.numeric(SynthOut$solution.w) %*% Outcomes[Periods == Pd & Clusters %in% ContClusts]
          }
        }
      }
    }
    ContrastVec <- ContrastFunc(Outcomes[Periods == Pd & StartTimes <= Pd],SCOutcomes)
    Contrasts <- data.frame(Period = append(Contrasts$Period, rep(Pd, NumInts)),
                            Cluster = append(Contrasts$Cluster, IntClusts),
                            StartTime = append(Contrasts$StartTime, StartTimeForDF),
                            Contrast = append(Contrasts$Contrast,ContrastVec),
                            MSPE = append(Contrasts$MSPE, MSPEs),
                            ClustMin = append(Contrasts$ClustMin, ClustMinVec))
  }
  Contrasts$PdMin <- rep(0,length(Contrasts$Contrast))
  Contrasts$SumMSPE <- rep(0,length(Contrasts$Contrast))
  for (i in 1:length(Contrasts$Contrast)) {
    Clusteri <- Contrasts$Cluster[i]
    ClustStart <- Contrasts$StartTime[i]
    Contrasts$SumMSPE[i] <- sum(1/Contrasts$MSPE[Contrasts$StartTime==ClustStart])
    Contrasts$PdMin[i] <- min(sum(ContrastPds >= ClustStart),
                                   length(unique(Periods[Periods < ClustStart])))
  }
  Contrasts$Wt2 <- (Contrasts$MSPE*Contrasts$SumMSPE)^(-1)
  
  TxEst1 <- mean(Contrasts$Contrast)
  TxEst2 <- sum(Contrasts$Contrast*Contrasts$Wt2)/sum(Contrasts$Wt2)
  return(list(Est.SCSWT1=TxEst1,Est.SCSWT2=TxEst2,ContrastMat=Contrasts))
}



## Crossover Synthetic Control Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with three elements:
### Est.COSC1: the estimated treatment effect for equal weighting (COSC-1 in the article)
### Est.COSC2: the estimated treatment effect for weighting by inverse-MSPE within equivalence groups (COSC-2 in the article)
### ContrastMat: the matrix of cluster-period-specific treatment effect estimates, as well as the weights used in SC-1 and SC-2.
####             Can be used for alternate weighting schemes

COSC.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc) {
  Pds <- sort(unique(Periods))
  Clusts <- sort(unique(Clusters))
  StartTimesOrd <- StartTimes[order(Clusters)]
  StartTimesByClust <- StartTimesOrd[seq(length(Pds),length(StartTimes),by=length(Pds))]
  names(StartTimesByClust) <- as.character(Clusts)
  COPds <- Pds[-1]
  Index <- Clusters*max(Periods)+Periods
  DataMat <- matrix(data=Outcomes[order(Index)], nrow=length(Pds), 
                    ncol=length(Clusts), byrow=FALSE,
                    dimnames=list(Pds,Clusts))
  ## Note: DataMat here only works if all clusters have same periods ##
  COMat <- ContrastFunc(DataMat[2:length(Pds),], DataMat[1:length(COPds),])
  TypeMat <- COMat
  for (i in 1:(length(COPds))) {
    Pdi <- COPds[i]
    TypeMat[i,] <- rep(0,dim(TypeMat)[2]) + (Pdi >= StartTimesByClust) + (Pdi > StartTimesByClust)
  }
  ContrastPds <- COPds[apply(TypeMat,1,function(x) sum(x==0)) > 0 & apply(TypeMat,1,function(x) sum(x==1)) > 0]
  Contrasts <- data.frame(Period = NULL, Cluster = NULL,Contrast = NULL, MSPE = NULL,
                          ClustMin = NULL)
  
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    COMatPd <- COMat[i,]
    TypeMatPd <- TypeMat[i,]
    n1 <- sum(TypeMatPd==1)
    n0 <- sum(TypeMatPd==0)
    IntClusts <- as.character(names(TypeMatPd[TypeMatPd==1]))
    ContClusts <- as.character(names(TypeMatPd[TypeMatPd==0]))
    ClustMinVec <- rep(min(n1,n0),n1)
    
    if (Pd == min(COPds)) {
      SCOutcomes <- rep(mean(COMatPd[TypeMatPd==0]),n1)
      MSPEs <- rep(10^(-8), n1)
    } else {
      PrevPds <- COPds[COPds < Pd]
      COMati <- COMat[as.character(PrevPds),,drop=FALSE]
      
      MSPEs <- rep(10^(-8), n1) ## Truncates MSPE at minimum of 10^(-8) to prevent inverse of zero issue
      if (n0 == 1) {
        SCOutcomes <- rep(COMatPd[TypeMatPd==0], n1)
        for (j in 1:n1) {
          MSPEs[j] <- max(MSPEs[j], 
                          mean((COMati[as.character(PrevPds),TypeMatPd==0]-COMati[as.character(PrevPds),IntClusts[j]])^2))
        }
      } else {
        SCOutcomes <- rep(0, n1)
        X0i <- COMati[,TypeMatPd==0,drop=FALSE]
        for (j in 1:n1) {
          X1ij <- COMati[,IntClusts[j],drop=FALSE]
          capture.output(SynthOut <- tryCatch({synth(X1=X1ij, X0=X0i, Z0=X0i, Z1=X1ij)},
                                              error=function(err) {1}), file="/dev/null")
          if (is.numeric(SynthOut)) {
            SCOutcomes[j] <- mean(COMatPd[TypeMatPd==0])
            MSPEs[j] <- max(MSPEs[j],mean((X1ij - 
                                             apply(X=X0i,MARGIN=1,FUN=mean))^2))
          } else {
            MSPEs[j] <- max(MSPEs[j],as.numeric(SynthOut$loss.v))
            SCOutcomes[j] <- as.numeric(SynthOut$solution.w) %*% COMatPd[TypeMatPd==0]
          }
        }
      }
    }
    ContrastVec <- as.numeric(COMatPd[TypeMatPd==1] - SCOutcomes)
    Contrasts <- data.frame(Period = append(Contrasts$Period, rep(Pd, n1)),
                            Cluster = append(Contrasts$Cluster, as.numeric(IntClusts)),
                            Contrast = append(Contrasts$Contrast, ContrastVec),
                            MSPE = append(Contrasts$MSPE, MSPEs),
                            ClustMin = append(Contrasts$ClustMin, ClustMinVec))
  }
  Contrasts$SumMSPE <- rep(0,length(Contrasts$Contrast))
  for (i in 1:length(Contrasts$Contrast)) {
    Periodi <- Contrasts$Period[i]
    Contrasts$SumMSPE[i] <- sum(1/Contrasts$MSPE[Contrasts$Period==Periodi])
  }
  Contrasts$Wt2 <- (Contrasts$MSPE*Contrasts$SumMSPE)^(-1)

  TxEst1 <- mean(Contrasts$Contrast)
  TxEst2 <- sum(Contrasts$Contrast*Contrasts$Wt2)/sum(Contrasts$Wt2)
  
  return(list(Est.COSC1=TxEst1,Est.COSC2=TxEst2,ContrastMat=Contrasts))
}


## Generic Ensemble Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Estimates: A named vector or list of treatment effect estimates that includes those to be used in the ensemble
## Prefix: The prefix that denotes the names of Estimates. For the above methods, "Est."
## Types: A vector of the names (excluding Prefix) of the Estimates to be used in the Method
## Weights: A vector of the weights to be given to the Types. Can sum to 1 or not; if not, weights will be scaled to sum to 1
# Outputs:
## The value of the ensemble estimate with given Weights to each of the given Types.

Ens.Effect.Est <- function(Estimates, Prefix, Types, Weights) {
  if (sum(Weights) != 1) {
    print("Warning: Weights do not sum to 1. Adjusting to maintain proportions given.")
    Weights = Weights/sum(Weights)
  }
  Vals <- Estimates[paste0(Prefix,Types)]
  out <- sum(Vals*Weights)
  return(out)
}

## Specific Ensemble Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Estimates: A named vector or list of treatment effect estimates that includes those to be used in the ensemble
###           Must include values for "SCSWT2" and "CO.CtWt"
## Prefix: The prefix that denotes the names of Estimates. For the above methods, "Est."
# Output:
## A list with one element:
### Est.ENS: the estimator for the ENS method described in the article,
####         an unweighted mean of SC-2 and CO-2
Ens1.Effect.Est <- function(Estimates, Prefix) {
  outEst <- Ens.Effect.Est(Estimates, Prefix, Types=c("SCSWT2","CO.CtWt"), Weights=c(1/2,1/2))
  return(list(Est.ENS=outEst))
}


### Methods for Inference for Above Estimators ###

## Permutated Effect Estimates ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
# Note: there is random permutation generation in this method. Set a seed before running for replicability.
# Outputs:
## A list with one element for each Type:
### PermEsts.[Type]: A vector of [NumPerms] effect estimates by the [Type] Method,
###                  one for each permutation conducted, with each permutation
###                  under the null hypothesis of no treatment effect.

Perm.Effects <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                         family, link, NumPerms, 
                         Type=c("MEM","CPI","NPWP",
                                "SCSWT1","SCSWT2",
                                "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                                "COSC1","COSC2","ENS")) {
  Starts <- as.numeric(tapply(StartTimes, Clusters, FUN=min))
  NumPds <- as.numeric(tapply(Periods, Clusters, FUN=length))
  
  RandStartsMat <- replicate(n=NumPerms, 
                             expr=rep(sample(x=Starts, size=length(Starts), replace=FALSE), 
                                      times=NumPds))
  
  outlist <- NULL
  
  if ("MEM" %in% Type) {
    PermEsts.MEM <- apply(X=RandStartsMat,
                          MARGIN=2,
                          FUN=function(x) MEM.Effect.Est(Periods, Outcomes, Clusters, x, 
                                                         Indivs,family,link)$Est.MEM)
    outlistint <- list(PermEsts.MEM = PermEsts.MEM)
    outlist <- append(outlist, outlistint)
  }
  
  if ("CPI" %in% Type) {
    PermEsts.CPI <- apply(X=RandStartsMat,
                          MARGIN=2,
                          FUN=function(x) CPI.Effect.Est(Periods, Outcomes, Clusters, x, 
                                                         Indivs,family,link)$Est.CPI)
    outlistint <- list(PermEsts.CPI = PermEsts.CPI)
    outlist <- append(outlist, outlistint)
  }
  
  if (link=="logit") {
    ContrastFunc <- logORFunc
  } else if (link=="log") {
    ContrastFunc <- logRRFunc
  } else if (link=="identity") {
    ContrastFunc <- RDFunc
  } else {
    print("Warning: Contrast function unknown for this family. Using difference for non-parametric methods.")
    ContrastFunc <- RDFunc
  }

  if ("NPWP" %in% Type) {
    PermEsts.NPWP <- apply(X=RandStartsMat,
                          MARGIN=2,
                          FUN=function(x) NPWP.Effect.Est(Periods, Outcomes, StartTimes=x, 
                                                          ContrastFunc=ContrastFunc)$Est.NPWP)
    outlistint <- list(PermEsts.NPWP = PermEsts.NPWP)
    outlist <- append(outlist, outlistint)
  }
  
  if (sum(c("SCSWT1","SCSWT2") %in% Type) > 0) {
    AllSCSWTs <- c("SCSWT1","SCSWT2")
    EstSCSWTs <- c("Est.SCSWT1","Est.SCSWT2")
    SCSWT.Vars <- EstSCSWTs[AllSCSWTs %in% Type]
    PermEsts.SCSWT <- apply(X=RandStartsMat,
                           MARGIN=2,
                           FUN=function(x) as.numeric(SCSWT.Effect.Est(Periods, Outcomes, Clusters, StartTimes=x, 
                                                           ContrastFunc=ContrastFunc)[SCSWT.Vars]))
    outlistint <- NULL
    namesprior <- NULL
    for (i in 1:length(SCSWT.Vars)) {
      namei <- paste0("PermEsts.",substr(SCSWT.Vars[i],5,10))
      outlistint <- append(outlistint, list(PermEsts.SCSWT[i,]))
      namesprior <- c(namesprior,namei)
    }
    names(outlistint) <- namesprior
    outlist <- append(outlist, outlistint)
  }
  
  if (sum(c("CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt") %in% Type) > 0) {
    AllCOs <- c("CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt")
    EstCOs <- c("Est.CO.Ctrl","Est.CO.Both","Est.CO.CtWt","Est.CO.BoWt")
    COTypes <- c("Ctrl","Both","CtWt","BoWt")
    CO.Vars <- EstCOs[AllCOs %in% Type]
    CtrlTypes <- COTypes[AllCOs %in% Type]
    PermEsts.CO <- apply(X=RandStartsMat,
                         MARGIN=2,
                         FUN=function(x) as.numeric(CO.Effect.Est(Periods, Outcomes, StartTimes=x,
                                                                  ContrastFunc=ContrastFunc,
                                                                  CtrlType=CtrlTypes)[CO.Vars]))
    outlistint <- NULL
    namesprior <- NULL
    if (length(CO.Vars)==1) {
      namei <- paste0("PermEsts.",substr(CO.Vars,5,11))
      outlistint <- append(outlistint, list(PermEsts.CO))
      namesprior <- c(namesprior,namei)
    } else {
      for (i in 1:length(CO.Vars)) {
        namei <- paste0("PermEsts.",substr(CO.Vars[i],5,11))
        outlistint <- append(outlistint, list(PermEsts.CO[i,]))
        namesprior <- c(namesprior,namei)
      }
    }
    names(outlistint) <- namesprior
    outlist <- append(outlist, outlistint)
  }
  
  if (sum(c("COSC1","COSC2") %in% Type) > 0) {
    AllCOSCs <- c("COSC1","COSC2")
    EstCOSCs <- c("Est.COSC1","Est.COSC2")
    COSC.Vars <- EstCOSCs[AllCOSCs %in% Type]
    PermEsts.COSC <- apply(X=RandStartsMat,
                            MARGIN=2,
                            FUN=function(x) as.numeric(COSC.Effect.Est(Periods, Outcomes, Clusters, StartTimes=x, 
                                                                        ContrastFunc=ContrastFunc)[COSC.Vars]))
    outlistint <- NULL
    namesprior <- NULL
    for (i in 1:length(COSC.Vars)) {
      namei <- paste0("PermEsts.",substr(COSC.Vars[i],5,9))
      outlistint <- append(outlistint, list(PermEsts.COSC[i,]))
      namesprior <- c(namesprior,namei)
    }
    names(outlistint) <- namesprior
    outlist <- append(outlist, outlistint)
  }
  
  if ("ENS" %in% Type) {
    estmatrix <- matrix(unlist(outlist), nrow=length(outlist), ncol=NumPerms, 
                        byrow=TRUE, dimnames <- list(names(outlist)))
    outlistint <- apply(X=estmatrix, MARGIN=2,
                        FUN=function(x) Ens1.Effect.Est(x, Prefix="PermEsts.")$Est.ENS)
    outlist <- append(outlist, list(PermEsts.ENS=outlistint))
  }
  
  return(outlist)
}

## Hypothesis Testing of H_0: beta = 0 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
## Estimate: a named vector or list of the observed (non-permuted) effect estimates,
###          with names in the form Est.[Type].
## Alternative: a character string specifying the type of alternative hypothesis:
###             "Both" to test H_a: beta != 0; "Greater" to test H_a: beta > 0; "Less" to test H_a: beta < 0
# Note: there is random permutation generation in this method (through calling Perm.Effects).
##      Set a seed before running for replicability.
# Outputs:
## A named vector with one element for each Type:
### PVal.[Type]: the p-value for the specified permutation-based hypothesis test
####             for the method [Type] using [NumPerms] permutations of the data
####             and the observed effect estimate [Estimate]

Null.Test <- function(Periods, Outcomes, Clusters, StartTimes, 
                      Indivs=NA, family, link, NumPerms, 
                      Type=c("MEM","CPI","NPWP",
                             "SCSWT1","SCSWT2",
                             "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                             "COSC1","COSC2","ENS"),
                      Estimate, Alternative="Both") {
  Perm.Ests <- Perm.Effects(Periods, Outcomes, Clusters, StartTimes, Indivs, 
                            family, link, NumPerms, Type)
  if (Alternative=="Both") {
    PValFunc <- function(Truth,Perms) mean(abs(Truth) <= abs(Perms), na.rm=TRUE)
  } else if (Alternative=="Greater") {
    PValFunc <- function(Truth,Perms) mean(Truth >= Perms, na.rm=TRUE)
  } else if (Alternative=="Less") {
    PValFunc <- function(Truth,Perms) mean(Perms <= Truth, na.rm=TRUE)
  } else {
    stop("Please enter 'Both', 'Greater', or 'Less' for the variable Alternative")
  }
  
  PVal <- NULL
  AllTypes <- c("MEM","CPI","NPWP",
                "SCSWT1","SCSWT2",
                "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                "COSC1","COSC2","ENS")
  for (i in 1:length(Type)) {
    outstr1 <- paste0("PVal.",Type[i])
    Perms <- unlist(Perm.Ests[paste0("PermEsts.",Type[i])])
    Truth <- unlist(Estimate[paste0("Est.",Type[i])])
    PVali <- PValFunc(Truth,Perms)
    namesprior1 <- names(PVal)
    PVal <- append(PVal, PVali)
    names(PVal) <- c(namesprior1, outstr1)
  }
  return(PVal)
}


### Full Analysis with Estimation ###
### and Hypothesis Tests for Specified Null Values ###
### for Selected Methods ###

# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Treatments: A vector of treatment indicators (0 for control, 1 for intervention)
###            in the same order as Periods.
#### Note: can specify either StartTimes or Treatments
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform for inference. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
## NullVals: a vector of null hypothesis values to test. By default, only tests 0 to get the p-value.
###          Can specify additional NullVals to test whether the 1-alpha CIs include those values.
## Alternative: a character string specifying the type of alternative hypothesis:
###             "Both" to test H_a: beta != 0; "Greater" to test H_a: beta > 0; "Less" to test H_a: beta < 0
## Alpha: the significance level at which to perform hypothesis tests and CI inclusion checks
# Note: there is random permutation generation in this method (through calling Perm.Effects).
##      Set a seed before running for replicability.
# Outputs:
## A data frame with one row for each component of NullVals. Each row has the elements:
### NullVal: the null hypothesis value being tested
### For each Type:
#### Est.[Type]: the treatment effect estimate for the [Type] method
#### PVal.[Type]: the p-value for H_0: beta = NullVal for the [Type] estimate
#### Res.[Type]: "Accept" if PVal.[Type] > alpha; otherwise, "Reject".
#####            "Accept" indicates failing to reject the null hypothesis in the test or 
#####             inclusion of NullVal in the (1-alpha) CI
#####            "Reject" indicates rejecting the null hypothesis in the test or 
#####             exclusion of NullVal from the (1-alpha) CI
#### AsyPVal.[Type]: only for MEM and CPI types; the same as PVal 
#####                but using asymptotic inference instead of permutation-based exact inference
#### AsyRes.[Type]: only for MEM and CPI types; the same as Res 
#####                but using asymptotic inference instead of permutation-based exact inference

SWT.Permutation.Analysis <- function(Periods, Outcomes, Clusters, 
                                     StartTimes=NA, Treatments=NA, 
                                     Indivs=NA,
                                     family, link, 
                                     NumPerms, 
                                     Type=c("MEM","CPI","NPWP",
                                            "SCSWT1","SCSWT2",
                                            "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                                            "COSC1","COSC2","ENS"), 
                                     NullVals=0,
                                     Alternative="Both", Alpha=.05) {
  if (sum(is.na(StartTimes)) > 0) {
    if (sum(is.na(Treatments)) > 0) {
      stop("Must Specify Either StartTimes or Vector of Treatment Indicators")
    } else {
      StartTimes <- StartTimeVec(Periods=Periods, Clusters=Clusters, Trts=Treatments)
    }
  }
  
  if (link=="logit") {
    ContrastFunc <- logORFunc
    ContrastInvApply <- logORInvApply
  } else if (link=="log") {
    ContrastFunc <- logRRFunc
    ContrastInvApply <- logRRInvApply
  } else if (link=="identity") {
    ContrastFunc <- RDFunc
    ContrastInvApply <- RDInvApply
  } else {
    print("Warning: Contrast function unknown for this family. Using difference for non-parametric methods.")
    ContrastFunc <- RDFunc
    ContrastInvApply <- RDInvApply
  }
  
  OutDF <- NULL
  for (m in 1:length(NullVals)) {
    if (NullVals[m] == 0) {
      AdjOutcomes <- Outcomes
    } else {
      AdjOutcomes <- ifelse(Periods >= StartTimes, pmax(0,ContrastInvApply(Outcomes,NullVals[m])), Outcomes)
    }
    
    Ests <- NULL
    if ("MEM" %in% Type) {
      MEM.res <- MEM.Effect.Est(Periods, AdjOutcomes, Clusters, StartTimes, Indivs, family, link)
      Ests <- append(Ests, list(Est.MEM=MEM.res$Est.MEM))
      Ests <- append(Ests, list(AsyPVal.MEM=MEM.res$PVal.MEM))
      }
    
    if ("CPI" %in% Type) {
      CPI.res <- CPI.Effect.Est(Periods, AdjOutcomes, Clusters, StartTimes, Indivs, family, link)
      Ests <- append(Ests, list(Est.CPI=CPI.res$Est.CPI))
      Ests <- append(Ests, list(AsyPVal.CPI=CPI.res$PVal.CPI))
      }
    
    if ("NPWP" %in% Type) {
      NPWP.res <- NPWP.Effect.Est(Periods, AdjOutcomes, StartTimes, ContrastFunc)
      Ests <- append(Ests, list(Est.NPWP=NPWP.res$Est.NPWP))
    }
    
    if (sum(c("SCSWT1","SCSWT2") %in% Type) > 0) {
      AllSCSWTs <- c("SCSWT1","SCSWT2")
      EstSCSWTs <- c("Est.SCSWT1","Est.SCSWT2")
      SCSWT.Vars <- EstSCSWTs[AllSCSWTs %in% Type]
      SCSWT.res <- SCSWT.Effect.Est(Periods, AdjOutcomes, Clusters, StartTimes, ContrastFunc)
      Ests <- append(Ests, SCSWT.res[SCSWT.Vars])
    }
    
    if (sum(c("CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt") %in% Type) > 0) {
      AllCOs <- c("CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt")
      EstCOs <- c("Est.CO.Ctrl","Est.CO.Both","Est.CO.CtWt","Est.CO.BoWt")
      COTypes <- c("Ctrl","Both","CtWt","BoWt")
      CO.Vars <- EstCOs[AllCOs %in% Type]
      CtrlTypes <- COTypes[AllCOs %in% Type]
      CO.res <- CO.Effect.Est(Periods, AdjOutcomes, StartTimes, ContrastFunc, 
                              CtrlType <- CtrlTypes)
      Ests <- append(Ests, CO.res[CO.Vars])
    }
    
    if (sum(c("COSC1","COSC2") %in% Type) > 0) {
      AllCOSCs <- c("COSC1","COSC2")
      EstCOSCs <- c("Est.COSC1","Est.COSC2")
      COSC.Vars <- EstCOSCs[AllCOSCs %in% Type]
      COSC.res <- COSC.Effect.Est(Periods, AdjOutcomes, Clusters, StartTimes, ContrastFunc)
      Ests <- append(Ests, COSC.res[COSC.Vars])
    }
    
    if ("ENS" %in% Type) {
      Ests <- append(Ests, Ens1.Effect.Est(unlist(Ests), Prefix="Est."))
    }
    
    if (NumPerms == 0) {
      OutDF = rbind(OutDF, data.frame(Ests))
    } else {
      Test.Res <- Null.Test(Periods, AdjOutcomes, Clusters, StartTimes, 
                            Indivs, family, link, NumPerms, 
                            Type, Estimate=Ests, Alternative)
      outrow <- data.frame(NullVal=NullVals[m])
      for (i in 1:length(Type)) {
        Estname <- paste0("Est.",Type[i])
        Pname <- paste0("PVal.",Type[i])
        Resname <- paste0("Res.",Type[i])
        appVals <- data.frame(a=unlist(Ests[Estname]),
                              b=Test.Res[Pname],
                              c=ifelse(Test.Res[Pname] <= Alpha, "Reject", "Accept"))
        names(appVals) <- c(Estname, Pname, Resname)
        outrow <- cbind(outrow, appVals)
        if (Type[i] %in% c("MEM","CPI")) {
          Pname2 <- paste0("AsyPVal.",Type[i])
          Resname2 <- paste0("AsyRes.",Type[i])
          appVals2 <- data.frame(a=unlist(Ests[Pname2]),
                                b=ifelse(unlist(Ests[Pname2]) <= Alpha, "Reject", "Accept"))
          names(appVals2) <- c(Pname2, Resname2)
          outrow <- cbind(outrow, appVals2)
        }
      }
      row.names(outrow) <- m
      OutDF = rbind(OutDF, outrow)
    }
  }
  return(OutDF)
}