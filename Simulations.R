######################################################
### Code to Replicate Simulations                  ###
### in Kennedy-Shaffer et al. 2019                 ###
### Update: 9/15/2019                              ###
### Contact: lee_kennedyshaffer@g.harvard.edu      ###
### See Latest Update at:                          ###
### https://github.com/leekshaffer/SW-CRT-analysis ###
######################################################

source("SW-CRT Analysis Methods.R")

#### Generating Random Simulated Data Sets ####

### Function to Create a Random Data Set with Risk Difference RD ###
MakeRandomData <- function(mu, thetas, prob1, NumClusts, tau, nu, NumInds, RD) {
  NumPds <- dim(thetas)[2]
  time <- seq(1,NumPds)
  switches <- sample(rep(seq(2,NumPds,by=1),NumClusts/(NumPds-1)), size=NumClusts, replace=FALSE)
  types <- 2 - rbinom(n=NumClusts, size=1, prob=prob1)
  DF <- data.frame(Cluster=NULL, Period=NULL, Trt=NULL, Outcome=NULL)
  clustFX <- rnorm(n=NumClusts, mean=mu, sd=tau)
  for (i in 1:NumClusts) {
    Cluster <- rep(i, NumPds)
    Period <- time
    Trt <- ifelse(Period >= switches[i], 1, 0)
    ControlProbs <- clustFX[i] + thetas[(types[i]),] + rnorm(n=NumPds, mean=0, sd=nu)
    TrueProbs <- ControlProbs + Trt*RD
    TrueProbsCut <- pmin(pmax(TrueProbs,0),1)
    Outcome <- rbinom(n=NumPds, size=NumInds, prob=TrueProbsCut)/NumInds
    DF <- rbind(DF, data.frame(Cluster, Period, Trt, Outcome))
  }
  return(DF)
}

### Function to Create NumSims Random Data Sets with Risk Difference RD, with set seed ###
MakeRandomDatas <- function(mu, thetas, prob1, NumClusts, tau, nu, NumInds, RD, NumSims, seed=NA) {
  if (is.numeric(seed)) {
    set.seed(seed)}
  Array <- replicate(n=NumSims, 
                     expr=as.matrix(MakeRandomData(mu, thetas, prob1, NumClusts, 
                                                   tau, nu, NumInds, RD)),
                     simplify="array")
  return(Array)
}

### Generating 1,000 Random Data Sets for each of 12 Scenarios in Simulation 1 ###
### for both 7 clusters and 21 clusters. Note: only 7-cluster scenarios are ###
### discussed in the article. These are scenarios 1,2,3,4,9,10,11,12,17,18,19,20. ###
mu_1 <- .3
thetas_a <- c(0, .08, .18, .29, .3, .27, .2, .13)
thetas_b <- c(0, .02, .03, .07, .13, .19, .27, .3)
thetas_1 <- matrix(c(thetas_a,thetas_b),nrow=2,ncol=length(thetas_a),byrow=TRUE)
tau_1 <- .06
K_1 <- 100
Varying <- data.frame(p1s = rep(c(1,1,.5,.5),6),
                      nus = rep(c(0,.01,0,.01),6),
                      Is = rep(c(rep(7,4),rep(21,4)),3),
                      betas = c(rep(-.2,8),rep(-.1,8),rep(0,8)))
save(Varying, file="Sim1_Data/Varying.Rda")

set.seed(11132019)
for (i in 1:(dim(Varying)[1])) {
  outnam <- paste0("SimData_Sim1_",i)
  assign(outnam, MakeRandomDatas(mu=mu_1, thetas=thetas_1, prob1=Varying$p1s[i], 
                                 NumClusts=Varying$Is[i], tau=tau_1,
                                 nu=Varying$nus[i], NumInds=K_1, RD=Varying$betas[i],
                                 NumSims=1000))
  
  save(list=c(outnam), file=paste0("Sim1_Data/",outnam,".Rda"))
}


###  Helper Functions for OR ###
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))


### Function to Create a Random Data Set with Log Odds Ratio logOR ###
MakeRandomData <- function(mu, thetas, prob1, NumClusts, tau, nu, NumInds, logOR) {
  NumPds <- dim(thetas)[2]
  time <- seq(1,NumPds)
  switches <- sample(rep(seq(2,NumPds,by=1),NumClusts/(NumPds-1)), size=NumClusts, replace=FALSE)
  types <- 2 - rbinom(n=NumClusts, size=1, prob=prob1)
  DF <- data.frame(Cluster=NULL, Period=NULL, Trt=NULL, Outcome=NULL)
  clustFX <- rnorm(n=NumClusts, mean=mu, sd=tau)
  for (i in 1:NumClusts) {
    Cluster <- rep(i, NumPds)
    Period <- time
    Trt <- ifelse(Period >= switches[i], 1, 0)
    ControlProbs <- clustFX[i] + thetas[(types[i]),] + rnorm(n=NumPds, mean=0, sd=nu)
    TrueProbs <- ControlProbs + Trt*logOR
    TrueProbsVal <- expit(TrueProbs)
    Outcome <- rbinom(n=NumPds, size=NumInds, prob=TrueProbsVal)/NumInds
    DF <- rbind(DF, data.frame(Cluster, Period, Trt, Outcome))
  }
  return(DF)
}

### Function to Create NumSims Random Data Sets with Log Odds Ratio logOR, with set seed ###
MakeRandomDatas <- function(mu, thetas, prob1, NumClusts, tau, nu, NumInds, logOR, NumSims, seed=NA) {
  if (is.numeric(seed)) {
    set.seed(seed)}
  Array <- replicate(n=NumSims, 
                     expr=as.matrix(MakeRandomData(mu, thetas, prob1, NumClusts, 
                                                   tau, nu, NumInds, logOR)),
                     simplify="array")
  return(Array)
}


### Generating 1,000 Random Data Sets for each of 12 Scenarios in Simulation 2 ###
### for both 7 clusters and 21 clusters. Note: only 7-cluster scenarios are ###
### discussed in the article. These are scenarios 1,2,3,4,9,10,11,12,17,18,19,20. ###
mu_1 <- logit(.3)
thetas_a <- log(c(1, 1.43, 2.15, 3.36, 3.5, 3.09, 2.33, 1.76))
thetas_b <- log(c(1, 1.1, 1.15, 1.37, 1.76, 2.24, 3.09, 3.5))
thetas_1 <- matrix(c(thetas_a,thetas_b),nrow=2,ncol=length(thetas_a),byrow=TRUE)
tau_1 <- .1
K_1 <- 100
Varying <- data.frame(p1s = rep(c(1,1,.5,.5),6),
                      nus = rep(c(0,.01,0,.01),6),
                      Is = rep(c(rep(7,4),rep(21,4)),3),
                      betas = c(rep(log(.4),8),rep(log(.66),8),rep(0,8)))
save(Varying, file="Sim2_Data/Varying.Rda")

set.seed(948315)
for (i in 1:(dim(Varying)[1])) {
  outnam <- paste0("SimData_Sim2_",i)
  assign(outnam, MakeRandomDatas(mu=mu_1, thetas=thetas_1, prob1=Varying$p1s[i], 
                                 NumClusts=Varying$Is[i], tau=tau_1,
                                 nu=Varying$nus[i], NumInds=K_1, logOR=Varying$betas[i],
                                 NumSims=1000))
  
  save(list=c(outnam), file=paste0("Sim2_Data/",outnam,".Rda"))
}


#### Data for Figures 2 and 7 ####
subgraphs <- c("a","b","c","d")
for (i in 17:20) {
  nameRD <- paste0("SimData_Sim1_",i)
  nameOR <- paste0("SimData_Sim2_",i)
  DF_RD <- data.frame(get(nameRD)[,,2])
  DF_OR <- data.frame(get(nameOR)[,,5])
  Data_Fig2 <- reshape(DF_RD[,c("Cluster","Period","Outcome")],
                   v.names="Outcome", idvar="Cluster", timevar="Period",
                   direction="wide")
  names(Data_Fig2) <- c("Cluster","1","2","3","4","5","6","7","8")
  assign(x=paste0("Data_Fig2",subgraphs[(i-16)]), value=Data_Fig2)
  save(list=c(paste0("Data_Fig2",subgraphs[(i-16)])),
       file=paste0("Fig_Data/Data_Fig2",subgraphs[(i-16)],".Rda"))
  Data_Fig7 <- reshape(DF_OR[,c("Cluster","Period","Outcome")],
                       v.names="Outcome", idvar="Cluster", timevar="Period",
                       direction="wide")
  names(Data_Fig7) <- c("Cluster","1","2","3","4","5","6","7","8")
  assign(x=paste0("Data_Fig7",subgraphs[(i-16)]), value=Data_Fig7)
  save(list=c(paste0("Data_Fig7",subgraphs[(i-16)])),
       file=paste0("Fig_Data/Data_Fig7",subgraphs[(i-16)],".Rda"))
}


#### Analysis of Simulation 1 ####

### Function to analyze Simulation 1 Data from a Matrix of Trial Results ###
AnalyzeMatrix <- function(TxEff, NumPerms, Matrix, 
                          Type=c("MEM","CPI","NPWP",
                                 "SCSWT1","SCSWT2",
                                 "CO.Ctrl","CO.Both","CO.CtWt",
                                 "COSC1","COSC2","ENS")) {
  if (TxEff==0) {
    Nulls <- 0
  } else {
    Nulls <- c(0,TxEff)
  }
  res <- SWT.Permutation.Analysis(Periods = Matrix[,'Period'],
                                  Outcomes = Matrix[,'Outcome'],
                                  Clusters = Matrix[,'Cluster'],
                                  StartTimes = NA,
                                  Treatments = Matrix[,'Trt'],
                                  Indivs = K_1,
                                  family=binomial,
                                  link="identity",
                                  NumPerms = NumPerms, 
                                  Type,
                                  NullVals = Nulls,
                                  Alternative="Both", Alpha=.05)
  
  if (NumPerms == 0) {
    outrow <- res[1,2:(length(res[1,])-1)]
  }
  else {
    outrow <- data.frame(a=1)
    for (i in 1:length(Type)) {
      Typei <- Type[i]
      Estname <- paste0("Est.",Typei)
      Pname <- paste0("PVal.",Typei)
      Resname <- paste0("Res.",Typei)
      Accname <- paste0("TestAcc.",Typei)
      CIname <- paste0("CICov.",Typei)
      Rowi <- data.frame(a=res[1,Estname], b=res[1,Pname], 
                         c=ifelse(res[1,Resname]=="Accept",1,0),
                         d=ifelse(res[length(Nulls),Resname]=="Accept",1,0))
      names(Rowi) <- c(Estname, Pname, Accname, CIname)
      outrow <- cbind(outrow, Rowi)
      if (Typei %in% c("MEM","CPI")) {
        Pname2 <- paste0("Asy",Pname)
        Resname2 <- paste0("Asy",Resname)
        Accname2 <- paste0("Asy",Accname)
        CIname2 <- paste0("Asy",CIname)
        Rowi2 <- data.frame(b=res[1,Pname2], 
                            c=ifelse(res[1,Resname2]=="Accept",1,0),
                            d=ifelse(res[length(Nulls),Resname2]=="Accept",1,0))
        names(Rowi2) <- c(Pname2, Accname2, CIname2)
        outrow <- cbind(outrow, Rowi2)
      }
    }
    outrow <- outrow[,2:(length(outrow[1,]))]
  }
  return(outrow)
}

### Analysis of Simulation 1 Data ###
load("Sim1_Data/Varying.Rda")

for (Scenario in c(1:4,9:12,17:20)) {
  DataNam <- paste0("SimData_Sim1_",Scenario)
  load(paste0("Sim1_Data/",DataNam,".Rda"))
  ResNam <- paste0("SimRes_Sim1_",Scenario)
  assign(ResNam, NULL)
  for (sim in 1:1000) {
    set.seed(1000*Scenario^2+sim)
    SimMat <- get(DataNam)[,,sim]
    Res <- AnalyzeMatrix(TxEff=Varying$betas[Scenario], NumPerms=500,
                         Matrix=SimMat,
                         Type=c("MEM","CPI","NPWP",
                                "SCSWT1","SCSWT2",
                                "CO.Ctrl","CO.Both","CO.CtWt",
                                "COSC1","COSC2","ENS"))
    Res2 <- cbind(data.frame(ArrayNum=sim), Res)
    assign(ResNam, rbind(get(ResNam), Res2))
    rm(Res,Res2)
  }
  save(list=c(ResNam), file=paste0("Sim1_Res/",ResNam,".Rda"))
  rm(list=c(ResNam,DataNam))
}

### Figures 3--6 Data ###
TypeNames <- c("MEM","CPI","NPWP",
               "SCSWT1","SCSWT2",
               "CO.Ctrl","CO.Both","CO.CtWt",
               "COSC1","COSC2","ENS",
               "MEM-a","CPI-a")
Data_Full_Sim1 <- data.frame(Scen=rep(c(1:4,9:12,17:20),each=length(TypeNames)),
                        Type=rep(TypeNames,1),
                        Mean=rep(NA,length(TypeNames)*12),
                        SD=rep(NA,length(TypeNames)*12),
                        Power=rep(NA,length(TypeNames)*12),
                        Coverage=rep(NA,length(TypeNames)*12))

Est.Var.Names <- c(paste0("Est.",TypeNames[1:11]),"Est.MEM","Est.CPI")
Test.Var.Names <- c(paste0("TestAcc.",TypeNames[1:11]),"AsyTestAcc.MEM","AsyTestAcc.CPI")
Cov.Var.Names <- c(paste0("CICov.",TypeNames[1:11]),"AsyCICov.MEM","AsyCICov.CPI")

for (Scen in c(1:4,9:12,17:20)) {
  ResNam <- paste0("SimRes_Sim1_",Scen)
  load(paste0("Sim1_Res/",ResNam,".Rda"))
  ResSub <- get(ResNam)[,Est.Var.Names]
  Data_Full_Sim1$Mean[Data_Full_Sim1$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
  Data_Full_Sim1$SD[Data_Full_Sim1$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) sd(x, na.rm=TRUE))
  ResSub <- get(ResNam)[,Test.Var.Names]
  Data_Full_Sim1$Power[Data_Full_Sim1$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) 1-mean(x, na.rm=TRUE))
  ResSub <- get(ResNam)[,Cov.Var.Names]
  Data_Full_Sim1$Coverage[Data_Full_Sim1$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
  rm(list=c(ResNam,"ResSub"))
}

## Figure 3 Data ##
Data_Fig3a <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% 9:12,
                             c("Scen","Type","Mean","SD")]
Data_Fig3a$MinEst <- Data_Fig3a$Mean - 0.5*Data_Fig3a$SD
Data_Fig3a$MaxEst <- Data_Fig3a$Mean + 0.5*Data_Fig3a$SD
Data_Fig3b <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% 17:20,
                             c("Scen","Type","Mean","SD")]
Data_Fig3b$MinEst <- Data_Fig3b$Mean - 0.5*Data_Fig3b$SD
Data_Fig3b$MaxEst <- Data_Fig3b$Mean + 0.5*Data_Fig3b$SD

## Figure 4 Data ##
Data_Fig4 <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% 17:20,
                            c("Scen","Type","Power")]
Data_Fig4$TIE <- Data_Fig4$Power

## Figure 5 Data ##
Data_Fig5a <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% 9:12,
                             c("Scen","Type","Coverage")]
Data_Fig5b <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% 17:20,
                            c("Scen","Type","Coverage")]

## Figure 6 Data ##
Data_Fig6 <- Data_Full_Sim1[Data_Full_Sim1$Scen %in% c(9:12,17:20),
                            c("Scen","Type","Power")]


## Figure 12 Data ##
load(paste0("Sim1_Res/SimRes_Sim1_20.Rda"))
Data_Fig12 <- cov(x=SimRes_Sim1_20[,Est.Var.Names[1:10]],use="pairwise.complete.obs")

rm(SimRes_Sim1_20)
rm(Varying)


save(Data_Fig3a, file="Fig_Data/Data_Fig3a.Rda")
save(Data_Fig3b, file="Fig_Data/Data_Fig3b.Rda")
save(Data_Fig4, file="Fig_Data/Data_Fig4.Rda")
save(Data_Fig5a, file="Fig_Data/Data_Fig5a.Rda")
save(Data_Fig5b, file="Fig_Data/Data_Fig5b.Rda")
save(Data_Fig6, file="Fig_Data/Data_Fig6.Rda")
save(Data_Fig12, file="Fig_Data/Data_Fig12.Rda")



### Analysis of Simulation 2 Data ###
load("Sim2_Data/Varying.Rda")

for (Scen in c(1:4,9:12,17:20)) {
  DataNam <- paste0("SimData_Sim2_",Scen)
  load(paste0("Sim2_Data/",DataNam,".Rda"))
  ResNam <- paste0("SimRes_Sim2_",Scen)
  assign(ResNam, NULL)
  for (sim in 1:1000) {
    set.seed(1000*Scen^2+sim)
    SimMat <- get(DataNam)[,,sim]
    Res <- AnalyzeMatrix(TxEff=Varying$betas[Scen], NumPerms=500,
                         Matrix=SimMat,
                         Type=c("MEM","CPI","NPWP",
                                "SCSWT1","SCSWT2",
                                "CO.Ctrl","CO.Both","CO.CtWt",
                                "COSC1","COSC2","ENS"))
    Res2 <- cbind(data.frame(ArrayNum=sim), Res)
    assign(ResNam, rbind(get(ResNam), Res2))
    rm(Res,Res2)
  }
  save(list=c(ResNam), file=paste0("Sim2_Res/",ResNam,".Rda"))
  rm(list=c(ResNam,DataNam))
}

### Figures 8--11 Data ###
TypeNames <- c("MEM","CPI","NPWP",
               "SCSWT1","SCSWT2",
               "CO.Ctrl","CO.Both","CO.CtWt",
               "COSC1","COSC2","ENS",
               "MEM-a","CPI-a")
Data_Full_Sim2 <- data.frame(Scen=rep(c(1:4,9:12,17:20),each=length(TypeNames)),
                             Type=rep(TypeNames,12),
                             Mean=rep(NA,length(TypeNames)*12),
                             SD=rep(NA,length(TypeNames)*12),
                             Power=rep(NA,length(TypeNames)*12),
                             Coverage=rep(NA,length(TypeNames)*12))

Est.Var.Names <- c(paste0("Est.",TypeNames[1:11]),"Est.MEM","Est.CPI")
Test.Var.Names <- c(paste0("TestAcc.",TypeNames[1:11]),"AsyTestAcc.MEM","AsyTestAcc.CPI")
Cov.Var.Names <- c(paste0("CICov.",TypeNames[1:11]),"AsyCICov.MEM","AsyCICov.CPI")

for (Scen in c(1:4,9:12,17:20)) {
  ResNam <- paste0("SimRes_Sim2_",Scen)
  load(paste0("Sim2_Res/",ResNam,".Rda"))
  ResSub <- get(ResNam)[,Est.Var.Names]
  Data_Full_Sim2$Mean[Data_Full_Sim2$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
  Data_Full_Sim2$SD[Data_Full_Sim2$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) sd(x, na.rm=TRUE))
  ResSub <- get(ResNam)[,Test.Var.Names]
  Data_Full_Sim2$Power[Data_Full_Sim2$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) 1-mean(x, na.rm=TRUE))
  ResSub <- get(ResNam)[,Cov.Var.Names]
  Data_Full_Sim2$Coverage[Data_Full_Sim2$Scen==Scen] <- apply(X=ResSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
  rm(list=c(ResNam,"ResSub"))
}

## Figure 8 Data ##
Data_Fig8a <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% 9:12,
                             c("Scen","Type","Mean","SD")]
Data_Fig8a$MinEst <- Data_Fig8a$Mean - 0.5*Data_Fig8a$SD
Data_Fig8a$MaxEst <- Data_Fig8a$Mean + 0.5*Data_Fig8a$SD
Data_Fig8b <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% 17:20,
                             c("Scen","Type","Mean","SD")]
Data_Fig8b$MinEst <- Data_Fig8b$Mean - 0.5*Data_Fig8b$SD
Data_Fig8b$MaxEst <- Data_Fig8b$Mean + 0.5*Data_Fig8b$SD

## Figure 9 Data ##
Data_Fig9 <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% 17:20,
                            c("Scen","Type","Power")]
Data_Fig9$TIE <- Data_Fig9$Power

## Figure 10 Data ##
Data_Fig10a <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% 9:12,
                             c("Scen","Type","Coverage")]
Data_Fig10b <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% 17:20,
                             c("Scen","Type","Coverage")]

## Figure 11 Data ##
Data_Fig11 <- Data_Full_Sim2[Data_Full_Sim2$Scen %in% c(9:12,17:20),
                            c("Scen","Type","Power")]

## Figure 13 Data ##
load(paste0("Sim2_Res/SimRes_Sim2_20.Rda"))
Data_Fig13 <- cov(x=SimRes_Sim2_20[,Est.Var.Names[1:10]],use="pairwise.complete.obs")

rm(SimRes_Sim2_20)
rm(Varying)

save(Data_Fig8a, file="Fig_Data/Data_Fig8a.Rda")
save(Data_Fig8b, file="Fig_Data/Data_Fig8b.Rda")
save(Data_Fig9, file="Fig_Data/Data_Fig9.Rda")
save(Data_Fig10a, file="Fig_Data/Data_Fig10a.Rda")
save(Data_Fig10b, file="Fig_Data/Data_Fig10b.Rda")
save(Data_Fig11, file="Fig_Data/Data_Fig11.Rda")
save(Data_Fig13, file="Fig_Data/Data_Fig13.Rda")






