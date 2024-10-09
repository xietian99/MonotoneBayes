# Sep 27 - Updated, rerun simulation
##
## This document corresponds to the draft simulation part 1 : Table 2, 3 and Figure 3.
## But we now fix N and N/L; instead previous of N and L.
## We consider N/L. = 2.5 , 5 , 10 and 25
## Then compare different prior choice for alpha's
library(MonotoneBayes)
library(Iso)
library(dplyr)
library(rstan)
source("TrueCurves.R")



N = as.numeric(Sys.getenv("N"))
distX = Sys.getenv("distX")
which_outcome = Sys.getenv("which_outcome")
arrayID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
NovL = as.numeric(Sys.getenv("NovL"))
# distX %in% c("Uniform", "Beta")
# which_outcome %in% c("fX1", "fX2", "fX3", "fX4")
# N = 50, maximun 20 min
# N = 500, 30 min
# N = 500
# distX = "Beta"
# which_outcome = "fX1"
# arrayID = 1
# NovL = 2.5
# Generate Datasets of Uniform Distribution
set.seed(108)
simSeeds = sample(.Machine$integer.max, 500)
set.seed(simSeeds[arrayID])
data = matrix(0, nrow = N, ncol = 6)
if (distX == "Uniform"){
  data[,1] = runif(N, 0, 1)
} else if (distX == "Beta") {
  data[,1] = rbeta(N, shape1 = 2, shape2 = 2)
}

minX = min(data[,1]); maxX = max(data[,1])
data[,1] = (data[,1] - minX) / (maxX - minX)
# Apr 10
# Try next time instead of hard coding
data[,2] = rbinom(N, 1, prob = fX1(data[,1]))
data[,3] = rbinom(N, 1, prob = fX2(data[,1]))
data[,4] = rbinom(N, 1, prob = fX9(data[,1]))
data[,5] = rbinom(N, 1, prob = fX5(data[,1]))
data[,6] = rbinom(N, 1, prob = fX7(data[,1]))
colnames(data) = c("x",paste("fX",c(1,2,9,5,7), sep = ""))
data = as.data.frame(data)
L = floor(N/NovL)
nodes.L = quantile(data$x,  seq(0,1,length.out = L+1))


dt.pred = NULL
basename = c("regHS.fix","regHS.invG","Laplacian","Gaussian")
TimeNames = c(basename, paste(basename, "eq", sep = "_"))
RunTime = rep(0, length = length(TimeNames))
recordTime = rep(0, length(RunTime)+1)
names(RunTime) = TimeNames
x.seq = seq(0,1,length.out = 100)
# Do Quantile Spaced Only for uniformly distributed
# # of Model
# Mine 2 * 2 + Comparison 2 * 2 + Logistic + Pava
c_sq = 1e5
c_alpha = 0.1
c_beta = 100
options(mc.cores = parallelly::availableCores())
t1 = Sys.time()
recordTime[1] =  Sys.time()

model.regHS.fix = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Regularized HS",
                                   fix = T, c_sq = c_sq, Eq.Space = F, nodes = nodes.L,
                                   seed = simSeeds[arrayID], iter = 4e3,
                                   L = L)

recordTime[2] =  Sys.time()
model.regHS.invG = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Regularized HS",
                                    c_alpha = c_alpha, c_beta = c_beta, Eq.Space = F, nodes = nodes.L,
                                    seed = simSeeds[arrayID],  iter = 4e3,
                                    L = L)
recordTime[3] =  Sys.time()
model.lap = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Laplacian",
                             Eq.Space = F, nodes = nodes.L,
                             seed = simSeeds[arrayID],  iter = 4e3,
                             L = L)
recordTime[4] =  Sys.time()
model.norm = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Gaussian",
                              Eq.Space = F, nodes = nodes.L,
                              seed = simSeeds[arrayID],  iter = 4e3,
                              L = L)

recordTime[5] =  Sys.time()
model.regHS.fix.eq = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Regularized HS",
                                      fix = T, c_sq = c_sq, Eq.Space = T,
                                      seed = simSeeds[arrayID],  iter = 4e3,
                                      L = L)

recordTime[6] =  Sys.time()
model.regHS.invG.eq = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Regularized HS",
                                       c_alpha = c_alpha, c_beta = c_beta, Eq.Space = T,
                                       #control = list(adapt_delta = 0.95),
                                       seed = simSeeds[arrayID],  iter = 4e3,
                                       L = L)
recordTime[7]  =  Sys.time()
model.lap.eq = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Laplacian",
                                seed = simSeeds[arrayID],  iter = 4e3,Eq.Space = T,
                                L = L)
recordTime[8]  =  Sys.time()
model.norm.eq = monotoneBayes(X = data[,"x"], Y = data[,which_outcome], prior = "Gaussian",
                                 seed = simSeeds[arrayID],  iter = 4e3,Eq.Space = T,
                                 L = L)
recordTime[9]  =  Sys.time()

RunTime[1:8] = recordTime[2:9] - recordTime[1:8]
rs.regHS.fix.L = cbind(newData,
                        as.data.frame(predictModel(model.regHS.fix, grid.eq = F, node.t = nodes.L, newdata = newData$x)),
                        "Method" = "regHS.fix", "L" = L, "Eq_Space" = 0, "distX" = distX)
rs.regHS.invG.L = cbind(newData,
                         as.data.frame(predictModel(model.regHS.invG, grid.eq = F, node.t = nodes.L, newdata = newData$x)),
                         "Method" = "regHS.invG", "L" = L, "Eq_Space" = 0, "distX" = distX)
rs.lap.L = cbind(newData,
                  as.data.frame(predictModel(model.lap, grid.eq = F, node.t = nodes.L, newdata = newData$x)),
                  "Method" = "Laplacian", "L" = L, "Eq_Space" = 0, "distX" = distX)
rs.norm.L = cbind(newData,
                   as.data.frame(predictModel(model.norm, grid.eq = F, node.t = nodes.L, newdata = newData$x)),
                   "Method" = "Gaussian", "L" = L, "Eq_Space" = 0, "distX" = distX)


rs.regHS.fix.L.eq = cbind(newData,
                           as.data.frame(predictModel(model.regHS.fix.eq, grid.eq = T, newdata = newData$x)),
                           "Method" = "regHS.fix", "L" = L, "Eq_Space" = 1, "distX" = distX)
rs.regHS.invG.L.eq = cbind(newData,
                            as.data.frame(predictModel(model.regHS.invG.eq, grid.eq = T, newdata = newData$x)),
                            "Method" = "regHS.invG", "L" = L, "Eq_Space" = 1, "distX" = distX)
rs.lap.L.eq = cbind(newData,
                     as.data.frame(predictModel(model.lap.eq,  grid.eq = T, newdata = newData$x)),
                     "Method" = "laplacian", "L" = L, "Eq_Space" = 1, "distX" = distX)
rs.norm.L.eq = cbind(newData,
                      as.data.frame(predictModel(model.norm.eq, grid.eq = T, newdata = newData$x)),
                      "Method" = "Gaussian", "L" = L, "Eq_Space" = 1, "distX" = distX)

model.list <- list(model.regHS.fix, model.regHS.invG, model.lap, model.norm,
                   model.regHS.fix.eq, model.regHS.invG.eq, model.lap.eq, model.norm.eq)
function.rhat <- function(model){return(max(summary(model)$summary[,"Rhat"]))}
function.Ndiv <- function(model){return(get_num_divergent(model))}
dt.diag = data.frame("Model" = TimeNames,
                     "Method" = rep(basename,2),
                     "N" = N,"L" = L, "Eq_Space" = c(rep(0,4),rep(1,4)), "distX" = distX,
                     "maxRhat" = unlist(lapply(model.list, function.rhat)),
                     "Num_Divergent" = unlist(lapply(model.list, function.Ndiv)),
                     "ArrayID" = arrayID,
                     "RunTime" = RunTime)


dt.pred =  bind_rows(rs.regHS.fix.L,  rs.regHS.invG.L,
                     rs.lap.L, rs.norm.L,
                     rs.regHS.fix.L.eq,  rs.regHS.invG.L.eq,
                     rs.lap.L.eq,  rs.norm.L.eq)

dt.pred$which_outcome = which_outcome
dt.pred$arrayID = arrayID
dt.pred$N = N


write.csv(dt.pred, paste("Outcome/Estimation/", distX, "_", which_outcome, "_N", N, "_Array", arrayID,".csv", sep = "")) #df
write.csv(dt.diag, paste("Outcome/RunTime/Time", distX, "_", which_outcome, "_N", N, "_Array",arrayID,".csv", sep = "")) #df
