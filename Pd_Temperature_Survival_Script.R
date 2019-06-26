#Load packages

library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(lme4)


##################### PERCENT DECLINE

# Load data
setwd("/home/dwalsh/Documents/Consulting/Lewis Campbell")


########################### RAW COUNT DATA

raw_data <-  read.csv("PD_Temp_RAW.csv")

raw_data$rounded_CFU <- ceiling(raw_data$corrected_CFU)



?ceiling

head(raw_data)

summary(raw_data)

sum_raw_data <- ddply(raw_data, c("Day", "Temperature", "Media", "Treatment"), summarise,
                      N = length(rounded_CFU),
                      mean = mean(rounded_CFU),
                      sd = sd(rounded_CFU),
                      se = sd / sqrt(N)
)



ggplot(sum_raw_data) + 
  geom_errorbar(aes(x= factor(Day), ymin =(mean - se), ymax = (mean + se), col = Treatment), width = 0.9) + 
  geom_point(aes(x= factor(Day), y= mean, fill = Media, shape = Media, col = Treatment), size = 4) + 
  scale_shape_manual(values = c(21, 22,23)) + 
  scale_colour_manual(values = c("red", "black"), name = "Treatment") + 
  facet_wrap(~Temperature)

# with lines and baselines  
ggplot(sum_raw_data) + 
  geom_errorbar(aes(x= Day, ymin =(mean - se), ymax = (mean + se)), width = 0.3) +
  geom_line(aes(x= Day, y= mean, col= Media), size = 1) + 
  geom_point(aes(x= Day, y= mean, fill = Media, shape = Media), size = 4) +
  scale_colour_manual(values = c("steelblue", "firebrick", "springgreen4"), name = "Media") +
  scale_fill_manual(values = c("steelblue", "firebrick", "springgreen4"), name = "Media") +
  scale_shape_manual(values = c(21, 22,23)) + 
  scale_x_continuous(limits=c(0, 150), breaks=c(1, 30, 60, 90, 120, 150)) +
  theme_bw(base_size = 16, 
           base_family = 'Helvetica') +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background =element_rect(fill="grey90")) +
  labs(x="Days post inoculation", y="CFU") +
  geom_hline(yintercept = 19, linetype = "dashed", col = "firebrick", size = 1) +
  geom_hline(yintercept = 40, linetype = "dashed", col = "steelblue", size = 1) +
  geom_hline(yintercept = 65, linetype = "dashed", col = "springgreen4", size = 1) + 
  facet_wrap(~Temperature)

# no lines
ggplot(sum_raw_data) + 
  geom_point(aes(x= Day, y= mean, fill = Media, shape = Media), size = 4) +
  scale_colour_manual(values = c("steelblue", "firebrick", "springgreen4"), name = "Media") +
  scale_fill_manual(values = c("steelblue", "firebrick", "springgreen4"), name = "Media") +
  scale_shape_manual(values = c(21, 22,23)) + 
  scale_x_continuous(limits=c(0, 150), breaks=c(1, 30, 60, 90, 120, 150)) +
  theme_bw(base_size = 16, 
           base_family = 'Helvetica') +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x="Days post inoculation", y="CFU") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red", size = 1) +
  facet_wrap(~Temperature)

###Estimation

## @knitr libraries
##Load necessary libraries
library(nimble)
library(Matrix)
library(coda)
library(knitr)
library(ggplot2)
library(xtable)

##@knitr rawdata
setwd("/home/dwalsh/Documents/Consulting/Lewis Campbell")



########################### RAW COUNT DATA

raw_data <-  read.csv("PD_Temp_RAW.csv")

raw_data$rounded_CFU <- round(raw_data$corrected_CFU)

raw_data$rounded_CFU2 <- if_else(raw_data$corrected_CFU > 0, ceiling(raw_data$corrected_CFU), raw_data$corrected_CFU)

raw_data$rounded_CFU <- raw_data$rounded_CFU2
?revalue

raw_data$Media <- revalue(raw_data$Media, c(SabDex = "SD"))

#Create dataset for use by Nimble
datain <- raw_data[raw_data$Day!=0,]

##total number of time points
nT<-length(unique(datain$Day))

##create num vector
num<-c(1,rep(2,nT-2),1)

##create adj vector
temp<-as.matrix(bandSparse(n=nT,k=c(1),symmetric = T))
temp2<-matrix(0,nT,nT)
for(i in 1:nrow(temp2)){
  temp2[i,] <-  temp[i,]*c(1:nT)
}
adj<-t(temp2)[which(t(temp2)!=0)]

nN <- length(adj)

##Create indicator matrix of media
ind <- matrix(0,nrow(datain),length(unique(datain$Media)))

for(i in 1:nrow(datain)){
  ind[i,as.numeric(datain$Media[i])] <- 1
}

##Create indicator matrix of media*temperature

ind2 <- matrix(0,nrow(datain),length(unique(datain$Media))*length(unique(datain$Temperature)))
temp<-rep(0,nrow(datain))
for(i in 1:nrow(datain)){
  temp[i] <- paste(datain$Temperature[i],datain$Media[i])
}
temp <- as.factor(temp)
for(i in 1:nrow(ind2)){
  ind2[i,as.numeric(temp[i])] <- 1
}

interactionind <- rep(0,nrow(datain))
for(i in 1:length(interactionind)){
  interactionind[i] <- which(ind2[i,]==1)
}



##Create needed datasets
analysisData <-list(adj = adj, weightsin = rep(1,nN), num = num, Y = datain$rounded_CFU)

##Create needed constants
analysisConst <-list(records = nrow(datain), offset = datain$baseline_CFU, 
                    temp = as.numeric(as.factor(datain$Temperature)), media = as.numeric(datain$Media), 
                    time = as.numeric(as.factor(datain$Day)),
                    levels = ncol(ind2),
                    nT = nT, nN = nN,
                    ind = ind2, inter = interactionind, ind2 = diag(rep(1,ncol(ind2)))
                    )

##Create initial values
##Function for creating initial values for each chain
initsgen <-function(){
      temp<-list(Binter = runif(ncol(ind2),-0.1,1), 
                 day = runif(9,0,0), day2 = runif(9,0,0), day3 = runif(9,0,0),
                 day4 = runif(9,0,0), day5 = runif(9,0,0), day6 = runif(9,0,0),
                 day7 = runif(9,0,0), day8 = runif(9,0,0), day9 = runif(9,0,0),
                 tau = runif(9,0.5,1))
    
  return(temp)
}

##@knitr mcmcspecs
##MCMC Parameters
nchains <- 3
reps <- 100000
nburnin <- 0.1*reps

##@knitr mcmcrun
analysisInits <- initsgen()


modelcode <- nimbleCode({
  
  ##Priors
    for(j in 1:(levels-1)){
       Binter[j] ~ dnorm(0, sd=sqrt(100))
    }
    Binter[levels] ~ dnorm(0, sd=sqrt(10))
  
  day[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                            num[1:nT], tau[1], zero_mean=1)
  day2[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                          num[1:nT], tau[2], zero_mean=1)
  day3[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                          num[1:nT], tau[3], zero_mean=1)
  day4[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[4], zero_mean=1)
  day5[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[5], zero_mean=1)
  day6[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[6], zero_mean=1)
  day7[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[7], zero_mean=1)
  day8[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[8], zero_mean=1)
  day9[1:nT] ~ dcar_normal(adj[1:nN], weights = weightsin[1:nN], 
                           num[1:nT], tau[9], zero_mean=1)
 
  for(i in 1:levels){
    tau[i] ~ dgamma(1, 1)
  }


  ##Likelihood
  for(i in 1:records){
    Y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- log(offset[i]) + Binter[inter[i]] +
                          ind[i,1]*day[time[i]] + ind[i,2]*day2[time[i]] + ind[i,3]*day3[time[i]]+
                          ind[i,4]*day4[time[i]] + ind[i,5]*day5[time[i]] + ind[i,6]*day6[time[i]] +
                          ind[i,7]*day7[time[i]] + ind[i,8]*day8[time[i]] + ind[i,9]*day9[time[i]] 
  }
  
  ##Derived Parameters
  ###Predictions
  for(k in 1:levels){
    for(l in 1:nT){
      rel.lambda[k,l] <- Binter[k] +  
        ind2[k,1]*day[l] + ind2[k,2]*day2[l] + ind2[k,3]*day3[l]+
        ind2[k,4]*day4[l] + ind2[k,5]*day5[l] + ind2[k,6]*day6[l] +
        ind2[k,7]*day7[l] + ind2[k,8]*day8[l] + ind2[k,9]*day9[l] 
    }
  }
  
  ###Estimated Differences
  ###Temp = 24c
  rel24.diff[1,1:nT] <- rel.lambda[2, 1:nT] - rel.lambda[1, 1:nT]
  rel24.diff[2,1:nT] <- rel.lambda[1, 1:nT] - rel.lambda[3, 1:nT]
  rel24.diff[3,1:nT] <- rel.lambda[2, 1:nT] - rel.lambda[3, 1:nT]
  
  ###Temp = 30c
  rel30.diff[1,1:nT] <- rel.lambda[5, 1:nT] - rel.lambda[4, 1:nT]
  rel30.diff[2,1:nT] <- rel.lambda[4, 1:nT] - rel.lambda[6, 1:nT]
  rel30.diff[3,1:nT] <- rel.lambda[5, 1:nT] - rel.lambda[6, 1:nT]
  
  ###Temp = 37c
  rel37.diff[1,1:nT] <- rel.lambda[8, 1:nT] - rel.lambda[7, 1:nT]
  rel37.diff[2,1:nT] <- rel.lambda[7, 1:nT] - rel.lambda[9, 1:nT]
  rel37.diff[3,1:nT] <- rel.lambda[8, 1:nT] - rel.lambda[9, 1:nT]
  
})



analysis1<-nimbleModel(code= modelcode, name="additive", constants = analysisConst,
                     data = analysisData, inits = analysisInits)

mcmcout1<-nimbleMCMC(model=analysis1, nchains=nchains, nburnin = nburnin,
                       niter=reps, summary=TRUE, WAIC=TRUE,
                       monitors = c("day","day2","day3","day4","day5","day6","day7",
                                    "day8","day9","tau","Binter","rel.lambda","rel24.diff",
                                    "rel30.diff", "rel37.diff"))



out<-mcmc.list(lapply(mcmcout1$samples, mcmc))
gr.statistic <- gelman.diag(out[,1:9], confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

saveRDS(out,"results.RDS")
saveRDS(mcmcout1$summary,"resultsummary.RDS")

##@knitr import
out<-readRDS("results.RDS")
out.summary<-readRDS("resultsummary.RDS")

## @knitr trace
###Trace plots of Interaction effects

for(i in 1:ncol(ind2)){
  traceplot(out[,paste("Binter[",i,"]",sep="")], main = paste(levels(temp)[i],sep=""))
}


for(i in 1:ncol(ind2)){
  traceplot(out[,paste("tau[",i,"]",sep="")], main = paste(levels(temp)[i],sep=""))
}

## @knitr table1
trt.effects <- out.summary$all.chains[grep("Binter",rownames(out.summary$all.chains)),]
rownames(trt.effects) <- levels(temp)

trt.effects <- trt.effects[order(trt.effects[,"Mean"],decreasing = TRUE),]

table.effects<- xtable(as.data.frame(trt.effects))
align(table.effects) <-c("l","c","c","c","c","c")
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c("Treatment & Mean & Median & SD & 95\\% LCL & 95\\% UCL \\\\\n")

print(table.effects,  table.placement="H",
      latex.environments = "center",
      add.to.row = addtorow,
      include.colnames = FALSE,
      include.rownames = TRUE,floating = F)

effect.sizes <- as.data.frame(trt.effects)

effect.sizes$Treatment <- row.names(effect.sizes)

summary(effect.sizes)

colnames(effect.sizes)[4]<- "CI_low"

colnames(effect.sizes)[5]<- "CI_high"
  

ggplot(effect.sizes) +
  geom_point(aes(x=Treatment, y = Mean)) +
  geom_errorbar(aes(x= Treatment, ymin =(CI_low), ymax = (CI_high)), width = 0.3)

## @knitr plots
###Plots of Time effects
plot.trtday <- function(name,name2,name3,dataset,timeid,title){
  mean.day <- dataset[grep(name,rownames(dataset)),"Mean"] 
  q2.5 <- dataset[grep(name,rownames(dataset)),"95%CI_low"] 
  q97.5 <-dataset[grep(name,rownames(dataset)),"95%CI_upp"] 
  
  mean.day2 <- dataset[grep(name2,rownames(dataset)),"Mean"] 
  q2.52 <- dataset[grep(name2,rownames(dataset)),"95%CI_low"] 
  q97.52 <-dataset[grep(name2,rownames(dataset)),"95%CI_upp"] 
  
  mean.day3 <- dataset[grep(name3,rownames(dataset)),"Mean"] 
  q2.53 <- dataset[grep(name3,rownames(dataset)),"95%CI_low"] 
  q97.53 <-dataset[grep(name3,rownames(dataset)),"95%CI_upp"] 
  
  group <- c(rep("BHI",length(mean.day)),rep("BHI+B",length(mean.day2)),rep("SD",length(mean.day3)))
  
  temp.data <- data.frame(Mean = c(mean.day,mean.day2,mean.day3),
                          x = rep(timeid,3),
                          LCL = c(q2.5, q2.52, q2.53),
                          UCL = c(q97.5, q97.52, q97.53),
                          Trt = group)
  
  p1 <-ggplot(data =temp.data,aes(x = x, y = Mean, group=Trt,color=Trt))+
    geom_line(aes(linetype = Trt))+
    geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Trt),alpha=.1,show.legend=NA,linetype=0)+
    theme_bw()+  scale_x_continuous(limits=c(0, 150), breaks=c(1, 30, 60, 90, 120, 150)) +
    ggtitle(paste("Time Effect at",title))+xlab("Time (Days)")+ylab("Log count")
  return(p1)
}

plot.diffday <- function(name,name2,name3,dataset,timeid,title){
  mean.day <- dataset[grep(name,rownames(dataset)),"Mean"] 
  q2.5 <- dataset[grep(name,rownames(dataset)),"95%CI_low"] 
  q97.5 <-dataset[grep(name,rownames(dataset)),"95%CI_upp"] 
  
  mean.day2 <- dataset[grep(name2,rownames(dataset)),"Mean"] 
  q2.52 <- dataset[grep(name2,rownames(dataset)),"95%CI_low"] 
  q97.52 <-dataset[grep(name2,rownames(dataset)),"95%CI_upp"] 
  
  mean.day3 <- dataset[grep(name3,rownames(dataset)),"Mean"] 
  q2.53 <- dataset[grep(name3,rownames(dataset)),"95%CI_low"] 
  q97.53 <-dataset[grep(name3,rownames(dataset)),"95%CI_upp"] 
  
  group <- c(rep("BHI+B - BHI",length(mean.day)),rep("BHI - SD",length(mean.day2)),rep("BHI+B - SD",length(mean.day3)))
  
  temp.data <- data.frame(Mean = c(mean.day,mean.day2,mean.day3),
                          x = rep(timeid,3),
                          LCL = c(q2.5, q2.52, q2.53),
                          UCL = c(q97.5, q97.52, q97.53),
                          Comparison = group)
  new_data <- return(temp.data)
  
  p1 <-ggplot(data =temp.data,aes(x = x, y = Mean, group=Comparison,color=Comparison))+
    geom_line(aes(linetype = Comparison))+
    geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Comparison),alpha=.1,show.legend=NA,linetype=0)+
    theme_bw(base_size = 16, 
             base_family = 'Helvetica') +  
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +  scale_x_continuous(limits=c(0, 150), breaks=c(1, 30, 60, 90, 120, 150)) +
    xlab("Days post incoculation")+ylab(expression(Delta*"Log CFU count")) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "black", size = 0.5) +
    facet_grid(~Comparison)
  return(p1)
}

## @knitr plot1
plot.trtday("rel.lambda\\[1","rel.lambda\\[2","rel.lambda\\[3",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 24C")

## @knitr plot2
plot.trtday("rel.lambda\\[4","rel.lambda\\[5","rel.lambda\\[6",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 30C")

## @knitr plot3
plot.trtday("rel.lambda\\[7","rel.lambda\\[8","rel.lambda\\[9",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 37C")


## @knitr plot4
plot1 <- plot.diffday("rel24.diff\\[1","rel24.diff\\[2","rel24.diff\\[3",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 24C")

## @knitr plot5
plot2 <- plot.diffday("rel30.diff\\[1","rel30.diff\\[2","rel30.diff\\[3",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 30C")

## @knitr plot6
plot3 <- plot.diffday("rel37.diff\\[1","rel37.diff\\[2","rel37.diff\\[3",out.summary$all.chains,as.numeric(unique(datain$Day)),"Temp - 37C")

library(grid)
library(gridExtra)

grid.arrange(plot1, plot2, plot3)


### Heatmap

df_heatmap_bat <- read.csv("Survivaldf.csv")

df_heatmap_bat$Significant <- as.factor(df_heatmap_bat$Significant)
df_heatmap_bat$Temp <- as.factor(df_heatmap_bat$Temp)
df_heatmap_bat$Day <- as.factor(df_heatmap_bat$Day)

df_heatmap_bat$Significant <- revalue(df_heatmap_bat$Significant, c("1" = "Significant", "0" = "Insignificant"))

expression(24~degree~C)

df_heatmap_bat_24 <- filter(df_heatmap_bat, Temp == 24)

ggplot(df_heatmap_bat, aes(Day, Comparison)) +
  geom_tile(aes(fill = Significant), color = "black", size = 0.7) +
  scale_fill_manual(values = c("white", "steelblue2")) +
  ylab("Comparison") +
  xlab("Days post inoculation") +
  theme_bw(base_size = 16, 
           base_family = 'Helvetica') +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(size = 14, angle = 180),
        strip.background =element_rect(fill="grey90")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  labs(fill = "Significance") + 
  facet_grid(Temp ~., switch = "y", scales = "free_y", space = "free_y")

theme(legend.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      plot.title = element_text(size=16),
      axis.title=element_text(size=16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      strip.text.y = element_text(size = 14, angle = 180),
      strip.background =element_rect(fill="grey90"))

  ### Bat_fur

fur_data <- read.csv("Bat_fur.csv")

fur_data$Temperature <- as.factor(fur_data$Temperature)

fur_data$Day <- as.factor(fur_data$Day)

fur_data$Pd <- as.factor(fur_data$Pd)

library(plyr)

fur_data$Pd <- revalue(fur_data$Pd, c("1" = "Viable", "0" = "Non-viable"))

library(ggplot2)
ggplot(fur_data, aes(Day, Temperature)) +
  geom_tile(aes(fill = Pd), color = "black", size = 0.7) +
  scale_fill_manual(values = c("white", "steelblue2")) +
  ylab("Temperature (degrees centigrade)") +
  xlab("Days post inoculation") +
  theme_bw(base_size = 16, 
           base_family = 'Helvetica') +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(size = 14, angle = 180),
        strip.background =element_rect(fill="grey90")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  labs(fill = "Pd Status") 

temps <- c(7, 24, 30, 37)
days <- c(240, 180, 60, 5)

bar_data <- cbind(temps, days)


colnames(bar_data) <- c("Temperature", "Survivalday")

bar_data <- as.data.frame(bar_data)

glimpse(bar_data)

bar_data$Temperature <- as.factor(bar_data$Temperature)
bar_data$Survivalday <- as.factor(bar_data$Survivalday)

ggplot(bar_data) + 
  geom_point(aes(x=Temperature, y=Survivalday))
  
ggplot(bar_data) + 
  geom_point(aes(x=Survivalday, y=Temperature), shape = c(21,22,23,24), size = 3, fill = c("red", "blue", "green", "black"))
  
ggplot(bar_data) + 
    geom_bar(aes(x=Temperature, y= Survivalday), stat = "identity", 
             fill = c("coral2", "firebrick", "steelblue", "springgreen"), col = "black", width = 0.6) + 
  coord_flip() + 
  scale_y_continuous(limits=c(0, 250), breaks=c(5, 15, 30, 60, 90, 120, 150, 180, 210, 240)) +
  theme_bw(base_size = 20, 
           base_family = 'Helvetica') +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  labs(x= "Incubation temperature (\u00b0C)", y = "Days post inoculation")
