# Oviposition preference analaysis for Euptychia westwoodi
# Chris Hamm and Jim Fordyce - May 2015 at Tirimbina
# Three major code chunks: 
# 	1) oviposition experiments
#	2) feeding experiments
# 	3) phylogenetic distance
#	
 
set.seed(15213748)

library("bayespref")
library("coda")
library("ggplot2")
library("picante")
library("phytools")
library("car")
library("dplyr")
library("multcomp")



setwd("~/Desktop/Projects/Euptychia_ovi_pref/")
# load("Data/Euptychia_Selaginella.RData")
# save.image("Data/Euptychia_Selaginella.RData")

(sessInf <- sessionInfo())


#####
##### Oviposition experiments
#####

Eupt_data <- read.csv("Data/Ovi_pref_data.csv", header = TRUE, row.names = 1)
head(Eupt_data)

Eupt_Sel <- Eupt_data[, 2:3]
head(Eupt_Sel)
Eupt_Sel <- Eupt_Sel[which(rowSums(Eupt_Sel) > 0), ]
Eupt_Sel
colSums(Eupt_Sel)

Eupt_run1 <- bayesPref(pData = Eupt_Sel, mcmcL = 1e4, pops = FALSE)
str(Eupt_run1)
names(Eupt_run1[[1]])

plot(Eupt_run1[[1]]$PopPref[1,], pch = 19, las = 1, xlab = "Iteration")

prefPlot(Eupt_run1[[1]], burn = 1e3, ymax = 60, pop = FALSE)
head(Eupt_run1[[1]]$PopPref[2, ])
length(Eupt_run1[[1]]$PopPref[2, ])

burnlow <- floor(0.1 *(length(Eupt_run1[[1]]$PopPref[2, ])))
postburn <- Eupt_run1[[1]]$PopPref[2, ][burnlow:length(Eupt_run1[[1]]$PopPref[2, ])]
plot(postburn, pch = 19, las = 1)

effectiveSize(postburn)
autocorr.plot(postburn, lag.max = 50)

# Running so the last state of a 1000 step run is saved for each of X runs
SpiderMonkey <- function(N){
	Name <- numeric(length = N)
		for(i in 1:N){
	temp <- bayesPref(pData = Eupt_Sel, mcmcL = 1e3, pops = FALSE)
	Name[i] <- temp[[1]]$PopPref[2, 1000]
	cat("\n", i, "of", N, "\n")
	}
	return(Name)
	}

Runt <- SpiderMonkey(N = 1e3)
autocorr.plot(Runt)
head(Runt)
plot(Runt, pch = 19, las = 1)
plot(density(Runt), lwd = 2, las = 1, ylim = c(0, 20), xlab = "", main = "")

source("PPlot.R")
# pdf(file = "Eupt_plot.pdf", bg = "white")
PPlot(Eupt_run1[[1]], burn = 1e3, ymax = 100, pop = FALSE)
# dev.off()


#####
##### Feeding experiments
####

Sel <- read.csv("Data/Sel1a.csv", header = TRUE)
Sel$Indiv <- as.factor(Sel$Indiv)
Sel$Time <- as.factor(Sel$Time)
#Sel <- na.omit(Sel)
str(Sel)
head(Sel)


# Treating response variable as the difference between start and end mass
M1 <- Sel %>% group_by(Indiv) %>% filter(Time == 1)
M2 <- Sel %>% group_by(Indiv) %>% filter(Time == 2)

Mdiff <- matrix(NA, nrow = 54, ncol = 2)
Mdiff[, 1] <- M2$Mass - M1$Mass
Mdiff[, 2] <- as.factor(M1$Trmt)
str(Mdiff)
summary(Mdiff)

a1 <- lm(Mdiff[, 1] ~ as.factor(Mdiff[, 2]))
a1
summary(a1)
Anova(a1)
TukeyHSD(aov(a1))

# pdf(file = "Images/Sel-box.pdf", bg = "white")
boxplot(Mdiff[, 1] ~ Mdiff[, 2], las = 1, col = c("light grey", "dark grey", "dark grey"), ylab = "Mass difference", names = c("G1", "S1", "S2"), pch = 19, ylim = c(-10, 30))
abline(h = 0, lwd = 2, lty = 2)
text(1, 2, "A")
text(2, 28, "B")
text(3, 21, "B")
# dev.off()

#####
##### Phylogenetic distance
#####

Plantae <- read.nexus("Data/Plantae_2.tre")
str(Plantae)

plot(Plantae, cex = 0.75)
edgelabels(Plantae$edge.length, cex = 0.5, frame = "none", adj = c(0, -0.5))

# create the data matrix with the species and traits
EA.data <- matrix(c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1), byrow = TRUE, nrow = 2, ncol = 12)
rownames(EA.data) <- c("Euptychia", "Adelpa")
colnames(EA.data) <- c("Physcomitrella", "Selaginella", "Poales", "Malpighiales", "Ericales", "Myrtales", "Rosales", "Asterales", "Lamiales", "Malvales", "Fagales", "Gentianales")
EA.data


# Do for Pierella too
EP.data <- matrix(c(1, 1, 1, 0, 0, 0, 1, 1), byrow = TRUE, nrow = 2)
rownames(EP.data) <- c("Euptychia", "Pieralla")
colnames(EP.data) <- c("Physcomitrella", "Selaginella", "Poales", "Zingiberales")
EP.data

pd(samp = EA.data, tree = Plantae, include.root = FALSE)
pd(samp = EP.data, tree = Plantae, include.root = FALSE)