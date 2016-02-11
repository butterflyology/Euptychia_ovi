# Oviposition preference analaysis for Euptychia westwoodi
# Chris Hamm and Jim Fordyce - May 2015 at Tirimbina
# Three major code chunks: 
# 	1) oviposition experiments
#	2) feeding experiments
#	3) ordinated diet breadth
# 	4) phylogenetic distance
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
library("ordiBreadth")

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
plot(postburn, pch = 19, las = 1, type = "l")

effectiveSize(postburn)
autocorr.plot(postburn, lag.max = 50)
1 / sqrt(effectiveSize(postburn)) # We have a problem with effective sampel size and the MCse, which we want below 0.05
geweke.diag(postburn) # values greater that 2.0 (on abs number scale) indicate a possbile problem with convergence. 

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

tnuR <- as.mcmc(Runt)
effectiveSize(tnuR) # very nice
1 / sqrt(effectiveSize(tnuR)) # problem solved
geweke.diag(tnuR) # suggests the chains have converged

head(Runt)
plot(Runt, pch = 19, las = 1, type = "l")
plot(density(Runt), lwd = 2, las = 1, ylim = c(0, 20), xlab = "", main = "")

source("PPlot.R")
# pdf(file = "Eupt_plot.pdf", bg = "white")
PPlot(Eupt_run1[[1]], burn = 1e3, ymax = 100, pop = FALSE)
# dev.off()


#####
##### Feeding experiments
####

# ANOVA
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

lma1 <- lm(Mdiff[, 1] ~ as.factor(Mdiff[, 2]))
lma1
summary(lma1)
Anova(lma1)
TukeyHSD(aov(lma1))

# pdf(file = "Images/Sel-box.pdf", bg = "white")
boxplot(Mdiff[, 1] ~ Mdiff[, 2], las = 1, col = c("light grey", "dark grey", "dark grey"), ylab = "Mass difference", names = c("G1", "S1", "S2"), pch = 19, ylim = c(-10, 30))
abline(h = 0, lwd = 2, lty = 2)
text(1, 2, "A")
text(2, 28, "B")
text(3, 21, "B")
# dev.off()


# ANCOVA
with(Sel, table(Time, Indiv))
# with(Sel, interaction.plot(Time, Indiv, Mass))

S1 <- Sel[Sel$Trmt == "S1", ]
S2 <- Sel[Sel$Trmt == "S2", ]
G1 <- Sel[Sel$Trmt == "G1", ]

par(mfrow = c(1, 2))
with(Sel[Sel$Time == "1", ], plot(Trmt, Mass, type = "b", pch = 19, las = 1, ylab = "Mass (mg)", main = "Time 1", col = c("grey", "dark green", "dark green")))
with(Sel[Sel$Time == "2", ], plot(Trmt, Mass, type = "b", pch = 19, yaxt = "n", main = "Time 2", col = c("grey", "dark green", "dark green")))
dev.off()

T1 <- Sel[Sel$Time == "1", ]
T2 <- Sel[Sel$Time == "2", ]

plot(T1$Mass, T2$Mass, col = T2$Trmt, pch = 19, las = 1, ylab = expression(paste("Mass (gm), time"[2])), xlab = expression(paste("Mass (gm), time"[1])), cex = 1.5, xlim = c(0, 50), ylim = c(0, 50))
abline(a = 0, b = 1, lwd = 3, lty = 2) # the 1 : 1 line
legend("bottomright", legend = c(expression(paste(italic("S. eurynota"))), expression(paste(italic("S. arthritica"))), expression(paste(italic("L. rusifolia")))), col = c("red", "green", "black"), pch = 19, pt.cex = 1.5, bty = "n")

# ANCOVA - the final mass of the individual is a function of its starting mass and the treatment group it was in. 
a1 <- lm(T2$Mass ~ T1$Mass * T2$Trmt) # interaction model
a2 <- lm(T2$Mass ~ T1$Mass + T2$Trmt) # additive model
anova(a1, a2) # NO differnce between models (present p value)
summary(a2)
par(mfrow = c(2, 2)); plot(a2)
dev.off()
# line 1 = intercept for G1, line2 = slope for T1 (they lost mass because we have an additive model)
# S1 is 11 mg heavier than G1 at the smallest size (which is 0 at T1), S2 is 12 mg higher than G1 - and we have the confidence intervals around the estimates

# estimate average mass gain for S1 and S2 at T2
cbind(as.vector(coef(a2)[-2]), c("G1", "S1", "S2"))
confint(a2)[-2, ]

S1delta <- T2$Mass[T2$Trmt == "S1"] - T1$Mass[T1$Trmt == "S1"]
S2delta <- T2$Mass[T2$Trmt == "S2"] - T1$Mass[T1$Trmt == "S2"]
G1delta <- T2$Mass[T2$Trmt == "G1"] - T1$Mass[T1$Trmt == "G1"]

t1 <- t.test(S1delta, S2delta)
t2 <- t.test(S1delta, G1delta)
t3 <- t.test(S2delta, G1delta)

p.adjust(p = c(t1[3], t2[3], t3[3]), method = "holm", n = 3)

#Just for kicks, compare only the final mass (not mass gain)
(t.test(T2$Mass[T2$Trmt == "S1"], T2$Mass[T2$Trmt == "S2"]))


pred.frame <- expand.grid(Mass1 = seq(min(na.omit(T1$Mass)), max(na.omit(T1$Mass)), length = 50), Trmt = unique(T2$Trmt))
mm <- model.matrix(~Mass1 + Trmt, data = pred.frame)
preds <- mm %*% coef(a2)
V <- vcov(a2)
pred.var <- diag(mm %*% V %*% t(mm))

predictions <- data.frame(pred.frame, preds,pred.var)
predictions <- with(predictions, data.frame(predictions,
	hi95 = preds + 1.96 * sqrt(pred.var),
	lo95 = preds - 1.96 * sqrt(pred.var)))

dat <- as.data.frame(cbind(Sel[Sel$Time == "1", ], Sel[Sel$Time == "2", 5]))
names(dat) <- c("Trmt", "Indiv", "Time", "HW", "Mass1", "Mass2")


dat$Trmt <- factor(dat$Trmt, levels = c("S1", "S2", "G1"))

P0 <- ggplot(data = dat) + xlab(expression(paste("Mass"[1], " (mg)"))) + ylab(expression(paste("Mass"[2], " (mg)"))) + facet_grid(~Trmt)

P1 <- P0 + geom_ribbon(data = predictions, aes(x = Mass1, ymin = lo95, ymax = hi95),fill = "grey") + geom_point(data = dat, aes(x = Mass1, y = Mass2),size = 4) + geom_line(data = predictions, aes(x = Mass1, y = preds)) + geom_abline(intercept = 0, slope = 1, lty = 2, lwd = 1.1) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# pdf(file = "Mass_plot_1.pdf", bg = "white")
P1
# dev.off()


#####
##### Ordinated diet breadth
#####

head(diet) # The diet breadth matrix

ord1 <- ordi.breadth(diet, dist.method = "jaccard")
summary.hbreadth(ord1, by = "Herbivore", round = 3)


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