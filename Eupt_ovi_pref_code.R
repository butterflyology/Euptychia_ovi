# Oviposition preference analaysis for Euptychia westwoodi
# Chris Hamm and Jim Fordyce - May 2015 at Tirimbina

set.seed(1234)
library("bayespref")
library("coda")

# load("Data/Eupt_ovi_pref.RData")

setwd("~/Desktop/Projects/Euptychia_ovi_pref")

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

plot(Eupt_run1[[1]]$PopPref[1,], pch = 19, las = 1, xlab = "Iteration", )

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
	cat("\n", i, "of", N, "\n")}
	return(Name)
	}

Runt <- SpiderMonkey(N = 1e3)
autocorr.plot(Runt)
head(Runt)
plot(Runt, pch = 19, las = 1)
plot(density(Runt), lwd = 2, las = 1, ylim = c(0, 20), xlab = "", main = "")



# save.image("Data/Eupt_ovi_pref.RData")