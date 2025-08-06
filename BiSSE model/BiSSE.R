## BiSSE likelihood inference on the phylogeny

#install.packages("hisse")
#install.packages("parallel")

setwd("")

library(hisse)
library(parallel)
set.seed(1234)

# Loading data
phyUltra <- ape::read.tree(file = "phylogeny.tre")
traits_data <- read.csv("data_traits.csv",nrows=17)
#traits_data <- read.csv("data_traits.csv")

# Preparing data
phyUltra_sub <- ape::drop.tip(phy = phyUltra, tip = "farlowella")
#phyUltra_sub <- ape::drop.tip(phy = phyUltra, tip = phyUltra$tip.label[-which(phyUltra$tip.label %in% traits_data$Species)])
traits_data <- traits_data[traits_data$Species %in% phyUltra_sub$tip.label,]

#phyUltra_sub <- ape::drop.tip(phy = phyUltra, tip = "farlowella")
traits_data <- as.matrix(traits_data[-1])
species <- ape::Ntip(phyUltra_sub)
sampling.f <- species/28

# Preparing model
trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse.red <- trans.rates.bisse
trans.rates.bisse.red[!is.na(trans.rates.bisse.red)] = 1
bisse.fit <- NA

# Likelihood inference
try(hisse.fit <- hisse(phy = phyUltra_sub, data = traits_data, f=rep(sampling.f, 2),
                       hidden.states=FALSE, turnover=c(1,2), eps=c(1,1),
                       trans.rate=trans.rates.bisse.red))
hisse.fit$solution[1]/(hisse.fit$solution[3]+1) # lambda 1
hisse.fit$solution[2]/(hisse.fit$solution[4]+1) # lambda 2
hisse.fit$solution[3] # turnover rate
hisse.fit$solution[5] # transition rate

# CI
supportHisse_fit <- hisse::SupportRegionHiSSE(hisse.fit)
quantile(supportHisse_fit$points.within.region[,2], 0.025)/(quantile(supportHisse_fit$points.within.region[,4], 0.025)+1) # CI low lambda 1
quantile(supportHisse_fit$points.within.region[,2], 0.975)/(quantile(supportHisse_fit$points.within.region[,4], 0.975)+1) # CI high lambda 1
quantile(supportHisse_fit$points.within.region[,3], 0.025)/(quantile(supportHisse_fit$points.within.region[,5], 0.025)+1) # CI low lambda 2
quantile(supportHisse_fit$points.within.region[,3], 0.975)/(quantile(supportHisse_fit$points.within.region[,5], 0.975)+1) # CI high lambda 2
quantile(supportHisse_fit$points.within.region[,4], 0.025) # CI low turnover rate
quantile(supportHisse_fit$points.within.region[,4], 0.975) # CI high turnover rate
quantile(supportHisse_fit$points.within.region[,6], 0.025) # CI low transition rate
quantile(supportHisse_fit$points.within.region[,6], 0.975) # CI high transition rate



##Plotting results


library(ggplot2)

# Data frame with point estimates and CIs
df <- data.frame(
  State = c("State 0", "State 1"),
  Lambda = c(0.032, 0.941),
  CI_low = c(0.016, 0.485),
  CI_high = c(0.141, 1.947)
)

# Plot
ggplot(df, aes(x = State, y = Lambda)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1) +
  theme_minimal() +
  ylab("Speciation Rate (Lambda)") +
  xlab("Trait State") +
  ggtitle("Speciation Rates by Trait State") +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

##Likelihood Ratio Test (LRT) to verify significance of data
# Full model (trait-dependent)
full_model <- hisse(phy = phyUltra_sub, data = traits_data,
                    turnover = c(1, 2), eps = c(1,1),
                    trans.rate = trans.rates.bisse.red,
                    hidden.states = FALSE, f = rep(sampling.f, 2))

# Null model (trait-independent)
null_model <- hisse(phy = phyUltra_sub, data = traits_data,
                    turnover = c(1, 1), eps = c(1,1),
                    trans.rate = trans.rates.bisse.red,
                    hidden.states = FALSE, f = rep(sampling.f, 2))

# LRT
LR = 2 * (full_model$loglik - null_model$loglik)
p_value = pchisq(LR, df=1, lower.tail=FALSE)

print(p_value)

