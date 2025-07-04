#Install INLA
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

#Install INLAjoint
#library(devtools)
#devtools::install_github('DenisRustand/INLAjoint', build_vignettes = TRUE)

library(INLA)
library(INLAjoint)
library(ggplot2)
library(boot)
library(JMbayes2)

head(pbc2) # primary biliary cholangitis dataset
pbc2$death <- ifelse(pbc2$status=="dead",1,0)
pbc2$trans <- ifelse(pbc2$status=="transplanted",1,0) # competing event
# extracting time to event data
SurvData <- pbc2[c(which(diff(as.numeric(pbc2[,which(colnames(pbc2)=="id")]))==1),
                   length(pbc2[,which(colnames(pbc2)=="id")])),c("id", "years", "death","drug", "sex")]

head(SurvData)


## Model 1 - Proportional hazards model with parametric baseline hazards
M1.exp <- inla(inla.surv(time = years, event = death) ~ drug + sex,
               data = SurvData, family="exponentialsurv")
summary(M1.exp)

M1.wei <- inla(inla.surv(years, death) ~ drug + sex,
               data = SurvData, family="weibullsurv")
summary(M1.wei)


# same model with INLAjoint
M1.exp_2 <- joint(formSurv = inla.surv(years, death) ~ drug + sex,
                  dataSurv = SurvData, basRisk="exponentialsurv")
summary(M1.exp_2)

summary(M1.exp_2, hr=TRUE) # hazard ratios instead of betas

# Model 2 - Weibull baseline risk
M1.wei_2 <- joint(formSurv = inla.surv(years, death) ~ drug + sex,
                  dataSurv = SurvData, basRisk="weibullsurv")
summary(M1.wei_2) # shape = 1 => exponential

# for weibull variants, see: inla.doc("weibullsurv")
# variant=0 is the commonly used weibull survival model
# variant=1 is a more stable, slightly differently parametrized but equivalent version
# => Can be formulated as AFT!


# Model 3 - smooth semiparametric baseline risk (random walk order 2)
M1_3 <- joint(formSurv = inla.surv(years, death) ~ drug + sex,
              dataSurv = SurvData, basRisk="rw2", NbasRisk = 30)
summary(M1_3, hr=T)
plot(M1_3) # posterior marginals and baseline risk

NewData <- SurvData[c(1,3),]
NewData$years=0
NewData$death=0
NewData
P <- predict(M1_3, NewData, horizon=14, survival=TRUE, id="id")
str(P)
P1 <- P$PredS[P$PredS$id==1,]
P2 <- P$PredS[P$PredS$id==3,]
# plot hazard conditional on sex
plot(P1$years, P1$Haz_quant0.5, type="l", ylim=c(0, 0.4), main="Hazard")
lines(P1$years, P1$Haz_quant0.025, type="l", lty=2)
lines(P1$years, P1$Haz_quant0.975, type="l", lty=2)
lines(P2$years, P2$Haz_quant0.5, col=2)
lines(P2$years, P2$Haz_quant0.025, col=2, lty=2)
lines(P2$years, P2$Haz_quant0.975, col=2, lty=2)
legend("topleft", c("Female", "Male"), lty=c(1,1), col=c(1,2))
# plot survival conditional on sex
plot(P1$years, P1$Surv_quant0.5, type="l", ylim=c(0,1), main="Survival probabilities")
lines(P1$years, P1$Surv_quant0.025, type="l", lty=2)
lines(P1$years, P1$Surv_quant0.975, type="l", lty=2)
lines(P2$years, P2$Surv_quant0.5, col=2)
lines(P2$years, P2$Surv_quant0.025, col=2, lty=2)
lines(P2$years, P2$Surv_quant0.975, col=2, lty=2)
legend("topright", c("Female", "Male"), lty=c(1,1), col=c(1,2))
summary(SurvData) # more data om females => less uncertainty

## Model 4 - Competing risks
SurvData2 <- pbc2[c(which(diff(as.numeric(pbc2[,which(colnames(pbc2)=="id")]))==1),
                   length(pbc2[,which(colnames(pbc2)=="id")])),c("id", "years", "death", "trans","drug", "sex")]
head(SurvData2)

M2 <- joint(formSurv = list(inla.surv(years, death) ~ drug + sex,
                            inla.surv(years, trans) ~ drug + sex),
            basRisk = c("rw2", "rw2"), dataSurv = SurvData2)
summary(M2)
summary(M2, hr=T)
plot(M2)

#Prediction
NewData2 <- SurvData2[c(2,3,5,14),]
NewData2$years <- 0
NewData2$death <- NewData2$trans <- 0
NewData2$id <- 1:4
NewData2

P2 <- predict(M2, NewData2, horizon=14, CIF=TRUE, id="id")
# CIF is preferred because it accounts for other competing events
ggplot(P2$PredS, aes(x=years, y=CIF_quant0.5, group=id)) +
  geom_line(aes(color=id)) +
  geom_line(aes(x=years, y=CIF_quant0.025, color=id), linetype="dashed")+
  geom_line(aes(x=years, y=CIF_quant0.975, color=id), linetype="dashed")+
  facet_wrap(~Outcome+id, ncol=4)

# Model 5 - Longitudinal model

LongData <- pbc2[, c("id", "year", "serBilir","drug")]
print(head(LongData, 10), row.names=F)
# mixed effects regression model
M1 <- joint(formLong = serBilir ~ year + drug  +  (1 + year|id),
            dataLong = LongData, id = "id", timeVar = "year",
            family = "lognormal")
summary(M1)
summary(M1, sdcor=TRUE)
inla.priors.used(M1)
plot(M1)
plot(M1, priors=TRUE, sdcor=TRUE)$Covariances


Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2


# empirical Bayes: same frequentist properties / cheaper computational cost

# Model 6 - Joint model with functions of time and current value + current slope association
f1 <- function(x) x^2
f2 <- function(x) x^3

M4 <- joint(formSurv =  inla.surv(time = years, event = death)  ~ drug,
            formLong = serBilir ~ (1 + year + f1(year) + f2(year))*drug +
              (f1(year) + f2(year) |id), family = "lognormal",
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", assoc = "CV_CS",
            basRisk = "rw2", NbasRisk=25, control=list(int.strategy="eb", safemode=F))
summary(M4)
plot(M4, sdcor=T)






# Model 7 - Joint model for longitudinal and competing risk
M6 <- joint(formSurv = list(inla.surv(time = years, event = death)  ~ sex + drug,
                            inla.surv(time = years, event = trans) ~ edema * sex),
            formLong = serBilir ~ year * (drug + sex) + (1+year|id),
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", family = "lognormal",
            basRisk = c("rw1", "rw1"), assoc = c("SRE", "SRE_ind"), control=list(int.strategy="eb", safemode=F))
summary(M6)

# Model 8 - Bivariate joint model with competing risks
M8 <- joint(formLong = list(serBilir ~ year * drug + sex + (1|id),
                            platelets ~ year + f1(year) + drug + sex + (1|id)),
            formSurv = list(inla.surv(time = years, event = death) ~ drug,
                            inla.surv(time = years, event = trans) ~ drug),
            dataLong = Longi, dataSurv=Surv, id = "id", corLong=TRUE, timeVar = "year",
            family = c("lognormal", "poisson"), basRisk = c("rw1", "rw1"),
            assoc = list(c("CV", "CV"), c("SRE", "")),
            control=list(int.strategy="eb"))
summary(M8)




# Model 9 - Multivariate joint model
#Takes around 24 minutes
Nsplines <- ns(Longi$year, knots=c(1,4))
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
f3 <- function(x) predict(Nsplines, x)[,3]

M16 <-joint(formSurv = list(inla.surv(time = years, event = death) ~ drug,
                            inla.surv(time = years, event = trans) ~ drug),
            formLong = list(serBilir ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            SGOT ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            albumin ~ (1 + year) * drug + (1 + year | id),
                            platelets ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                              (1 + f1(year) + f2(year) + f3(year) | id),
                            spiders ~ (1 + year) * drug + (1 + year | id)),
            dataLong = Longi, dataSurv = Surv, id = "id", timeVar = "year", corLong = F,
            family = c("lognormal", "lognormal", "gaussian", "poisson", "binomial"),
            basRisk = c("rw2", "rw1"), NbasRisk = 15, assoc = list(c("CV_CS", "CV"),
                                                                   c("CV", ""), c("CV", "CV"), c("CV", "CV"), c("CV", "")),
            control=list(int.strategy="eb"))
summary(M16)

plot(M16)


#Extra examples

## Multi-state model (e.g., illness-death)
data(SurvMS) # package INLAjoint
E12 <- inla.surv(time = SurvMS[[1]]$Tstop, event = SurvMS[[1]]$status) # transition 1->2
E13 <- inla.surv(time = SurvMS[[2]]$Tstop, event = SurvMS[[2]]$status) # transition 1->3
E23 <- inla.surv(time = SurvMS[[3]]$Tstop, truncation=SurvMS[[3]]$Tstart,
                 event =SurvMS[[3]]$status) # transition 2->3

M3 <- joint(formSurv=list(E12 ~ X, E13 ~ X, E23 ~ X),
            basRisk = c("rw2", "rw1", "exponentialsurv"),
            dataSurv = SurvMS)
summary(M3)

## Mixture cure model
# library(smcure)
# ?bmt
#data(bmt) # package smcure
surv.obj.cure <- inla.surv(time = bmt$Time, event = bmt$Status,
                           cure = cbind("Int"=1, "TRT"=bmt$TRT))
M4 <- joint(formSurv = surv.obj.cure ~ TRT, basRisk = "weibullsurv",
            dataSurv = bmt, control = list(variant=0))
summary(M4)
smp.cure <- inla.hyperpar.sample(500, M4)[, 2:3]
quantile(inv.logit(smp.cure[,1]), c(0.025, 0.5, 0.975)) # allogeneic
quantile(inv.logit(rowSums(smp.cure)), c(0.025, 0.5, 0.975)) # autologous
# The proportion of the population that is ``cured'' and will not experience the event
# is 28% [16%, 43%] for the treatment group allogeneic and
# 20% [10%, 35%] for the treatment group autologous.

# For the remaining patients who are not considered cured, we can estimate the hazard ratio of treatment:
summary(M4, hr=T)
# In the non-cured fraction of patients, those who received autologous treatment
# have a hazard of experiencing the event 2.08 [1.23, 3.33] times higher than
# those who received allogeneic treatment, indicating a significantly greater risk
# associated with the autologous treatment option.





## Shared frailty model for recurrent events
# data(readmission) # package frailtypack
# ?readmission
M5 <- joint(formSurv=inla.surv(t.stop, event) ~ sex + (1|id), id="id", dataSurv = readmission)
summary(M5)

summary(M5, hr=T)
# This suggests that females are associated with a 30\% [12\%, 46\%] reduced susceptibility to rehospitalization
# when accounting for unobserved or latent individual-specific factors captured by the frailty term.






## Shared frailty model for recurrent and terminal event
terminalData <- readmission[readmission$event==0,]
M6 <- joint(formSurv=list(inla.surv(t.stop, event) ~ sex + (1|id), # risk of rehospitalization
                           inla.surv(t.stop, death) ~ sex), id="id", # risk of death
             basRisk=c("rw2", "exponentialsurv"), assocSurv=TRUE,
             dataSurv = list(readmission,terminalData))
summary(M6)




# Joint model for longitudinal and multi-state
data(SurvMS) # package INLAjoint
data(LongMS) # package INLAjoint
E12 <- inla.surv(time = SurvMS[[1]]$Tstop, event = SurvMS[[1]]$status) # transition 1->2
E13 <- inla.surv(time = SurvMS[[2]]$Tstop, event = SurvMS[[2]]$status) # transition 1->3
E23 <- inla.surv(time = SurvMS[[3]]$Tstop, truncation=SurvMS[[3]]$Tstart,
                 event =SurvMS[[3]]$status) # transition 2->3
M7 <- joint(formSurv=list(E12 ~ X, E13 ~ X, E23 ~ X),
            formLong=list(y ~ time + X + (1+time|id)),
            basRisk = c("rw2", "rw1", "exponentialsurv"), timeVar = "time",
            assoc = list(c("CV", "CV", "CV")), id="id",
            dataSurv = SurvMS, dataLong = LongMS,
            cutpoints=seq(0, max(SurvMS[[3]]$Tstop), len=15))
summary(M7)









