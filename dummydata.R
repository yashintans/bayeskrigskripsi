# Data asli
dta.ori <- read.csv("CHJabar.csv",sep=";")

# Generate dari distribusi normal
djan <- abs(rnorm(214, mean(dta.ori$JAN), sd(dta.ori$JAN)))
dfeb <- abs(rnorm(214, mean(dta.ori$FEB), sd(dta.ori$FEB)))
dmar <- abs(rnorm(214, mean(dta.ori$MAR), sd(dta.ori$MAR)))
dapr <- abs(rnorm(214, mean(dta.ori$APR), sd(dta.ori$APR)))
dmei <- abs(rnorm(214, mean(dta.ori$MAY), sd(dta.ori$MAY)))
djun <- abs(rnorm(214, mean(dta.ori$JUN), sd(dta.ori$JUN)))
djul <- abs(rnorm(214, mean(dta.ori$JUL), sd(dta.ori$JUL)))
daug <- abs(rnorm(214, mean(dta.ori$AUG), sd(dta.ori$AUG)))
dsep <- abs(rnorm(214, mean(dta.ori$SEP), sd(dta.ori$SEP)))
dokt <- abs(rnorm(214, mean(dta.ori$OKT), sd(dta.ori$OKT)))
dnov <- abs(rnorm(214, mean(dta.ori$NOV), sd(dta.ori$NOV)))
ddes <- abs(rnorm(214, mean(dta.ori$DES), sd(dta.ori$DES)))

# Dummy data
pos.id <- 1:214
pos.name <- paste0("Pos ", pos.id)
dta0 <- data.frame(pos.id,pos.name,dta.ori$CURRENT.LATITUDE,dta.ori$CURRENT.LONGITUDE,
		dta.ori$CURRENT.ELEVATION.M, djan, dfeb, dmar, dapr, dmei, djun, djul, 
		daug, dsep, dokt, dnov, ddes)
dta0$ann <- rowSums(dta0[,6:17])
colnames(dta0)[1:18] <- colnames(dta.ori)
head(dta0)

# Export to csv
write.csv(dta0,"CHJabar.csv", row.names=FALSE)

