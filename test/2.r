library("gamlss")

cgid <- "!chr10:100011340"
y <- c(12,26,15,13,9,10,10,9,18,20,10,10,18,16,NaN,4,8,10,5,12,10,14,0,10,12,8,6,30,21,16,4,4,12,14,14,10,14)
total_reads <- c(14,37,31,21,9,16,12,9,28,30,12,16,18,21,NaN,24,16,12,11,22,12,20,18,24,29,25,14,36,25,27,6,10,21,30,21,17,30)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100011341"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(21,54,26,50,13,18,10,NaN,32,29,24,16,NaN,22,12,26,25,26,31,20,27,22,26,19,19,28,42,34,15,45,5,26,31,44,26,46,44)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100011387"
y <- c(0,0,0,0,0,0,0,0,0,0,0,2,0,0,NaN,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(14,40,32,20,9,17,13,10,28,30,13,17,18,22,NaN,25,16,13,12,26,12,21,20,24,32,25,14,37,18,28,6,10,22,32,22,18,30)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100011388"
y <- c(18,47,23,42,11,15,7,NaN,24,27,20,12,NaN,20,5,23,20,19,28,16,25,17,19,18,16,25,42,33,15,35,4,24,21,40,19,41,36)
total_reads <- c(18,54,26,49,13,19,8,NaN,32,29,22,16,NaN,22,12,26,24,25,32,18,30,22,23,19,17,28,42,35,15,44,5,24,29,45,26,46,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100026933"
y <- c(16,28,16,24,1,12,19,8,12,5,18,6,NaN,29,0,11,10,14,11,16,7,30,7,13,12,8,10,21,10,18,13,14,12,18,32,15,29)
total_reads <- c(30,55,46,28,28,15,22,8,27,13,23,11,NaN,33,6,21,28,22,17,35,21,40,40,19,19,30,21,55,11,34,22,24,23,35,38,30,39)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100026991"
y <- c(0,0,0,0,0,0,0,0,0,0,0,3,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0)
total_reads <- c(30,48,48,28,28,14,21,8,27,16,23,13,NaN,33,6,21,28,23,19,35,22,40,40,19,17,32,22,55,12,33,22,24,24,37,38,30,39)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027909"
y <- c(2,1,2,2,0,0,0,0,2,0,0,1,0,2,0,1,0,1,0,2,2,0,2,0,0,0,2,0,0,4,0,2,2,0,3,0,2)
total_reads <- c(77,66,62,59,49,25,31,26,57,53,81,63,18,42,40,56,48,32,56,55,59,65,64,46,48,57,49,74,26,85,43,30,53,46,57,59,87)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027910"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(34,49,42,36,36,15,36,12,44,44,42,35,13,26,NaN,48,34,21,35,50,35,44,49,40,29,27,33,38,13,71,24,34,29,45,33,45,44)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027919"
y <- c(0,0,0,0,0,0,0,0,3,0,0,2,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,0,0,1,0,2,0,0)
total_reads <- c(76,66,62,62,39,24,30,28,49,50,69,60,18,42,20,54,45,27,53,50,59,63,58,44,44,54,46,70,24,86,39,29,53,44,51,52,82)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027920"
y <- c(4,2,0,2,0,0,0,0,0,4,0,2,6,0,NaN,0,0,0,0,0,0,2,0,2,0,2,1,0,2,2,2,0,0,0,2,2,0)
total_reads <- c(71,100,86,74,72,30,73,24,89,90,86,72,26,54,NaN,97,68,43,74,100,70,89,100,82,58,56,65,75,26,142,48,68,59,90,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027923"
y <- c(0,1,0,2,0,0,0,0,0,0,2,0,0,3,0,0,0,0,2,2,0,0,0,2,0,1,0,0,0,0,0,0,2,0,2,2,2)
total_reads <- c(74,67,62,63,39,24,30,30,49,49,71,60,18,42,21,54,46,28,54,51,58,63,58,45,44,53,47,71,24,87,39,29,55,45,50,53,84)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027924"
y <- c(4,4,2,0,0,0,2,0,0,2,0,0,0,0,NaN,0,2,0,0,0,0,3,0,0,0,0,4,2,4,4,2,2,0,0,0,2,0)
total_reads <- c(71,100,86,73,72,30,73,24,89,90,86,68,26,54,NaN,97,68,42,75,99,69,89,100,81,58,56,65,75,26,142,48,68,60,89,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027926"
y <- c(2,0,0,0,0,0,2,0,2,0,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0)
total_reads <- c(76,67,62,63,39,24,30,30,51,50,71,60,18,42,21,54,46,28,54,52,59,63,58,45,44,54,47,72,24,87,39,29,55,45,51,53,84)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027927"
y <- c(0,0,2,0,0,0,0,0,0,0,2,0,0,0,NaN,2,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,2,0,2)
total_reads <- c(71,100,87,74,72,30,73,24,89,90,86,72,26,54,NaN,97,68,43,75,100,70,89,100,82,58,56,65,75,26,142,48,68,60,90,66,90,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027928"
y <- c(1,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0)
total_reads <- c(76,67,61,63,39,24,30,30,51,50,71,60,18,42,21,54,46,28,54,51,59,63,58,45,44,54,47,72,24,87,36,29,55,45,51,53,84)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027929"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,4,0,0,6,0,0,0,0,0,0,2,0,0,0,0,0)
total_reads <- c(71,100,87,74,72,30,73,24,89,90,86,72,26,54,NaN,97,68,43,75,100,70,89,100,82,58,56,65,75,26,142,48,68,60,90,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027931"
y <- c(3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,4,2,0,0,0,0,0,2,0,0,2,0,4,0,0)
total_reads <- c(76,67,62,63,39,24,30,30,51,50,71,60,18,42,21,54,46,28,54,52,59,62,59,45,44,54,47,72,24,87,39,29,55,45,51,53,84)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027932"
y <- c(2,4,0,2,0,0,2,0,0,0,0,2,0,2,NaN,2,0,0,0,0,0,2,0,0,0,0,6,1,0,0,0,2,0,2,0,0,0)
total_reads <- c(71,100,87,74,72,30,74,24,89,90,86,72,26,54,NaN,97,68,43,75,100,70,89,100,82,58,56,65,75,26,142,48,68,60,90,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027939"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0)
total_reads <- c(78,66,63,63,39,24,30,30,52,49,71,60,18,43,21,54,46,28,54,52,59,63,59,46,44,55,47,72,24,87,40,29,55,45,51,53,85)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027940"
y <- c(0,0,0,1,0,0,2,0,4,0,0,0,0,0,NaN,4,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(70,99,87,73,72,30,74,24,89,90,86,72,26,54,NaN,97,68,43,75,100,70,89,100,81,58,56,64,75,26,142,48,68,60,90,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027945"
y <- c(6,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0)
total_reads <- c(78,65,64,65,42,24,30,30,52,49,72,60,18,43,21,56,46,28,55,52,61,63,61,46,44,55,47,72,24,88,40,29,57,45,52,53,86)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027946"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0)
total_reads <- c(71,99,86,72,72,30,74,24,89,90,86,70,26,54,NaN,97,68,43,75,100,70,89,100,81,58,54,64,75,26,142,48,68,60,90,66,91,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027958"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(108,73,79,77,69,37,44,27,73,76,97,74,22,61,30,68,62,42,80,80,81,85,91,50,61,76,63,105,27,106,52,40,70,55,68,65,101)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027959"
y <- c(0,0,1,0,0,0,0,0,1,0,0,0,0,1,NaN,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,1,0)
total_reads <- c(36,42,42,31,33,15,37,12,45,44,41,33,13,26,NaN,46,31,18,35,50,34,41,49,40,29,28,31,36,13,70,24,34,29,40,32,45,44)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027962"
y <- c(0,0,0,0,0,0,0,NaN,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,14,18,14,27,13,13,NaN,22,28,25,15,6,19,14,12,18,14,27,40,21,24,30,6,17,23,19,38,6,20,13,11,19,12,18,17,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027963"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
total_reads <- c(42,30,13,29,23,16,22,NaN,26,25,24,19,8,25,10,17,13,15,16,20,14,29,30,14,16,31,26,40,6,30,8,17,24,20,16,33,40)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027982"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,14,18,14,27,13,14,NaN,22,28,25,16,6,19,14,12,18,14,27,40,22,24,30,6,17,23,19,38,6,20,13,11,19,13,18,17,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027983"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(42,30,12,29,23,16,22,NaN,26,25,25,19,8,25,10,17,13,15,16,20,14,29,30,14,16,31,27,40,6,30,8,17,24,21,16,33,40)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027989"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,14,18,14,27,13,14,NaN,22,28,25,16,6,19,14,11,18,14,27,40,22,24,30,6,17,23,19,38,6,20,13,11,19,13,18,17,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027990"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(42,30,12,27,23,16,22,NaN,26,25,25,19,8,25,10,17,14,15,16,20,14,29,30,14,16,31,27,40,6,30,8,17,24,21,16,34,41)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027993"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,14,18,14,27,13,14,NaN,22,28,25,16,6,19,14,11,18,14,27,40,23,24,30,6,17,23,19,38,6,20,13,11,17,13,18,17,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027994"
y <- c(0,1,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0)
total_reads <- c(42,30,12,27,23,16,22,NaN,26,25,25,19,8,25,10,17,14,14,16,20,13,29,30,14,16,31,27,40,6,30,8,17,24,20,16,34,41)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027995"
y <- c(0,0,0,0,0,1,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,14,18,14,27,13,14,NaN,22,29,25,16,6,19,14,11,18,14,27,40,23,24,30,6,17,23,18,38,6,20,13,11,17,12,18,17,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100027996"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(42,30,12,27,23,16,22,NaN,26,25,25,19,8,25,10,17,14,15,16,20,13,29,30,14,16,31,27,40,6,30,8,17,24,21,16,34,41)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028009"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
total_reads <- c(35,20,18,15,27,14,14,NaN,22,27,23,17,6,20,14,12,18,14,27,39,26,23,31,6,19,23,24,35,6,22,11,11,19,14,18,18,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028010"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(41,28,14,26,22,14,19,NaN,26,24,25,18,8,25,6,18,14,14,19,22,14,30,31,14,14,31,26,34,6,29,8,18,25,20,15,34,40)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028018"
y <- c(0,0,0,0,0,0,0,NaN,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(35,20,18,17,27,14,14,NaN,22,27,23,17,6,20,14,14,18,14,27,39,25,23,31,6,19,23,23,35,6,22,11,11,19,14,18,18,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028019"
y <- c(0,0,0,0,0,0,0,NaN,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0)
total_reads <- c(41,29,14,27,22,14,19,NaN,26,25,23,18,8,25,6,18,15,14,19,22,14,29,31,14,14,31,26,34,6,29,8,18,25,19,15,34,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028020"
y <- c(0,0,0,1,0,0,0,NaN,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(35,20,18,17,27,14,14,NaN,22,27,23,17,6,20,14,14,18,14,27,39,25,23,31,6,19,22,23,35,6,22,11,11,19,14,18,18,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028021"
y <- c(0,0,0,0,0,0,0,NaN,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(42,29,14,27,22,14,20,NaN,26,25,24,18,8,25,6,18,14,14,19,22,14,30,31,14,15,31,25,35,6,29,8,18,25,20,16,34,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028032"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,20,18,17,27,15,14,20,22,29,25,16,6,20,15,14,18,14,27,40,25,26,31,6,19,24,26,38,6,22,13,11,20,14,18,18,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028033"
y <- c(1,0,0,0,0,0,0,NaN,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(43,30,14,27,23,14,20,NaN,26,25,25,19,8,25,10,18,14,15,19,22,14,30,32,16,16,31,25,40,6,32,8,19,26,21,16,33,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028041"
y <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(39,20,18,17,27,15,14,20,22,29,25,16,6,20,15,14,18,14,28,40,24,26,31,6,19,24,26,38,6,22,13,11,20,14,18,18,18)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028042"
y <- c(0,0,0,1,0,0,0,NaN,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
total_reads <- c(42,30,13,27,23,14,20,NaN,26,24,25,19,7,25,10,18,14,15,19,22,14,27,32,16,15,31,25,39,6,31,8,19,26,21,16,33,41)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028046"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(54,40,48,28,42,26,27,34,41,41,51,39,19,39,26,38,46,23,52,56,42,47,47,23,37,34,47,62,11,58,16,16,46,28,27,42,43)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028050"
y <- c(0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(30,45,51,24,29,21,25,35,35,25,50,45,19,38,21,42,52,18,48,34,33,46,29,33,34,17,47,47,10,65,6,10,55,27,16,44,51)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028051"
y <- c(0,0,0,0,0,0,0,0,0,0,0,2,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(20,39,35,46,33,40,30,32,29,26,40,48,22,30,NaN,38,14,32,26,23,18,62,19,49,22,13,27,48,12,49,19,17,19,38,30,14,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028063"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(30,46,51,24,29,22,25,35,35,25,50,45,19,38,21,42,52,18,48,34,33,46,29,33,34,17,44,47,10,66,6,10,55,27,16,44,51)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028064"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(20,39,35,46,33,40,30,32,29,26,40,48,22,30,NaN,38,14,32,26,23,18,63,19,49,22,13,27,48,12,49,19,16,19,38,30,14,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028065"
y <- c(0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(30,46,51,24,29,22,25,35,35,25,50,45,19,38,21,42,52,18,48,34,33,46,29,33,34,17,44,47,10,66,6,10,55,27,16,44,51)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028066"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,2,0,NaN,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(20,39,35,46,33,40,30,32,29,26,40,48,22,30,NaN,38,14,32,26,23,18,63,19,49,22,13,27,48,12,49,19,16,19,38,30,14,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028069"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0)
total_reads <- c(30,46,51,24,29,22,25,35,35,25,50,45,19,38,21,42,52,18,48,34,33,46,29,33,34,17,44,47,10,66,6,10,55,27,16,44,51)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028070"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0)
total_reads <- c(20,38,35,46,33,40,30,32,29,26,40,48,22,30,NaN,38,14,32,26,23,18,63,19,49,22,13,27,48,12,49,19,16,19,38,30,14,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028074"
y <- c(4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2)
total_reads <- c(30,46,51,24,29,22,25,25,35,25,50,45,19,38,21,42,52,18,48,34,33,46,29,33,34,17,44,47,10,66,6,10,55,27,16,44,51)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028075"
y <- c(0,2,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(20,38,35,46,33,39,30,32,29,26,40,48,22,30,NaN,38,14,32,26,23,18,63,19,48,22,13,27,48,12,49,19,15,19,38,30,14,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028078"
y <- c(0,0,2,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
total_reads <- c(92,89,77,49,56,55,45,36,63,41,71,57,56,71,27,69,79,40,74,51,63,82,62,52,51,50,68,100,31,84,36,31,75,56,57,74,88)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028086"
y <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(77,66,56,37,42,44,34,25,47,30,47,37,50,52,17,51,55,31,51,34,48,58,49,37,36,43,44,77,26,54,33,26,50,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028087"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0)
total_reads <- c(73,70,63,52,57,38,31,NaN,65,55,56,60,23,39,10,58,56,27,56,49,35,71,71,48,48,65,58,66,14,57,35,33,57,43,36,36,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028088"
y <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,2,0)
total_reads <- c(77,66,54,37,42,44,34,25,47,30,47,37,49,52,17,51,55,30,51,34,48,58,49,37,36,43,44,77,26,54,33,26,50,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028089"
y <- c(0,0,0,0,0,1,0,NaN,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
total_reads <- c(73,69,63,52,57,38,31,NaN,65,55,56,60,23,39,10,58,56,27,56,49,35,71,71,48,48,65,58,66,14,57,35,34,57,43,36,36,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028090"
y <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
total_reads <- c(77,66,54,37,42,44,34,25,47,30,47,37,49,52,17,51,55,31,51,34,48,58,49,37,36,43,44,77,26,54,33,26,51,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028091"
y <- c(0,0,2,0,0,0,0,NaN,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,1,0,0,0,0)
total_reads <- c(73,70,63,52,57,38,31,NaN,65,55,56,60,23,39,10,58,56,27,56,49,35,71,71,48,48,65,58,66,14,57,35,34,57,43,36,36,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028093"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(149,128,100,71,83,86,66,49,92,54,92,71,97,98,33,89,107,57,97,64,91,110,95,71,69,83,83,147,50,100,64,51,90,85,95,101,122)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028094"
y <- c(0,0,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0)
total_reads <- c(132,125,111,95,102,69,58,5,112,95,96,105,39,71,17,98,99,47,97,91,58,128,132,82,86,119,104,118,24,102,59,55,101,71,61,64,133)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028100"
y <- c(0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,2,0,0,0)
total_reads <- c(147,129,96,73,82,86,65,37,90,54,89,72,92,100,32,89,105,60,97,66,90,115,92,70,67,80,83,145,50,98,64,52,91,84,96,100,124)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028101"
y <- c(0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,3,0,0,0,0,0,0,0,0)
total_reads <- c(134,132,122,95,86,68,60,6,126,97,102,112,43,74,19,114,111,51,105,96,65,137,135,87,90,121,112,130,27,112,69,65,109,80,66,70,145)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028103"
y <- c(0,6,0,0,8,0,0,0,0,1,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
total_reads <- c(149,132,100,73,83,87,66,37,91,57,89,73,92,102,32,94,105,62,99,66,91,116,94,73,69,82,86,148,52,101,63,52,94,84,98,102,124)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028104"
y <- c(0,0,0,2,0,0,0,2,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0)
total_reads <- c(137,135,122,100,87,68,61,6,129,105,109,116,46,74,20,116,112,52,105,97,70,139,139,93,91,123,116,130,27,114,69,67,113,85,69,72,148)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028122"
y <- c(2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,4,0,0,1,0,0,0,0,0,2,0,0,0,0,2,4,0,0,0)
total_reads <- c(148,132,101,71,83,88,66,37,90,57,87,71,92,103,31,92,103,62,98,66,89,115,95,72,69,81,87,148,52,101,62,51,91,82,98,104,121)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028123"
y <- c(0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0)
total_reads <- c(141,133,123,100,84,68,61,6,129,106,109,116,46,72,19,116,110,53,103,97,68,138,140,95,90,125,117,131,27,114,69,67,113,85,74,72,148)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028128"
y <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)
total_reads <- c(76,66,72,37,42,56,44,25,47,30,47,51,50,52,17,51,55,31,49,34,46,59,49,46,47,43,44,77,26,54,33,26,48,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028129"
y <- c(0,0,1,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1)
total_reads <- c(72,68,62,52,53,34,31,NaN,65,55,54,57,24,37,10,58,55,26,53,49,35,71,71,48,46,64,61,65,14,57,36,35,57,43,40,36,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028136"
y <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,3,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0)
total_reads <- c(77,66,58,37,42,46,34,25,47,30,47,37,50,52,17,51,55,31,49,34,46,59,49,37,36,43,44,77,26,54,33,26,48,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028137"
y <- c(0,0,0,0,0,0,2,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
total_reads <- c(73,68,63,52,53,34,31,NaN,65,55,54,58,23,37,10,58,56,27,54,49,35,71,71,48,46,65,60,66,14,57,36,34,57,43,40,36,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028142"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)
total_reads <- c(77,66,58,36,42,46,34,25,47,30,47,37,50,52,17,51,55,31,49,34,46,59,49,37,36,43,44,77,26,54,33,26,47,43,50,54,63)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028369"
y <- c(0,1,0,1,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,NaN,0,0,2,0,0,0)
total_reads <- c(21,29,36,36,NaN,17,15,7,35,15,30,21,8,36,13,48,35,13,19,37,18,48,24,20,8,29,20,46,9,42,NaN,15,33,26,27,27,38)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028384"
y <- c(0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,4,0,0,0,4,2,2,0,0,3,0,2,0,2,0,2,0,0,0,0,2)
total_reads <- c(36,59,72,68,7,32,27,14,65,28,57,40,13,70,21,94,68,22,35,70,35,90,42,40,15,52,39,89,15,82,5,27,64,50,47,51,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028385"
y <- c(0,0,0,0,0,0,0,0,0,2,0,4,0,0,0,4,0,2,0,0,3,5,2,0,0,0,0,0,4,0,0,0,0,0,8,0,0)
total_reads <- c(46,42,42,56,7,20,35,24,35,55,44,47,14,34,16,69,53,27,44,54,42,81,44,24,27,31,46,50,35,89,18,27,36,23,63,26,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028388"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,1,0,1,0,0,0,0,0,0,0,0,0)
total_reads <- c(36,59,72,68,7,32,27,14,65,28,57,40,13,70,21,94,68,22,35,70,35,90,42,40,15,52,39,89,15,82,5,27,64,50,47,51,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028389"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0)
total_reads <- c(46,42,42,56,7,20,35,23,35,55,44,47,14,34,16,69,53,27,44,54,42,81,44,24,27,31,46,50,35,89,18,27,35,23,63,26,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028394"
y <- c(0,0,0,2,0,0,0,0,0,0,4,2,0,2,0,4,0,0,0,5,2,6,0,0,0,0,0,2,0,0,0,0,0,4,0,0,2)
total_reads <- c(36,59,72,68,7,32,27,14,65,28,57,40,13,70,21,94,68,22,35,70,35,90,42,40,15,52,39,89,15,82,5,27,64,50,47,51,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028395"
y <- c(0,0,0,0,3,0,0,0,0,2,0,3,0,0,0,2,0,2,4,0,0,2,0,1,0,0,0,0,0,0,0,0,0,2,4,0,0)
total_reads <- c(45,42,42,56,7,20,35,23,34,54,44,47,14,34,16,67,53,27,44,54,42,81,44,24,27,31,46,50,35,89,18,26,36,23,62,26,75)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028398"
y <- c(0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,2,2,5,0,1,0,1,0,1,0,0,0,0,0,1,0,1,4,0,0)
total_reads <- c(100,87,91,77,28,71,64,9,90,70,81,37,23,73,43,91,91,49,53,121,53,107,64,73,66,74,56,116,17,119,47,45,76,75,94,69,106)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028414"
y <- c(3,0,0,0,0,1,1,NaN,2,1,0,0,0,0,0,1,3,0,0,4,2,1,0,0,0,0,1,0,3,2,0,0,0,1,3,2,0)
total_reads <- c(84,61,55,43,25,57,53,NaN,59,57,54,18,17,40,35,48,58,37,36,87,36,65,46,54,59,51,37,73,10,79,45,33,45,53,75,45,70)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028415"
y <- c(4,0,1,1,0,4,0,NaN,1,4,1,0,0,1,0,1,1,3,1,4,0,3,0,3,0,0,1,0,0,0,0,0,2,0,3,0,1)
total_reads <- c(79,54,60,51,32,39,32,NaN,46,54,64,31,18,56,50,50,49,49,41,61,51,59,49,39,66,68,57,81,6,54,49,49,50,46,62,69,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028423"
y <- c(0,0,6,2,0,4,0,NaN,6,9,0,2,2,2,0,2,0,0,0,2,2,2,2,2,0,7,2,0,0,0,0,0,0,2,4,0,0)
total_reads <- c(162,116,109,86,49,109,104,NaN,114,109,104,36,33,79,65,88,113,74,70,170,70,128,86,106,117,98,74,142,18,153,89,64,87,104,142,85,136)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028424"
y <- c(2,0,0,0,0,0,0,0,1,0,4,2,0,0,4,0,0,0,4,0,0,4,0,0,8,6,2,0,0,2,0,0,2,2,4,2,0)
total_reads <- c(146,108,116,102,60,78,58,6,89,106,122,54,36,108,88,101,95,97,80,121,100,114,97,76,131,135,111,160,12,105,98,97,95,87,121,136,146)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028429"
y <- c(6,4,0,0,0,2,3,NaN,4,7,0,0,2,0,0,2,6,0,0,6,2,5,4,4,0,5,0,0,0,0,0,0,1,4,6,2,0)
total_reads <- c(163,116,110,86,48,110,104,NaN,114,111,105,36,33,79,65,88,113,74,70,170,71,128,86,106,117,98,74,142,19,155,89,65,88,104,142,86,138)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028430"
y <- c(2,2,2,0,0,2,0,0,0,0,0,0,0,0,6,0,2,3,2,0,2,8,4,2,14,8,4,4,0,0,0,0,0,2,6,0,0)
total_reads <- c(146,108,116,102,60,78,58,6,89,106,122,55,36,108,88,101,95,97,80,121,100,115,97,76,131,135,111,160,12,105,98,97,95,87,121,136,147)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028451"
y <- c(5,0,1,0,1,1,1,NaN,1,6,0,0,0,0,0,0,2,1,0,5,0,3,2,5,5,2,0,1,0,3,2,2,0,3,3,0,3)
total_reads <- c(84,61,55,43,25,56,53,NaN,60,57,54,18,17,41,35,44,58,37,37,87,36,65,46,52,59,51,37,73,10,79,45,33,45,53,73,45,70)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028452"
y <- c(2,1,0,2,0,1,0,NaN,3,5,0,1,0,1,2,2,0,3,1,1,0,2,1,2,5,2,2,1,0,3,0,1,1,4,2,3,1)
total_reads <- c(79,54,58,51,32,40,32,NaN,44,54,62,28,18,56,48,51,49,48,40,59,48,57,49,38,66,68,56,79,6,53,46,49,49,44,62,68,74)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100028468"
y <- c(0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(84,61,55,43,25,55,53,NaN,61,57,54,18,17,41,37,42,58,38,37,87,36,65,46,54,60,51,37,74,10,80,45,33,48,53,71,46,72)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100035482"
y <- c(47,52,44,38,20,43,12,30,65,49,46,25,27,44,36,45,39,53,52,45,48,40,38,50,57,49,53,101,17,93,24,36,40,32,55,41,61)
total_reads <- c(53,58,46,46,21,50,15,30,74,53,56,28,27,50,36,52,43,53,54,58,51,55,42,53,73,53,55,118,18,103,27,43,40,37,66,51,69)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100035542"
y <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0)
total_reads <- c(53,54,47,47,28,50,15,30,74,55,57,29,28,52,36,52,44,54,54,60,50,58,43,52,64,55,55,119,17,106,27,43,40,40,66,51,68)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100050815"
y <- c(17,26,16,21,14,17,7,12,22,12,14,9,8,13,NaN,18,12,11,11,25,17,7,8,18,11,6,28,21,NaN,15,9,14,18,22,16,26,25)
total_reads <- c(19,30,19,25,18,18,13,12,26,19,20,10,13,13,NaN,21,12,14,17,28,18,11,9,20,21,7,32,23,NaN,19,9,15,20,28,23,30,31)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100050872"
y <- c(0,2,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0,0,1,0,0,0,NaN,0,0,0,0,1,0,2,0)
total_reads <- c(19,32,19,25,18,18,13,12,26,20,20,9,13,14,NaN,21,12,14,17,28,18,10,10,19,20,7,33,22,NaN,17,9,14,18,28,23,31,30)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100061625"
y <- c(8,24,11,18,0,22,10,0,24,9,20,12,2,11,6,16,4,12,20,19,7,35,20,21,18,24,22,54,15,13,10,8,16,13,20,29,26)
total_reads <- c(50,74,32,39,8,45,26,32,39,35,41,21,13,23,33,34,51,50,52,44,22,65,50,53,44,51,53,83,19,36,26,31,33,33,66,73,42)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100061626"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,0,0,0,1,0,0,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,0,0)
total_reads <- c(12,28,23,30,21,18,11,5,29,17,24,9,7,17,NaN,30,26,26,20,28,7,36,30,17,18,21,33,29,NaN,25,25,23,23,28,15,44,28)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")

cgid <- "!chr10:100061674"
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
total_reads <- c(47,69,26,36,8,44,23,30,37,34,36,18,12,23,34,29,51,50,46,39,19,62,40,52,40,47,52,79,19,36,26,30,30,37,56,72,39)
survival <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
Sex <- c(1,1,2,1,2,2,2,1,2,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1,2)
capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ survival+Sex, family=ZIBB)),file="NUL")
capture.output(s <- summary(fit, save=TRUE), file="NUL")
pval <- s$coef.table[,4]
coef <- s$coef.table[,1]
cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\t")
cat("\n")
