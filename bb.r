# fit beta-binomial models
# alternative
calc_logl_alt<-function(param) {
    phi=exp(param[1])/(exp(param[1])+1)
    mu=param[2]
    gamma=param[3]
    theta=exp(c*gamma+mu)/(exp(c*gamma+mu)+1)

    alpha=theta*(1-phi)/phi
    beta=(1-theta)*(1-phi)/phi
    logl=sum(lbeta(alpha+x, beta+y-x))-sum(lbeta(alpha, beta))
    return(logl)
}

# null
calc_logl_null<-function(param) {
    phi=exp(param[1])/(exp(param[1])+1)
    mu=param[2]
    gamma=0
    theta=exp(mu)/(exp(mu)+1)

    alpha=theta*(1-phi)/phi
    beta=(1-theta)*(1-phi)/phi

    logl=sum(lbeta(alpha+x, beta+y-x))-length(x)*lbeta(alpha, beta)
    return(logl)
}

# fit the beta-binomial for each site, where mcounts is a vector of methylated counts, tcounts is a vector of total counts and cvt is a vector of covariates
bbfit<-function (mcounts, tcounts, cvt) {
    x=mcounts; y=tcounts; c=cvt
    param=c(-1,0.5,0)
    phi=exp(param[1])/(exp(param[1])+1)
    mu=param[2]
    gamma=param[3]
    theta=exp(c*gamma+mu)/(exp(c*gamma+mu)+1)

    alpha=theta*(1-phi)/phi
    beta=(1-theta)*(1-phi)/phi

    phi=1/(alpha+beta+1); mu=log(alpha/beta); gamma=0;
    par=c(-2, mu)

    for (i_iter in 1:100) {
        if (i_iter==1) {
            null=optim(par,calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
        } else {
            null=optim(alt$par[1:2],calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
        }
        alt=optim(c(null$par,0),calc_logl_alt, control=list(fnscale=-1),  method = c("CG"))

        while (null$value==0 | alt$value==0 | alt$value<(null$value-0.01)) {
            alpha=runif(1); beta=runif(1)
            phi=1/(alpha+beta+1); mu=log(alpha/beta); gamma=0;
            par=c(log(phi), mu)

            null=optim(par,calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
            alt=optim(c(null$par,0),calc_logl_alt, control=list(fnscale=-1),  method = c("CG"))
        }

        lr=2*(alt$value-null$value)
        if (lr<0) {lr=0}

        if (i_iter==1) {
            lr_old=lr
            lr_new=lr
        } else {
            lr_new=lr
            if (abs(lr_new-lr_old)<0.01 & lr!=0) {break}
            lr_old=lr
        }
    }

    return(1-pchisq(lr,df=1))
}

