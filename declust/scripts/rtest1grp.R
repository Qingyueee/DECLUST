rtest <- function(y, group, single.only = FALSE, robust.var = TRUE, maxiter = 5, err.limit = 0.01, espXWX=1e-16, upreg = TRUE) {
  # constant matrices
  X <- model.matrix(~group, contrasts = list(group = "contr.sum"))
  gsize <- table(group)
  ng <- length(gsize)
  a <- as.factor(1:ng)
  X0 <- model.matrix(~a, contrasts = list(a = "contr.sum"))
  
  # # new group label
  # unique_values <- data.frame(group = unique(group), a = as.numeric(as.character(unique(a))))
  # print(unique_values)
  
  # starting estim
  m <- tapply(y, group, median)
  b <- solve(t(X0) %*% X0, t(X0) %*% m)
  if (sum(b < 0.1) == ng) {
    return(list(rcode = 0, b = b))
  } ## NOTE: no need to continue
  err <- exp(25)
  
  ## IWLS
  iter <- 1
  while (max(err) > err.limit & iter <= maxiter) {
    oldb <- b
    iter <- iter + 1
    yfit <- c(X %*% b)
    res <- y - yfit
    phi <- mean(res^2 / ifelse(yfit < 0.1, 0.1, yfit))
    phi <- ifelse(phi > 10, 10, phi)
    v <- phi * yfit
    v <- ifelse(v > 0.1, v, 0.1) ## limit lower bound: avoid blowing up 1/v
    
    ## robust weight: need to control within each group!
    q3 <- tapply(abs(res), group, quantile, 0.75)
    q3 <- rep(q3, gsize) # long version
    rw <- ifelse(abs(res) <= q3, 1, q3 / (abs(res) + 0.0001))
    ## full weight:
    w <- rw / (v)
    
    ## estimate
    XWX <- t(X) %*% (X * w)
    XWy <- t(X) %*% (y * w)
    if (abs(det(XWX)) <= 1e-20) {
      XWX <- XWX + diag(espXWX, ng) #it was 1e-16 by Tian
    } # if XWX is singular or computationally singular,
    # add a small value on diagonal to make it invertible - Tian-18062019
    b <- solve(XWX, XWy, tol = 1e-25) # use a smaller tolerance - Tian-18062019
    err <- (oldb - b) / oldb
    # print(c(b))
  } # for iter
  
  ## estimates, including last group: called '5' here.
  b5 <- -sum(b[-1])
  n4 <- ng - 1 ##
  
  ## robust variance formula (see page 394 of In All Likelihood)
  XWXi <- solve(XWX, tol = 1e-20) # use a smaller tolerance - Tian-18062019
  yfit <- c(X %*% b)
  res <- y - yfit
  w2 <- w^2 * res^2
  XWX2 <- t(X) %*% (X * w2)
  V <- XWXi %*% XWX2 %*% XWXi
  
  if (all(diag(V) >= 0)) { # Tian-18062019
    
    if (robust.var) {
      se <- sqrt(diag(V))[-1]
      se5 <- sqrt(rep(1, n4) %*% V[-1, -1] %*% rep(1, n4))
    }
    
    if (!robust.var) {
      se <- sqrt(diag(XWXi))[-1]
      se5 <- sqrt(rep(1, n4) %*% XWXi[-1, -1] %*% rep(1, n4))
    }
  } else {
    cat(paste0("i = ", i, ". Caught an error, the estimated variance is negative. Try non robust variance.\n"))
    se <- sqrt(diag(XWXi))[-1]
    se5 <- sqrt(rep(1, n4) %*% XWXi[-1, -1] %*% rep(1, n4))
  } # End of if else  # Tian-18062019
  
  ## collect results: drop the intercept!!
  b.all <- c(b[-1], b5)
  se.all <- c(se, se5)
  
  ## choose max-group
  stat <- b.all / se.all
  
  if (upreg) {
    imax <- which((stat) == max((stat))) ## looking only on over-expression
  } else {
    imax <- which((abs(stat)) == max((abs(stat))))
  } ## looking for over and low-expression
  
  group_levels <- levels(group)
  
  if (length(imax) > 1) {
    return(list(rcode = 0, b = b))
  }
  
  ## second-stage: test if the rest is equal to each other
  ## last group requires special handling, since it is outside the IWLS
  
  if (single.only) {
    ## collect output
    out <- c(maxgroup = group_levels[imax], tstat = stat[imax])
    return(list(rcode = 1, b = b.all, se = se.all, stat = stat, out = out))
  } else {
    ## contrast matrix in the second stage test
    cmat <- (diag(ng - 1) - contr.sum(ng)[-1, ])[-1, ]
    cmat.last <- (diag(ng - 1) - contr.sum(ng)[-1, ])[c(-ng + 1), ]
    if (imax < ng) {
      irest <- c(1:n4)[-imax]
      V2 <- (V[-1, -1])[c(imax, irest), c(imax, irest)] ## rearrange to match the order of imax and irest
      cvar <- cmat %*% V2 %*% t(cmat)
      b2 <- cmat %*% c(b.all[imax], b.all[irest])
    }
    if (imax == ng) { ##  irest = 1:n4
      V2 <- V[-1, -1]
      cvar <- cmat.last %*% V2 %*% t(cmat.last)
      b2 <- cmat.last %*% b.all[1:n4]
    }
    
    ## test singularity cvar
    # emin= min(eigen(cvar)$value)
    # if (emin>0.1) cstat = c( t(b2)%*%solve(cvar)%*% b2) else cstat=  c(t(b2) %*% solve(cvar+ 0.1*diag(ncol(cvar))) %*% b2) # Setia
    
    # emin= min(eigen(cvar)$value)/max(eigen(cvar)$value)
    # sfactor=min(max(eigen(cvar)$value),1)
    # if (emin>0.1) cstat = c( t(b2)%*%solve(cvar)%*% b2) else cstat=c(t(b2) %*% solve(cvar+ sfactor*diag(ncol(cvar))) %*% b2)#nghia-25052015
    
    # Tian-18062019
    emin <- min(Re(eigen(cvar)$value)) / max(Re(eigen(cvar)$value)) # keep only the real part of eigen values
    sfactor <- min(max(Re(eigen(cvar)$value)), 1)
    if (abs(det(cvar)) <= 1e-16) {
      cvar <- cvar + diag(1e-16, ncol(cvar))
    } # if cvar is singular or computationally singular,
    # add a small value on diagonal to make it invertible
    if (emin > 0.1) cstat <- c(t(b2) %*% solve(cvar) %*% b2) else cstat <- c(t(b2) %*% solve(cvar + sfactor * diag(ncol(cvar))) %*% b2)
    
    ## collect output
    out <- c(maxgroup = group_levels[imax], tstat = stat[imax], chisqstat = cstat)
    
    return(list(rcode = 1, b = b.all, se = se.all, stat = stat, out = out))
  }
} ## end function
