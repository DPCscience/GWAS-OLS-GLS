###  GWAS-OLS-GLS


OLS=Ordinary Least Squares
GLS=Generalized Least Squares
GWAS=Genome Wide Association Study


*Contact*: Gustavo de los Campos (gdeloscampos@gmail.com)

The following funcitons can be used to carry out GWAS using OLS and GLS for correlated errors. 
There are multiple software that can perform GWAS using OLS and GLS more efficiently than the code presented below. The only  objective of this repository is to illustrate how some of the required computations can be perofrmed and to demonstrate a few ways of improving computational speed.

### 1: GWAS-OLS using the lm() function

```R
GWAS.OLS.lm<-function(y,W=NULL,X,verbose=F){
     n=length(y)
    p=ncol(X)
    if(is.null(W)){
        W=matrix(nrow=n,ncol=1,1)
    }else{
        W=cbind(1,W)
    }
    OUT=matrix(nrow=p,ncol=4,NA)
    colnames(OUT)<-c('estimate','SE','z-score','pValue')
    for(i in 1:p){
         Z=cbind(X[,i],W)
        fm=lm(y~Z-1)
        OUT[i,]=summary(fm)$coef[1,]
        if(verbose){ print(i) }
    }
    return(OUT)
}
```

### 2: GWAS-OLS using explicit matrix operations

The following code improves speed relative to GWAS.OLS.lm (see above)

```R
GWAS.OLS.fast<-function(y,W=NULL,X,verbose=F){
    ## Brute force GLS
    n=length(y)
    p=ncol(X)
    if(is.null(W)){ W=matrix(nrow=n,ncol=1,1) }
    OUT=matrix(nrow=p,ncol=4,NA)
    colnames(OUT)<-c('estimate','SE','z-score','pValue')
    for(i in 1:p){
        Z=cbind(X[,i],W)
         C=crossprod(Z)     # this is X'X
        rhs=crossprod(Z,y)  # this is X'y
        CInv=chol2inv(chol(C)) # gets solve(X'X)
        sol<-crossprod(CInv,rhs) # this gets the OLS estimate
        OUT[i,1]=sol[1]
        error=y-Z%*%sol
        RSS=sum(error^2)
        vE=RSS/(n-ncol(Z))
        OUT[i,2]=sqrt(CInv[1,1]*vE)
        if(verbose){ print(i) }
    }
    OUT[,3]=OUT[,1]/OUT[,2]
    OUT[,4]=2*(1-pnorm(abs(OUT[,3])))
     return(OUT)
}
```

### 3) GLS using standard matrix operations


```R
GWAS.GLS<-function(y,W=NULL,X,G,vU,vE,verbose=F){
    ## Brute force GLS
    n=length(y)
    p=ncol(X)
    if(is.null(W)){ W=matrix(nrow=n,ncol=1,1) }
    R=G*vU
    diag(R)<-diag(R)+vE
    RInv=chol2inv(chol(R))
    OUT=matrix(nrow=p,ncol=4,NA)
    colnames(OUT)<-c('estimate','SE','z-score','pValue')
    Z=cbind(0,W)
    for(i in 1:p){
        Z[,1]=X[,i]
         C=crossprod(Z,crossprod(RInv,Z))
        rhs=crossprod(Z,RInv)%*%y
        CInv=chol2inv(chol(C))
        sol<-crossprod(CInv,rhs)
        OUT[i,1]=sol[1]
        error=y-Z%*%sol
        RSS=sum(error^2)
        vE=RSS/(n-ncol(Z))
        OUT[i,2]=sqrt(CInv[1,1])
        if(verbose){ print(i) }
    }
    OUT[,3]=OUT[,1]/OUT[,2]
    OUT[,4]=2*(1-pnorm(abs(OUT[,3])))
    return(OUT)

}

```

### 4) GLS using the eigenvalue decomposition of G

```R
 GWAS.GLS.eig<-function(y,W=NULL,X,G=NULL,V=NULL,d=NULL,vU,vE,verbose=F){
    ##  GLS eith eigenvectors
    #X'[VDV'D*vG+IvE]^-1X=X'V[D*k+I]^-1V'X/vE
    # Z=SV'X  where S=diag(1/sqrt(d*vG+vE))
    # Z'Z= X'VS'SV'X
    if(is.null(V)){ EVD=eigen(G) }
    V=EVD$vectors
    d=EVD$values
    n=length(y)
    p=ncol(X)
  
    if(is.null(W)){ 
    	W=matrix(nrow=n,ncol=1,1) 
    }else{
    	W=cbind(1,W)
    }
    dStar=d*vU+vE
    OUT=matrix(nrow=p,ncol=6,NA)
    colnames(OUT)<-c('estimate','SE','z-score','pValue-z','RSS-Dif','pValue-Chisq')
    for(i in 1:ncol(V)){ V[,i]=V[,i]/sqrt(dStar[i]) }  # V%*%diag(1/sqrt(d*vG+vE))
    W=crossprod(V,W)
    y=crossprod(V,y)
    C0=crossprod(W)
    rhs0<-crossprod(W,y)
    CInv0=chol2inv(chol(C0))
    sol0<-crossprod(CInv0,rhs0)
    error0<-y-W%*%sol0
    WRSS0<-sum(error0^2)
    for(i in 1:p){
        Z=crossprod(V,X[,i])
        Z=cbind(Z,W) 
        C=crossprod(Z)
        rhs=crossprod(Z,y)
        CInv=chol2inv(chol(C))
        sol<-crossprod(CInv,rhs)
        OUT[i,1]=sol[1]  # estimate
        OUT[i,2]<-sqrt(CInv[1,1]) # SE
        error=y-Z%*%sol
        WRSS1=sum(error^2)      
        OUT[i,5]<-WRSS0-WRSS1 # LRT-statistic
        if(verbose){ print(i) }
    }
    
    OUT[,3]=OUT[,1]/OUT[,2] # z-stat
    OUT[,4]=2*(1-pnorm(abs(OUT[,3]))) # p-value normal test
    OUT[,6]=1-pchisq(df=1,q=OUT[,5])
    return(OUT)
}

```

### 5) GLS Using eigenvectors and a general interface for contrasts 

```R
GWAS.GLS.eig2<-function(y,W=NULL,X,G=NULL,V=NULL,d=NULL,vU,vE,verbose=F,getContrasts,bruteForce=T){
    ##  GLS eith eigenvectors
    #X'[VDV'D*vG+IvE]^-1X=X'V[D*k+I]^-1V'X/vE
    # Z=SV'X  where S=diag(1/sqrt(d*vG+vE))
    # Z'Z= X'VS'SV'X
    if(is.null(V)){ EVD=eigen(G) }
    V=EVD$vectors
    d=EVD$values
    n=length(y)
    p=ncol(X)
    if(is.null(W)){ 
    	W=matrix(nrow=n,ncol=1,1) 
    }else{
    	W=cbind(1,W)
    }
    dStar=d*vU+vE
    OUT=matrix(nrow=p,ncol=3,NA)
    colnames(OUT)<-c('RSS-Dif','df','pValue')
    for(i in 1:ncol(V)){ V[,i]=V[,i]/sqrt(dStar[i]) }  # V%*%diag(1/sqrt(d*vG+vE))
    W=crossprod(V,W)
    y=crossprod(V,y)
    C0=crossprod(W)
    rhs0<-crossprod(W,y)
    CInv0=chol2inv(chol(C0))
    sol0<-crossprod(CInv0,rhs0)
    error0<-y-W%*%sol0
    WRSS0<-sum(error0^2)
    TMP<-getContrasts(x=X[,i])
    dimZ<-ncol(TMP)
    dimC<-ncol(C0)+dimZ
    for(i in 1:p){
    	TMP<-getContrasts(x=X[,i])
        Z=crossprod(V,TMP)
        if(bruteForce){
          Z=cbind(Z,W) 
          C<-crossprod(Z)
          rhs=crossprod(Z,y)
        }else{
         C=matrix(nrow=dimC,ncol=dimC,0)
         C[(dimZ+1):dimC,(dimZ+1):dimC]<-C0
         C[1:dimZ,1:dimZ]=crossprod(Z)
         C[1:dimZ,(dimZ+1):dimC]<-crossprod(Z,W)
         C[(dimZ+1):dimC,1:dimZ]<-C[1:dimZ,(dimZ+1):dimC]
         rhs<-c(crossprod(Z,y),rhs0)
        }
               
        CInv=chol2inv(chol(C))
        sol<-crossprod(CInv,rhs)
        OUT[i,2]<-dimZ
        if(bruteForce){
	        error=y-Z%*%sol
    	}else{
    	    error=y-Z%*%sol[1:dimZ]-W%*%sol[(dimZ+1):dimC]
        }
        WRSS1=sum(error^2)      
        OUT[i,1]<-WRSS0-WRSS1 # LRT-statistic
        
        if(verbose){ print(i) }
    }
    OUT[,3]=1-pchisq(df=1,q=OUT[,1])
    return(OUT)
}
```

### 6) A few examples/tests of the code

```
library(BGLR)
data(wheat)
X=wheat.X
y=wheat.Y[,1]

system.time( OUT.OLS<-GWAS.OLS.lm(y=y,X=X,verbose=F) )
system.time( OUT.OLS.fast<-GWAS.OLS.fast(y=y,X=X,verbose=F) )
G=tcrossprod(scale(X))/ncol(X)
system.time( OUT.GLS<-GWAS.GLS(y=y,X=X,G=G,vU=.5,vE=.5,verbose=F) )
EVD<-eigen(G)
system.time( OUT.GLS2<-GWAS.GLS.fast(y=y,X=X,V=EVD$vectors,d=EVD$values,vU=.5,vE=.5,verbose=F) )

 system.time( OUT.GLS3<-GWAS.GLS.fast2(y=y,X=X,V=EVD$vectors,d=EVD$values,vU=.5,vE=.5,verbose=F,
                                            getContrasts=function(x){ return(x) } 
                                       )
 ```
