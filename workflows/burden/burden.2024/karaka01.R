#library("RSpectra",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")
#library("SPAtest",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")
#library("SKAT",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")

library("RSpectra")
library("SPAtest")
library("SKAT")


##library("KBAC",lib.loc="/commun/data/users/karaka/BURDEN_JVARKIT/lib")
#besoin de : vcf.title,genotypes, population, variants

prepareData<-function(genotypes,population,variants,assoctest) {
  data <- matrix(genotypes,nrow = nrow(population), byrow = FALSE)

  n <- nrow(population)
  phenot <- population[,3] 
  genot <- data
  genot[genot==-9] <- 0
  perm <- 100
  #liste de test a faire, les resultats sortent dans le meme ordre
  
  #********************
  #necessaire pour CAST
  #********************
  somme<- apply(genot,1,FUN = function(data){sum(as.numeric(data))})
  SV1 <- rep(0,nrow(genot))
  SV1[somme>= 1] <- 1
  matSV1 <- table(phenot,factor(SV1,c(0,1)))

  #*****************************
  #necessaire pour SKAT et SKATO
  #*****************************
  if (assoctest=="SKAT" || assoctest=="SKATO"){
    obj <- SKAT_Null_Model(phenot~1, out_type="D", Adjustment=F)
  }  else {
    obj <- NULL
  }
  
  if (assoctest=="SKAT_adjusted" || assoctest =="SKATO_adjusted"){
    obj2 <- SKAT_Null_Model(phenot~1, out_type="D", Adjustment=T)
  } else {
    obj2 <- NULL
    }
  
  list(genot=genot,matSV1=matSV1,MAFs=variants$AF,obj=obj,obj2=obj2,n=n)
  }

applyCAST <- function(genotypes,population,variants) {
    data <- prepareData(genotypes,population,variants,"CAST")
    rez <- fisher.test(data$matSV1)
    return (list(p.value=rez$p.value,
        round(rez$estimate,4),
        matSV1=data$matSV1,
        CASE_ALT= data$matSV1[2,2],
        CASE_REF= data$matSV1[2,1],
        CTRL_ALT= data$matSV1[1,2],
        CTRL_REF= data$matSV1[1,1]
        ))
    }

applySKAT <-  function(genotypes,population,variants) {
    data <- prepareData(genotypes,population,variants,"SKAT")
    MAFs<- data$MAFs
    
    rez <- SKAT(Z=as.matrix(data$genot),weights=1/sqrt(data$n*MAFs*(1-MAFs)), obj=data$obj, kernel="linear.weighted", method="davies")

    list(p.value= rez$p.value)
    }

applySKATO <-  function(genotypes,population,variants) {
    data <- prepareData(genotypes,population,variants,"SKATO")
    MAFs<- data$MAFs
    
    rez <- SKAT(Z=as.matrix(data$genot),weights=1/sqrt(data$n*MAFs*(1-MAFs)), obj=data$obj, kernel="linear.weighted", method="optimal")

    list(p.value= rez$p.value)
    }

applySKATAdjusted <-  function(genotypes,population,variants) {
    data <- prepareData(genotypes,population,variants,"SKAT_adjusted")
    MAFs<- data$MAFs
    
    rez <- SKAT(Z=as.matrix(data$genot),weights=1/sqrt(data$n*MAFs*(1-MAFs)), obj=data$obj2, kernel="linear.weighted", method="davies")
    list(p.value= rez$p.value)
    }

applySKATOAdjusted <-  function(genotypes,population,variants) {
    data <- prepareData(genotypes,population,variants,"SKATO_adjusted")
    MAFs<- data$MAFs
    
    rez <- SKAT(Z=as.matrix(data$genot),weights=1/sqrt(data$n*MAFs*(1-MAFs)), obj=data$obj2, kernel="linear.weighted", method="optimal")
    list(p.value= rez$p.value)
    }

