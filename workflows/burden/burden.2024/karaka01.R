library("RSpectra",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")
library("SPAtest",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")
library("SKAT",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")
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
    data <- prepareData(genotypes,population,variants,"SKAT_adjusted")
    MAFs<- data$MAFs
    
    rez <- SKAT(Z=as.matrix(data$genot),weights=1/sqrt(data$n*MAFs*(1-MAFs)), obj=data$obj2, kernel="linear.weighted", method="optimal")
    list(p.value= rez$p.value)
    }


matildeFunction=function(){
  data <- matrix(genotypes,nrow = nrow(population), byrow = FALSE)

  MAFs<- variants$AF
  n <- nrow(population)
  phenot <- population[,3] 
  genot <- data
  genot[genot==-9] <- 0
  perm <- 100
  #liste de test a faire, les resultats sortent dans le meme ordre
  assotests <- c("CAST", "SKAT", "SKATO", "SKAT_adjusted", "SKATO_adjusted","KBAC")
  
  #********************
  #necessaire pour CAST
  #********************
  somme<- apply(genot,1,FUN = function(data){sum(as.numeric(data))})
  SV1 <- rep(0,nrow(genot))
  SV1[somme>= 1] <- 1
  matSV1 <- table(phenot,factor(SV1,c(0,1)))
  tesSV1 <- fisher.test(matSV1)

  #*****************************
  #necessaire pour SKAT et SKATO
  #*****************************
  if (sum(c("SKAT","SKATO")%in%assotests)>=1){
    obj=SKAT_Null_Model(phenot~1, out_type="D", Adjustment=F)
  }
  
  if (sum(c("SKAT_adjusted","SKATO_adjusted")%in%assotests)>=1){
    obj2=SKAT_Null_Model(phenot~1, out_type="D", Adjustment=T)
  }
  
  
  resultats_tests=NULL
  for (test in assotests){
    res=switch(test,
		CAST=fisher.test(matSV1)$p.value,
               	SKAT=SKAT(Z=as.matrix(genot),weights= 1/sqrt(n*MAFs*(1-MAFs)), obj=obj, kernel="linear.weighted", method="davies")$p.value,
               	SKATO=SKAT(Z=as.matrix(genot),weights=1/sqrt(n*MAFs*(1-MAFs)), obj=obj, kernel="linear.weighted", method="optimal")$p.value,
               	SKAT_adjusted=SKAT(Z=as.matrix(genot),weights=1/sqrt(n*MAFs*(1-MAFs)), obj=obj2, kernel="linear.weighted", method="davies")$p.value,
               	SKATO_adjusted=SKAT(Z=as.matrix(genot),weights=1/sqrt(n*MAFs*(1-MAFs)),obj=obj2, kernel="linear.weighted", method="optimal")$p.value,
               	)    
    resultats_tests=c(resultats_tests, res)
  }
  print(vcf.title)

  fail_flag=-9999.999

  writeLines(paste0("UPDATE VCF set ",
                 "KAST=",ifelse(is.nan(resultats_tests[1]),fail_flag,resultats_tests[1]),
		 ",SKAT=" ,ifelse(is.nan(resultats_tests[2]),fail_flag,resultats_tests[2]),
		 ",SKATO=" ,ifelse(is.nan(resultats_tests[3]),fail_flag,resultats_tests[3]),
		 ",SKAT_ADJUSTED=" ,ifelse(is.nan(resultats_tests[4]),fail_flag,resultats_tests[4]),
		 ",SKATO_ADJUSTED=" ,ifelse(is.nan(resultats_tests[5]),fail_flag,resultats_tests[5]),
                 ",KBAC=" ,ifelse(is.nan(resultats_tests[6]),fail_flag,resultats_tests[6]),		 
                 ",SVCAS=",ifelse(is.nan(matSV1[2,2]),fail_flag,matSV1[2,2]),
                 ",SVCTRL=",ifelse(is.nan(matSV1[1,2]),fail_flag,matSV1[1,2]),
                 ",ODDR=",ifelse(is.infinite(tesSV1$estimate),fail_flag,round(tesSV1$estimate,4)),
                 " WHERE ID=\'",vcf.title,"\';"),stdout())
}



##matildeFunction()
