    {
    SEX=$5;
    if(SEX=="male") SEX="1";
    else if(SEX=="female") SEX="2"; 
    
    PHENO="0";
    if($6=="affected" || $6=="case") PHENO="2";
    else if($6=="unaffected" || $6=="control") PHENO="1";
        
    printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,SEX,PHENO);
    }