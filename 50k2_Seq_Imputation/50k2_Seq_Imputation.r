# The following codes are for setting up an imputatioin pipeline for 50k SNP data upto sequence level using EAGLE and Minimac3
# You need your 50k data either in PLINK format or in vcf format, you also need your reference sequence in vcf format
# The steps are as follows: 50k - 700k - SEQ
    # if you have 50k data in PLINK binary or txt format, you need to convert them into vcf format each chromosome at a time
    # for each chromosome use EAGLE to phase the 50k vcf file
    # Use 700k data as reference; to do that you need 700k data phased as well one chromosome at a time
    # Once you habe 50k and 700k phased data use Minimac3 to impute 50k data to 700k. The minimac3 output is already phased so you do not need to phase it again
    # Now you can merge imputed 700k and reference 700k and then impute them upto sequence level using the reference sequence or you can seperately impute them upto sequence level


#########################################################################################################################################################################################

#Phase the 50k PLINK file one chromosome at a time using EAGLE

#PLINK will create .VCF files from binary plink file 
#EAGLE with phase and impute the .vcf files
#INPUT files for plink .bed,.bim and .fam
#            for EAGLE .vcf and genetic_map_1cMperMb.txt (this file is required if the .bim file does not have and CM distance and this file is avaiable in EAGLE software)

#Output will be phased and missing value imputed 50k .vcf file one chromosome per file


plinkFileName="mergedHANWOO"  #50k PLINK file name


#read the .bim file to get the number of chromosomes (if the data in PLINK binary format)
bim=read.table(paste(plinkFileName,".bim",sep=""),header=F)
length(unique(bim[,1]))
noOFchr=unique(bim[,1])

#read the .bim file to get the number of chromosomes (if the data in PLINK txt format)
#map=read.table(paste(plinkFileName,".map",sep=""),header=F)
#length(unique(map[,1]))
#noOFchr=unique(map[,1])

for(i in max(noOFchr):1){
#for(i in 29:28){    
  #save the EAGLE phased data into .vcf format using .vcf files
  plinkCommand=paste("./plink --bfile ",plinkFileName," --chr-set 30 --recode vcf --chr ",i," --out tmp",sep="")
  system(plinkCommand)
  if(!file.exists("phased50k")){
    dir.create("phased50k")
  }
  eagleCommand=paste("./eagle --vcf=tmp.vcf --chromX=30 --geneticMapFile=genetic_map_1cMperMb.txt  --outPrefix=phased50k/phased_chr",i," --vcfOutFormat=v --numThreads=8",sep="")
  system(eagleCommand)
  
  print("Chromosome ")
  print(i)
  print("Done")
  print("\n\n\n\n")
  
}

###
#Now do the same thing with 700k to create phased 700k if the 700k data is in PLINK binary format 
#You will have phased 700k data into a folder "700kPhased"

#########################

#The following code segment is for phasing 700k data using PLINK and EAGLE if you have 700k unphased data in PLINK txt format

for(i in max(noOFchr):1){
#for(i in 28:29){
  
  
  if(!file.exists("700kPhased")){
    dir.create("700kPhased")
  }
  #plinkCommand=paste("./plink --bfile 700kUnphased/700k_chr",i," --cow --recode vcf --out 700kUnphased/tmp",sep="")
  #if plink files are in txt format
  plinkCommand=paste("./plink --file 700kUnphased/chr",i," --cow --recode vcf --out 700kUnphased/tmp",sep="")
  system(plinkCommand)
  
  eagleCommand=paste("./eagle --vcf=700kUnphased/tmp.vcf --chromX=30 --geneticMapFile=genetic_map_1cMperMb.txt  --outPrefix=700kPhased/700k_phased_chr",i," --vcfOutFormat=v --numThreads=8",sep="")
  system(eagleCommand)
  print(i)
}

#######################################################

#Now you have 50k phased data as target dataset and 700k phased data as reference data set
#Use Minimac3 to impute 50k phased data into 700k  
#This is only for autosomes and not for sex chromosomes. Minimac3 needs different options for sex chromosome imputation
#Minimac3 needs --refHaps i.e reference data .vcf format and ofcourse phased data
#               --haps  i.e. the data to be imputed .vcf format and ofcourse phased data
#               --rsid  this will save the output .vcf files with SNPs' rsID, if not supplied this option ID column of .vcf will be chr:position. This is helpfull if your sequence .vcf files contain lots of missing rs ids. In this case ingnore this --rsid option
#               --prefix  imputed file name prefix
#               --cpus  number of cupus for parallel execution
#The output will be 700k imputed data each chromosome as a seperate file
for(i in max(noOFchr):1){

  #Minimacomp_compiled (actual Minimac3 has 1-22 and X chromosome limitation so the source code of Minimac3 needs to be recompiled)
  # Otherwise you can change the first column of your vcf file to a number between 1 - 22 to work with minimac3 but I found this option is not very convenient 

  
  if(!file.exists("Imputed700k")){
    dir.create("Imputed700k")
  }

  miniMacCommand=paste("./minimac3 --refHaps 700kPhased/700k_phased_chr",i,".vcf --haps phased50k/phased_chr",i,".vcf  --rsid --prefix Imputed700k/imputed700k_chr",i," --cpus 8 --log",sep="")
  system(miniMacCommand)
  unzipCommand=paste("gunzip Imputed700k/imputed700k_chr",i,".dose.vcf.gz",sep="")
  system(unzipCommand)
  print(i)
}


## Merging imputed 700k and Reference 700k
## Merging can be done by followong codes
for(i in max(noOFchr):1){
    if(!file.exists("mergedImputed_N_Ref700k")){
        dir.create("mergedImputed_N_Ref700k")
    }
    system(paste("bgzip Imputed700k/imputed700k_chr",i,".dose.vcf",sep=""))
    system(paste("tabix -p vcf Imputed700k/imputed700k_chr",i,".dose.vcf.gz",sep=""))
    system(paste("bgzip 700kPhased/700k_phased_chr",i,".vcf",sep=""))
    system(paste("tabix -p vcf 700kPhased/700k_phased_chr",i,".vcf.gz",sep=""))
    system(paste("vcf-merge Imputed700k/imputed700k_chr",i,".dose.vcf.gz 700kPhased/700k_phased_chr",i,".vcf.gz >mergedImputed_N_Ref700k/merged700k_chr",i,".vcf",sep=""))
}

#####################################################
#Now you have all the samples imputed upto  700k
#And you are ready to impute them to seq level
#But you need to phase the reference sequence data first. 
#The following code is for making reference sequence each chromosome seperately from the whole genomoe Reference sequence
#And you need VCFTOOLS to do that
#Using VCFTOOLS Devide the Reference VCF files into each chromosome. If you have your reference sequences already splited into chromosome you do not need to do the following
for(i in max(noOFchr):1){
  VCFcommand=paste("vcftools --vcf CattleWholeGenomeRefSeq/WholeRefSeq.vcf --recode --chr ",i," --out CattleGenomeRefSeq/RefSeq_Chr",i,sep="")
  print(VCFcommand)
  system(VCFcommand)
}

#the output of the VCFTOOLS contains .recode.vcf extension
#But the output of VCFTOOLS are not phased
#You need to phase the data before the actual imputatioun

#phase each Chr RefSeq using EAGLE

for(i in max(noOFchr):1){
  if(!file.exists("PhasedCattleGenomeRefSeq")){
    dir.create("PhasedCattleGenomeRefSeq")
  }
  eagleCommand=paste("./eagle --vcf=CattleGenomeRefSeq/RefSeq_Chr",i,".recode.vcf --chromX=30 --geneticMapFile=genetic_map_1cMperMb.txt  --outPrefix=PhasedCattleGenomeRefSeq/RefSeq_phased_chr",i," --vcfOutFormat=v --numThreads=8",sep="")
  print(eagleCommand)
  system(eagleCommand)
  
}

#Now you have phased Reference sequence and phased 700k target sequence
#And you are good to go for imputation to sequence level

#imputation of merged 50k imputated and Ref 700k together 


for(i in max(noOFchr):1){
   if(!file.exists("mergedImputedSequence")){
    dir.create("mergedImputedSequence")
  }
  miniMacCommand=paste("./minimac3 --refHaps PhasedCattleGenomeRefSeq/RefSeq_phased_chr",i,".vcf --haps mergedImputed_N_Ref700k/merged700k_chr",i,".vcf --rsid --prefix mergedImputedSequence/merged_imputedTOseq_chr",i," --cpus 8 --log",sep="")
  system(miniMacCommand)
  print(i)
}




#imputation of Reference 700k to sequence 


for(i in max(noOFchr):1){
   if(!file.exists("700K_Imputed2Sequence")){
    dir.create("700K_Imputed2Sequence")
  }
  miniMacCommand=paste("./minimac3 --refHaps PhasedCattleGenomeRefSeq/RefSeq_phased_chr",i,".vcf --haps 700kPhased/700k_phased_chr",i,".vcf.gz --rsid --prefix 700K_Imputed2Sequence/700k_Imputed2Seq_chr",i," --cpus 8 --log",sep="")
  system(miniMacCommand)
  print(i)
}



#imputation of imputated 700k to sequence 


for(i in max(noOFchr):1){
   if(!file.exists("50K_Imputed2Sequence")){
    dir.create("50K_Imputed2Sequence")
  }
  miniMacCommand=paste("./minimac3 --refHaps PhasedCattleGenomeRefSeq/RefSeq_phased_chr",i,".vcf --haps Imputed700k/imputed700k_chr",i,".dose.vcf.gz --rsid --prefix 50K_Imputed2Sequence/50k_Imputed2Seq_chr",i," --cpus 8 --log",sep="")
  system(miniMacCommand)
  print(i)
}
