###calculate the ASE##################################

###read coverage files and only polymorphic site

###nucleotied reverse
nuc_rev <- function(x){
a = as.vector(x)
    if(a=="A") {
    a = "T"
    } else if (a=="C") {
    a = "G"
    } else if (a=="G") {
    a = "C"
    } else if (a=="T") {
    a = "A"
    }
return(a)
}

###To perform ASE analysis, three types of files are necessary to obtain:
#1. genome snp file which shows snp between species, the genome snp skips highly divergent regions
#2. depth mapping to hybrid genome (for those that passed RNA Integrity test-please, refer to https://github.com/icalic/RNA-Integrity-)
#3. vcf mapping to hybrid genome, here showing the wrong assigned reads (for those that passed RNA Integrity test-please, refer to https://github.com/icalic/RNA-Integrity-)

###set the working directory

setwd("/home/icalic/Documents/ASE_analysis")

###load the snp information of orthlogs between ath and aly (ortholog alignment), data looks like 
#     gene          Algene                           tha lyr chr    Atpos      scaf         Alpos    Atdir Aldir
# [1,] "AT1G01010.1" "fgenesh1_pg.C_scaffold_1000119" "C" "A" "Chr1" "3769"	"scaffold_1" "538505" "+"   "-"  


load("/home/icalic/Documents/ASE_analysis/AthAly_orth_genome.snp.RData")


###get the depth of mapping files by "samtools depth ..."
dpfiles = list.files(pattern="dp.txt")

##Ath SNP cover, only focus on polymorphic site
chr = unique(genome.snp[,5])
for(k in 1:length(dpfiles))
{
snp.cover = cbind(genome.snp, NA)
for(i in 1:length(chr))
{
cmd = paste("grep", chr[i], dpfiles[k], "> temp.txt")
system(cmd)
data = read.table("temp.txt")
m = which(genome.snp[,5]==chr[i])
rpos = genome.snp[m,6]
pos = data[,2]
cover= data[,3]
n = match(rpos,pos)
snp.cover[m,11] = cover[n]
cat(i, "\n")
} 
snp.cover[is.na(snp.cover)] = 0
save(snp.cover, file=sub("dp.txt", "snpcover.RData",dpfiles[k]), compress=T)
cat(k, "\n")
}

##Aly cover for snp, only focus on polymorphic site
scfile = list.files(pattern="snpcover")
scaf = unique(genome.snp[,7])
for(k in 1:length(dpfiles))
{
load(scfile[k])
snp.cover = cbind(snp.cover, NA)
for(i in 1:length(scaf))
{
cmd = paste("grep", scaf[i], dpfiles[k], "> temp.txt")
system(cmd)
data = read.table("temp.txt")
m = which(genome.snp[,7]==scaf[i])
rpos = genome.snp[m,8]
pos = data[,2]
cover= data[,3]
n = match(rpos,pos)
snp.cover[m,12] = cover[n]
cat(i, "\n")
} 
snp.cover[is.na(snp.cover)] = 0
save(snp.cover, file=sub("dp.txt", "snpcover.RData",dpfiles[k]), compress=T)
cat(k, "\n")
}

###adjust coverage by wrong mapped reads basing on snp information
###if the read have a snp which can map to another species, it will be assigned to another species
###get the snp information by samtools, only focus on orthologous snp

varfiles = list.files(path="/home/icalic/Documents/ASE_analysis/vcf",full.names=T)

varfiles = list.files(pattern="vcf")
snpcoverfiles = list.files(pattern="snpcover")
Atpos = paste(genome.snp[,5], genome.snp[,6], sep=":")
Alpos = paste(genome.snp[,7], genome.snp[,8], sep=":")

for(k in 1:length(snpcoverfiles))
{
load(snpcoverfiles[k])
vcf = read.table(varfiles[k])
info = as.vector(vcf[,8])
a = unlist(strsplit(info, ";"))
dp4 = a[grep("DP4=", a)]
dp4 = sub("DP4=", "", dp4)
dp4 = matrix(as.numeric(unlist(strsplit(dp4, ","))), nc=4, byrow=T)
ad = cbind(dp4[,1]+dp4[,2], dp4[,3]+dp4[,4])
vcf1 = cbind(as.matrix(vcf[,4:5]),ad)
vpos = paste(as.vector(vcf[,1]), vcf[,2], sep=":")
vpos = sub("om_", "", vpos)

m = match(Atpos, vpos)
n = match(Alpos, vpos)
snp.cover.var = cbind(snp.cover, vcf1[m,], vcf1[n,])
rownames(snp.cover.var) = 1:nrow(snp.cover.var)

snp.cover.adjust = snp.cover.var
n = which(!is.na(snp.cover.var[,13]))
for(i in n)
{
if(snp.cover.var[i,9] == snp.cover.var[i,10])
	{
	if(snp.cover.var[i,3]==snp.cover.var[i,13] & snp.cover.var[i,4]==snp.cover.var[i,14])
		{
		snp.cover.adjust[i,12] = as.numeric(snp.cover.var[i,12])+as.numeric(snp.cover.var[i,16])
		snp.cover.adjust[i,11] = as.numeric(snp.cover.var[i,11]) - as.numeric(snp.cover.var[i,16])
		}
	}
if(snp.cover.var[i,9] != snp.cover.var[i,10])
	{
	if(snp.cover.var[i,3]==snp.cover.var[i,13] & snp.cover.var[i,4]==nuc_rev(snp.cover.var[i,14]))
		{
		snp.cover.adjust[i,12] = as.numeric(snp.cover.var[i,12])+as.numeric(snp.cover.var[i,16])
		snp.cover.adjust[i,11] = as.numeric(snp.cover.var[i,11]) - as.numeric(snp.cover.var[i,16])
		}
	}
}

n = which(!is.na(snp.cover.var[,17]))
for(i in n)
{
if(snp.cover.var[i,9] == snp.cover.var[i,10])
	{
	if(snp.cover.var[i,4]==snp.cover.var[i,17] & snp.cover.var[i,3]==snp.cover.var[i,18])
		{
		snp.cover.adjust[i,12] = as.numeric(snp.cover.var[i,12])-as.numeric(snp.cover.var[i,20])
		snp.cover.adjust[i,11] = as.numeric(snp.cover.var[i,11]) + as.numeric(snp.cover.var[i,20])
		}
	}
if(snp.cover.var[i,9] != snp.cover.var[i,10])
	{
	if(snp.cover.var[i,4]==snp.cover.var[i,17] & snp.cover.var[i,3]==nuc_rev(snp.cover.var[i,18]))
		{
		snp.cover.adjust[i,12] = as.numeric(snp.cover.var[i,12])-as.numeric(snp.cover.var[i,20])
		snp.cover.adjust[i,11] = as.numeric(snp.cover.var[i,11]) + as.numeric(snp.cover.var[i,20])
		}
	}
}
save(snp.cover.adjust, file=sub("snpcover","AD", snpcoverfiles[k]))
}

###make a table for all samples of snp coverage

load("/home/icalic/Documents/ASE_analysis/AthAly_orth_genome.snp.RData")
adfiles = list.files(pattern="AD")
snpcover.all = genome.snp
for(i in 1:length(adfiles))
{
load(adfiles[i])
snpcover.all = cbind(snpcover.all, snp.cover.adjust[,11:12])
}

name = sub(".AD.RData","",adfiles)
name = rep(name,each=2)
names = paste(name,c("At","Al"),sep="_")
colnames(snpcover.all) = c(colnames(genome.snp),names)
hybrid_dehydration_ts_allele = snpcover.all
save(hybrid_dehydration_ts_allele, file="hybrid_dehydration_ts_allele.RData")

###filter SNP with DNA samples 
###filter SNP with 3 steps in DNA samples: 
1. SNP away from intron > 50 bp 
2. low divergent regions between Ath and Aly 
3. Median of the allele ratio.

load("/home/icalic/Documents/ASE_analysis/AthAly_orth_genome_snpwithdistTposDs200.RData")
load("/home/icalic/Documents/ASE_analysis/Alyr.allele.DNA.RData")
hybrid.allele.dts = cbind(hybrid_dehydration_ts_allele, Alyr.allele.DNA)

##gene away intron greater than 50 bp  (we begin with 1.34 million ortholg SNPs but with filter this number goes to 948704 SNPs)

genome.snp<-as.data.frame(genome.snp)
n = which(as.numeric(genome.snp[,11])>50)
genome.snp = genome.snp[n,]
hybrid.allele.dts = hybrid.allele.dts[n,]

###snp coverage higher than 5 in DNA sample (the number of leftover SNPs goes down to 932212 SNPs)

hybrid.allele.dts<-as.data.frame(hybrid.allele.dts)
hybrid.allele.dts$cmdna1At<-as.numeric(hybrid.allele.dts$cmdna1At)
hybrid.allele.dts$cmdna1Al<-as.numeric(hybrid.allele.dts$cmdna1Al)
hybrid.allele.dts$cmdna2At<-as.numeric(hybrid.allele.dts$cmdna2At)
hybrid.allele.dts$cmdna2Al<-as.numeric(hybrid.allele.dts$cmdna2Al)

n = which(hybrid.allele.dts[,33]> 5 & hybrid.allele.dts[,34]>5 & hybrid.allele.dts[,35]> 5 & hybrid.allele.dts[,36]>5)
genome.snp = genome.snp[n,]
hybrid.allele.dts = hybrid.allele.dts[n,]

###snp with the density less than 10 SNPs 200 bp (351936 SNPs kept)

n = which(as.numeric(genome.snp[,15])<11)
genome.snp = genome.snp[n,]
hybrid.allele.dts = hybrid.allele.dts[n,]

###the median ratio coverage of all of samples

hybrid.allele.dts1<-matrix(as.numeric(as.matrix(hybrid.allele.dts[,11:36])),nr=nrow(hybrid.allele.dts))
rownames(hybrid.allele.dts1) = as.character(hybrid.allele.dts[,1])
head(hybrid.allele.dts1)
gene = unique(genome.snp[,1])
lowsnp.mr = matrix(nr=length(gene), nc=ncol(hybrid.allele.dts1))

for(i in 1:length(gene))
 {
     unit = hybrid.allele.dts1[rownames(hybrid.allele.dts1)==gene[i],]

     if(length(unit)==ncol(hybrid.allele.dts1)) lowsnp.mr[i,] = unit
     else {
         for(j in 1:(ncol(hybrid.allele.dts1)/2))
         {
             x = unit[,(2*j-1):(2*j)]
             ratio = (as.numeric(x[,2])+0.25)/(as.numeric(x[,1])+0.25)
             n = order(ratio)[trunc(length(ratio)/2)+1]
             lowsnp.mr[i, (2*j-1):(2*j)] = x[n,]
         }
     }
     if(i/100 == trunc(i/100))cat(i,'\n')
 }


#The outcome of the loop is a list of numbers from 100 until 17400 (if runs for a single stage for instance seedling)

rownames(lowsnp.mr) = gene
colnames(lowsnp.mr) = colnames(hybrid.allele.dts[,11:36])
hybrid.allele.lowsnp = lowsnp.mr
save(hybrid.allele.lowsnp, file='hybrid.allele.lowsnp.RData')

#plot allele ratio distribution for all sample stages

jpeg("hybrid.allele.dts.lowsnp.jpeg",width=1200,height=1200,pointsize=30)
b = hybrid.allele.lowsnp
plot(density(log2(b[,24]/b[,23])),lwd=2,main='',xlab='log2(A.lyrata/A.thaliana)',ylim=c(0,2.5),xlim=c(-4,4))
lines(density(log2(b[,26]/b[,25])),lwd=2)

color =c("red", "red", "red", "blue", "blue", "blue", "blue", "green", "green", "green", "green")

for(i in 1:11)
{
lines(density(log2((b[,2*i]+0.2)/(b[,2*i-1]+0.2))),lwd=2,col=color[i])
}
legend('topright',c('AthxAly-DNA','AthxAly-FloweringRNA', 'AthxAly-SeedlingRNA', 'AthxAly-VegetativeRNA'),col=c('black','red', 'blue', 'green'),lty=1,lwd=2,bty='n')
dev.off()

###scaling to 1M reads 
exp = cbind(hybrid.allele.lowsnp[,1]+hybrid.allele.lowsnp[,2])
for(i in 2:11)
{
exp = cbind(exp, hybrid.allele.lowsnp[,2*i-1]+hybrid.allele.lowsnp[,2*i])
}

tcover = colSums(exp)
sizefactor = tcover/1000000
sizefactor = c(rbind(sizefactor,sizefactor))
allele.norm = round(hybrid.allele.lowsnp/matrix(sizefactor, nr=nrow(hybrid.allele.lowsnp), nc=26,byrow=T))
hybrid.allele.lowsnp.ms = allele.norm+1
save(hybrid.allele.lowsnp.ms,file="hybrid.allele.lowsnp.ms.RData")

###plot scaled allele cover
jpeg("hybrid.allele.lowsnp.ms.coverage.jpeg",width=1200,height=1200,pointsize=30)
b = hybrid.allele.lowsnp.ms
plot(density(log2(b[,24]/b[,23])),lwd=2,main='',xlab='log2(A.lyrata/A.thaliana)',ylim=c(0,2.5),xlim=c(-4,4))
lines(density(log2(b[,26]/b[,25])),lwd=2)

color =c("red", "red", "red", "blue", "blue", "blue", "blue", "green", "green", "green", "green")

for(i in 1:11)
{
lines(density(log2((b[,2*i]+0.2)/(b[,2*i-1]+0.2))),lwd=2,col=color[i])
}
legend('topright',c('AthxAly-DNA','AthxAly-FloweringRNA', 'AthxAly-SeedlingRNA', 'AthxAly-VegetativeRNA'),col=c('black','red', 'blue', 'green'),lty=1,lwd=2,bty='n')
dev.off()

####ASE analysis

load("hybrid.allele.lowsnp.ms.RData")
length(colnames(hybrid.allele.lowsnp.ms))
readcount<-hybrid.allele.lowsnp.ms

#####RNA vs DNA samples

#rename columns named "seedling" with "0-basal"

colnames(readcount)[7]<-"0-basal-93326_At"
colnames(readcount)[8]<-"0-basal-93326_Al"
colnames(readcount)[9]<-"0-basal-93328_At"
colnames(readcount)[10]<-"0-basal-93328_Al"
colnames(readcount)[11]<-"0-basal-93330_At"
colnames(readcount)[12]<-"0-basal-93330_Al"
colnames(readcount)[13]<-"0-basal-93378_At"
colnames(readcount)[14]<-"0-basal-93378_Al"

treat <- c(("flowering", "flowering", "flowering", "0-basal","0-basal","0-basal", "0-basal","vegetative", "vegetative", "vegetative", "vegetative"), rep("1DNA",2))

n<- nrow(readcount_new)
pval<-matrix(nrow=n, ncol=4)
esti<-matrix(nrow=n, ncol=4)
std<-matrix(nrow=n, ncol=4)
dimnames(pval)<-list(c(1:n), c("intercept", "0-basal-dna", "0-basal-flowering", "0-basal-vegetative"))
colnames(pval)
for (i in 1:n){
  x<-readcount_new[i,seq(1,26,2)]
  y<-readcount_new[i,seq(2,26,2)]
  x<-as.vector(as.numeric(x))
  y<-as.vector(as.numeric(y))
  count<-cbind(x,y)
  mod<-glm(count~treat, family=quasibinomial)
  pval[i,]<- summary(mod)$coef[, 4] 
  esti[i,]<- summary(mod)$coef[, 1]
  std[i,]<- summary(mod)$coef[, 2]
}

#get the summary og glm model

summary(mod)

#####################################################
#Call:
#glm(formula = count ~ treat, family = quasibinomial)

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.3680  -0.3116  -0.1011   0.1643   0.6277  

#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      0.78016    0.14692   5.310 0.000487 ***
#treat1DNA       -0.72609    0.19799  -3.667 0.005176 ** 
#treatflowering   0.76029    0.29578   2.570 0.030164 *  
#treatvegetative -0.08701    0.20925  -0.416 0.687279    
---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for quasibinomial family taken to be 0.1628092)

#Null deviance: 7.1364  on 12  degrees of freedom
#Residual deviance: 1.4851  on  9  degrees of freedom
#AIC: NA

#Number of Fisher Scoring iterations: 4

#get the overview for pval

summary(pval)

#########################################################################
intercept        0-basal-dna      0-basal-flowering 0-basal-vegetative
 Min.   :0.00000   Min.   :0.00000   Min.   :0.00000   Min.   :0.0000    
 1st Qu.:0.00078   1st Qu.:0.02351   1st Qu.:0.03709   1st Qu.:0.1908    
 Median :0.04189   Median :0.17927   Median :0.25441   Median :0.4656    
 Mean   :0.26955   Mean   :0.29707   Mean   :0.37599   Mean   :0.4973    
 3rd Qu.:0.46531   3rd Qu.:0.52395   3rd Qu.:0.69045   3rd Qu.:0.8039    
 Max.   :1.00000   Max.   :1.00000   Max.   :1.00000   Max.   :1.0000    
 NA's   :143       NA's   :143       NA's   :143       NA's   :143


#get the structure for pval

head(pval)

#######################################################################
 intercept 0-basal-dna 0-basal-flowering 0-basal-vegetative
1 1.0000000000 0.890122824      1.000000e+00         0.78503698
2 0.0011789866 0.026989899      7.665843e-01         0.48691508
3 0.0009710315 0.002695769      4.199793e-01         0.91254364
4 0.0885823387 0.368400572      5.718339e-05         0.01373735
5 0.0050259614 0.230912598      1.215997e-01         0.72252571
6 0.0082802342 0.021454889      6.273456e-01         0.69644377

#get warnings


warnings()
######################################################################

Warning message:
In matrix(sizefactor, nr = nrow(hybrid.allele.lowsnp),  ... :
  data length [22] is not a sub-multiple or multiple of the number of rows [16782]

#pval=c()
#esti=c()
#std=c()

###Benjamini hochberg correction
pval.adj<-matrix(nrow=n, ncol=4)
row.names(pval.adj)<-rownames(readcount)
colnames (pval.adj)<-c("intercept", "0-basal-dna","0-basal-flowering", "0-basal-vegetative")

for (i in 1:4){
  pval.adj[,i]<-p.adjust(pval[,i], method="fdr")
}


#get the structure of pval.adjusted

head(pval.adj)

############################################################################################
  intercept 0-basal-dna 0-basal-flowering 0-basal-vegetative
AT1G01010.1 1.000000000  0.94800958       1.000000000          1.0000000
AT1G01020.1 0.004353564  0.10267145       0.973976519          0.9378811
AT1G01030.1 0.003721952  0.02127342       0.685439826          1.0000000
AT1G01040.1 0.156500481  0.55793658       0.002571552          0.4258838
AT1G01050.1 0.014199400  0.42031632       0.318868601          1.0000000
AT1G01060.2 0.021460252  0.08866231       0.868925611          1.0000000

#create histogram distribution for "0-basal-dna" 

hist(pval.adj[,2])

#create QQplot for distribution of pval and padj for "0-basal-dna"

qqplot(pval[,2],pval.adj[,2])

#create histogram distribution for "0-basal-flowering"

hist(pval.adj[,3])

#create QQplot for distribution of pval and padj for "0-basal-flowering"

qqplot(pval[,3],pval.adj[,3])

#create histogram distribution for "0-basal-vegetative"

hist(pval.adj[,4])

#create QQplot for distribution of pval and padj for "0-basal-vegetative"

qqplot(pval[,4],pval.adj[,4])

#save the histogram of pval seedling
jpeg("Histogram_pval_seedling", width=1200, height=1200, pointsize=30)
plot(hist(pval.adj[,2]))
dev.off()

#save the QQplot for seedling
jpeg("QQplot_pvalpadj_seedling.jpeg", width=1200, height=1200, pointsize=30)
qqplot(pval[,2],pval.adj[,2])
dev.off()

#save the histogram of pval for vegetative
jpeg("Histogram_pval_vegetative", width=1200, height=1200, pointsize=30)
plot(hist(pval.adj[,3]))
dev.off()

#save the QQplot for vegetative
jpeg("QQplot_pvalpadj_vegetative.jpeg", width=1200, height=1200, pointsize=30)
qqplot(pval[,3],pval.adj[,3])
dev.off()

###dna ratio for correction

###AthxAly samples, this is to be run if single stage (e.g., seedling alone)
dna.p = rowMeans(cbind(readcount[,10]/(readcount[,9]+readcount[,10]), readcount[,12]/(readcount[,11]+readcount[,12])))

###DNA-Seq columns indicated in the file when run across stages (seedling and vegetative)
dna.p = rowMeans(cbind(readcount[,18]/(readcount[,17]+readcount[,18]), readcount[,20]/(readcount[,19]+readcount[,20])))

###DNA-Seq columns indicated in the file when run across all three stages 
dna.p = rowMeans(cbind(readcount[,24]/(readcount[,23]+readcount[,24]), readcount[,26]/(readcount[,25]+readcount[,26])))

###AthxAly samples indicated in the file 

###colxAly-seedling
p = pval[,2]
q = pval.adj[,2]
sd = std[,2] 
estimate = -esti[,2]
cover = rowSums(readcount[,1:8])/4
Athcover = round(rowMeans(readcount[,seq(1,8,2)]))
Alycover = round(rowMeans(readcount[,seq(2,8,2)])/dna.p*(1-dna.p))
log2ratio = log2(Alycover/Athcover)
seedling.ase = cbind(cover,Athcover,Alycover,log2ratio,p,q, estimate, sd)
seedling.Aly = rownames(seedling.ase[seedling.ase[,6]<0.01 & seedling.ase[,7]> 0 & seedling.ase[,1] > 9,])
seedling.Ath = rownames(seedling.ase[seedling.ase[,6]<0.01 & seedling.ase[,7]< 0 & seedling.ase[,1] > 9,])

jpeg("seedling.jpeg",width=1200,height=1200,pointsize=30)
n = which(seedling.ase[,1] < 20)
plot(seedling.ase[-n,7],-log10(seedling.ase[-n,6]), xlab="log2(A.lyrata/A.thaliana)", ylab="-log10(q)",cex.lab=1.4,ylim=c(0,5))
abline(h=-log10(0.01))
abline(v=0)
dev.off()

plot(seedling.ase[-n,4],seedling.ase[-n,7],xlab='log2ratio',ylab='estimate')

###colxAly-vegetative
p = pval[,3]
q = pval.adj[,3]
sd = std[,3] 
estimate = -esti[,3]
cover = rowSums(readcount[,9:16])/4
Athcover = round(rowMeans(readcount[,seq(9,16,2)]))
Alycover = round(rowMeans(readcount[,seq(10,16,2)])/dna.p*(1-dna.p))
log2ratio = log2(Alycover/Athcover)
vegetative.ase = cbind(cover,Athcover,Alycover,log2ratio,p,q, estimate, sd)
vegetative.Aly = rownames(vegetative.ase[vegetative.ase[,6]<0.01 & vegetative.ase[,7]> 0 & vegetative.ase[,1] > 9,])
vegetative.Ath = rownames(vegetative.ase[vegetative.ase[,6]<0.01 & vegetative.ase[,7]< 0 & vegetative.ase[,1] > 9,])

jpeg("vegetative.jpeg",width=1200,height=1200,pointsize=30)
n = which(vegetative.ase[,1] < 20)
plot(vegetative.ase[-n,7],-log10(vegetative.ase[-n,6]), xlab="log2(A.lyrata/A.thaliana)", ylab="-log10(q)",cex.lab=1.4,ylim=c(0,15))
abline(h=-log10(0.01))
abline(v=0)
dev.off()
plot(vegetative.ase[-n,4],vegetative.ase[-n,7],xlab='log2ratio',ylab='estimate')


###colxAly-flowering vs seedling
p_floseed = pval[,3]
q_floseed = pval.adj[,3]
sd_floseed = std[,3] 
estimate_floseed = -esti[,3]
cover_flowering = rowSums(readcount_flowering[,1:6])/3
cover_seedling = rowSums(readcount_seedling[,1:8])/4
cover_flowering_seedling = (rowSums(readcount_flowering[,1:6])+rowSums(readcount_seedling[,1:8]))

readcount_flowering_seedling<-readcount[,1:14]
Athchange_fs = round(rowMeans(readcount_flowering_seedling[,seq(1,14,2)]))
Alychange_fs = round(rowMeans(readcount_flowering_seedling[,seq(2,14,2)])/dna.p*(1-dna.p))
log2ratio_fs = log2(Alychange_fs/Athchange_fs)

#log2ratio_flowering_versus_seedling = Alychange_floseed - Athchange_floseed
flowering_seedling_final.ase = cbind(cover_flowering_seedling, Athchange_fs, Alychange_fs,log2ratio_fs,p_floseed,q_floseed, estimate_floseed, sd_floseed)
flowering_seedling_final.ase.Aly = rownames(flowering_seedling_final.ase[flowering_seedling_final.ase[,6]<0.01 & flowering_seedling_final.ase[,7]> 0 & flowering_seedling_final.ase[,1] > 9,])
flowering_seedling_final.Ath = rownames(flowering_seedling_final.ase[flowering_seedling_final.ase[,6]<0.01 & flowering_seedling_final.ase[,7]< 0 & flowering_seedling_final.ase[,1] > 9,])


jpeg("Flower_X_Seed.jpeg",width=1200,height=1200,pointsize=30)
n = which(flowering_seedling.ase[,1] < 70)
plot(flowering_seedling.ase[-n,7],-log10(flowering_seedling.ase[-n,6]), xlab="log2(A.lyrata/A.thaliana)", ylab="-log10(q)",cex.lab=1.4,ylim=c(0,5))
abline(h=-log10(0.01))
abline(v=0)
dev.off()


###colxAly-vegetative vs seedling
p_vegseed = pval[,4]
q_vegseed = pval.adj[,4]
sd_vegseed = std[,4] 
estimate_vegseed = -esti[,4]
cover_vs_new = rowSums(readcount_new[,1:16])/8
Athcover_vs_new = round(rowMeans(readcount_new[,seq(1,16,2)]))
Alycover_vs_new = round(rowMeans(readcount_new[,seq(2,16,2)])/dna.p*(1-dna.p))
log2ratio_vs_new = log2(Alycover_vs_new/Athcover_vs_new)
seed_veg_test_final.ase = cbind(cover_vs_new,Athcover_vs_new,Alycover_vs_new,log2ratio_vs_new,p_vegseed,q_vegseed, estimate_vegseed, sd_vegseed)
seed_veg_test_final.ase.Aly = rownames(seed_veg_test_final.ase[seed_veg_test_final.ase[,6]<0.01 & seed_veg_test_final.ase[,7]> 0 & seed_veg_test_final.ase[,1] > 9,])
seed_veg_test_final.ase.Ath = rownames(seed_veg_test_final.ase[seed_veg_test_final.ase[,6]<0.01 & seed_veg_test_final.ase[,7]< 0 & seed_veg_test_final.ase[,1] > 9,])

jpeg("Veg_Y_Seed.jpeg",width=1200,height=1200,pointsize=30)
n = which(seed_veg_test_final.ase[,1] < 80)
plot(seed_veg_test_final.ase[-n,7],-log10(seed_veg_test_final.ase[-n,6]), xlab="log2(A.lyrata/A.thaliana)", ylab="-log10(q)",cex.lab=1.4,ylim=c(0,5))
abline(h=-log10(0.01))
abline(v=0)
dev.off()

######save data as list- i change names
###save flowering vs seedling ASE results

head(flowering_seedling_final.ase)

write.table(flowering_seedling_final.ase, "/home/icalic/Documents/ASE_analysis/Flowering_vs_SeedlingASE.txt")

###save vegetative vs seedling ASE results

head(seed_veg_test_final.ase)

write.table(seed_veg_test_final.ase, "/home/icalic/Documents/ASE_analysis/Vegetative_vs_Seedling.txt")








