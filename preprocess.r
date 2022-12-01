library(data.table)
library(R.utils)
library(bedr)

unknown=sample[which(sample[,5]=="unknown"),]
sample=sample[which(sample[,5]!="unknown"),]
sample[which(sample[,5]=="embryonic"),5]="adult"

state=list()
for (group in unique(sample$group)){

bed=c()
for(id in unique(sample[which(sample$group==group),2])){
file=tmp[grep(id,tmp)]
file=fread(file,data.table=F)
file=file[which(file[,4] %in% c("TssA","TssFlnk","TssFlnkU","TssFlnkD","Tx","EnhG1","EnhG2","EnhA1","EnhA2")),]
bed=rbind(bed,file)
}
bed=bed[,1:3]
fwrite(bed,"int.bed",collapse="")
system("bedtools merge -i int.bed -c 1 -o count > int1.bed")
bed=fread("int.bed",data.table=F)
bed=bed[which(bed[,4]>ceiling(length(unique(sample[which(sample$group==group),2]))/2)),]
state[[group]]=bed[,1:3]
}
save.image()

fwrite(state,"state.bed",col.names=F,sep="\t")



CELL=fread("col.txt",data.table=F,head=F)
mat=fread("matrix.tsv.gz",data.table=F)
mat=mat[,1:2]
colnames(mat)=c()
CRE=fread("cCREs.bed.gz",data.table=F)

ATAC=list()
for(i in 1:222){
cell=CELL[i,1]
id=mat[which(mat[,2]==i),1]
bed=CRE[id,]
bed=bed[which(bed[,1] %in% paste("chr",1:22,sep="")),]
bed=bedr.merge.region(bed)
ATAC[[cell]]=bed
}

short_arm - short arm gaps (count: 5; size range: 5,000,000 - 16,990,000 bases)
heterochromatin - heterochromatin gaps (count: 11; size range: 20,000 - 30,000,000 bases)
telomere - telomere gaps (count: 48; all of size 10,000 bases)
contig - gaps between contigs in scaffolds (count: 285; size range: 100 - 400,000 bases)
scaffold - gaps between scaffolds in chromosome assemblies (count: 470; size range: 10 - 624,000 bases)

dat=fread("H_C.bed",data.table=F)
cen=fread("centromere.txt",data.table=F)
cen=bedr.sort.region(cen[,1:3])

gap=fread("gap.txt",data.table=F)
gap=bedr.sort.region(gap[,1:3])

HSD=fread("HSD.csv",data.table=F)
HSD=bedr.sort.region(HSD[1:218,1:3])

gap=rbind(rbind(cen,gap),HSD)
gap=bedr.merge.region(gap)
a.sub1 <- bedr.subtract.region(dat[,1:3],gap)

for(i in names(clean)){
int=dat[[i]]
int=int[which(int[,3] %in% clean[[i]][,3]),]
clean[[i]]=data.frame(clean[[i]],int[,4])
}
save(clean,file="clean.RData")

for i in {9..22};
do
./table_annovar.pl ~/HCG/chr$i.vcf humandb/ -buildver hg38 -out chr$i -remove -protocol refGene -operation g -nastring . -vcfinput -polish -xref example/gene_xref.txt
echo $i
done


total=2657302487
load("/lustre/home/acct-bmelgn/bmelgn-3/epimap/state/.RData")
dat=fread("clean.bed",data.table=F)
l=names(state)
stat=c()
mat=data.frame()
N=list()

for(i in l){
activate=state[[i]]
fwrite(activate,"int/activate.bed",sep="\t",col.names=F)
system("bedtools subtract -a int/activate.bed -b gap.txt > int/int.bed")
activate=fread("int/int.bed",data.table=F)
n=sum(activate[,3]-activate[,2])
N[[i]]=n
system("bedtools subtract -a clean.bed -b int/int.bed > int/int1.bed")
system("bedtools subtract -a clean.bed -b int/int1.bed > int/int2.bed")
system("sed -i "s/$/\t1/" int/int2.bed")
system("sed -i "s/$/\t0/" int/int1.bed")
system("cat int/int1.bed int/int2.bed > int/int.bed")
system("bedtools sort -i int/int.bed > int/activate.bed")

activate=fread("int/activate.bed",data.table=F)
rowstat=c(n,34537810-nrow(activate))
stat=rbind(stat,rowstat)
col=activate[,6]
col=as.integer(col)
mat=cbind(mat,col)
}
fwrite(mat,"state.txt",sep="\t",col.names=F)

srun -p debug -n 1 --pty /bin/bash
conda activate deepsea
module load cuda/10.0.130-gcc-4.8.5

bedtools intersect -a clean.bed -b asmc.bed -wa -wb >int.bed

sh run_pipeline.sh test.vcf hg38 res/

join -e "0" -a1 -a2 -o 0 1.2 2.2 1.3 1.4\
   <(awk '{print $1"-"$2"-"$3"-"$4" "$6" "$8" "$9}' gnomad.txt | sort -k 1) \
   <(tail -n +2 topmed.txt | awk '{print $1"-"$2"-"$3"-"$4" "$5}' - | sort -k 1) \
 | sed -e 's/-/ /g' > out.txt
 
 for i in {20..22};
 do
 awk '{print $1 " "  $2 " . " $3 " " $4 " . PASS " $5 "|" $6 "|" $7 "|" $8 "|" }' chr$i.all | gzip > chr$i.vcf.gz
 done
 
type=c("AG","AL","DG","DL")
f=function(x){
	prob=max(x[2:5])
	if(prob<0.4){out=c(0,0,0,0)}else{
		index=which.max(x[2:5])
 		t=type[index]
 		pos=x[index+5]
 		out=c(x[1],prob,t,pos)
 	}
 	return(out)
}

d=fread("spliceAI.vcf",data.table=F)
d=d[,3:17]
d[,c2:9]=apply(d[,2:9],2,as.numeric)
d$out=apply(d[,1:9],1,f)
d$count=ifelse(d$out==0,0,1)
l=unique(d[,1])
d=split(d,f=d[,1])
p=0.0007
res=c()
for(gene in l){
int=d[[gene]]
x=length(which(int$count==1))
n=nrow(int)
model=binom.test(x=x,n=n,p=p)
res=rbind(res,c(gene=gene,p=model$p.value,count=x,all=n))
}

library(data.table)
library(R.utils)

for (i in 1:11){
hcg=fread(paste(c("chr",i,".vcf"),collapse=""),data.table=F)
hcg=hcg[,1:5]
all=fread(paste(c("../gnomad/chr",i,".vcf"),collapse=""),data.table=F,head=F)
l=match(hcg[,2],all[,2])
all=all[,1:9]
all[is.na(all)]=0
all[,8:9]=apply(all[,8:9],2,as.numeric)
hcg$overlap=all[l,8]+all[l,9]
L=ifelse(hcg[,5]==all[l,5],l,NA)
hcg$reverse=all[L,8]+all[L,9]
d=data.frame(hcg,gid=l)
fwrite(d,paste(c("chr",i,".vcf"),collapse=""),sep="\t")
}

am=function(x){
d=abs(x)
val=max(d)
index=ifelse(val>0.272,name[which.max(d)],NA)
return(c(val,index))
}

thred=function(x){
ifelse(x[3:length(x)]>x[2],1,ifelse(x[3:length(x)] < x[1],-1,0))
}

labels=function(x){
up=paste(name[which(x==1)],collapse="|")
down=paste(name[which(x==-1)],collapse="|")
}
name=sub(".hg19","",name)
name=sub("CLIP","",name)
name=sub("rep","",name)

cat "/lustre/home/acct-bmelgn/bmelgn-3/deepsea/models/seqweaver/data.hg19.colnames" | while read name; do name="${name/hg19/}"; name="${name/rep/}"; name="${name/CLIP/}"; num=$(grep -o -i $name all.seq | wc -l);  echo ${name}_${num} >> count_217.txt ; done
#19985388




for(i in 0:12){
setwd(as.character(i))
tmp=list.files(pattern="annot*")
int=fread(tmp,data.table=F)
setwd("seqweaver_human")
tmp=list.files(pattern="*txt")
row=fread(tmp,data.table=F)
l=!duplicated(row[,3])
row=row[l,1:3]
row$gene=int[match(row[,3],int[,3]),7]
data=h5read("annot_diffs.h5","data")
data=data[,l]
row=row[,3:4]
L=apply(data,2,am)
L=t(L)
row$seq_max=L[,1]
row$max_var=L[,2]
data=cbind(bound[-1,],data)
data=apply(data,1,thred)
data=apply(data,1,labels)
row$up=data[,1]
row$down=data[,2]
fwrite(row,paste(c("~/deepsea/",i,".seq"),collapse=""),sep="\t")
rm(data)
rm(row)
gc()
setwd("~/deepsea")
}

GC=GC[which(GC[,1]>0),]
d=data.frame(count=res$count,all=res$all,length=GC[match(res[,4],GC$name),1],GC=GC[match(res[,4],GC$name),2],name=GC[match(res[,4],GC$name),3],id=GC[match(res[,4],GC$name),4])
d[,1]=d[,1]+1
d=na.omit(d)
model=glm(count~all+GC+length,family = poisson(link = "log"), data = d)
pred=predict(model,d[,2:4],type="link", se.fit=TRUE)
d$pred=pred$fit
d$se=pred$se
d$fold=exp(log(d[,1])-d$pred)
d[,1]=d[,1]-1
d=d[order(d$fold,decreasing=T),]

data:  table(tran[, 3:4])
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.841721 0.878883
sample estimates:
odds ratio
 0.8601072

bgrbp=ASMC[which(ASMC$rbp==1 & ASMC[,2]<0.3763),1]
47989264/190457853

NM=c()
for(chr in 1:22){
int=rbp[which(rbp[,1]==paste("chr",chr,sep="")),]
gnomad=fread(paste(c("../gnomad/chr",chr,".vcf"),collapse=""),data.table=F)
gnomad[,8:9]=apply(gnomad[,8:9],2,as.numeric)
nm=c()
for(i in 1:nrow(int)){
l=which(gnomad[,2]>int[i,3]-100 & gnomad[,2]<int[i,3]+100)
nm=rbind(nm,c(int[i,6],length(l),length(which(gnomad[l,8]+gnomad[l,9]==1))))
}
NM=rbind(NM,nm)
rm(gnomad)
gc()
}
save(NM,file="0117_BGatRBP.RData")

sei=sei[-which(sei[,4] %in% c(0,1,4,8,11,15,18,20,22:24,29,32:35,39)),]
#0117: Rscript test.r
#null distribution not yet generated

sei=fread("int.bed",data.table=F)
sei$type=ifelse(sei[,3]-sei[,2]>4096,NA,1)
sei=split(sei,f=sei[,1])
newbed=c()
for(chr in 2:22){
int=sei[[paste("chr",chr,sep="")]]
removed=c()
for(index in 1:nrow(int)){
if(is.na(int[index,4])){next}
target=index:(index+3)
low=int[index,2]
up=low+4048
covered=target[which(int[target,3]<up)]
if(length(covered)==0){covered=index}
removed=c(removed,covered)
int[covered,4]=NA
coverlength=int[max(covered),3]-low
flank=floor((4096-coverlength)/2)
low=int[index,2]-flank
bed=c(int[index,1],low,low+4096)
newbed=rbind(newbed,bed)
}
int=int[-removed,]
colnames(int)=c("chr","start","end")
newbed=data.frame(newbed)
colnames(newbed)=c("chr","start","end")
newbed=rbind(newbed,int[,1:3])
}
save(newbed,file="newbed.RData")

#0117: Rscript test.r
library(seqinr)
library(data.table)
command=paste("samtools faidx ~/deepsea/resources/hg38_UCSC.fa",range, sep=" ")
system("module load samtools")
nt=c("A","C","T","G")
nvar=fread("RBPnvar.txt",data.table=F)
fasta=read.fasta("ref.fa")
saturateList=c()
for(i in 1:nrow(nvar)){
id=nvar[i,1]
range=paste(c(nvar[i,10],":",nvar[i,11]-100,"-",nvar[i,11]+100),collapse="")
ref=as.vector(toupper(fasta[[range]]))
for(j in 0:200){
row=data.frame(chr=nvar[i,10],pos=nvar[i,11]-100+j,cid=id,ref=ref[j+1],alt=setdiff(nt,ref[j+1]))
saturateList=rbind(saturateList,row)
}
}
fwrite(saturateList,"../deepsea/RBPsat.vcf",sep="\t",col.names=F)

for(i in 2:17){
h5file=paste(c(i,"/seqweaver_human/annot.",i,"_diffs.h5"),collapse="")
rowfile=paste(c(i,"/seqweaver_human/annot.",i,"_row_labels.txt"),collapse="")

d=h5read(h5file,"data")
d=t(d)
rowlabel=fread(rowfile,data.table=F)
rbp=rbind(rbp,d)
ROW=rbind(ROW,rowlabel)
}

row$dot=NA
for(ID in id){
l=which(row[,3]==ID)
int=rbp[l,]
INT=highlight[which(var[,3]==ID)[1],]
dot=apply(int,1,function(x){sum(x*INT)})
row[l,"dot"]=dot
}

for(i in name){
down=fread(paste(i,"down",sep="."),data.table=F)
up=fread(paste(i,"up",sep="."),data.table=F)
all=rbind(down,up)
if(nrow(all)==0){next}
l=table(all[,2])
l=l[D$id]

d=data.frame(D[,c("all","length","GC")],count=as.vector(l))
d[which(is.na(d[,4])),4]=0
d[,4]=d[,4]+1
d=na.omit(d)
model=glm(count~all+GC+length,family = poisson(link = "log"), data = d)
pred=predict(model,d[,1:3],type="link", se.fit=TRUE)
d$pred=pred$fit
d$se=pred$se
fold=exp(log(d[,4])-d$pred)
D=data.frame(D,int=fold)
colnames(D)[ncol(D)]=i
}

res_fisher=c()
for(i in name){
int=data.frame(a=ifelse(D[,i]>2.5,1,0),b=D$loeuf)
int=table(int)
model=fisher.test(int)

or=model$estimate
low=model$conf.int[1]
up=model$conf.int[2]
p=model$p.value
res_fisher=rbind(res_fisher,c(i,"LOEUF",p,low,or,up))

int=data.frame(a=ifelse(D[,i]>2.5,1,0),b=D$SCG)
int=table(int)
model=fisher.test(int)

or=model$estimate
low=model$conf.int[1]
up=model$conf.int[2]
p=model$p.value
res_fisher=rbind(res_fisher,c(i,"SCG",p,low,or,up))
}

for chr in {1..22};do
awk -v s=1 'BEGIN{OFS="\t"} {print $1,$2-s,$2,$3,$6}' $chr/seqweaver_human/annot.int_row_labels.txt | tail -n +2 > chr.bed
bedtools intersect -a block.bed -b chr.bed -wa -wb | awk '{if($10 == $4) print $9,$5,sqrt(($3-$7-500)^2)}' > snp2seq/$chr.txt
done

cat *txt | sort -k 3n | sort -k 1,1 -u > snp2seq

library(rhdf5)
library(data.table)
load("name.RData")
name=name[-1]
ha=fread("ancestor/hg_anc_diff.txt.gz",data.table=F)
ha=round(ha,2)
rha=fread("ancestor/row.txt.gz",data.table=F)
rha[which(is.na(rha[,5])),5]=0
s2b=fread("snp2seq/snp2seq",data.table=F)
for(i in 1:21){
	file=paste(i,"/seqweaver_human/annot.int_diffs.h5",sep="")
	d=h5read(file,"data")
	d=t(d)
	d=round(d,2)
	file=paste(i,"/seqweaver_human/annot.int_row_labels.txt",sep="")
	row=fread(file,data.table=F)	
	bim=fread(paste(c("bim/1000G.EUR.QC.",i,".bim"),collapse=""),data.table=F)
	l=match(bim[,2],row[,3])
	bim=data.frame(snp=bim[,2],A1=ifelse(is.na(l),bim[,5],row[l,5]),A2=ifelse(is.na(l),bim[,6],row[l,4]))
	d=d[l,]
	d[is.na(d)]=0
	l=match(bim[,1],s2b[,1])
	l=match(s2b[l,3],rha[,5])
	d1=ha[l,]
	d1[is.na(d1)]=0
	d=d * d1
	colnames(d)=name
	d$all=rowSums(d)
	d=d*1000
	
	bim$x=0
	for(j in colnames(d)){
		bim[,4]=d[,j]
		colnames(bim)=c("SNP","A1","A2",j)
		fwrite(bim,paste(c("seqannot/",j,".",i,".sannot.gz"),collapse=""),sep="\t",scipen=10)
	}
}
