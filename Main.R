x=c(6,4,3,8,1,3,10,5)
length(x)
unique(x)
which(x==3)
rev(x)
sort(X)
sum(x)
mean(x)
median(x)
quantile(x)
summary(x)
paste0('Ali',' ','madad')
paste0(seq(1,3),' mississippi')
x[1]
x[c(1,2,3)]
vec=c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE)
x[vec]
x<8
x[x<8]
x %in% c(2,3,5)
x[x %in% c(2,3,5)]
?length
chr_vec=c("chr3","chrX","chr2")
is.character(chr_vec)
sort(chr_vec)
Ordered_vec=paste0("chr",c(1:22,'X','Y','M'))
sort(chr_vec)
chr_fac=factor(chr_vec,
               levels=Ordered_vec,
               ordered = TRUE)
sort(chr_fac)
my_First_Matrix=matrix(1:12,nrow=4,ncol=3)
matrix(1,nrow = 10,ncol=10)
dim(my_First_Matrix)
as.vector(my_First_Matrix)

gene=list(name='GAPDH',
          Protein_coding=TRUE,
          chromosome=13)
gene$name
gene[[1]]
gene[['name']]

my_First_Matrix[3:4,1:2]

df=data.frame(gene_symbol=c("G1","G2","G3"))
df

#This will be a code section ----
install.packages("BiocManager")
install.packages('Biobase')
install('Biobase') 
install(c('multtest','tidyverse'))


library(BiocManager)

library(multtest)
library(ggplot2)

data(golub)


ccnd3_exp=golub[1042,]



golub_factor=factor(golub.cl,levels = 0:1,labels=c("ALL","AML"))
df=data.frame("Exp"=ccnd3_exp,
              "Class"=golub_factor)
g=ggplot(data=df,aes(x=Class,y=Exp)) 
g +  geom_violin(aes(fill=golub_factor))
g +   geom_jitter(width=0.2)
g +geom_boxplot()+geom_jitter(width=0.2)
g+geom_boxplot()+geom_point()

ggplot(data=df,aes(x=Class,y=Exp)) +
  geom_boxplot(aes(fill=Class)) +
  geom_jitter()+
  theme_bw()+
  ggtitle("Some Title")


golub.gnames[c(829,1042),]
df$CST3_exp=golub[828,]
df

## This wont plot well in gg; its not standard --> Tidying-->>

golub.gnames[c(829,1042),]
gene_expression=golub[c(829,1042),]
gene_expression=as.vector(t(gene_expression))
rep(1:5,times=3)
rep(1:5,each=3)
Cancer_Type=rep(golub_factor,times=2)
gene_label=rep(c('CST3','CCND3'),each=38)

df2=data.frame(gene_expression,
               gene_label,
               Cancer_Type)
df2
ggplot(df2,aes(x=Cancer_Type,y=gene_expression))+
  geom_boxplot(aes(color=gene_label))
  
ggplot(df2,aes(x=Cancer_Type,y=gene_expression,colour=Cancer_Type))+
  geom_violin(aes(fill=Cancer_Type),alpha=0.1)+
  geom_jitter(aes(color=Cancer_Type))+
  facet_grid(cols=vars(gene_label))+
  theme_bw()

getwd()



library(readr)
library(dplyr)
library(magrittr)
library(tximport)

Sample_Table=read_csv("SraRunTable.txt") %>% select('Sample Name',source_name,treatment,
                                                    Cell_Line,Cell_type, time_point) %>%
  slice(seq(1,48,by=4))
File_Names=paste0(pull(Sample_Table,'Sample Name'),'/quant.sf')
names(File_Names)=pull(Sample_Table,'Sample Name')
gene_map=read_csv("TG.txt",col_names =c("enstid","enstid") )

Count_Data=tximport(files=File_Names,
         type = "salmon",
         tx2gene=gene_map ,
         ignoreTxVersion = TRUE)
Count_Data['counts']
