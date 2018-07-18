
library(data.table)
library(dplyr)
library(magrittr)

count_kmers<-function(target, db,
                      adj = "fdr",
                      signLev=0.01,
                      target_preCalculated=F,
                      db_precalculated=F,
                      with.gaps=F,
                      groovy= "/raid/users/Artem/.sdkman/candidates/groovy/current/bin/groovy",
                      v_kmer="/raid/users/Artem/tools/vdj_kmer/kmerv_new.groovy",
                      cdr_mapping="/raid/users/Artem/tools/kmer/cdr1cdr2_1subst.txt") {
  

    
  # arguments
  # target - cloneset or precalculated repertoire k-mers set for which we want to find enriched k-mers
  # db - cloneset or precalculated repertoire k-mers set, serving as background for enrichment calculation
  # adj - multiple comparisons adjustments, possible values:  "holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"
  # threshold - adjusted Pvalue threshold
  # with.gaps - should we also investigate k-mer accounts with gaps 
  
  if (with.gaps){
    with.gaps<-"true"
  } else {
    with.gaps<-"false"
  }
  
    
 if (target_preCalculated){
    
    tg_kmers<-fread(target,skip=1)
    colnames(tg_kmers)<-c("kmer","v","count")
    len<- scan(file = target,n = 1)
    
  } else {

    len<-as.integer(system(paste0("cat ", target, " | wc -l"),intern=T))-1
    tg_kmers<-fread(input = paste(groovy,v_kmer,target, cdr_mapping, with.gaps,"false"),skip=1) # calculate kmers
    colnames(tg_kmers)<-c("kmer","v","count")
  }

  
  
if (db_precalculated){
      kmer_db <-fread(db,skip = 1)
      db_len<- scan(file = d,n = 1)
      db_name<- strsplit(tail(strsplit(d,split = "/")[[1]],1),
                         split = ".",fixed = T)[[1]][1]
      colnames(kmer_db)<-c("kmer","v","count")
      
    } else {
      
      db_len<-as.integer(system(paste0("cat ",db, " | wc -l"),intern=T))-1
      kmer_db<-fread(input = paste(groovy,v_kmer,db, cdr_mapping,with.gaps,"false"),skip=1) # calculate kmers
      colnames(kmer_db)<-c("kmer","v","count")
    }
    
tg_kmers<-merge(x = tg_kmers,y = kmer_db,by = c("kmer","v"),suffixes = c(".target",".db"), all.x =T)

tg_kmers %<>% mutate(norm.count.db=count.db/db_len, norm.count.target=count.target/len)%>% 
  mutate(log2fc=log2(norm.count.target/norm.count.db))


# tg_kmers$p.hyper<-phyper( q=tg_kmers$count.target,
#                                    m=tg_kmers$count.db,
#                                    n = db_len - tg_kmers$count.db,
#                                    k =len,
#                                    lower.tail = F)
#     
#     
# tg_kmers$p.hyper.adj<-p.adjust(tg_kmers$p.hyper,method=adj)


tg_kmers$p.binom<-pbinom( q=tg_kmers$count.target,
                    size=len,
                    prob=tg_kmers$count.db/db_len,
                    lower.tail = F)


tg_kmers$p.binom.adj<-p.adjust(tg_kmers$p.binom,method=adj)

tg_kmers %<>% mutate(`minus_log10p`=-log10(p.binom.adj))

  return(tg_kmers)
}



