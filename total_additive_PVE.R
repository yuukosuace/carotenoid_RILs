loci = read.table('sign_loci.txt',stringsAsFactors=F)
phe = read.table('6pop_phe_JLM.txt',header=T,row=-1,stringsAsFactors=F)
geo = read.table('Pop7_Bin_ParentId.csv',header=T,row=1,stringsAsFactors=F,)

inter = intersect(rownames(geo),rownames(phe))
geo = geo[inter,]
phe = phe[inter,]

d = data.frame()

for(i in 1:ncol(phe)){
  n1 = which(loci[,2] == colnames(phe)[i])
  n2 = which(loci[n1,1] %in% colnames(geo))
  
  tmp_geo = geo[,n2]
  
  dd = data.frame(y=phe[,i],gep[,n2])
  
  l = lm(y~.,dd)
  res = summary(l)
  r = res$adj.r.squared
  tmp_d = data.frame(phe=colnames(phe)[i],r2=r)
  d = rbind(d,tmp_d)
}

write.table(d,'R2.txt',sep='\t',col.names = T,quote=F,row.names = F)
