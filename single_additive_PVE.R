loci = read.table('sign_loci.txt',stringsAsFactors=F)
phe = read.table('6pop_phe_JLM.txt',header=T,row=-1,stringsAsFactors=F)
geo = read.table('Pop7_Bin_ParentId.csv',header=T,row=1,stringsAsFactors=F,)

inter = intersect(rownames(geo),rownames(phe))
geo = geo[inter,]
phe = phe[inter,]

d = data.frame()
for(i in 1:nrow(loci)){
  n1 = which(colnames(geo) == loci[i,1])
  n2 = which(colnames(phe) == loci[i,2])
  
  rd = data.frame(g=geo[,n1],p=phe[,n2])
  
  l = lm(rd$p ~ rd$g)
  res = summary(l)
  r = res$adj.r.squared
  tmp_d = data.frame(loci=loci[i,1],phe=loci[i,2],r2=r)
  d = rbind(d,tmp_d)
}

write.table(d,'R2.txt',sep='\t',col.names = T,quote=F,row.names = F)
