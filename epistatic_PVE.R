loci = read.table('sign_loci.txt',stringsAsFactors=F)
phe = read.table('6pop_phe_JLM.txt',header=T,row=-1,stringsAsFactors=F)
geo = read.table('Pop7_Bin_ParentId.csv',header=T,row=1,stringsAsFactors=F,)

inter = intersect(rownames(geo),rownames(phe))
geo = geo[inter,]
phe = phe[inter,]

res_d = data.frame()
for(i in 1:ncol(phe)){
  n = which(loci$V2 == colnames(phe)[i])
  tmp_sig = loci[n,]
  m = which(colnames(geo) %in% tmp_sig$V1)
  tmp_gen = geo[,m]
  tmp_phe = phe[,i]
  names(tmp_phe) = rownames(phe)
  intt = intersect(rownames(tmp_gen),names(tmp_phe))
  tmp_phe = tmp_phe[intt]
  tmp_gen = tmp_gen[intt,]
  tmp_d = data.frame(trait = tmp_phe)
  tmp_d = cbind(tmp_d,tmp_gen)
  
  cmb = t(combn(tmp_sig$V1,2))
  for(mkn in 1:nrow(cmb)){
    mk1 = cmb[mkn,1]
    mk2 = cmb[mkn,2]
    
    m1 = tmp_d[,mk1]
    m2 = tmp_d[,mk2]
    y = tmp_d$trait
    
    lm1 = lm(y~m1+m2)
    lm2 = lm(y~m1+m2+m1:m2)
    pvalue = anova(lm1,lm2)[2,6]
    
    fl = aov(y~m1+m2+m1:m2)
    aov.res=summary(fl)[[1]]
    total.var = sum(aov.res[,2])
    add.var = sum(aov.res[1,2]+aov.res[2,2])
    epi.var = aov.res[3,2]
    add.R2 = add.var/total.var
    epi.R2 = epi.var/total.var
    
    tmp_dd = data.frame(trait=trait[i],m1=mk1,m2=mk2,
                        pvalue=pvalue,add.R2=add.R2,epi.R2=epi.R2)
    res_d = rbind(res_d,tmp_dd)
  }
}
res_d

write.table(res_d,'epi.txt',row.names = F,quote=F,sep='\t')


write.table(d,'R2.txt',sep='\t',col.names = T,quote=F,row.names = F)
