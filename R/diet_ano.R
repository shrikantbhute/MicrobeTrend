diet_ano<-function(pseq, treat, t_col, day=c(y),...)
{
  
  otu_tab<-otu_table(pseq)
  old_map<-data.frame(sample_data(pseq))
  tree <-phy_tree(pseq)
  tax<-tax_table(pseq)
  category<-old_map[,t_col]
  old_map$category<-category
  
  Treat<-old_map[,treat]
  old_map$Treat<-Treat
  
  new_map<-sample_data(old_map)
  
  #Creating new phyloseq object
  n_pseq<- phyloseq(otu_tab,tax,new_map,tree)
  
  
  y<-as.list(seq(0:200))
  
  cat("Performing ANOSIM on all day",day,"samples\n")
  
  
  d<-prune_samples(sample_data(n_pseq)$category==day,n_pseq)
  
  ano.matrix<-data.frame(Matrix::t(otu_table(d)),check.names = FALSE) # check.names is important,else data.frame will add 'X' to every column head if they are not syntactically valid
  
  ano = with(data.frame(sample_data(d)), anosim(ano.matrix,Treat))
  
  fout<-capture.output(summary(ano))
  
  write.csv(fout,paste(day,"_anosim_summary.csv"),row.names=FALSE)
  
  cat("ANOSIM analysis on",day,"data is completed. Results are written to the file",day,"_anosim_summary.csv\n")
  
}