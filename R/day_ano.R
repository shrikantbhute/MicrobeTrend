day_ano<-function(pseq, time_col, treat, diet=c(x),...)
{
  
  otu_tab<-otu_table(pseq)
  old_map<-data.frame(sample_data(pseq))
  tree <-phy_tree(pseq)
  tax<-tax_table(pseq)
  category<-old_map[,treat]
  old_map$category<-category
  
  Time<-old_map[,time_col]
  old_map$Time<-Time
  
  new_map<-sample_data(old_map)
  
  #Creating new phyloseq object
  n_pseq<- phyloseq(otu_tab,tax,new_map,tree)
  
  
  
  cat("Performing ANOSIM on all",diet,"samples\n")
  
  
  d<-prune_samples(sample_data(n_pseq)$category==diet,n_pseq)
  
  ano.matrix<-data.frame(Matrix::t(otu_table(d)),check.names = FALSE) # check.names is important,else data.frame will add 'X' to every column head if they are not syntactically valid
  
  ano = with(data.frame(sample_data(d)), anosim(ano.matrix,Time))
  
  fout<-capture.output(summary(ano))
  
  write.csv(fout,paste(diet,"_anosim_summary.csv"),row.names=FALSE)
  
  cat("ANOSIM analysis on",diet,"data is completed. Results are written to the file",diet,"_anosim_summary.csv\n")
  
}