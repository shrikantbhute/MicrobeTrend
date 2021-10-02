dietsim<-function(pseq, treat, diet= c(x), t_col)
{
  otu_tab<-otu_table(pseq)
  old_map<-data.frame(sample_data(pseq))
  tree <-phy_tree(pseq)
  tax<-tax_table(pseq)
  category<-old_map[,treat]
  old_map$category<-category
  Time <- old_map[,var_col]
  old_map$Time<-Time
  new_map<-sample_data(old_map)

  #Creating new phyloseq object
  n_pseq<- phyloseq(otu_tab,tax,new_map,tree)

  x<-as.character(unique(sample_data(n_pseq)$category))

  cat("Performing SIMPER on all",diet,"samples\n")

  d<-prune_samples(sample_data(n_pseq)$category==diet,n_pseq)

  simper.matrix<-data.frame(Matrix::t(otu_table(d)),check.names = FALSE) # check.names is important,else data.frame will add 'X' to every column head if they are not syntactically valid

  sim=with(data.frame(sample_data(d)),simper(simper.matrix,Time))
  fout = tidy(capture.output(summary(sim)))
  fout %>% map(write.table,file= paste(diet,"_simper_summary.csv"),append= T, sep=',')

  cat("SIMPER analysis on",diet,"data is completed. Results are written to the file Diet_simper_summary.csv\n")

  }
