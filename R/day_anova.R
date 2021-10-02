
day_anova<- function(pseq,tlevel,taxa, t_col, day,treat,...)
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
  
  
  
  
  tlevel<-match.arg(tlevel,c("Phylum","Class","Order","Family","Genus","Species"))
  pglom<-tax_glom(n_pseq,taxrank = tlevel)
  
  
  glomed.df<-psmelt(pglom)
  
  taxa==taxa
  
  if ("Phylum" %in% tlevel) filtereddf<- glomed.df %>% filter(Phylum %in% paste(taxa)) %>% droplevels()
  else if ("Class" %in% tlevel) filtereddf<- glomed.df %>% filter(Class %in% paste(taxa)) %>% droplevels()
  else if ("Order" %in% tlevel) filtereddf<- glomed.df %>% filter(Order %in% paste(taxa)) %>% droplevels()
  else if ("Family" %in% tlevel) filtereddf<- glomed.df %>% filter(Family %in% paste(taxa)) %>% droplevels()
  else if ("Genus" %in% tlevel) filtereddf<- glomed.df %>% filter(Genus %in% paste(taxa)) %>% droplevels()
  else filtereddf<- glomed.df %>% filter(Species %in% paste(taxa)) %>% droplevels()
  
  
  anovadf<- filter(filtereddf,category==day)
  
  test<-aov(Abundance~Treat,data= anovadf)
  
  results<- tidy(capture.output(TukeyHSD(test)))
  
  
  write.csv(results,paste(day,"_",taxa,"_ANOVA_summary.csv"),row.names=FALSE)
  
  
}