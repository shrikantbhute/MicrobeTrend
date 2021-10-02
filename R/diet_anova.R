
diet_anova<- function(pseq,tlevel,taxa, treat, diet,t_col,...)
{
  
  
  otu_tab<-otu_table(pseq)
  old_map<-data.frame(sample_data(pseq))
  tree <-phy_tree(pseq)
  tax<-tax_table(pseq)
  category<-old_map[,treat]
  old_map$category<-category
  
  Time<-old_map[,t_col]
  old_map$Time<-Time
  
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
  
  
  anovadf<- filter(filtereddf,category==diet)
  anovadf$Time<-as.factor(anovadf$Time)
  test<-aov(Abundance~Time,data= anovadf)
  
  results<- tidy(capture.output(TukeyHSD(test)))
  
  
  write.csv(results,paste(diet,"_",taxa,"_ANOVA_summary.csv"),row.names=FALSE)
  
  
}