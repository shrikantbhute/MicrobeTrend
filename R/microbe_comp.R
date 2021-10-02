microbe_comp<-function(pseq, tlevel, treat, diet, cutoff, x_var, cols=NULL)
{
  
  # Filtering the diet data
  cat("Filtering the phyloseq object for",diet," diet samples\n")
  

  sampids<- as.character(get_variable(pseq,paste(treat)))==diet
  df<-prune_samples(sampids,pseq)
  
  # Aglomerating phyloseq object at the given taxonomic level
  cat("Aglomerating the phyloseq object at the",tlevel,"level\n")
  
  tlevel<-match.arg(tlevel,c("Phylum","Class","Order","Family","Genus","Species"))
  
  pglom<-tax_glom(df,taxrank = tlevel) %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% psmelt() %>% 
    filter(Abundance>cutoff) %>% arrange_(tlevel)
  
  
  Taxa<-pglom[,ncol(pglom)]
  
  plt = ggplot(pglom, aes(x = Sample, y = Abundance, fill = Taxa)) +
    facet_grid(paste(x_var), scales="free",space = "free")+
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     size=14,face = "bold", color = "black"),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold")) +
    theme(axis.text.x = element_text(color = "black",size=12,face = "bold"),
          axis.text.y = element_text(color = "black",size=12,face = "bold")) +
    ylab("Relative Abundance") +
    theme(legend.title = element_text(color = "black",size=14,face = "bold")) +
    theme(legend.text = element_text(color = "black",size=12,face = "bold")) 
  
  print (plt)
  
  
}