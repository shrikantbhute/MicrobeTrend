ttplot<-function(pseq, tlevel, taxa, x_var,treat, cols)
{
  # sample_data(pseq)$Day<-as.character(sample_data(pseq)$Day)
  
  
  # otu_tab<-otu_table(pseq)
  # old_map<-data.frame(sample_data(pseq))
  # tree <-phy_tree(pseq)
  # tax<-tax_table(pseq)
  # category<-old_map[,treat]
  # old_map$category<-category
  # new_map<-sample_data(old_map)
  # 
  # n_pseq<- phyloseq(otu_tab,tax,new_map,tree)
  
  
  tlevel<-match.arg(tlevel,c("Phylum","Class","Order","Family","Genus","Species"))
  
  pglom<-tax_glom(pseq,taxrank = tlevel)
  
  pglom<- transform_sample_counts(pglom, function(x) x/sum(x))
  
  glomed.df<-psmelt(pglom)
  
  taxa==taxa
  
  if ("Phylum" %in% tlevel) filtereddf<- glomed.df %>% filter(Phylum %in% paste(taxa)) %>% droplevels()
  else if ("Class" %in% tlevel) filtereddf<- glomed.df %>% filter(Class %in% paste(taxa)) %>% droplevels()
  else if ("Order" %in% tlevel) filtereddf<- glomed.df %>% filter(Order %in% paste(taxa)) %>% droplevels()
  else if ("Family" %in% tlevel) filtereddf<- glomed.df %>% filter(Family %in% paste(taxa)) %>% droplevels()
  else if ("Genus" %in% tlevel) filtereddf<- glomed.df %>% filter(Genus %in% paste(taxa)) %>% droplevels()
  else filtereddf<- glomed.df %>% filter(Species %in% paste(taxa)) %>% droplevels()
  
  
  
  
  # Grouping the table
  cat("Grouping the samples for the Day and the",tlevel,"level\n")
  
  
  x<-filtereddf %>%  group_by_(x_var, treat)
  
  avgx <-x%>%summarise(Taxa=mean(Abundance))
  # avgx$Taxa<-log(avgx$Taxa) # log transforming Clostridioides count
  names(avgx)[3]<-paste(taxa)
  
  ttp<-ggplot(avgx,aes_string(x=paste(x_var),y=taxa,color=paste (treat),group=paste (treat)))+
    geom_line(size=1)+geom_point(size=4,alpha=0.6) +
    scale_color_manual(values=cols) + theme_minimal() +
    # scale_x_discrete(limits=c("0","3","10","13","16","17","18","19","20","21","22","30","47")) +
    theme(axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.text.x = element_text(color = "black",size=12,face = "bold"),
          axis.text.y = element_text(color = "black",size=12,face = "bold")) +
    theme(legend.title = element_text(color = "black",size=14,face = "bold")) +
    theme(legend.text = element_text(color = "black",size=12,face = "bold"))
  ttp 
  
  
}