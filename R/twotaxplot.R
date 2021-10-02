twotaxplot<-function(pseq,treat, diet,tlevel, taxa1, taxa2, x_var, cols,...)
{

  # sample_data(pseq)$Day<-factor(sample_data(pseq)$Day)


  otu_tab<-otu_table(pseq)
  old_map<-data.frame(sample_data(pseq))
  tree <-phy_tree(pseq)
  tax<-tax_table(pseq)
  category<-old_map[,treat]
  old_map$category<-category
  new_map<-sample_data(old_map)
  
  #Creating new phyloseq object
  n_pseq<- phyloseq(otu_tab,tax,new_map,tree)




  # Filtering the diet data
  cat("Filtering the phyloseq object for",diet," diet samples\n")

  x<-as.character(unique(sample_data(n_pseq)$category))
  df<-prune_samples(sample_data(n_pseq)$category==diet,n_pseq)



  # Aglomerating phyloseq object at the given taxonomic level
  cat("Aglomerating the phyloseq object at the",tlevel,"level\n")

  tlevel<-match.arg(tlevel,c("Phylum","Class","Order","Family","Genus","Species"))
  pglom<-tax_glom(df,taxrank = tlevel)
  
  pglom<- transform_sample_counts(pglom, function(x) x/sum(x))
  glomed.df<-psmelt(pglom)
  taxatofilter <-c(taxa1,taxa2)

  # Filtering the Otus for the taxa provided
  cat("Filtering",taxa1,"and", taxa2, "from the table\n")

  if ("Phylum" %in% tlevel) filtereddf<- glomed.df %>% filter(Phylum %in% taxatofilter) %>% droplevels()
  else if ("Class" %in% tlevel) filtereddf<- glomed.df %>% filter(Class %in% taxatofilter) %>% droplevels()
  else if ("Order" %in% tlevel) filtereddf<- glomed.df %>% filter(Order %in% taxatofilter) %>% droplevels()
  else if ("Family" %in% tlevel) filtereddf<- glomed.df %>% filter(Family %in% taxatofilter) %>% droplevels()
  else if ("Genus" %in% tlevel) filtereddf<- glomed.df %>% filter(Genus %in% taxatofilter) %>% droplevels()
  else filtereddf<- glomed.df %>% filter(Species %in% taxatofilter) %>% droplevels()

  # Grouping the table
  cat("Grouping the samples for the Day and the",tlevel,"level\n")

  if ("Phylum"%in%tlevel) x<- filtereddf %>% group_by_(x_var,tlevel)
  else if ("Class"%in%tlevel) x<- filtereddf %>% group_by_(x_var,tlevel)
  else if ("Order"%in%tlevel) x<- filtereddf %>% group_by_(x_var,tlevel)
  else if ("Family"%in%tlevel) x<- filtereddf %>% group_by_(x_var,tlevel)
  else if ("Genus"%in%tlevel) x<- filtereddf %>% group_by_(x_var,tlevel)
  else x<- filtereddf %>% group_by_(x_var,tlevel)

  # averaging and taking log
  cat("Averaging and taking log of",tlevel," level abundance\n")

  avgx <-x%>%summarise(Abundance=mean(Abundance))
  # avgx$Abundance<-log(avgx$Abundance)
  tl<- colnames(avgx[2])
  # avgx[,1]<-factor(avgx[,1], levels =c(unique(avgx[,1])))

  # x_lab <-factor(avgx[1])
  # avgx[,1]<-factor(avgx[,1], levels=x_lab)

  # avgx[1]<-as.array(avgx[1])

  # Plotting the taxa
  cat("Plotting the Day vs log abundance plot for the",diet,"samples at ", tlevel," level using the",taxa1, "and", taxa2, "\n")

  plt<-ggplot(avgx,aes_string(x=paste(x_var),y="Abundance",color=tl,group=tl))+geom_line(size=1)+geom_point(size=4,alpha=0.6) +
    theme_minimal() +
    scale_color_manual(values=cols)+
    # scale_x_discrete(labels=x_lab) +
    theme(axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.text.x = element_text(color = "black",size=12,face = "bold"),
          axis.text.y = element_text(color = "black",size=12,face = "bold")) +
    theme(legend.title = element_text(color = "black",size=14,face = "bold")) +
    theme(legend.text = element_text(color = "black",size=12,face = "bold"))
  plt

}
