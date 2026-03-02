install.packages("tidyverse")
library(tidyverse)

OTU <- read.table("FAB2_2022_ITS_mumu_curated.txt", header = T)

Tax_raw <- read.table("Taxonomy_sintax_ITS2.txt", fill = T)

colnames(Tax_raw) <- c("OTUId", "Taxonomy", "Strand", "Tax_above_threshold")

Tax_raw1 <- Tax_raw %>% separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "[,]")

Tax_raw1 <- Tax_raw1 %>% mutate(Kingdom = gsub("k:", '', Kingdom),
                                Phylum = gsub("p:", '', Phylum),
                                Class = gsub("c:", '', Class),
                                Order = gsub("o:", '', Order),
                                Family = gsub("f:", '', Family),
                                Genus = gsub("g:", '', Genus),
                                Species = gsub("s:", '', Species))
Tax_raw2 <- Tax_raw1 %>% separate(Kingdom, c("Kingdom", "Kingdom_support"), sep = "[(]") %>% 
  separate(Phylum, c("Phylum", "Phylum_support"), sep = "[(]") %>%
  separate(Class, c("Class", "Class_support"), sep = "[(]") %>%
  separate(Order, c("Order", "Order_support"), sep = "[(]") %>%
  separate(Family, c("Family", "Family_support"), sep = "[(]") %>%
  separate(Genus, c("Genus", "Genus_support"), sep = "[(]") %>%
  separate(Species, c("Species", "Species_support"), sep = "[(]") 

Tax_raw2 <- Tax_raw2 %>% mutate(Kingdom_support = gsub("[)]", '', Kingdom_support),
                                Phylum_support = gsub("[)]", '', Phylum_support),
                                Class_support = gsub("[)]", '', Class_support),
                                Order_support = gsub("[)]", '', Order_support),
                                Family_support = gsub("[)]", '', Family_support),
                                Genus_support = gsub("[)]", '', Genus_support),
                                Species_support = gsub("[)]", '', Species_support))

Tax1 <- Tax_raw2 %>% select(c("OTUId","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
Tax2 <- Tax_raw2 %>% select(c("Kingdom_support", "Phylum_support", "Class_support", "Order_support", "Family_support", "Genus_support", "Species_support"))
Tax3 <- Tax_raw2 %>% select(c("Strand", "Tax_above_threshold"))

Tax <- cbind(Tax1, Tax2, Tax3)

OTU_with_taxonomy <- left_join(OTU, Tax, by = "OTUId")

write.table(OTU_with_taxonomy, "FAB2_ITS_OTU_with_taxonomy_ITS2.txt", sep = "\t", quote = F, row.names = F)


