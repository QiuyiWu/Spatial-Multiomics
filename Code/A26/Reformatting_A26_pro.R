library(tidyr)

a26.pro <- read_tsv("../../../data/A26 Proteomics with QC variables.txt")
a26.pro <- data.frame(a26.pro)
a26.pro <- a26.pro[, c("R.Condition", "PG.ProteinAccessions", "Y", "Plate", "ID.of.Bonnie.metadata", "PG.Qvalue", "PG.NrOfPrecursorsUsedForQuantification")]
# Remove SPQC_tech 
a26.pro <- a26.pro[a26.pro$R.Condition != "SPQC_tech",]

a26.pro <- a26.pro %>%
  mutate(pro.sampletype = ifelse(grepl("SPQC", R.Condition),"SPQC", "Non-SPQC"))

a26.pro.selected <- a26.pro %>%
  select(ID.of.Bonnie.metadata, pro.sampletype, Plate, PG.ProteinAccessions, Y)




# Convert long data to wide data
wide_data <- a26.pro.selected %>%
  pivot_wider(
    names_from = PG.ProteinAccessions,
    values_from = Y
  ) %>%
  distinct()

wide_data <- wide_data %>%
  select(ID.of.Bonnie.metadata, pro.sampletype, Plate, everything()) %>% data.frame()

pro.a26.adj <- wide_data 
colnames(pro.a26.adj)[1:3] <- c("pro.sample", "pro.sampletype", "pro.plate")
pro.a26.adj$pro.plate <- ifelse(pro.a26.adj$pro.plate == 1, "Plate1", 
                                ifelse(pro.a26.adj$pro.plate == 2,"Plate2", 
                                       ifelse(pro.a26.adj$pro.plate == 3,"Plate3", 
                                              ifelse(pro.a26.adj$pro.plate == 4,"Plate4", 
                                                     ifelse(pro.a26.adj$pro.plate == 5,"Plate5", "Plate6")))))



save(pro.a26.adj, file = "Adjusted_Protein_A26.RData")


