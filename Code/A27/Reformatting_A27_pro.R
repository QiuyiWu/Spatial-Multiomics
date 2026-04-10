library(tidyr)
a27.pro <- read_tsv("../../../data/10788_SNE_combine_final_012026_Report_protein_export (Normal).tsv")
# SAME DukeID: "ID124351_01_SPD60_OA10222_10788_040325" "ID124351_01_SPD60_OA10222_10788_050725", REMOVE 050725 ONE. 
a27.pro <- a27.pro[!a27.pro$R.FileName == "ID124351_01_SPD60_OA10222_10788_050725",]

a27.pro <- data.frame(a27.pro)
a27.pro.select1 <- a27.pro[, c("R.FileName", "PG.ProteinAccessions", "PG.Quantity")]
# add the plate to the data
plate.df <- read_tsv("../data/20260120_122736_Identifications.tsv")
plate.df <- data.frame(plate.df)
plate.df$FileName <- sub("\\.htrms$", "", plate.df$FileName)
plate.df <- plate.df[, c("FileName", "Plate.", "SPQC.type")]
colnames(plate.df) <- c("R.FileName", "Plate", "SampleType")
a27.pro.select1 <- a27.pro.select1 %>% left_join(plate.df, by = "R.FileName")



colnames(a27.pro.select1)[c(3)] <- c("Y")
a27.pro.select1$Y <- log2(a27.pro.select1$Y)
a27.pro.select1$ID.metadata <- substr(a27.pro.select1$R.FileName, 1,8)
a27.pro.select1$pro.sampletype <- ifelse(is.na(a27.pro.select1$SampleType), "Non-SPQC", "SPQC")   
a27.pro.selected <- a27.pro.select1[, -1]
# Convert long data to wide data
wide_data <- a27.pro.selected %>%
  pivot_wider(
    id_cols = c(ID.metadata, Plate, pro.sampletype, SampleType),
    names_from = PG.ProteinAccessions,
    values_from = Y,
    values_fn = mean,
    values_fill = NA
  )

wide_data <- wide_data %>%
  select(ID.metadata, Plate, pro.sampletype, SampleType, everything()) %>% data.frame()
pro.a27.adj <- wide_data 
colnames(pro.a27.adj)[1:3] <- c("pro.sample", "pro.plate", "pro.sampletype")


save(pro.a27.adj, file = "Adjusted_Protein_A27.RData") #  566 x 10253


