## Script to get every TCGA barcode from folder with RNAseq downloads
# Useful for checking if every sample was processed per file

library(TCGAutils)

# Give path to gdc download
path_downloads <- commandArgs(trailingOnly = TRUE)[1]

# Name of the output file
output_file <- paste(path_downloads, "TCGA_barcodes.txt", sep = "")

# Check if output file is already present, if not create txt file with barcode and corresponding uuid
if (file.exists(output_file)) {
  cat(output_file)
} else { 
  uuids <- list.dirs(path_downloads, full.names = FALSE, recursive = FALSE) # Reading in each file UUID
  
  barcode <- c() # Getting barcode of each file UUID shortened and unshortened
  full_bar <- c() 
  for (i in uuids){
    new <- UUIDhistory(i)
    bar <- UUIDtoBarcode(new[which.max(new$data_release), "uuid"], from_type = "file_id")
    barcode <- c(barcode, substr(bar[1, 2], 1, 12))
    full_bar <- c(full_bar, substr(bar[1, 2], 1, 28))
  }
  
  table <- data.frame(Barcode = c(full_bar),
                      Shortened = c(barcode),
                      UUID = c(uuids)) # Create table with barcodes and uuids
  
  write.table(table, file = output_file, quote = F, row.names = F) # Write the table
  cat(output_file) 
}










