# Creates histogram of cufflinks results
# Modified from Mae Rose Gott's script
# Allison Nau 2021/05/07

#Reads table
cuff <- read.table("/projectnb2/bf528/users/dachshund/project5_anau/P0_1_cufflinks/genes.fpkm_tracking")
#Set column names
colnames(cuff) <- cuff[1,] 
#delete row of column names
cuff <- cuff [-1,] 

#Prepare the data
cuff$FPKM <- as.numeric(cuff$FPKM) #Change to numeric

# Save old par:
oldpar<-par(no.readonly = T)


# Create histograms
png(filename="cufflink.density_scaling.png", width=600, height=600)
par(mfrow=c(2,1))
hist(log(cuff[cuff$FPKM>1,"FPKM"]), main = "A) FPKM greater than 1", 
     xlab="log(FPKM)", col='lightblue', breaks=25)
hist(10**(cuff[cuff$FPKM<=1 & cuff$FPKM>0, "FPKM"]), main="B) FPKM between 0 and 1", 
     xlab="10^FPKM", col='lightblue', breaks=25)
dev.off()

# Reset to oldpar:
# par(oldpar)