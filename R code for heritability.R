#Directory
setwd("C:/Users/g429833/Desktop/The2017_row dat/Curated_data2017")

# Package and libraries ( first install the package and load the library)
library(agricolae)
library(BGLR)
library(coda)
library(tidyverse)
library(ggplot2)
library(BMTME)
library(lme4)
library(pheatmap)
library(plyr)

# importing data
The2017=read.csv("data_curated2017.csv", header = TRUE)
head(The2017)
names(The2017)
str(The2017)
# Data transformation from numeric to class or from class to numeric
data3 = transform (The2017,  Year = factor(Year), Loc = factor(Loc),  Taxa = as.character(Taxa), Block  =factor(Block  ),   
                   Row = factor(Row), Column  =factor(Row ), HD = as.numeric(HD),  HT = as.numeric(HT))
str(data3) # check if the format is changed 

Data2017=data3[,-c( 1:7)] # avoid columns that are not needed in the analysis
head(Data2017) #check the columns


############################################################################################
require(ggpubr)
require(tidyverse)
require(Hmisc)
require(corrplot)
library(data.table) 
library(readr)
library(plyr)
library(plotrix)
library(tidyr)
library(ggcorrplot)
library(lavaan)
library(semPlot)
library(OpenMx)
library(knitr)
library(kableExtra)
library(GGally)
library(MplusAutomation)

#write.csv(data_long, "data_long.csv")# save the result
Data_CV_2017=Data2017[,-c(4:5)] # avoid columns that are not needed in the analysis
# change from wide to long format
data_long = gather(Data_CV_2017, Trait, value, Yield:Oil, factor_key=TRUE)

#Trait and location aggregate
aggregat_data <- ddply(data_long, c("Trait" , "Loc"), summarise,
                       Mean = mean(value, na.rm=TRUE),
                       Max =  max(value, na.rm=TRUE),
                       Min =  min(value, na.rm=TRUE),
                       CV = sd(value, na.rm=TRUE)/(mean(value, na.rm=TRUE)),
                       SE=std.error(value,na.rm=TRUE))
write.csv(aggregat_data, "Trait2017_Location_Mean_max_min_cv_SE.csv")

#Trait and location aggregate
aggregat_data1 <- ddply(data_long, c("Trait"), summarise,
                        Mean = mean(value, na.rm=TRUE),
                        Max =  max(value, na.rm=TRUE),
                        Min =  min(value, na.rm=TRUE),
                        CV = sd(value, na.rm=TRUE)/(mean(value, na.rm=TRUE)),
                        SE=std.error(value,na.rm=TRUE))
write.csv(aggregat_data1, "Trait2017_Mean_max_min_cv_SE.csv")

#####################################################################################

# Correlation
Data_Corr_2017=na.omit(Data2017[,-c(1:5)]) # avoid columns that are not needed in the analysis
M<-cor(na.omit(Data_Corr_2017))
corrplot(M, type = "upper",  method = "number")
write.csv(M, "Corr2017.csv")

cor_5 <- rcorr(as.matrix(Data_Corr_2017))
M <- cor_5$r
p_mat <- cor_5$P
corrplot(M, type = "lower", order = "hclust", 
         p.mat = p_mat, sig.level = 0.05)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
#jpeg("Corr_plot2017.jpg", width = 600, height = 400)
corrplot(M, method = "color", col = col(200),  
         type = "lower", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
         # Combine with significance level
         p.mat = p_mat, sig.level = 0.05,  
         # hide correlation coefficient on the principal diagonal
         diag = T, 
         number.cex= 14/ncol(p_mat),
         tl.cex=1.5,
         cl.cex=1.5)
#dev.off()

# an alternative plot
model.matrix(~0+., data=Data_Corr_2017) %>% 
  cor(use="pairwise.complete.obs") %>% 
  
  ggcorrplot(show.diag = F, type="lower",  lab=TRUE, lab_size=2, col = col(200),
             # Add coefficient of correlation
             tl.srt = 45, #Text label color and rotation
             # Combine with significance level
             p.mat = p_mat, sig.level = 0.05 )

#########################################################################################




#BLUP analysis
##################################################################################

Data2017=data3[,-c( 1:7)] # avoid columns that are not needed in the analysis
head(Data2017) #check the columns

# Paste extra column (i.e Env) that have year and location combined. Eg. Volga  and 2017 converted into Volga2017
Data2017$Env <- as.factor(paste0(Data2017$Loc,"_",Data2017$Year))
head(Data2017) # check the columns
colnames(Data2017)
nc=ncol(Data2017)

Data2017=Data2017[,c(1:5,7)]
# varible where the loop will store the result
BLUPs_byloc <- as.data.frame(unique(Data2017$Taxa))
BLUPs_comb <- as.data.frame(unique(Data2017$Taxa))
colnames(BLUPs_byloc) <- "Taxa" # giving column names
colnames(BLUPs_comb) <- "Taxa" # giving column names

dd=data.frame(matrix(NA,1, 6))
names(dd)=c("trt_byloc",        "grp",        "var1", "var2",       "vcov",     "sdcor")


# Trait and environment to work with
traits <- colnames(Data2017[,-c(1:5,nc)]) # names of traits that will be analyzed
unique_Env=as.vector(unique(Data2017$Env)) # unique list of location with year
###############################################################################################
#BLUPS by location
for (i in unique_Env){
  df <- Data2017[Data2017$Env==i, ]
  
  ### To get BLUP value, run the the following <for> loops
  for (j in 1:length(traits)){
    trait <- traits[j] # iterating each trait
    trait_df <- df[,trait] # extracting column data of each trait
    if (all(is.na(trait_df)) == TRUE) {# avoids inclusion of missing columns in some locations
      next # jump the empty column and continue to the <next> iteration of the remaining traits
    }
    set.seed(1)
    fm <- lmer(trait_df ~ (1|Taxa) + (1|Row) , data=df)  #  <lmer> R package BLUP model
    rr <- ranef(fm) # extracting random effect i.e intercept of each taxa
    blup <- as.data.frame(round(rr$Taxa, digits = 5)) # BLUP value of the trait
    taxa <- rownames(blup)
    trt <- paste0(trait,"_", i)
    blup_df <- as.data.frame(cbind(taxa,blup$`(Intercept)`))
    colnames(blup_df) <- c("Taxa",trt)
    BLUPs_byloc <- BLUPs_byloc %>% left_join(blup_df, by = "Taxa")
    
    
    
    # Source of variation 
    var_res=summary(fm)
    var_trait= as.data.frame(var_res$varcor)
    trt_byloc=  rep(paste0(trait,"_", i), nrow(var_trait))
    blup_df <- as.data.frame(cbind(trt_byloc,var_trait))
    dd <- rbind(dd, blup_df)
    dd1=dd
  }
}
head(BLUPs_byloc)
nrow(BLUPs_byloc)
head(dd)



###################################################################################
#BLUPS combined
for (j in 1:length(traits)){
  trait <- traits[j] # iterating each trait
  trait_df <- Data2017[,trait] # extracting column data of each trait
  if (all(is.na(trait_df)) == TRUE) {# avoids inclusion of missing columns in some locations
    next # jump the empty column and continue to the <next> iteration of the remaining traits
  }
  set.seed(1)
  #Explicit nesting
  #fm1 <- lmer(Yield ~  (1|Taxa) + (1|Location)+ (1|Location:Taxa)+ (1|Location:Row) + (1|Location:ColumnPYT), data=Data2017)  #  <lmer> R package BLUP model
  #Implicit nesting
  fm_comb <- lmer(trait_df ~ (1|Loc) + (1|Taxa) + (1|Loc:Row) , data=Data2017)  #  <lmer> R package BLUP model
  rr <- ranef(fm_comb) # extracting random effect i.e intercept of each taxa
  blup <- as.data.frame(round(rr$Taxa, digits = 5)) # BLUP value of the trait
  taxa <- rownames(blup)
  trt <- paste0(trait,"_", "Comb_2017")
  blup_df <- as.data.frame(cbind(taxa,blup$`(Intercept)`))
  colnames(blup_df) <- c("Taxa",trt)
  BLUPs_comb <- BLUPs_comb %>% left_join(blup_df, by = "Taxa")
  
  # Source of variation 
  var_res=summary(fm_comb)
  var_trait= as.data.frame(var_res$varcor)
  trt_byloc=  rep(paste0(trait,"_", "comb"), nrow(var_trait))
  blup_df <- as.data.frame(cbind(trt_byloc,var_trait))
  dd <- rbind(dd, blup_df)
  dd2=dd
}
head(dd)
head(BLUPs_comb)
nrow(BLUPs_comb)

#################################################################################
BLUP2017=data.frame(BLUPs_byloc, BLUPs_comb)
write.csv(BLUP2017, "BLUPS_2017.csv", row.names=FALSE) # Save BLUP result

LMER_var_dat2017=rbind(dd1, dd2)
write.csv(LMER_var_dat2017, "LMER_var_dat2017.csv", row.names=FALSE) # Save BLUP result

#####################################################################################





# Reshaping the result 
##############################################
library(reshape2)
library(plyr)

Blup_2017_data=read.csv("BLUPS_2017.csv", header=TRUE)
head(Blup_2017_data)
dm <- melt(Blup_2017_data)
dm=dm[,-2]
var <- as.data.frame(str_split_fixed(dm$variable,'_',3))
df2 <- cbind(dm[,1],var,paste0(var$V3,var$V2),dm[,3])
colnames(df2) <- c('Taxa','Trait','Location','Year',"env", 'value')
head(df2)
dataead <- df2 %>% spread(Trait,value)
write.csv(dataead, 'dataead.csv', row.names=FALSE)


#############################################################################################
#######################################################################################
# Genomic based heritability


# Create directory (i.e.folder). Put the input data into directory

# First load the packages
library(heritability)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(plyr)


######## Make a data frame for the list of genotypes with phenotype  data
# data import. The data have wide format ( Taxa, Location, year....Trait1...Trait2.... )
BLUP2017_spr <- read.csv("BLUP2017_spread.csv",header=T) 
BLUP2017_spr=data.frame(na.omit(BLUP2017_spr))
#BLUP2017=BLUP2017[order(BLUP2017$Taxa), ]
#nrow(BLUP2017)


pg=read.table("Transpose of Geno_GWAS2017.txt", header = TRUE, sep = "", row.names=1)# Marker data
pg=pg[unique(BLUP2017_spr$Taxa),]

rownames(pg)==unique(BLUP2017_spr$Taxa) # check if <FALSE> 
#
pg[1:5,1:5]
##change -1,0,1 to 0,1,2
pg <- pg + 1
pg[1:5,1:5]




taxa <- as.vector(BLUP2017_spr$Taxa)

M=tcrossprod(scale(pg))  # centered and scaled markers'
G=M/mean(diag(M)) # GBLUP calculated from <M>
G = as.matrix(G)
ncol(G)
G[1:10,1:10]

# Method2
M=pg

############################################

## Inner joining only Taxa that are found between the BLUPs phenotype data and genotype data
Geno_names=as.data.frame(row.names(pg))
colnames(Geno_names)="Taxa"
data <- BLUP2017_spr%>% inner_join(Geno_names, by="Taxa") # data to work with after joining
data[1:10,1:10]
nrow(data)


# The row  and column names of phenotype and genotype are the some
rownames(G) == unique(data$Taxa) # G from method 2
colnames(G) == unique(data$Taxa)

# Name of location and combined data over location
unique_Env=as.vector(unique(data$Location )) # unique list of location with year

# Create new matrix upon which the result from each <Location> are attached
Her1<- matrix(nrow=1, ncol=7,NA)
colnames(Her1) <- c("Environment", "Trait", "h2", 'va', 've', 'lower_CI', 'upper_CI')

for (i in unique_Env){
  # Subletting data by its location. Please, check <colnames(df)> and ncol(df)
  df <- data[data$Location==i, ] 
  traits <- colnames(df[,-c(1:4)]) # Trait names that will be analyzed
  
  # Matrix upon which the result from each <traits> within a <Location> are attached
  Her <- matrix(nrow=length(traits), ncol=7,NA) 
  colnames(Her) <- c("Environment", "Trait", "h2", 'va', 've', 'lower_CI', 'upper_CI')
  
  for (j in 1:length(traits)){
    trait <- traits[j] # iterating each trait
    trait_df <- df[,trait] # extracting column data of each trait
    
    if (all(is.na(trait_df)) == TRUE) {# avoids inclusion of missing columns in some locations
      next # jump the empty column and continue to the <next> iteration of the remaining traits
    }
    if(all(na.omit(trait_df==0))==TRUE){ # to avoid traits which have all zero BLUPs due to singularity
      next
    }
    set.seed(1)
    # The heritability model
    Model_out <- marker_h2(data.vector = trait_df , geno.vector=df$Taxa, max.iter = 1000, K = G)
    Her[j,1] <- i
    Her[j,2] <- trait
    Her[j,3] <- Model_out$h2
    Her[j,4] <- Model_out$va
    Her[j,5] <- Model_out$ve
    Her[j,6] <- Model_out$conf.int1[1]
    Her[j,7] <- Model_out$conf.int1[2]
    
  }
  # binding the result from all the loops
  Her1= rbind(Her1, Her)
}

write.csv(Her1, "Heritability2017_from_BLUP1.csv")
#remove(Her1,Her)# Remove from the environment if rerunning the code is needed.

#remove(Her1,Her)# Remove from the environment if rerunning the code is needed.

#########################################################################################
jpeg("Heritability_plot2017.jpg", width = 800, height = 400)

# plot genomic heritability
#import <gh=read.csv('Heritability2017_from_BLUP.csv',header=T)>
# or use the  <Her1> from the memory from analysis above
gh=as.data.frame(na.omit(Her1))# removes <NA> rows 
gh$Environment = as.factor( gh$Environment)
gh$h2 = as.numeric( gh$h2)
gh$lower_CI = as.numeric( gh$lower_CI)
gh$upper_CI = as.numeric( gh$upper_CI)
gh$SE =  (gh$upper_CI - gh$lower_CI)/2
dim(gh)
str(gh)
ns <- ggplot(gh,aes(x=Trait, y=h2))  + geom_bar(aes(fill=Environment),width=0.7, position = "dodge", stat='identity', size = 5)
ns + theme_bw() + labs(title =  "PYT 2017", x = "Traits") + theme( text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.major.y = element_line(size=0.1,color = "gray"), axis.line.x = element_blank())
dev.off()



####################
setwd("C:/Users/girma/Desktop/SPARSE WRITE UP/Heritability")
library(dataRetrieval)
library(dplyr) # for `rename`
library(tidyr) # for `gather`
library(ggplot2)
library(cowplot)



Herta_plot=read.csv("Heritability_curated all.csv", head = TRUE)

g1 <- ggplot(data = Herta_plot[,c(1:4)], aes(fill=Environment, x =Trait, y = h2), ylab("Heritability (h2)")) +
  geom_bar(position='dodge', stat='identity', width = 0.6)+
  theme_bw() +   scale_fill_manual(values = c('coral2','cyan3', 'navy', 'seagreen4','maroon4') ) +
  xlab("Traits") + ylab("Heritability (h2)")
g1 + facet_grid(Year ~ .) + theme(text = element_text(size = 23)) 


Herta_plot=read.csv("Heritability_curated all.csv", head = TRUE)

library(dplyr)
Herta_plot <- read.csv("Heritability_curated_all.csv", header = TRUE)

Herta_plot <- Herta_plot %>%
  mutate(Trait = recode(Trait,
                        "Glucan" = "BG",
                        "Groat_pct" = "GP",
                        "HD" = "HD",
                        "HT" = "HT",
                        "Mid" = "MK",
                        "Oil" = "FA",
                        "Plump" = "PK",
                        "Protein" = "PR",
                        "Thin" = "TK",
                        "TKW" = "TKW",
                        "TW" = "TW",
                        "Yield" = "YD"))

trait_order <- c("FA", "BG", "PR", "PK", "MK", "TK", "GP", "TKW", "TW", "YD", "HT", "HD")
Herta_plot$Trait <- factor(Herta_plot$Trait, levels = trait_order)

library(ggplot2)
g1 <- ggplot(data = Herta_plot[, c(1:4)], aes(fill = Environment, x = Trait, y = h2)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.6) +
  theme_bw() +
  scale_fill_manual(values = c('coral2', 'cyan3', 'navy', 'seagreen4', 'maroon4')) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  xlab("Traits") +
  ylab(expression(paste("Heritability  (",  h^2, ")")))

g1 + facet_grid(Year ~ .) + theme(text = element_text(size = 16))

##
#Heritbality for combined BLUP

library(ggplot2)
comb=Herta_plot[which(Herta_plot$Environment=='Combined'),]
# Calculate mean values
mean_values <- aggregate(h2 ~ Year + Trait, data = comb[, c(1:4)], FUN = mean)

g1 <- ggplot(data = comb[, c(1:4)], aes(fill = Year, x = Trait, y = h2)) +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.6) +
  geom_text(data = mean_values, aes(label = round(h2, 2), y = h2 + 0.05), position = position_dodge(width = 0.6), size = 4) +  # Add mean labels on top
  theme_bw() +
  scale_fill_manual(values = c('cyan3', 'navy', 'seagreen4')) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.1)) +
  xlab("Traits") +
  ylab(expression(paste("Heritability  (", h^2, ")"))) +
  facet_grid(Year ~ .) +
  theme(text = element_text(size = 16))

g1


