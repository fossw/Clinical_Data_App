library(lme4)
library(doBy)
library(multcomp)
library(stringr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(tidyr)
library(ggsignif)
library(shiny)
library(plotly)


update_model <- function(formula, model_number) {
  paste0("model", model_number, " <- lmer(\"", formula, "\", data=filtered)")
}

update_script <- function(updated_code, model_name) {
  # Path to your source script
  script_path <- "PFRpackage.R"  # Update this path
  # Read the source script
  script_lines <- readLines(script_path)
  # Find the line to replace
  line_to_replace <- grep(paste0(model_name, " <-"), script_lines)
  original_line <- script_lines[line_to_replace]
  indentation <- sub("^(\\s*).*", "\\1", original_line)  # Get leading whitespace
  updated_code <- paste0(indentation, updated_code)
  script_lines[line_to_replace] <- updated_code
  writeLines(script_lines, script_path)
  
}

MLMapply <- function(data) { # Applies mixed linear models to the data
  data <- data[data$Intensity > 0,] # Filters out 0 values
  
  data$Log2Y <- log2(data$Intensity) # Creates a new column for log change of quantifier data
  data
  data<- ddply(data, c("PFR"), transform, zScore = scale(Log2Y)) # Creates a new column for zScore, calculated from Log2Y
  
  # Calculate the missing values of the data set, filter out incomplete data
  counts <- data %>%  group_by(PFR, Treatment)  %>%   dplyr::summarize(counts = n()) 
  maxCounts <- counts %>%  group_by(Treatment)  %>%   dplyr::summarize(max = max(counts))
  counts <- merge(counts, maxCounts)
  data <- merge(data, counts)
  
  total <- nrow(data)
  
  data <-  subset(data, counts>=max/2)
  tooFew <- data %>%  group_by(PFR) %>% dplyr::summarize(groupCount = n_distinct(Treatment)) # count how many distinct groups
  data <- merge(data, tooFew)
  data <-  subset(data, groupCount > 1 )
  
  
  missing <- (1 -(nrow(data)/total))*100
  sprintf("Missing Values: %4f %% ", missing)
  
  
  data<- ddply(data, c("PFR"), transform, zScore = scale(Log2Y))
  # Save data subset specifically formatted for boxplot
  data4box <- data
  # Create new data frame for p value between each treatment group with std. error
  pValues <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(pValues) <- c("PFR","Treatment1","Treatment2","p.value","estimate","std.error","fold")
  # analyze data points by treatment group
  TreatmentLevels <- unique(data$Treatment)
  # data subset for heatmap data
  heatmap <- data.frame(matrix(ncol = length(TreatmentLevels)*2 + 1, nrow = 0)) # 2 columns per treatment + PFR number/any label we want.
  PFRs <- unique(data$PFR)
  CovAll <- data.frame(Name = numeric(),  grp = numeric(), percent = numeric())
  
  n_iter <- length(PFRs)
  pb <- txtProgressBar(min = 0, # Progress bar that increases after each data row iterated 
                       max = n_iter, 
                       style = 3,    
                       width = 50,   
                       char = "=")
  # Only suppressing warnings to view progress bar. This can be removed.  
  suppressWarnings({
    for (i in seq_along(PFRs)) {
      pfr <- PFRs[i]
      filtered <- data[data$PFR == pfr,]
      suppressMessages({
        # Mixed linear model comparing treatment to replicates as a random effect
        model1 <- lmer("zScore ~ -1 + (1|Treatment) + (1|BioRep) + (1|TechRep)", data=filtered)
      })
      covData <- as.data.frame(VarCorr(model1))
      covData["PFR"] = pfr
      covData$percent <- 100*covData$vcov / sum(covData$vcov)
      CovAll <- rbind(CovAll,covData[,c("PFR","grp","percent")])
      suppressMessages({
        # Mixed linear model comparing treatment to replicates as a fixed effect
        model2 <- lmer("zScore ~ Treatment -1 + (1|BioRep) + (1|BioRep:TechRep)", data=filtered)
      })
      lsm <- LSmeans(model2, effect="Treatment")$coef
      row.names(lsm) <- LSmeans(model2, effect="Treatment")$grid$Treatment
      
      # we need to get this into a format that works with the output dataframe, and to make sure that they line up.
      estimates <- list(pfr)
      for(treatment in TreatmentLevels){
        estimates <- append(estimates,lsm[treatment,"estimate"])
        estimates <- append(estimates,lsm[treatment,"std.error"])
      }
      
      
      
      g2 <- glht(model2, mcp(Treatment="Tukey")) # get p values from this step
      L2 <- g2$linfct
      blah <- linest(model2, L=L2)$coef[c("p.value")] # This gives the Diffs of LSMs, including p value
      blah[c('Treatment1', 'Treatment2')] <- str_split_fixed(rownames(blah), ' - ', 2)
      blah["PFR"] = pfr
      suppressMessages({
        # Mixed linear model comparing treatment to replicates as a fixed effect with log intensity value
        model3 <- lmer("Log2Y ~ Treatment -1 + (1|BioRep) + (1|BioRep:TechRep) + (1|BioRep:TechRep)", data=filtered)
      })
      #Extract p values from model
      g3 <- glht(model3, mcp(Treatment="Tukey"))
      L3 <- g3$linfct
      blah2 <- linest(model3, L=L3)$coef[c("estimate","std.error")] # This gives the Diffs of LSMs
      blah2[c('Treatment1', 'Treatment2')] <- str_split_fixed(rownames(blah2), ' - ', 2)
      blah2["PFR"] = pfr
      lsm2 <- LSmeans(model3, effect="Treatment")$coef
      row.names(lsm2) <- LSmeans(model3, effect="Treatment")$grid$Treatment
      
      for(treatment in TreatmentLevels){
        estimates <- append(estimates,lsm2[treatment,"estimate"])
        estimates <- append(estimates,lsm2[treatment,"std.error"])
      }
      
      heatmap <- rbind(heatmap,estimates) # add to the output df
      
      merged <- merge(blah,blah2,c("PFR","Treatment1","Treatment2"))
      
      merged["fold"] <- 2**merged$estimate
      
      pValues <- rbind(pValues,merged)
  
      setTxtProgressBar(pb, i)
      
    }
  }) 
  
  titles <- list("PFR")
  for(treatment in TreatmentLevels)
  {
    titles <- append(titles,sprintf("%s_estimateZ",treatment))
    titles <- append(titles,sprintf("%s_Std.ErrorZ",treatment))
  }
  
  colnames(heatmap) <- titles
  
  pValues <- pValues[order(pValues$p.value),]
  # p values adjusted by Benjamini Hochberg
  pValues["qValue"] <- p.adjust(pValues$p.value, method = "BH")
  
  colnames(pValues) <- c("PFR","Treatment1","Treatment2","p.value","estimate","std.error","fold","qValue")
  
  metaData <- unique(data %>% dplyr::select(c("PFR":"TargetMass","Accession":"Search_Qvalue")))
  
  finalData <- merge(pValues,metaData, c("PFR"))
  finalData <- merge(finalData,heatmap, c("PFR"))
  
  returning <- list("data4box" = data4box, "finalData" = finalData,
                    "pb" = pb, "missing" = missing)
  return(returning)
}

MLMBoxPlot <- function(data4box, finalData) {
  # Function to format data for boxplot generation
  finalData <- finalData[!duplicated(finalData), ]
  data4box <- data4box[!duplicated(data4box), ]
  finalData<- finalData[, c('PFR', 'qValue')]
  sub <- data4box[, c('PFR', 'Treatment', "Log2Y")]
  sub$qValue <- finalData$qValue[match(sub$PFR,finalData$PFR)]
  PFRs <- unique(sub$PFR)
  
  to_return <- list("sub" = sub, "PFRs" = PFRs)
  return(to_return)
  
}

MLMvolcanoPlot <- function(A,B, data, pfr) {
  # Function to generate volcano plot based on two treatment groups
  AvB <- data[data$Treatment1== A & data$Treatment2 == B ,]
  
  mode <- 1
  
  if (nrow(AvB) == 0)
  {
    AvB <- data[data$Treatment1== B & data$Treatment2 == A ,]
    mode <- -1
    
  }
  
  AvB$diffexpressed <- "NO"
  AvB$delabel <- NA
  
  AvB$diffexpressed[mode*log2(AvB$fold) >0.5 & AvB$qValue < 0.05 ] <- "Up"
  AvB$diffexpressed[mode*log2(AvB$fold) < -0.5 & AvB$qValue < 0.05 ] <- "Down"
  
  mycolors <- c("blue", "red", "gray")
  names(mycolors) <- c("Down", "Up", "NO")
  maxFold <- max(abs(log2(AvB$fold)))
  
  filtered <- AvB[AvB$PFR == pfr,]
  ggplot(data=AvB, aes(x=mode*log2(fold), y=-log10(qValue), col=diffexpressed)) +
    geom_text(data = filtered, aes(label = PFR)) +
    geom_point(size=1)+
    theme_minimal()+
    ggtitle(sprintf("%s vs %s",A,B)) +
    xlab(expression(log[2] ~ Fold ~ Change)) + ylab(expression(-log[10] ~ FDR)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(size=12, face="bold", hjust=0.5)) +
    theme(axis.title.x = element_text( size=10, face="bold", color="black")) +
    theme(axis.title.y = element_text( size=10, face="bold", color="black") ) +
    theme(axis.text.x = element_text( size=8,  color="black")) +
    theme(axis.text.y = element_text( size=8,  color="black") ) +
    theme(legend.position = "none") +
    geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype = "longdash") +
    geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") +
    scale_colour_manual(values = mycolors)+
    xlim((0-maxFold),(0+maxFold))

}

MLMHeatMap <- function(data) {
  # Function to create heatmap
  plot_ly(data, type = "heatmap")
}
QCdist <- function(massdata) {
  # Function to create a histogram of data distribution for QC purposes
  ggplot(aes(x=TargetMass), data=massdata) +
  geom_histogram(position = position_nudge(x = 0, y = 0),bins = 10,
                 color="#e9ecef", fill = "darkslateblue") +

  geom_boxplot(position = position_nudge(x = 0, y = 180),
               width=10, alpha = 0.7,
               color="black", fill = "darkslateblue")+
  theme(legend.position="none")+
  labs(title ="", x="Target Mass")+
  theme_bw()
}
