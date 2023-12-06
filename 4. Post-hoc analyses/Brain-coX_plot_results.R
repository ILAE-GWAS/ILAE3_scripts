# R code to summarise and plot brain-coX results

library(tidyverse)
library(readxl)
library(patchwork)

setwd("/brain-coX_results")

gge.data <- read_excel("GGE_cand_resultsv2.xlsx")
all.data <- read_excel("allEpi_cand_resultsv2.xlsx")
jme.data <- read_excel("JME_cand_results.xlsx")

gge.data %>% 
    group_by(Locus) %>% 
    summarise(mean = mean(`Number of Times Prioritised`))

gge.long <- gge.data %>% 
    pivot_longer(cols = c("Hawrylycz", "Miller", "Kang" ,"Colantuoni","Hernandez", "Trabzuni", "Zhang"),
                 names_to = "Resource") %>% 
    filter(value == 1) %>% 
    select(-(value))

gge.data$Prioritised <- factor(gge.data$Prioritised,                 # Relevel group factor
                         levels = c("Yes", "No"))

gge.data$Locus <- factor(gge.data$Locus,                 # Relevel group factor
                         levels = c("1q43", "2p16.1", "2q12.1", "2q24.3",
                                    "2q32.2", "3p21.31", "3p22.3", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1",
                                    "9q21.32", "10q24.32", "12q13.13", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32"))

p = gge.data %>% 
    ggplot(aes(x=Locus, fill = Prioritised)) +
    geom_bar() +
    theme_classic() + ylab("Candidate gene count") + 
    scale_fill_manual(values=c('#E69F00','#999999')) +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6)) +
    ggtitle("Generalised epilepsy candidates")

p

all.data$Prioritised <- factor(all.data$Prioritised,                 # Relevel group factor
                               levels = c("Yes", "No"))

all.data$Locus <- factor(all.data$Locus,                 # Relevel group factor
                         levels = c("2p16.1", "2q24.3", "9q21.13", "10q24.32"))

p1 = all.data %>% 
    ggplot(aes(x=Locus, fill = Prioritised)) +
    geom_bar() +
    theme_classic() + ylab("Candidate gene count") + 
    scale_fill_manual(values=c('#e65800','#999999')) +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6)) +
    ggtitle("`All epilepsy` candidates")

p1

jme.data$Prioritised <- factor(jme.data$Prioritised,                 # Relevel group factor
                               levels = c("Yes", "No"))

jme.data$Locus <- factor(jme.data$Locus,                 # Relevel group factor
                         levels = c("4p12", "8q23.1", "16p11.2"))

p2 = jme.data %>% 
    ggplot(aes(x=Locus, fill = Prioritised)) +
    geom_bar() +
    theme_classic() + ylab("Candidate gene count") + 
    scale_fill_manual(values=c('yellowgreen','#999999')) +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6)) +
    ggtitle("JME candidates")

p2

p1 / p / p2
