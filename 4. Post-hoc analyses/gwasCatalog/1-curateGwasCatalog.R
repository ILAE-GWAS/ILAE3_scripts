## ILAE GWAS 3.0 - curate gwascatalog hits

library(tidyverse)
library(corrplot)

catalog.dir <- " "

jme <- read_tsv(file.path(catalog.dir, "JME/gwascatalog.txt")) %>%
  filter(Region %in% c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32"))

count(jme)

all <- read_tsv(file.path(catalog.dir, "All-Epilepsy/gwascatalog_more.txt")) %>%
  filter(Region %in% c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32")) 

count(all)

gge <- read_tsv(file.path(catalog.dir, "GGE/gwascatalog_more.txt")) %>%
  filter(Region %in% c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32")) 

count(gge)

cae <- read_tsv(file.path(catalog.dir, "CAE/gwascatalog.txt")) %>%
  filter(Region %in% c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32"))

count(cae)

gwascataloghits <- rbind(all, gge, jme, cae)

addmargins(t(table(gwascataloghits$Region, useNA = "ifany")))


## Filtering

#GWAS Catalog traits filtered for:
#   * Pval < 5*10-8
#   * European ancestry (predominantly)
#   * GWAS result (exc. MTAG, exome-wide)
#   * Sumstats available (e.g. OR, direction of effect)
#   * When similar/same trait – keep largest &/or most recent entry only
#   * Exclude epilepsy (not FS)

gwascataloghits1 <- gwascataloghits %>% 
    filter(P <= 5*10^-8) %>%
  unique()

addmargins(t(table(gwascataloghits1$Region, useNA = "ifany")))

gwascataloghits2 <- gwascataloghits1 %>% 
  filter(!(grepl('epilepsy|Epilepsy', Trait))) %>%
  filter(!(grepl('epilepsy|Epilepsy', Study))) %>% 
  filter(!(grepl('MTAG', Study))) %>% 
  filter(!(grepl('MTAG', Trait))) %>% 
  separate(Strongest, into = c("strongest.rsid","EA"), sep = "-", remove = TRUE) %>% 
  filter(!(EA == "?")) %>% 
  group_by(Study, Region) %>%
  arrange(P) %>% 
  slice(1) %>% # select only lowest p-value trait per study per locus
  ungroup()

addmargins(t(table(gwascataloghits2$Region, useNA = "ifany"))) 

## High level trait groupings

gwascataloghits3 <- gwascataloghits2 %>% 
  mutate(Trait_Group = case_when(grepl('pressure|atrial|cholesterol|A1|triglyceride|diuretics|infarction|Electrocardiogram', Trait, ignore.case = TRUE) ~ "coronary", 
                                 grepl('hematocrit|cell|hemoglobin|corpuscular|lymphocyte|neutrophil|platelet|eosinophil|leukemia', Trait, ignore.case = TRUE) ~ "blood cell",
                                 grepl('smoking', Trait, ignore.case = TRUE) ~ "smoking",
                                 grepl('neurotic|schizophrenia|depress|irritable|nervous|stress|neuropsychiatric|well-being|miserable|worry|Neurociticism', Trait, ignore.case = TRUE) ~ "psychiatric",
                                 grepl('sleep|daytime|morning|chronotype|insomnia|snoring', Trait, ignore.case = TRUE) ~ "sleep",
                                 grepl('cogniti|ability|intelligence|education|math', Trait, ignore.case = TRUE) ~ "cognition",
                                 grepl('waist|body|height|hip circumference|joint space', Trait, ignore.case = TRUE) ~ "misc",
                                 grepl('sexual|alcohol|adventurousness|caffeine|walking|sun-seeking|vegetable|participation', Trait, ignore.case = TRUE) ~ "behavioural",
                                 grepl('graves|psoriasis|diabetes|lupus', Trait, ignore.case = TRUE) ~ "autoimmune",
                                 grepl('seizure|Parkinson|Alzheimer|PHF-tau', Trait, ignore.case = TRUE) ~ "neurological",
                                 grepl('sex hormone|biological sex|baldness|testosterone', Trait, ignore.case = TRUE) ~ "misc",
                                 grepl('osteoporosis|grip strength|heel bone', Trait, ignore.case = TRUE) ~ "misc",
                                 grepl('glomerular|urinary|non-albumin|phosphate', Trait, ignore.case = TRUE) ~ "misc",
                                 grepl('hair', Trait, ignore.case = TRUE) ~ "misc",
                                 TRUE ~ "misc"))

addmargins(t(table(gwascataloghits3$Region, gwascataloghits3$Trait_Group, useNA = "ifany"))) 

## Misc group

misc.group <- gwascataloghits3 %>% 
    filter(Trait_Group == "misc") # %>% 
    select(Region, Trait, PMID)

## Plots

gwascataloghits3$Region <- factor(gwascataloghits3$Region,                 # Relevel group factor
                         levels = c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32"))

dummy.data <- read_tsv(file.path(catalog.dir, "no_hits/gwascatalog_dummy.txt"))

gwascataloghits4 <- gwascataloghits3 %>%
  mutate(count = 1) %>% 
  rbind(dummy.data)
  
gwascataloghits4 %>% 
  filter(Trait_Group != "misc") %>% 
  ggplot(aes(x = Region, y = count, fill = Trait_Group)) +
  geom_bar(stat = "identity") + 
  theme_classic() + ylab("Publication count") +
  theme(legend.position="bottom", legend.title = element_blank(), axis.title.x = element_blank()) +
    ggtitle("GWAS Catalog entries") +
    scale_fill_brewer(palette="Set1") +
    theme(axis.text.x = element_text(size=16, angle=45, hjust=1),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 16),
          title = element_text(size = 20))

write_tsv(gwascataloghits3, file.path(catalog.dir, "catalog_filtered_more.tsv"))

## Plot matrix

cor.mat <- read_tsv(file.path(catalog.dir, "../matrix.txt"))

names <- c("1q43", "2p16.1", "2q12.1", 
                                    "2q32.2", "2q24.3", "3p21.31", "3p22.3", "4p12", "4p15.1",
                                    "5q22.3", "5q31.2", "6q22.33", "7p14.1", "8q23.1", "9q21.13",
                                    "9q21.32", "10q24.32", "12q13.13", "16p11.2", "16p13.3",
                                    "17p13.1", "17q21.32", "19p13.3", "21q21.1", "21q22.1", "22q13.32")

row.names(cor.mat) <- names
as.matrix(cor.mat)

cor.mat[is.na(cor.mat)] <- 0

corrplot(as.matrix(cor.mat), method = 'color', is.corr = FALSE, col = COL1("Purples"))
