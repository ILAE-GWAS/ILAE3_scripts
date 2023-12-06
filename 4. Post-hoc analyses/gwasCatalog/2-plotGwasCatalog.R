## ILAE GWAS 3.0 - Plot GWAS Catalog data

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
