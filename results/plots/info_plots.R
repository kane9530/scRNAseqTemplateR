library(tidyverse)
library(readxl)
library(scales)
library(showtext)
library(ggrepel)

# Add Roboto fonts from google using showtext
font_add_google("Roboto")
showtext_auto()

# Load tidy data from Svensson et al., (2020)
d <- read_excel("../data/scRNAseq_studies_db.xlsx") %>%
  as_tibble() %>%
  rename(total_cells = `Reported cells total`) %>%
  filter(total_cells < 2e6)
     
# Convert Date column to datetime format and obtain year
extractYear <- function(x) {
  x2 <- x %>%
    as.POSIXct() %>%
    format(format="%Y")
}

d <- d %>%
  mutate(year = extractYear(Date))

# Plot 1: change in techniques used over years

chromium <- grepl(".*[C|c]hromium.*", d$Technique)
smart_seq <- grepl(".*[S|s]mart-seq2?.*", d$Technique)
#drop_seq <- grepl(".*[D|d]rop-seq.*", d$Technique)
#smarter <- grepl(".*SMARTer.*", d$Technique)
#indrops <- grepl(".*[I|i]n[D|d]rops*", d$Technique)
#labels <- c(chromium, smart_seq, drop_seq, smarter, indrops)

p_all <- ggplot(data = d, aes(x=year, y =total_cells,
                                   ymin=1e2, ymax=1.5e6)) +
  geom_boxplot(aes(group=year), lwd=1)+
  geom_jitter(alpha=0.5, size =3, aes(color=Technique))+
  scale_y_continuous(trans = log10_trans())+
  labs(x="Date", y="Total cells profiled")+
  ggtitle(label = "Landscape of Single-cell RNA sequencing experiments",
          subtitle = "Data from Svensson et al., 2020")+
  theme_minimal()+
  theme(legend.position="none",
        text = element_text(family = "Roboto"))

p_chromium <- ggplot(data = d, aes(x=year, y =total_cells,
                           ymin=1e2, ymax=1.5e6)) +
  geom_boxplot(aes(group=year), lwd=1)+
  geom_jitter(data = subset(d, chromium),alpha=0.5, size =3, color="grey3")+
  scale_y_continuous(trans = log10_trans())+
  labs(x="Date", y="Total cells profiled")+
  ggtitle(label = "Landscape of Single-cell RNA sequencing experiments",
          subtitle = "Data from Svensson et al., 2020")+
  theme_minimal()+
  theme(legend.position="none",
        text = element_text(family = "Roboto"))

p_smartseq <- ggplot(data = d, aes(x=year, y =total_cells,
                                   ymin=1e2, ymax=1.5e6)) +
  geom_boxplot(aes(group=year), lwd=1)+
  geom_jitter(data = subset(d, smart_seq),alpha=0.5, size =3, color="grey3")+
  scale_y_continuous(trans = log10_trans())+
  labs(x="Date", y="Total cells profiled")+
  ggtitle(label = "Landscape of Single-cell RNA sequencing experiments",
          subtitle = "Data from Svensson et al., 2020")+
  theme_minimal()+
  theme(legend.position="none",
        text = element_text(family = "Roboto"))

# Plot 2: Tissue types against total number of experiments

tissue_counts <- d %>%
  separate_rows(Tissue, sep=",", convert=TRUE) %>%
  mutate(Tissue_trim = stringr::str_trim(Tissue, side = "both")) %>% 
  group_by(Tissue_trim) %>%
  count() %>%
  arrange(desc(n)) %>%
  na.omit() %>%
  ungroup()%>%
  top_n(20)

p_counts <- ggplot(data=tissue_counts,
                   aes(x = reorder(Tissue_trim,-n), y = n))+
  geom_col(aes(fill=Tissue_trim))+
  coord_flip()+
  labs(x="Number of experiments", y="Tissue Type")+
  ggtitle(label = "scRNA-seq experiments by tissue type",
          subtitle = "Data from Svensson et al., 2020")+
  theme_minimal()+
  theme(legend.position="none",
        text = element_text(family = "Roboto"))