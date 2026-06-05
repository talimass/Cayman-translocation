library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)
library(tidygraph)


#save.image(file='kinship5.RData')
#load('kinship.RData')
#### all clones/no relatives ####
setwd("~/haifa/cayman/rna/connect/admixture/new/revision/")
# Read KING table

king <- read.csv("idb.kin0", stringsAsFactors=FALSE, sep = "\t") # all


# Define thresholds
clone_cutoff <- 0.354      # definite clones
firstdeg_cutoff <- 0.177   # first-degree relatives

# Create a new column for relationship type
king$rel_type <- ifelse(king$KINSHIP > clone_cutoff, "Clone",
                        ifelse(king$KINSHIP > firstdeg_cutoff, "1st-degree", NA))

# Keep only edges above 0.177
edges <- subset(king, !is.na(rel_type), select=c(ID1, ID2, KINSHIP, rel_type))

# Build graph
g <- graph_from_data_frame(edges, directed=FALSE)
# Short label: remove everything until first '-' (including '-')
V(g)$label <- sub("-[^-]*$", "", V(g)$name)

# Site from full name (keeps working even after shortening)
V(g)$site <- str_extract(V(g)$name, "(CC|MF)")

# Depth class from full name
V(g)$depth <- dplyr::case_when(
  grepl("(40|DS|DD)", V(g)$name) ~ "Deep",
  grepl("(10|SS|SD)", V(g)$name) ~ "Shallow",
  TRUE ~ NA_character_
)

# Plot
kinship <- ggraph(g, layout="fr") +
  geom_edge_link(aes(color=rel_type, width=KINSHIP), alpha=0.6) +
  
  # nodes: shape = site, color = depth
  geom_node_point(aes(shape=site, color=depth), size=4) +
  
  # labels: shortened
  geom_node_text(aes(label=label), repel=TRUE, size=3) +
  
  #ggtitle("Pairwise kinship relationships among samples") +
  
  # edge colors 
  scale_edge_color_manual(values=c("Clone"="#6F2DBD", "1st-degree"="#F4A261")) +
  
  # node colors
  scale_color_manual(values=c("Shallow"="#f75f55", "Deep"="#00A9FF"), na.value="grey70") +
  
  # node shapes by site
  scale_shape_manual(values=c("CC"=16, "MF"=17)) +
  
  scale_edge_width_continuous(name="Kinship", range=c(0.3, 2)) +
  theme_void() +
  labs(edge_color="Relationship", color="Origin", shape="Site")

kinship
#saveRDS(kinship, "fig_kinship.rds")
ggsave("kinship.noartclones.jpg", kinship, width = 10, height = 10)

#### removing artificial clones ####
king <- read.csv("idb_no_art.kin0", stringsAsFactors=FALSE, sep = "\t") # no_art_clones


# Define thresholds
clone_cutoff <- 0.354      # definite clones
firstdeg_cutoff <- 0.177   # first-degree relatives

# Create a new column for relationship type
king$rel_type <- ifelse(king$KINSHIP > clone_cutoff, "Clone",
                        ifelse(king$KINSHIP > firstdeg_cutoff, "1st-degree", NA))

# Keep only edges above 0.177
edges <- subset(king, !is.na(rel_type), select=c(ID1, ID2, KINSHIP, rel_type))

# Build graph
g <- graph_from_data_frame(edges, directed=FALSE)
# Short label: remove everything until first '-' (including '-')
V(g)$label <- sub("-[^-]*$", "", V(g)$name)

# Site from full name (keeps working even after shortening)
V(g)$site <- stringr::str_extract(V(g)$name, "(CC|MF)")

# Depth class from full name
V(g)$depth <- dplyr::case_when(
  grepl("(40|DS|DD)", V(g)$name) ~ "D",
  grepl("(10|SS|SD)", V(g)$name) ~ "S",
  TRUE ~ NA_character_
)


# This leaves only the numbers (e.g., "./117-MF-SD" becomes "117")
V(g)$label <- gsub("^\\./|\\-.*$", "", V(g)$name)

# Create a combined group for coloring 
V(g)$group <- paste(V(g)$depth, V(g)$site, sep=".")

pt_to_mm <- function(pt) pt / 2.845

kinship <- ggraph(g, layout="fr") +
  # Edges
  geom_edge_link(aes(color=rel_type, width=KINSHIP), alpha=0.6) +
  
  # Nodes
  geom_node_point(aes(color=group), shape=16, size=2) +
  
  # Labels
  geom_node_text(aes(label=label), repel=TRUE, size=pt_to_mm(6)) +
  
  # Edge colors
  scale_edge_color_manual(values=c("Clone"="#3E2DBD", "1st-degree"="#F4A261")) +
  
  scale_color_discrete(name="Group") +
  
  scale_edge_width_continuous(name="Kinship", range=c(0.3, 0.8)) +
  
  theme_void() +
  labs(edge_color="Relationship") +
  theme(legend.position = "right")

kinship

# for revision - inch to mm
pt_to_mm <- function(pt) pt / 2.845

kinship <- ggraph(g, layout = "fr") +
  
  geom_edge_link(
    aes(color = rel_type, width = KINSHIP),
    alpha = 0.7
  ) +
  
  geom_node_point(
    aes(color = group),
    shape = 16,
    size = 2.0
  ) +
  
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = pt_to_mm(6)
  ) +
  
  scale_edge_color_manual(
    values = c(
      "Clone" = "#3E2DBD",
      "1st-degree" = "#F4A261"
    ),
    name = "Relationship"
  ) +
  
  scale_color_discrete(name = "Group") +
  
  scale_edge_width_continuous(
    range = c(0.28, 1),
    guide = "none"   # removes the 0.4 linewidth legend
  ) +
  
  theme_void(base_size = 6) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.tag = element_text(size = 9, face = "bold"),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )
kinship

saveRDS(kinship, "fig_kinship.rds")
ggsave("kinship.noartclones.jpg", kinship, width = 6, height = 6)
kinship <- readRDS("fig_kinship.rds")


