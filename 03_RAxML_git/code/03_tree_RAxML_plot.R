# Packages needed ---------------------------------------------------------

library(treedataverse)
library(ggbreak)

# 0. Import data ---------------------------------------------------------

# RAxML analysis output with the best tree with bootstrap support values
tree <- read.raxml("03_RAxML/result/RAxML_bipartitionsBranchLabels.best_bootstrap_tree_1000_plus")

ggtree(tree) +
  geom_nodelab(aes(label=node))


ggtree(tree_rooted) +
  geom_nodelab(aes(label=node))

nodes_bs <- tree@data

# 1. Root the tree -----------------------------------------------------
tree_rooted <- root(tree, outgroup = "UFG3804", resolve.root = T)

# Assign equal branch lenghts to the rooted tree
tree_rooted@phylo$edge.length <- rep(0.01, nrow(tree_rooted@phylo$edge))

# Plot tree with equal branch lenght and names of the nodes
ggtree(tree_rooted) +
  geom_tiplab(aes(label=label),size = 2) +
  geom_nodelab(aes(label=label))

# Assign original branch lenghts ------------------------------------------
# Creo un data frame con la informacion sobre cada eje usado
edges_tree <- as.data.frame(tree@phylo$edge)|> 
  mutate(edge.length = tree@phylo$edge.length) |> 
  rename(parent_node = V1,
         child_node = V2)

# Cargo la info de los nuevos ejes generados en el rooted tree, y le asocio 
# el lenghts
edges_tree_rooted <- as.data.frame(tree_rooted@phylo$edge)|>
  rename(parent_node = V1,
         child_node = V2) |> 
  left_join(edges_tree, by = join_by(child_node)) |> 
  select(!3)

# Asigno un largo al nuevo eje
edges_tree_rooted[166,3] <- 0.001

# Assign original branch lenghts in the rooted tree
tree_rooted@phylo$edge.length <- edges_tree_rooted$edge.length

# Reduce branch lenght of the outgroup
tree_rooted@phylo$edge.length[[83]] <- 1.78356e-02


# Assign tip labels -------------------------------------------------------
tree_rooted@phylo$tip.label <- tree@phylo$tip.label


ggtree(tree_rooted) +
  geom_tiplab(aes(label=label),size = 2) 
  # geom_nodelab(aes(label=label))



# Assign bootstrap values -------------------------------------------------------
tree_rooted@data$node <- tree@data$node
tree_rooted@data$bootstrap <- tree@data$bootstrap
# Data frames with populations ID, and subsps
pops_df <- read.csv("00_Sampling_sites/Sampling_sites83.csv", sep = ";")|>
  select(!X) |>
  select(!Populations) |> 
  rename(label = Sample,
         Populations = Populations2)# Rename to match the "label" column of tree_tibble

# Data frame with phenotypes data
phenotypes_df <- read.csv("03_RAxML/data/Phenotypes.csv", sep = ";") |> 
  select(!2:6) |> 
  rename(label = Sample.ID)

# Merge data
pops_df2 <- full_join(pops_df, phenotypes_df, by = "label")

# 1. Uno data del arbol con otros datos -----------------------------------
tree_pop_tibble <- full_join(tree_rooted_tibble, pops_df2, by = "label")


# # Nombro al outgroup en las columnas de pop y sub especie
# tree_pop_tibble[83, 17:18] <- "Outgroup"
# 
# # Copio el elemnto para luego generar un tibble con las subespecies ordenadas tb
# # Se usa en el paso 3b!
# tree_subsp_tibble <- tree_pop_tibble
# tree_pheno_tibble <- tree_pop_tibble
# 
# # Ordeno poblaciones de Norte a Sur para luego graficar
# tree_pop_tibble$Populations <- factor(tree_pop_tibble$Populations, 
#                                       levels = c("Guyana", 
#                                                  "Para",
#                                                  "Amapa",
#                                                  "Marajo",
#                                                  "Caatinga",
#                                                  "East",
#                                                  "West",
#                                                  "Outgroup"))

# Lo vuelvo a convertir en treedata
tree_pop <- as.treedata(tree_pop_tibble)
tree_pop

# Agrupo en un clado al ingroup
tree_pop1 <- groupClade(.data = tree_pop,
                    .node = 85,
                    group_name = "Outgroup")

ggtree(tree_pop1, branch.length = "branch.lenght")


# 2. Manipulating tree data -----------------------------------------------

# Reroot tree
tree_pop_rooted <- root(tree_pop, outgroup = "UFG3804")
tree_pop_rooted <- rename_taxa(tree = tree_pop_rooted, 
                               data = tree_pop_tibble, 
                               key = node, 
                               value = label)
tree_pop_rooted@phylo$edge.length <- tree_pop@phylo$edge.length
tree_pop_rooted

# Shorten tree branch to visualize
tree_pop_rooted@phylo$edge.length[[83]] <- mean(tree_pop_rooted@phylo$edge.length)
tree_pop_rooted@phylo$edge.length[[83]] <- 1.078356e-02

# 3. Tree visualization ------------------------------------------------------

## 3.a. Plot by populations -----------------------------------------------------

tree_plot_pop <- 
  ggtree(tree_pop_rooted) +
  geom_tippoint(size = 3, shape = 15, aes(color = Populations)) +
  geom_tiplab(aes(label=label),size = 2) +
  # geom_hilight() +
  scale_colour_manual(values = c("#d4b9da", # Guyana
                                 "#df65b0", # Para
                                 "#e7298a", # Amapa
                                 "#a50f15", # Marajo 
                                 "#54278f", # Caatinga
                                 "#dab600", # East
                                 "#6baed6", # West
                                 "#6e6e6e")) + # Outgroup
  geom_nodepoint(aes(subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 50)) +
  #geom_label2(aes(label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 50)) +
  geom_treescale(x=0, y=45)
  

tree_plot_pop

ggsave(filename = "03_RAxML/plot/03_RAxMLtree.pdf",
       plot = tree_plot_pop,
       width = 7,
       height = 7.5)
ggsave(filename = "03_RAxML/plot/03_RAxMLtree.png",
       plot = tree_plot_pop,
       width = 13,
       height = 7.5)

## 3.b. Plot by subsps -----------------------------------------------------
tree_subsp_tibble$Subespecie <- factor(tree_subsp_tibble$Subespecie, 
                                       levels = c("cayana", "huberi", "flava", 
                                                  "sincipitalis", "margaritae",
                                                  "chloroptera", "Outgroup"))
tree_subsp <- as.treedata(tree_subsp_tibble)
tree_subsp_rooted <- root(tree_subsp, outgroup = "UFG3804")
tree_subsp_rooted <- rename_taxa(tree = tree_subsp_rooted, 
                               data = tree_subsp_tibble, 
                               key = node, 
                               value = label)

tree_plot_subsp <- 
  ggtree(tree_subsp_rooted) + 
  geom_tippoint(size = 3, shape = 15, aes(color = Subespecie)) +
  geom_tiplab(size = 2) +
  scale_colour_manual(values = c("#df65b0",  # cayana
                                 "#a50f15",  # huberi
                                 "#54278f",  # flava
                                 "#dab600",  # sincipitalis
                                 "#6baed6",  # margaritae
                                 "#007e19",  # chloroptera
                                 "#6e6e6e")) + # Outgroup
  geom_label2(aes(label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 70))
tree_plot_subsp

# Save the plot
ggsave(filename = "03_RAxML/plot/03_RAxMLtree_subsp.pdf",
       plot = tree_plot_subsp,
       width = 13,
       height = 7.5)
ggsave(filename = "03_RAxML/plot/03_RAxMLtree_subsp.png",
       plot = tree_plot_subsp,
       width = 13,
       height = 7.5)

## 3.c. Plot by phenotype -----------------------------------------------------
tree_pheno_tibble$Phenotype <- factor(tree_subsp_tibble$Phenotype, 
                                       levels = c("cayana", 
                                                  "flava", 
                                                  "DZ ?"))
tree_pheno <- as.treedata(tree_pheno_tibble)
tree_pheno_rooted <- root(tree_pheno, outgroup = "UFG3804")
tree_pheno_rooted <- rename_taxa(tree = tree_pheno_rooted, 
                                 data = tree_pheno_tibble, 
                                 key = node, 
                                 value = label)

tree_plot_pheno <- 
  ggtree(tree_pheno_rooted) + 
  geom_tippoint(size = 3, shape = 15, aes(color = Phenotype)) +
  geom_tiplab(size = 2) +
  scale_colour_manual(values = c("#f0f0f0",  # cayana
                                 "#252525",  # flava
                                 "#cb181d")) + # DZ?
  geom_label2(aes(label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 70))
tree_plot_pheno

# Save the plot
ggsave(filename = "03_RAxML/plot/03_RAxMLtree_pheno.pdf",
       plot = tree_plot_pheno,
       width = 13,
       height = 7.5)
ggsave(filename = "03_RAxML/plot/03_RAxMLtree_subsp.png",
       plot = tree_plot_subsp,
       width = 13,
       height = 7.5)

ggtree(tree) +
  geom_tiplab(aes(label=label),size = 2) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 75)) +
  geom_label2(aes(label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) > 50), 
              size) 
