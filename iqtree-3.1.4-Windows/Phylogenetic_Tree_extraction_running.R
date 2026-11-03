##### THIS IS THE CORRECT COMMAND LINE TO RUN, THIS WILL ETRACT THE SPECIFIC SPECIES AND SAVE IN NEWICK FORMAT

PrimateFamily_Tree <- read.tree("C:/Users/adesi/Downloads/primate_species.nwk")
PrimateFamily_Tree
plot(PrimateFamily_Tree, show.node.label = TRUE, node.pos = TRUE)
#plotTree(PrimateFamily_Tree, fsize= 1.5, lwd = 2)
PrimateFamily_Tree$tip.label

# Specify a list of species
Pr.genus <-  c('Nycticebus_coucang',
              'Carlito_syrichta',
       'Otolemur_garnettii',
                'Propithecus_coquereli',
              'Eulemur_rufifrons',
             'Lemur_catta',
               'Microcebus_murinus',
               'Cebus_imitator',
               'Sapajus_apella',
             'Piliocolobus_tephrosceles',
             'Aotus_nancymaae',
             'Saimiri_boliviensis',
             'Pongo_abelii',
             'Symphalangus_syndactylus',
               'Nomascus_leucogenys',
              'Mandrillus_leucophaeus',
             'Papio_anubis',
               'Trachypithecus_francoisi',
              'Rhinopithecus_roxellana',
               'Pongo_pygmaeus',
              'Rhinopithecus_bieti',
             'Gorilla_gorilla_gorilla',
             'Callithrix_jacchus',
              'Colobus_angolensis_palliatus',
             'Chlorocebus_sabaeus',
              'Macaca_mulatta',
             'Cercocebus_atys',
             'Macaca_fascicularis',
            'Macaca_thibetana_thibetana',
             'Macaca_nemestrina',
             'Pan_paniscus',
             'Homo_sapiens',
             'Pan_troglodytes')

Pr.genus
# CHECK FOR SPECIES IN TREE VS IN DATA
chk <- name.check(PrimateFamily_Tree, Pr.genus)
summary(chk)
print(chk)
chk$data_not_tree
chk$tree_not_data


#DROP SPECIES NOT IN MY DATA
keep <- drop.tip(PrimateFamily_Tree, chk$tree_not_data)
keep 



# Number of tips in your tree
Ntip(PrimateFamily_Tree)

# Number of unmatched taxa
length(chk$tree_not_data)

# Number of taxa matched between tree and data
length(setdiff(
  PrimateFamily_Tree$tip.label,
  chk$tree_not_data
))

##IDENTIFY SPEICIES THAT EXIST IN BOTH DATA AND TREE
#First 20 tree species
head(PrimateFamily_Tree$tip.label, 20)
head(Pr.genus, 20)

intersect(PrimateFamily_Tree$tip.label, Pr.genus)

#IDENTIFY SPECIES TO REMOVE AND THOSE TO RETAIN 
species_to_Drop <- setdiff(PrimateFamily_Tree$tip.label, Pr.genus)
keep <- drop.tip(PrimateFamily_Tree, species_to_Drop)
species_to_Drop
keep
plot(keep, show.node.label = TRUE, use.edge.length = TRUE)

setwd("C:/Users/adesi/Downloads")
getwd()


#SAVE OUTPUT AS NEWICK FILE
write.tree(keep, file= "example_1.Newick")



#EXTRACT NEEDED SPECIES FROM TREE, PLOT   
Primate_genus_nodes_1 <- sapply(Pr.genus, grep, x= PrimateFamily_Tree$tip.label)
Primate_genus_nodes_1
plotTree(PrimateFamily_Tree, type= "fan", fsize = 0.6, lwd=1, ftype="i")
add.arrow(PrimateFamily_Tree, tip= Primate_genus_nodes_1, arr1= 0.15, col="red")


#GET MOST RECENT COMMON ANCESTOR OF SPECIES OF INTEREST FROM TREE
#xx= getMRCA(PrimateFamily_Tree, Primate_genus_nodes_1)

#first convert from list to integer
bbv <- unlist(Primate_genus_nodes_1)
bb <- getMRCA(PrimateFamily_Tree, bbv)
bb


#EXTRACT SPECIES AND COMMON ANCESTORS, PLOT
Pr.clade <- extract.clade(PrimateFamily_Tree, bb)
Pr.clade
plot(Pr.clade, show.node.label = TRUE, use.edge.length = FALSE)
plotTree(Pr.clade, type= "fan", fsize = 0.6, lwd=1, ftype="i")

#SAVE OUTPUT AS NEWICK FILE
write.tree(Pr.clade, file= "example.Newick")




