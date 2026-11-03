install.packages("ape")
install.packages("phangorn")
install.packages("geiger")

library(ape)
library(phangorn)
library(geiger)

packageVersion("ape")
packageVersion("phangorn")
packageVersion("phytools")
packageVersion("geiger")

ZEB2.string <- "(Nycticebus_coucang:0.008123,((Carlito_syrichta:0.011805,((((((((((Cebus_imitator:0.011608,Sapajus_apella:0.000000)N13:0.003237,Callithrix_jacchus:0.000810)N12:0.000000,(Aotus_nancymaae:0.001616,Saimiri_boliviensis:0.002428)N14:0.000807)N11:0.000807,((Pongo_abelii:0.002426,Pongo_pygmaeus:0.000000)N16:0.000808,(Trachypithecus_francoisi:0.000808,(Rhinopithecus_roxellana:0.000000,Rhinopithecus_bieti:0.000000)N18:0.000000)N17:0.000808)N15:0.000000)N10:0.000000,(Piliocolobus_tephrosceles:0.003233,(Symphalangus_syndactylus:0.000000,Nomascus_leucogenys:0.000808)N20:0.001616)N19:0.000000)N9:0.000000,Colobus_angolensis_palliatus:0.000809)N8:0.000000,(((Gorilla_gorilla_gorilla:0.001619,Homo_sapiens:0.000808)N23:0.000000,Pan_troglodytes:0.000000)N22:0.000000,Pan_paniscus:0.000000)N21:0.000790)N7:0.000788,Macaca_thibetana_thibetana:0.000000)N6:0.000001,Macaca_nemestrina:0.000000)N5:0.000001,(((((Mandrillus_leucophaeus:0.002952,Papio_anubis:0.015602)N28:0.011136,Cercocebus_atys:0.003231)N27:0.000001,Chlorocebus_sabaeus:0.000808)N26:0.000001,Macaca_mulatta:0.000808)N25:0.000001,Macaca_fascicularis:0.000000)N24:0.000001)N4:0.000812)N3:0.001658,(Propithecus_coquereli:0.001661,((Eulemur_rufifrons:0.000932,Lemur_catta:0.002304)N31:0.000808,Microcebus_murinus:0.002303)N30:0.000891)N29:0.000811)N2:0.006448,Otolemur_garnettii:0.005678)N1;
"
ZEB2_Tree <- read.tree(text= ZEB2.string)
ZEB2_Tree

plot(ZEB2_Tree, no.margin = TRUE)
help("plot.phylo")
args(plot.phylo)




##########Primates_Family
PrimateFamily.string <-  "(((Lorisidae:36.27257000,Galagidae:36.27257000)'13':21.37560000,((((Indriidae:26.28885000,Palaeopropithecidae:26.28885000)'12':7.85196000,(Lepilemuridae:32.29407000,Cheirogaleidae:32.29407000)'10':1.84674000)'18':1.32791000,Lemuridae:35.46872000)'17':16.44565000,Daubentoniidae:51.91437000)'9':5.73380000)'28':14.00024000,((((Atelidae:20.31717000,(Cebidae:18.67159000,Aotidae:18.67159000)'26':1.64558000)'24':1.14495000,Pitheciidae:21.46212000)'22':20.96798000,((Hylobatidae:19.59604000,Hominidae:19.59604000)'8':9.30190000,Cercopithecidae:28.89794000)'43':13.53216000)'51':22.48917000,Tarsiidae:64.91927000)'49':6.72914000);"
PrimateFamily_Tree <- read.tree(text = PrimateFamily.string)





##### THIS IS THE CORRECT COMMAND LINE TO RUN, THIS WILL ETRACT THE SPECIFIC SPECIES AND SAVE IN NEWICK FORMAT
PrimateFamily_Tree <- read.tree("C:/Users/adesi/Downloads/primate_species.nwk")
PrimateFamily_Tree
plot(PrimateFamily_Tree, show.node.label = TRUE, node.pos = TRUE)
#plotTree(PrimateFamily_Tree, fsize= 1.5, lwd = 2)
PrimateFamily_Tree$tip.label

# Specify a list of family
#Pr.genus <-  c("Pitheciidae", "Aotidae", "Cebidae", "Atelidae")

Pr.genus <-  c("Arctocebus_calabarensis", "Arctocebus_aureus", "Perodicticus_potto_ibeanus", "Perodicticus_potto_subsp._potto", "Perodicticus_potto_edwarsi", "Loris_tardigradus_nordicus")

Primate_genus_nodes_1 <- sapply(Pr.genus, grep, x= PrimateFamily_Tree$tip.label)
Primate_genus_nodes_1
plotTree(PrimateFamily_Tree, type= "fan", fsize = 0.6, lwd=1, ftype="i")
add.arrow(PrimateFamily_Tree, tip= Primate_genus_nodes_1, arr1= 0.15, col="red")

xx= getMRCA(PrimateFamily_Tree, Primate_genus_nodes_1)
xx

Pr.clade <- extract.clade(PrimateFamily_Tree, xx)
Pr.clade
plot(Pr.clade, show.node.label = TRUE, use.edge.length = FALSE)

write.tree(Pr.clade, file= "example.Newick")









class(Primate_genus_nodes_1)
class(PrimateFamily_Tree)
class(Pr.genus)





########### Primates_Species
Primates_string <- read.tree("/Users/adesi/Downloads/primate_species.nwk") 
plot(Primates_string, no.margin = TRUE, font = 2, tip.color = "Blue")
Primates_string
#plotTree(Primates.string, fsize= 0.4, lwd = 2)


#Let's examine the components of the tree
str(Primates_string)
Primates_string$edge
Primates_string$tip.label
Primates_string$edge.length


# Let's specify a list of specific Species
Pr.species <- c("Sapajus_apella", "Callithrix_jacchus", "Aotus_nancymaae", "Saimiri_boliviensis", "Pongo_abelii", "Pongo_pygmaeus")
Pr.species
Primate_nodes_1 <- sapply(Pr.species, grep, x= Primates_string$tip.label)
Primate_nodes_1


#view these species indicated by arrows
plotTree(Primates.string, type= "fan", fsize=0.6, lwd=1, ftype="i")
add.arrow(Primates.string, tip= Primate_nodes_1, arr1= 0.15, col="red")
nodes

# Drop these species from the tree
Primates.woNodes <- drop.tip(Primates_string, Pr.species)
Primates.woNodes

#Extract only these species
Primate_nodes_1a <- getMRCA(Primates_string, Pr.species[
  -which(Pr.species%in%c("Callithrix_jacchus"))])

class(Primates_string)
class(PrimateFamily_Tree)
class(Primate_nodes_1)
class(Pr.species)

vv= getMRCA(Primates_string, Primate_nodes_1)
Primate_nodes_1

Primates_string$tip.label
Pr.species



