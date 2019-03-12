# -----------------------------------------------------------------------------
# C. Elegans Connectome Investigation
# -----------------------------------------------------------------------------
# 10 March 2019
# 1. Reviewing exercises by Sergio Peignier.
# -----------------------------------------------------------------------------

library(here)           # here
library(tidyverse)      # readr --> read_csv().

neurons <- read_csv(file=here("data", "peignier", "connectome.csv"))
muscles <- read_csv(file=here("data", "peignier", "neurons_to_muscles.csv"))
sensory <- read_csv(file=here("data", "peignier", "sensory.csv"))

cells <- muscles %>% transmute(CELL=Muscle, TYPE="MUSCLE")
cells <- cells %>% distinct(CELL, TYPE)
sensors <- sensory %>% transmute(CELL=Neuron, TYPE=Function)

cells <- full_join(cells,sensors)
cells <- full_join(neurons %>% transmute(CELL=Neuron,TYPE="*Neuron"),cells)
cells <- cells %>% distinct(CELL, TYPE)
cells <- cells %>% arrange(CELL, TYPE) %>% group_by(CELL) %>% filter(row_number()==n()) %>% ungroup()

cells$TYPE[cells$TYPE=="*Neuron"]="neuron"

connections <- neurons %>% spread(Target, `Number of Connections`)

aggregate(x=connections, by=list(Neuron=connections$Neuron), FUN=sum)

c2 <- connections %>% select(-X1)
c2[is.na(c2)]<-0
c3 <- aggregate(c2[,3:287], c2[,1:2], FUN=sum)

neuro <- neurons %>% select(-X1,-Neurotransmitter)
test <- neuro %>% filter(Neuron=="ADAL"|Neuron=="ADAR")

c4<-full_join(cells, c3, by=c("CELL"="Neuron"))



# A bit of fiddling with igraph.
library(igraph)
# Undirected star graph with 10 nodes.
g <- graph.star(n=10, mode = "undirected")
plot(g)
# Add edges between selected pairs of nodes.
g <- add_edges(g, c(8,5, 6,3, 6,10, 5,2, 4,8, 2,7, 6,5, 7,9))
plot(g)
# Replot as a circle with nodes ordered by vertex number.
l <- layout_in_circle(g, order = V(g))
plot(g, layout = l)
# Verify that the graph is connected or bipartite.
is.connected(g)
is.bipartite(g)
avg_degree <- 2*ecount(g)/vcount(g)
diameter(g)
sapply(maximal.cliques(g), length)
# New plot with node size proportional to closeness.
plot(g, vertex.size = closeness(g)*500)
# Colour the nodes.
V(g)$color <- ifelse(V(g)%%2 == 0, "blue", "red")
plot(g, vertex.size = closeness(g)*500)
# A bipartite graph.
matches <- data.frame(name = rep(c("Jerry", "Lilly", "Karl", "Jenny"), each = 4),
                      subject = rep(c("Maths", "English", "Biology", "French"), 4),
                      weight = c(81, 78, 24, 58, 76, 60, 62, 83, 35, 59, 50, 56, 72, 90, 88, 86))
g <- graph_from_data_frame(matches, directed = FALSE)

V(g)$type <- V(g)$name %in% matches[,1]
g %>%
    add_layout_(as_bipartite()) %>%
    plot()
# Add some colour
V(g)$color <- ifelse(V(g)$type==TRUE, "red", "blue")
plot(g)
#
max_bipartite_match(g, weights = g$weight)
# Add some vertices and edges ...
g <- add_vertices(g, 3, name = c("Becky", "Ben", "Mark"))
g <- add_edges(g, c("Becky", "Maths", "Ben", "French", "Mark", "Biology"))
# Spectra
spectra <- spectrum(g)['vectors']
plot(x = 1:length(spectrum(g)$vectors), y = spectrum(g)$vectors, xlab = "Vector Index", ylab = "Value")
# Largest independent set of nodes.
independence.number(g)
#
min_cut(g)
#
edge_connectivity(g)


# For the C Elegans neurons.
n <- graph_from_data_frame(neuro, directed=FALSE)
V(n)$type <- V(n)$Neuron %in% neuro[,1]
n %>% add_layout_(as_bipartite()) %>% plot()

t <- graph_from_data_frame(test, directed=FALSE)
V(t)$type <- V(t)$name %in% test[,1]
t %>% add_layout_(as_bipartite()) %>% plot()
l <- layout_in_circle(t,order=V(t))
plot(t,layout=layout_nicely(t),edge.width=E(t)$`Number Of Connections`)
V(t)$color <- ifelse(V(t)%%2 == 0, "green", "pink")

# For petri nets
library(petrinetR)
