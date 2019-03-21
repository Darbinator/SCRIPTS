library(ggplot2)
library(tidyverse)
library(scales)

dir = setwd("~/Desktop/")

nbReads = read.csv("~/Desktop/tableaux/raw_nbReads.txt", sep="", head=T)

point <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)
ReadsGr <- ggplot(nbReads, aes(Echantillon,Nombre, fill=Condition)) + geom_text(aes(label = nbReads$Echantillon),position = position_stack(vjust = 1.05)) + geom_bar(stat='identity') + scale_y_continuous(labels = point,breaks = seq(0, 100000000, by = 10000000)) + ylab("Nombre") + ggtitle("Nombre de reads initial")


AS_events = read.csv("~/Desktop/tableaux/AS_Events.txt", sep="\t", head=T)

p <- ggplot(AS_events, aes(Comparaison,Nombre, fill=Evenement, label = Label)) +  geom_bar(stat='identity',position='fill') + geom_text(size=3,position = position_fill(vjust = 0.5)) 

p + geom_text(aes(label = AS_events$Pourcentage), position = 'fill',vjust = 3)
qplot(AS_events$Comparaison, AS_events$Nombre, fill = AS_events$Evenement)
p + geom_text(aes(label = AS_events$Nombre), position = 'stack')

ggsave("~/Desktop/Images/nb_AS_events.png", width = 5, height = 5)
