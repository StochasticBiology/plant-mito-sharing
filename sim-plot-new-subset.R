library(ggplot2)
library(dplyr)

pdf = data.frame(mean=rep(2:10, 10), size=rep(1:10, each=9))
pdf$prob = (1-dpois(0, pdf$mean))**pdf$size
ggplot(pdf, aes(x=mean, y=size, fill=prob)) + geom_tile()

df = read.csv("sim-out-new.csv")

# target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
# swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA

target.labels = c("Random", "DNA-bearing", "No nucleoprotein", "DNA-bearing,\nno nucleoprotein", "Random by\nbatch")
swap.labels = c("None", "Subunits", "Subunits+DNA", "Nucleoproteins", "Complexes")
social.labels = c("All mitos", "DNA-bearing")
df$target = target.labels[df$target+1]
df$swap = swap.labels[df$swap+1]
df$social = social.labels[df$social+1]
mean.100 = mean(df$avprot[df$expression==100])
mean.1000 = mean(df$avprot[df$expression==1000])
df$meanprot = 0
df$meanprot[df$expression==100] = round(mean.100, digits=0)
df$meanprot[df$expression==1000] = round(mean.1000, digits=0)
hist(df$avprot[df$expression==1000])
hist(df$avprot[df$expression==100])
hist(df$avdna[df$expression==100])
#df = df[df$expression == 100 & df$import == 100 & df$t == 100,]

# next question -- without fusion, how does targetting influence completeness 
df_mean <- df[df$swap == "None" & df$social == "All mitos",] %>%
  group_by(across(-c(rep, completes, completes2, 
                     avprot, avdna, avprotempty, avprotfull))) %>%   # group by all other columns
  summarise(
    completes = mean(completes, na.rm = TRUE),
    completes2 = mean(completes2, na.rm = TRUE),
    avprot = mean(avprot, na.rm=TRUE),
    avdna = mean(avdna, na.rm=TRUE),
    avprotempty = mean(avprotempty, na.rm=TRUE),
    avprotfull = mean(avprotfull, na.rm=TRUE),
    .groups = "drop"
  )
ggplot(df_mean, aes(x = target, y = completes, color=factor(meanprot))) + 
  geom_point() + theme(axis.text.x = element_text(angle=90))
ggplot(df_mean, aes(x = target, y = avprotempty, color=factor(meanprot))) + 
  geom_point() + theme(axis.text.x = element_text(angle=90))
ggplot(df_mean, aes(x = target, y = avprotfull, color=factor(meanprot))) + 
  geom_point() + theme(axis.text.x = element_text(angle=90))
# ^ targetting 1 and 3 usually best -- only those with DNA
# this one should still be zero as DNAs can never meet
ggplot(df_mean, aes(x = target, y = completes2, color=factor(expression))) + 
  geom_point() 

# next question -- with random targetting, how does different exchange influence completeness 
df_mean <- df[df$target == "Random" & df$nfuse == 100,] %>%
  group_by(across(-c(rep, completes, completes2, 
                     avprot, avdna, avprotempty, avprotfull))) %>%   # group by all other columns
  summarise(
    completes = mean(completes, na.rm = TRUE),
    completes2 = mean(completes2, na.rm = TRUE),
    avprot = mean(avprot, na.rm=TRUE),
    avdna = mean(avdna, na.rm=TRUE),
    avprotempty = mean(avprotempty, na.rm=TRUE),
    avprotfull = mean(avprotfull, na.rm=TRUE),
    .groups = "drop"
  )
ggplot(df_mean, aes(x = swap, y = completes, color=factor(meanprot))) + 
  geom_point() + facet_wrap(~ social) + theme(axis.text.x = element_text(angle=90))
ggplot(df_mean, aes(x = swap, y = completes2, color=factor(meanprot))) + 
  geom_point() + facet_wrap(~ social) + theme(axis.text.x = element_text(angle=90))
# ^ swapping only between mtDNA-holding mitos helps a lot


# next question -- targetting vs swapping
df_mean <- df[df$nfuse == 100,] %>%
  group_by(across(-c(rep, completes, completes2, avprot, avdna))) %>%   # group by all other columns
  summarise(
    completes = mean(completes, na.rm = TRUE),
    completes2 = mean(completes2, na.rm = TRUE),
    avprot = mean(avprot, na.rm=TRUE),
    avdna = mean(avdna, na.rm=TRUE),
    .groups = "drop"
  )
ggplot(df_mean, aes(x=target, y=swap, fill=completes)) + 
  geom_tile() + facet_grid(expression ~ social) + theme(axis.text.x = element_text(angle = 90)) 
ggplot(df_mean, aes(x=target, y=swap, fill=completes2)) + 
  geom_tile() + facet_grid(expression ~ social) + theme(axis.text.x = element_text(angle = 90)) 

my.aov = aov(completes ~ target +swap, data=df_mean[df_mean$expression==100,])
summary(my.aov)
my.aov = aov(completes2 ~ target +swap, data=df_mean[df_mean$expression==100,])
summary(my.aov)

TukeyHSD(my.aov)
interaction.plot(df$target, df$swap, df$completes)

ggplot(df_mean, aes(x = import, y = completes, fill = factor(target), color = factor(swap))) +
  geom_line() + facet_wrap(~ "Social"+social)
ggplot(df_mean, aes(x = import, y = completes2, fill = factor(target), color = factor(swap))) +
  geom_line() + facet_wrap(~ "Social"+social)
# ^ swapping only between mtDNA-holding mitos helps a lot

