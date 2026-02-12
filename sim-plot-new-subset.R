library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(viridis)

pdf = data.frame(mean=rep(2:10, 10), size=rep(1:10, each=9))
pdf$prob = (1-dpois(0, pdf$mean))**pdf$size
ggplot(pdf, aes(x=mean, y=size, fill=prob)) + geom_tile()

# number of dnas in my compartment is m = 1 + Poisson(0.5)
# if i have m dnas i need at least m copies of every subunit = (1-P(k < m))**n
df.res = data.frame()
for(meanprot in 1:10) {
  for(meandna in (1:4)/4) {
    for(n in c(1,5,10)) {
      p = 0
      for(m in 1:10) {
        x = 1-sum(dpois(0:(m-1), meanprot))
        p = p + dpois(m-1, meandna) * x**n
      }
      p1 = (1-dpois(0, meandna))*(1-dpois(0, meanprot))**n
      df.res = rbind(df.res, data.frame(meanprot = meanprot, meandna = meandna, n=n, p=p, p1=p1))
    }
  }
}
df.res[df.res$meanprot==3 & df.res$meandna==0.5 & df.res$n == 5,]
ggarrange(
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    labs(x = "Mean proteins per mito", y = "Mean DNA per mito", fill = "MtDNAs\nin\nDNA-\ncomplex"),
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p1)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    labs(x = "Mean proteins per mito", y = "Mean DNA per mito", fill = "MtDNAs\nwith\npartner\nand\ncomplex"),
  nrow=2
)

df = read.csv("sim-out-new-5.csv")

# target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
# swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA

target.labels = c("Random", "DNA-bearing", "No DNA-complex", "DNA-bearing,\nno DNA-complex", "Random by\nbatch")
swap.labels = c("None", "Subunits", "DNA + Subunits", "DNA-Complexes", "Complexes", "DNA-complexes +\nComplexes")
social.labels = c("DNA-bearing", "All mitos")
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

#### should demo that t=900 and t=1000 are comparbale

df = df[df$poisson == 1 & df$t == 1000,]
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

if(FALSE) {
  df_sub <- df[df$swap == "None" & df$social == "All mitos",]
  ggplot(df_sub, aes(x = target, y = completes, fill=factor(t))) + 
    geom_boxplot() + theme(axis.text.x = element_text(angle=90)) + 
    facet_grid(poisson ~ expression)
  
  ggplot(df, aes(x = target, y = completes, group=factor(t),
                 shape = factor(poisson), color=factor(meanprot))) + 
    geom_box() + theme(axis.text.x = element_text(angle=90))
  
  ggplot(df_mean, aes(x = target, y = completes, 
                      shape = factor(poisson), color=factor(meanprot))) + 
    geom_point() + theme(axis.text.x = element_text(angle=90))
}

df_sub <- df[df$swap == "None" & df$social == "All mitos",]
ggarrange(
  ggplot(df_sub, aes(x = target, y = completes, color=factor(meanprot))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule (no exchange)", y = "MtDNAs in nucleoprotein", 
         color = "Mean\nproteins\nper\nmito") +
    ylim(0,100),
  ggplot(df_sub, aes(x = target, y = completes2, color=factor(meanprot))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule", y = "MtDNAs with partners and nucleoprotein", 
         color = "Mean\nproteins\nper\nmito") + 
    ylim(0,100),
  ncol = 2
)

ggplot(df_sub, aes(x = target, y = avprotempty, color=factor(meanprot))) + geom_boxplot() +
  geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
  labs(x = "Targetting rule (no exchange)", y = "Mean proteins in DNA-empty mitos", color = "Mean\nproteins\nper\nmito")
ggplot(df_sub, aes(x = target, y = avprotfull, color=factor(meanprot))) + geom_boxplot() +
  geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
  labs(x = "Targetting rule (no exchange)", y = "Mean proteins in DNA-bearing mitos", color = "Mean\nproteins\nper\nmito")


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

df_sub <- df[df$target == "Random" & df$expression == 100,]
ggarrange(
  ggplot(df_sub, aes(x = swap, y = completes, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "MtDNAs in nucleoprotein", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100),
  ggplot(df_sub, aes(x = swap, y = completes2, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "MtDNAs with partners and nucleoprotein", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100),
  ncol=2
)

ggarrange(
  ggplot(df_sub, aes(x = swap, y = avprotempty, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "Mean proteins in mtDNA-empty mitos", 
         color = "Mitos\nallowed\nto fuse") ,
  ggplot(df_sub, aes(x = swap, y = avprotfull, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "Mean proteins in mtDNA-bearing mitos", 
         color = "Mitos\nallowed\nto fuse") ,
  ggplot(df_sub, aes(x = swap, y = maxdna, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "Max DNAs in a mito", 
         color = "Mitos\nallowed\nto fuse") ,
  ggplot(df_sub, aes(x = swap, y = onedna, color=factor(social))) + geom_boxplot() +
    geom_beeswarm(dodge.width=0.75) +   theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule (random targetting)", y = "Proportion of mtDNAs alone in mito", 
         color = "Mitos\nallowed\nto fuse") ,
  ncol=2,nrow=2
)

# next question -- targetting vs swapping
df_sub = df[df$expression == 100,]
df_mean <- df_sub %>%
  group_by(across(-c(rep, completes, completes2, avprot, avdna))) %>%   # group by all other columns
  summarise(
    completes = mean(completes, na.rm = TRUE),
    completes2 = mean(completes2, na.rm = TRUE),
    avprot = mean(avprot, na.rm=TRUE),
    avdna = mean(avdna, na.rm=TRUE),
    .groups = "drop"
  )
ggarrange(
  ggplot(df_mean, aes(x=target, y=swap, fill=completes)) + 
    geom_tile() + facet_wrap(~ "Fusion rule:"+social) + 
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Targetting rule", y = "Exchange rule", fill = "MtDNAs\nwith\nnucleo-\nprotein"),
  
  ggplot(df_mean, aes(x=target, y=swap, fill=completes2)) + 
    geom_tile() + facet_wrap(~ "Fusion rule:"+social) + 
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Targetting rule", y = "Exchange rule", fill = "MtDNAs\nwith\npartner\nand\nnucleo-\nprotein"),
  nrow = 2)

my.aov = aov(completes ~ target *swap, data=df_sub)
summary(my.aov)
interaction.plot(df_sub$target, df_sub$swap, df_sub$completes)
my.aov = aov(completes ~ target + swap, data=df_sub)
TukeyHSD(my.aov)

my.aov = aov(completes2 ~ target *swap, data=df_sub)
summary(my.aov)
interaction.plot(df_sub$target, df_sub$swap, df_sub$completes2)
my.aov = aov(completes2 ~ target + swap, data=df_sub)
TukeyHSD(my.aov)

ggplot(df_mean, aes(x=completes, y=completes2, 
                    shape=factor(target), color=factor(swap))) + 
  geom_point() + facet_wrap(~ "Fusion rule:"+social) + 
  labs(x = "MtDNAs in DNA-complex", y = "MtDNAs with partner and complex",
       shape = "Targetting\nrule", color = "Exchange\nrule")
