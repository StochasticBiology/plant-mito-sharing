library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(viridis)
library(ggpubr)

sf = 2
bs.size = 0.25
dw.size = 0.5

pdf = data.frame(mean=rep(2:10, 10), size=rep(1:10, each=9))
pdf$prob = (1-dpois(0, pdf$mean))**pdf$size
ggplot(pdf, aes(x=mean, y=size, fill=prob)) + geom_tile()

# number of dnas in my compartment is m = 1 + Poisson(0.5)
# if i have m dnas i need at least m copies of every subunit = (1-P(k < m))**n
df.res = data.frame()
for(meanprot in 1:10) {
  for(meandna in (1:6)/6) {
    for(n in c(1,5,9)) {
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
df.res[df.res$meanprot==4 & abs(df.res$meandna-0.33) < 0.05 & df.res$n == 5,]
fig.1 = ggarrange(
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    labs(x = "Mean copies per mito\nof each protein type", y = "Mean DNA per mito", fill = "MtDNAs\nin\nDNA-\ncomplex"),
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p1)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    labs(x = "Mean copies per mito\nof each protein type", y = "Mean DNA per mito", fill = "MtDNAs\nwith\npartner\nand\ncomplex"),
  labels=c("A", "B"),
  nrow=2
)

png("fig-1.png", width=600*sf, height=400*sf, res=72*sf)
print(fig.1)
dev.off()

###### now use the data from simulations

df.0 = read.csv(paste0("sim-damage-0.csv", collapse=""))
df.1 = read.csv(paste0("sim-damage-1.csv", collapse=""))
df.2 = read.csv(paste0("sim-damage-2.csv", collapse=""))

df.0$expt = "BER"
df.1$expt = "Template"
df.2$expt = "BER-complex"

df = rbind(df.0[df.0$expression==min(df.0$expression),], 
           df.1[df.1$expression==min(df.1$expression),],
           df.2[df.2$expression==min(df.2$expression),])
df$meanprot[df$expt=="Template"] = round(mean(df$avprot[df$expt == "Template"]))
df$meanprot[df$expt=="BER"] = round(mean(df$avprot[df$expt == "BER"]))
df$meanprot[df$expt=="BER-complex"] = round(mean(df$avprot[df$expt == "BER-complex"]))

# target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
# swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA

target.labels = c("Random", "DNA-bearing", "No DNA-complex", "DNA-bearing,\nno DNA-complex", "Random by\nbatch", "Damaged")
swap.labels = c("None", "Subunits", "DNA + Subunits", "DNA-Complexes", "Complexes", "DNA-complexes +\nComplexes")
social.labels = c("DNA-bearing", "All mitos")
df$target = target.labels[df$target+1]
df$swap = swap.labels[df$swap+1]
df$social = social.labels[df$social+1]

#### demo that t=900 and t=1000 are comparbale
set.1 = df[df$t==900,]
set.2 = df[df$t==1000,]
comp.df = data.frame(at.900 = set.1$completes, at.1000 = set.2$completes,
                     at.900.1 = set.1$completes2, at.1000.1 = set.2$completes2,
                     at.900.2 = set.1$avdamage, at.1000.2 = set.2$avdamage)
fig.s0 = ggarrange(
  ggplot(comp.df, aes(x=at.900, y=at.1000)) + 
    geom_hex(fill = "#0000FF", aes(alpha=after_stat(sqrt(count)))) + geom_abline() +
    scale_fill_continuous() + 
    labs(x= "mtDNA with complex\nat t=900", y="mtDNA with complex\nat t=1000",
         alpha="sqrt(\nnumber\nof sims)"),
  ggplot(comp.df, aes(x=at.900.1, y=at.1000.1)) + 
    geom_hex(fill = "#0000FF", aes(alpha=after_stat(sqrt(count)))) + geom_abline() +
    scale_fill_continuous() + 
    labs(x= "mtDNA with complex\nand partner\nat t=900", y="mtDNA with complex\nand partner\nat t=1000",
         alpha="sqrt(\nnumber\nof sims)"),
  ggplot(comp.df, aes(x=at.900.2, y=at.1000.2)) + 
    geom_hex(fill = "#0000FF", aes(alpha=after_stat(sqrt(count)))) + geom_abline() +
    scale_fill_continuous() + 
    labs(x= "Proportion\ndamaged mtDNA\nat t=900", y="Proportion\ndamaged mtDNA\nat t=1000",
         alpha="sqrt(\nnumber\nof sims)"),
  nrow = 1, labels=c("A", "B", "C")
)
fig.s0

# looking at dependence on mutation and scan rates
df_big = df[df$t == 1000, c("scanrate", "mut", "target", "swap", "social", "rep", "avdamage", "expt")]
df_mean <- df_big %>%
  group_by(across(-c(rep, avdamage))) %>%   # group by all other columns
  summarise(
    avdamage = mean(avdamage, na.rm=TRUE),
    .groups = "drop"
  )
fig.s3 = ggarrange(
  ggplot(df_mean[df_mean$expt=="BER" & df_mean$social=="DNA-bearing",], 
         aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
    geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
    theme(axis.text.x = element_text(angle=90)) +
    labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA") +
    facet_grid("scan="+scanrate ~ "mu="+mut),
ggplot(df_mean[df_mean$expt=="Template" & df_mean$social=="DNA-bearing",], 
       aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA") +
  facet_grid("scan="+scanrate ~ "mu="+mut),
ggplot(df_mean[df_mean$expt=="BER-complex" & df_mean$social=="DNA-bearing",], 
       aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA") +
  facet_grid("scan="+scanrate ~ "mu="+mut),
nrow=3,
labels=c("A", "B", "C")
)
fig.s3

######## from now on, zoom in to a particular mutation and scan rate
df = df[df$t == 1000 & df$mut == 0.04 & df$scanrate == 1,]

# switch off exchange and ask questions about targetting
df_sub <- df[df$swap == "None" & df$social == "All mitos",]
fig.2 = ggarrange(
  ggplot(df_sub, aes(x = target, y = completes, color=factor(meanprot))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    #geom_jitter(position="dodge",size=bs.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "MtDNAs with\nprotein", 
         color = "Mean\nproteins\nper\nmito") +
    ylim(0,100) + facet_wrap(~expt),
  ggplot(df_sub, aes(x = target, y = completes2, color=factor(meanprot))) + 
    geom_boxplot() +
    # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "MtDNAs with\npartner and\nprotein", 
         color = "Mean\nproteins\nper\nmito") + 
    ylim(0,100) + facet_wrap(~expt),
  ncol = 2,
  labels=c("A", "B")
)
fig.2

fig.2.new = 
  ggarrange(
    ggplot(df_sub[df_sub$expt == "BER",], 
           aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
      # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
      theme(axis.text.x = element_text(angle=90), legend.position="none") +
      labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
           color = "Mean\nproteins\nper\nmito") ,
    ggarrange(
      ggplot(df_sub[df_sub$expt == "Template",], 
             aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
        # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
        theme(axis.text.x = element_text(angle=90), legend.position="none") +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
             color = "Mean\nproteins\nper\nmito") 
    ),
    ggarrange(
      ggplot(df_sub[df_sub$expt == "BER-complex",], 
             aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
        # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
        theme(axis.text.x = element_text(angle=90), legend.position="none") +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
             color = "Mean\nproteins\nper\nmito") 
    ),
    labels=c("A", "B", "C"),
    nrow=1
  )
fig.2.new

# next question -- with random targetting, how does different exchange influence completeness 
df_sub <- df[df$target == "Random",]
fig.3 = ggarrange(
  ggplot(df_sub, aes(x = swap, y = completes, color=factor(social))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "MtDNAs in\nDNA-complex", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100) + facet_wrap(~expt),
  ggplot(df_sub, aes(x = swap, y = completes2, color=factor(social))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +  
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "MtDNAs with\npartners and\ncomplex", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100) + facet_wrap(~expt),
  nrow=2,
  labels=c("A", "B")
)
fig.3

fig.3.new = ggarrange(
  ggplot(df_sub[df_sub$expt == "BER",], aes(x = swap, y = avdamage, color=factor(social))) + 
    geom_boxplot() +
    #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA", 
         color = "Mitos\nallowed\nto fuse") ,
  ggplot(df_sub[df_sub$expt == "Template",], aes(x = swap, y = avdamage, color=factor(social))) + 
    geom_boxplot() +
    #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA", 
         color = "Mitos\nallowed\nto fuse") ,
  ggplot(df_sub[df_sub$expt == "BER-complex",], aes(x = swap, y = avdamage, color=factor(social))) + 
    geom_boxplot() +
    #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA", 
         color = "Mitos\nallowed\nto fuse") ,
  labels=c("A", "B", "C"),
  nrow=1
)
fig.3.new


# next question -- targetting vs swapping
df_sub_strip <- df[df$scanrate ==1,c("expt", "target", "swap", "social", "avdamage")]
df_mean <- df_sub_strip %>%
  group_by(across(-c(avdamage))) %>%   # group by all other columns
  summarise(
    avdamage = mean(avdamage, na.rm=TRUE),
    .groups = "drop"
  )

fig.4.new = 
  ggplot(df_mean, aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(color="#FFFFFF44") + facet_grid(expt ~ "Fusion rule:"+social) + 
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Targetting rule", y = "Exchange rule", fill = "Proportion\ndamaged\nmtDNA")
fig.4.new


df_sub = df_sub_strip
my.aov = aov(avdamage ~ target *swap, data=df_sub)

#my.aov = aov(completes ~ target *swap, data=df_sub)
summary(my.aov)
interaction.plot(df_sub$target, df_sub$swap, df_sub$avdamage)
my.aov = aov(completes ~ target + swap, data=df_sub)
TukeyHSD(my.aov)

my.aov = aov(completes2 ~ target *swap, data=df_sub)
summary(my.aov)
interaction.plot(df_sub$target, df_sub$swap, df_sub$completes2)
my.aov = aov(completes2 ~ target + swap, data=df_sub)
TukeyHSD(my.aov)

sf = 2

png(paste0("fig-s0-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.s0)
dev.off()
png(paste0("fig-s1-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2)
dev.off()
png(paste0("fig-s2-new.png", collapse=""), width=600*sf, height=750*sf, res=72*sf)
print(fig.3)
dev.off()
png(paste0("fig-s3-new.png", collapse=""), width=600*sf, height=750*sf, res=72*sf)
print(fig.s3)
dev.off()


png(paste0("fig-2-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2.new)
dev.off()
png(paste0("fig-3-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.3.new)
dev.off()
png(paste0("fig-4-new.png", collapse=""), width=650*sf, height=450*sf, res=72*sf)
print(fig.4.new)
dev.off()

