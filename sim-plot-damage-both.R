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

nsubs = 5

df.5 = read.csv(paste0("sim-damage-5.csv", collapse=""))
df.1 = read.csv(paste0("sim-damage-1.csv", collapse=""))

df.5$expt = 5
df.1$expt = 1

# XXX TO DO merge dfs from 1 and 5, look over first-step delays, just think about min copy number

df = rbind(df.5[df.5$expression==min(df.5$expression),], 
           df.1[df.1$expression==min(df.1$expression),])
df$meanprot[df$expt==1] = round(mean(df$avprot[df$expt == 1]))
df$meanprot[df$expt==5] = round(mean(df$avprot[df$expt == 5]))

# target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
# swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA

target.labels = c("Random", "DNA-bearing", "No DNA-complex", "DNA-bearing,\nno DNA-complex", "Random by\nbatch", "Damaged")
swap.labels = c("None", "Subunits", "DNA + Subunits", "DNA-Complexes", "Complexes", "DNA-complexes +\nComplexes")
social.labels = c("DNA-bearing", "All mitos")
df$target = target.labels[df$target+1]
df$swap = swap.labels[df$swap+1]
df$social = social.labels[df$social+1]
if(FALSE){
  mean.min = mean(df$avprot[df$expression==min(df$expression)])
mean.max = mean(df$avprot[df$expression==max(df$expression)])
df$meanprot = 0
df$meanprot[df$expression==min(df$expression)] = round(mean.min, digits=0)
df$meanprot[df$expression==max(df$expression)] = round(mean.max, digits=0)
hist(df$avprot[df$expression==min(df$expression)])
hist(df$avprot[df$expression==max(df$expression)])
}
#df = df[df$expression == 100 & df$import == 100 & df$t == 100,]

#### should demo that t=900 and t=1000 are comparbale
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

df = df[df$t == 1000 & df$mut == 0.04,]
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

df_sub <- df[df$swap == "None" & df$social == "All mitos" & df$scanrate==1,]
fig.2 = ggarrange(
  ggplot(df_sub, aes(x = target, y = completes, color=factor(meanprot))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    #geom_jitter(position="dodge",size=bs.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "MtDNAs with\nprotein", 
         color = "Mean\nproteins\nper\nmito") +
    ylim(0,100),
  ggplot(df_sub, aes(x = target, y = completes2, color=factor(meanprot))) + 
    geom_boxplot() +
   # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "MtDNAs with\npartner and\nprotein", 
         color = "Mean\nproteins\nper\nmito") + 
    ylim(0,100),
  ncol = 2,
  labels=c("A", "B")
)
fig.2

fig.2.dmg = 
  ggplot(df_sub, aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
   # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
  theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "Proportion damaged mtDNA", 
         color = "Mean\nproteins\nper\nmito") +
     facet_wrap(~scanrate)
fig.2.dmg

fig.2.new = 
  ggarrange(
    ggplot(df_sub[df_sub$expt == 5,], 
           aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
      # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
      theme(axis.text.x = element_text(angle=90), legend.position="none") +
      labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
           color = "Mean\nproteins\nper\nmito") ,
    ggarrange(
      ggplot(df_sub[df_sub$expt == 1,], 
             aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
        # geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
        theme(axis.text.x = element_text(angle=90), legend.position="none") +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
             color = "Mean\nproteins\nper\nmito") 
  ),
  labels=c("A", "B")
  )
fig.2.new

ggplot(df_sub, aes(x = target, y = avprotempty, color=factor(meanprot))) + geom_boxplot() +
  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
  labs(x = "Targetting rule (no exchange)", y = "Mean proteins in DNA-empty mitos", color = "Mean\nproteins\nper\nmito")
ggplot(df_sub, aes(x = target, y = avprotfull, color=factor(meanprot))) + geom_boxplot() +
  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
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

df_sub <- df[df$target == "Random" & df$scanrate ==1,]
fig.3 = ggarrange(
  ggplot(df_sub, aes(x = swap, y = completes, color=factor(social))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "MtDNAs in\nDNA-complex", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100),
  ggplot(df_sub, aes(x = swap, y = completes2, color=factor(social))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +  
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "MtDNAs with\npartners and\ncomplex", 
         color = "Mitos\nallowed\nto fuse") + ylim(0,100),
  ncol=2,
  labels=c("A", "B")
)
fig.3

fig.3.dmg = 
  ggplot(df_sub, aes(x = swap, y = avdamage, color=factor(social))) + 
  geom_boxplot() +
  #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
  theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion damaged mtDNA", 
         color = "Mitos\nallowed\nto fuse") +  facet_grid(~mut)
  fig.3.dmg
  
  fig.3.new = ggarrange(
    ggplot(df_sub[df_sub$expt == 5,], aes(x = swap, y = avdamage, color=factor(social))) + 
      geom_boxplot() +
      #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
      theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA", 
           color = "Mitos\nallowed\nto fuse") ,
    ggplot(df_sub[df_sub$expt == 1,], aes(x = swap, y = avdamage, color=factor(social))) + 
      geom_boxplot() +
      #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
      theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA", 
           color = "Mitos\nallowed\nto fuse") ,
    labels=c("A", "B")
  )
  fig.3.new
  
if(nsubs > 1) {
  fig.s1 = ggarrange(
    ggplot(df_sub, aes(x = swap, y = maxdna, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "Max DNAs in a mito", 
           color = "Mitos\nallowed\nto fuse") ,
    ggplot(df_sub, aes(x = swap, y = avprotempty, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "Mean proteins in\nmtDNA-empty mitos", 
           color = "Mitos\nallowed\nto fuse") ,
    #ggplot(df_sub, aes(x = swap, y = avprotfull, color=factor(social))) + geom_boxplot() +
    #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
    #  labs(x = "Exchange rule (random targetting)", y = "Mean proteins in\nmtDNA-bearing mitos", 
    #       color = "Mitos\nallowed\nto fuse") ,
    ggplot(df_sub, aes(x = swap, y = onedna, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "Proportion of mtDNAs\nalone in mito", 
           color = "Mitos\nallowed\nto fuse") ,
    ggplot(df_sub, aes(x = swap, y = freesubs, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "Mean free subunits\nper mito", 
           color = "Mitos\nallowed\nto fuse") ,
    ncol=2,nrow=2,
    labels=c("A", "B", "C", "D")
  )
  
  
  ggarrange(
    ggplot(df_sub, aes(x = swap, y = completes, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "MtDNAs in DNA-complex", 
           color = "Mitos\nallowed\nto fuse") + ylim(0,100),
    ggplot(df_sub, aes(x = swap, y = freesubs, color=factor(social))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) +   theme(axis.text.x = element_text(angle=90)) +
      labs(x = "Exchange rule (random targetting)", y = "Mean free subunits per mito", 
           color = "Mitos\nallowed\nto fuse")
  )
}


# next question -- targetting vs swapping
  df_sub_strip <- df[df$scanrate ==1,c("expt", "target", "swap", "social", "avdamage")]
df_mean <- df_sub_strip %>%
  group_by(across(-c(avdamage))) %>%   # group by all other columns
  summarise(
 #   completes = mean(completes, na.rm = TRUE),
#    completes2 = mean(completes2, na.rm = TRUE),
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
 
fig.s3 = 
  ggplot(df_mean, aes(x=target, y=swap, fill=100-avdamage)) + 
    geom_tile() + facet_wrap(~ "Fusion rule:"+social) + 
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Targetting rule", y = "Exchange rule", fill = "Proportion\nintact\nDNA") +
  facet_wrap(~ expt)
fig.s3

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

fig.4 = ggplot(df_mean, aes(x=completes, y=completes2, 
                            shape=factor(target), color=factor(swap))) + 
  geom_point() + facet_wrap(~ "Fusion rule:"+social) + 
  labs(x = "MtDNAs in DNA-complex", y = "MtDNAs with partner and complex",
       shape = "Targetting\nrule", color = "Exchange\nrule") +
  facet_wrap(~ expt)
fig.4

sf = 2

  png(paste0("fig-s1-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2)
dev.off()
png(paste0("fig-s2-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.3)
dev.off()
png(paste0("fig-s0-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.s0)
dev.off()

if(FALSE)
{
png(paste0("fig-4-delay-", nsubs, ".png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.4)
dev.off()

png(paste0("fig-s1-delay-", nsubs, ".png", collapse=""), width=600*sf, height=500*sf, res=72*sf)
print(fig.s1)
dev.off()
png(paste0("fig-s2-delay-", nsubs, ".png", collapse=""), width=600*sf, height=500*sf, res=72*sf)
print(fig.s2)
dev.off()
png(paste0("fig-s3-delay-", nsubs, ".png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.s3)
dev.off()
png(paste0("fig-2-delay-dmg-", nsubs, ".png", collapse=""), width=300*sf, height=250*sf, res=72*sf)
print(fig.2.dmg)
dev.off()
png(paste0("fig-3-delay-dmg-", nsubs, ".png", collapse=""), width=300*sf, height=250*sf, res=72*sf)
print(fig.3.dmg)
dev.off()
png(paste0("fig-all-delay-", nsubs, ".png", collapse=""), width=600*sf, height=750*sf, res=72*sf)
print(ggarrange(fig.2 + theme(plot.margin = margin(20, 5, 5, 5)), 
                fig.3 + theme(plot.margin = margin(20, 5, 5, 5)), 
                fig.4 + theme(plot.margin = margin(20, 5, 5, 5)), 
                labels=c("i", "ii", "iii"), 
                label.x = 0.02,
                label.y = 1.0,
                hjust = 0,
                vjust = 1,
                nrow=3))
dev.off()
}

png(paste0("fig-2-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2.new)
dev.off()
png(paste0("fig-3-new.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.3.new)
dev.off()
png(paste0("fig-4-new.png", collapse=""), width=650*sf, height=450*sf, res=72*sf)
print(fig.4.new)
dev.off()
