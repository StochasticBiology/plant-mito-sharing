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

if(FALSE) {
  # previous model, where proteins get "used up" by mtDNAs
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
}
# number of dnas in my compartment is m = 1 + Poisson(0.5)
# i need at least m copies of every subunit = (1-P(k < m))**n. BER m=1 n=5, TR m=2 n=1
df.res = data.frame()
for(meanprot in 1:10) {
  for(meandna in (1:6)/6) {
    # n = number of protein types
    for(n in c(1,5,9)) {
      p = 0
      # m = count of protein
      # this 
      for(m in 1:10) {
        x = 1-sum(dpois(0:(m-1), meanprot))
        p = p + dpois(m-1, meandna) * x**n
      }
      # P(some DNA)*P(at least one copy of each n protein)
      p = (1-dpois(0:1, meanprot))**n
      p1 = (1-dpois(0, meandna))*(1-dpois(0, meanprot))**n
      df.res = rbind(df.res, data.frame(meanprot = meanprot, meandna = meandna, n=n, p=p, p1=p1))
    }
  }
}

df.res[df.res$meanprot==4 & abs(df.res$meandna-0.33) < 0.05 & df.res$n == 5,]
fig.1 = ggarrange(
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    theme_minimal() +
    labs(x = "Mean copies per mito\nof each protein type", y = "Mean DNA per mito", fill = "MtDNAs\nin\nDNA-\ncomplex"),
  ggplot(df.res, aes(x=meanprot, y=meandna, fill=p1)) + geom_tile() + 
    scale_fill_viridis() + facet_wrap(~ paste("n =",n)) + 
    theme_minimal() +
    labs(x = "Mean copies per mito\nof each protein type", y = "Mean DNA per mito", fill = "MtDNAs\nwith\npartner\nand\ndimer"),
  labels=c("A", "B"),
  nrow=2
)

png("fig-1.png", width=600*sf, height=400*sf, res=72*sf)
print(fig.1)
dev.off()

###### now use the data from simulations

if(FALSE) {
  # use this to look at behaviour in unbiologically low copy number circumstances
  df.0 = read.csv(paste0("sim-damage-0.csv", collapse=""))
  df.1 = read.csv(paste0("sim-damage-3.csv", collapse=""))
  df.2 = read.csv(paste0("sim-damage-4.csv", collapse=""))

  expt.labels = c("BER", "BER-small", "BER-small-complex")
} else {
  df.0 = read.csv(paste0("sim-damage-update-0-0.csv", collapse=""))
  df.1 = read.csv(paste0("sim-damage-update-5-0.csv", collapse=""))
  df.2 = read.csv(paste0("sim-damage-update-2-0.csv", collapse=""))
  
  expt.labels = c("BER", "Template", "BER-complex")
}
expt.set = c("BER", "Template")

df.0$expt = expt.labels[1]
df.1$expt = expt.labels[2]
df.2$expt = expt.labels[3]

df = rbind(df.0[df.0$expression==min(df.0$expression),], 
           df.1[df.1$expression==min(df.1$expression),],
           df.2[df.2$expression==min(df.2$expression),])
df$meanprot[df$expt==expt.labels[2]] = round(mean(df$avprot[df$expt == expt.labels[2]]))
df$meanprot[df$expt==expt.labels[1]] = round(mean(df$avprot[df$expt == expt.labels[1]]))
df$meanprot[df$expt==expt.labels[3]] = round(mean(df$avprot[df$expt == expt.labels[3]]))

# target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
# swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA

target.labels = c("Random", "DNA-bearing", "Incomplete", "DNA-bearing,\nincomplete", "Random by\nbatch", "Damaged")
swap.labels = c("None", "Subunits", "DNA + Subunits", "DNA-Metabolons", "Metabolons", "DNA-Metabolons +\nMetabolons")
social.labels = c("DNA-bearing", "All mitos")
df$target = target.labels[df$target+1]
df$swap = swap.labels[df$swap+1]
df$social = social.labels[df$social+1]

# exclude random by batch
df = df[df$target != "Random by\nbatch",]

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
fig.s3.list = list(
  ggplot(df_mean[df_mean$expt==expt.labels[1] & df_mean$social=="DNA-bearing",], 
         aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
    geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
    theme(axis.text.x = element_text(angle=90)) +
    labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA (%)") +
    facet_grid("scan="+scanrate ~ "mu="+mut),
ggplot(df_mean[df_mean$expt==expt.labels[2] & df_mean$social=="DNA-bearing",], 
       aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA (%)") +
  facet_grid("scan="+scanrate ~ "mu="+mut),
ggplot(df_mean[df_mean$expt==expt.labels[3] & df_mean$social=="DNA-bearing",], 
       aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(size=2,color="#FFFFFF88") + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(x="Targetting rule", y="Exchange rule\n(DNA-bearing mitos fuse)", fill="Proportion\ndamaged\nmtDNA (%)") +
  facet_grid("scan="+scanrate ~ "mu="+mut) )
fig.s3 = ggarrange(plotlist=fig.s3.list[1:2], nrow=2, labels=c("BER pathway", "TR pathway"), label.y = 1.05 )
#fig.s3 = ggarrange(plotlist=fig.s3.list[1:3], nrow=3, labels=c("A", "B", "C"))

fig.s3

######## from now on, zoom in to a particular mutation and scan rate
df = df[df$t == 1000 & df$mut %in% c(0.004, 0.04) & df$scanrate == 1,]

# switch off exchange and ask questions about targetting
df_sub <- df[df$swap == "None" & df$social == "All mitos" & df$expt %in% expt.set,]
fig.2 = ggarrange(
  ggplot(df_sub[df_sub$expt == "BER",], aes(x = target, y = completes, color=factor(meanprot))) + 
    geom_boxplot() +
    #geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    #geom_jitter(position="dodge",size=bs.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Targetting rule\n(no exchange)", y = "MtDNAs with\nprotein", 
         color = "Mean\nproteins\nper\nmito") +
    ylim(0,100) + facet_wrap(~expt),
  ggplot(df_sub[df_sub$expt == "Template",], aes(x = target, y = completes2, color=factor(meanprot))) + 
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

bs.size = 2
fig.2.new.list = 
  list(
    ggplot(df_sub[df_sub$expt == expt.labels[1],], 
           aes(x = target, y = avdamage, fill=target)) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) + #theme_minimal()+   
      theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
      scale_fill_viridis_d() +
      labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
           color = "Mean\nproteins\nper\nmito") ,
    ggarrange(
      ggplot(df_sub[df_sub$expt == expt.labels[2],], 
             aes(x = target, y = avdamage, fill=target)) + geom_boxplot() +
        geom_beeswarm(size=bs.size, dodge.width=dw.size) + #theme_minimal()+   
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
        scale_fill_viridis_d() +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
             color = "Mean\nproteins\nper\nmito")  
    ),
    ggarrange(
      ggplot(df_sub[df_sub$expt == expt.labels[3],], 
             aes(x = target, y = avdamage)) + geom_boxplot() +
        geom_beeswarm(size=bs.size, dodge.width=dw.size) + #theme_minimal() +   
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA", 
             color = "Mean\nproteins\nper\nmito") 
    ))

fig.2.new = ggarrange(plotlist=fig.2.new.list[1:2], nrow=1, 
                      labels=c(" BER pathway", "  TR pathway"))
#fig.2.new = ggarrange(plotlist=fig.2.new.list[1:3], nrow=1, labels=c("A", "B", "C"))
fig.2.new

bs.size = 2
fig.2.max.list = 
  list(
    ggplot(df_sub[df_sub$expt == expt.labels[1] & df_sub$target != "Incomplete",], 
           aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
      geom_beeswarm(size=bs.size, dodge.width=dw.size) + theme_minimal()+   
      theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
      #scale_fill_viridis_d() +
      annotate("text", x = 1, y = 10, label = "BER") +
      labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA (%)") ,
    ggarrange(
      ggplot(df_sub[df_sub$expt == expt.labels[2] & df_sub$target != "Incomplete",], 
             aes(x = target, y = avdamage, color=factor(meanprot))) + geom_boxplot() +
        geom_beeswarm(size=bs.size, dodge.width=dw.size) + theme_minimal()+   
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
        #scale_fill_viridis_d() +
        annotate("text", x = 1, y = 87, label = "TR") +
        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA (%)", 
             color = "Mean\nproteins\nper\nmito")  
    ),
    ggarrange(
      ggplot(df_sub[df_sub$expt == expt.labels[3],], 
             aes(x = target, y = avdamage)) + geom_boxplot() +
        geom_beeswarm(size=bs.size, dodge.width=dw.size) + #theme_minimal() +   
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +

        labs(x = "Targetting rule\n(no exchange)", y = "Proportion\ndamaged mtDNA (%)", 
             color = "Mean\nproteins\nper\nmito") 
    ))

fig.2.new.max = ggarrange(plotlist=fig.2.max.list[1:2], nrow=1, 
                          labels=c("A", "B"))
#fig.2.new = ggarrange(plotlist=fig.2.new.list[1:3], nrow=1, labels=c("A", "B", "C"))
fig.2.new.max

# why is "Damaged" worse than "DNA-bearing"?
# if we send unmatched proteins to a random mito, this change disappears -> perhaps just can't find damaged DNA?
# yes, confirmed

# next question -- with random targetting, how does different exchange influence completeness 
df_sub <- df[df$target == "Random" & df$expt %in% expt.set,]
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

bs.size = 1
dw.size = 0.75
fig.3.new.list = list(
  ggplot(df_sub[df_sub$expt == expt.labels[1],], aes(x = swap, y = avdamage, 
                                                     color=factor(social))) + 
    geom_boxplot() +
    #scale_fill_viridis_d() + 
    geom_beeswarm(size=bs.size, dodge.width = dw.size, cex=0.25, alpha=0.5) +   
    theme_minimal() + 
    theme(axis.text.x = element_text(angle=45,hjust=1)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA (%)", 
         color = "Mitos\nallowed\nto fuse") +
    annotate("text", x = 1, y = 9, label = "BER"),
  ggplot(df_sub[df_sub$expt == expt.labels[2],], aes(x = swap, y = avdamage, 
                                                     color=factor(social))) + 
    geom_boxplot() +
    geom_beeswarm(size=bs.size, dodge.width = dw.size, cex=0.25, alpha=0.5) +   
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45,hjust=1)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA (%)", 
         color = "Mitos\nallowed\nto fuse") +
    annotate("text", x = 1, y = 75, label = "TR"),
  ggplot(df_sub[df_sub$expt == expt.labels[3],], aes(x = swap, y = avdamage, 
                                                     color=factor(social))) + 
    geom_boxplot(color="black") +
    #scale_fill_viridis_d() + 
    #  geom_beeswarm(size=bs.size, dodge.width=dw.size) +   
    theme(axis.text.x = element_text(angle=90)) +
    labs(x = "Exchange rule\n(random targetting)", y = "Proportion\ndamaged mtDNA (%)", 
         color = "Mitos\nallowed\nto fuse"))

fig.3.new = ggarrange(plotlist=fig.3.new.list[1:2], nrow=2, 
                      labels=c("A", "B"))
fig.3.new
#fig.3.new = ggarrange(plotlist=fig.3.new.list[1:3], nrow=1, labels=c("A", "B", "C"))

fig.3.new


# next question -- targetting vs swapping
df_sub_strip <- df[df$scanrate ==1 & df$expt %in% expt.set,c("expt", "target", "swap", "social", "avdamage")]
df_mean <- df_sub_strip %>%
  group_by(across(-c(avdamage))) %>%   # group by all other columns
  summarise(
    avdamage = mean(avdamage, na.rm=TRUE),
    .groups = "drop"
  )

fig.4.new = 
  ggplot(df_mean[df_mean$target != "Incomplete",], aes(x=target, y=swap, label=round(avdamage, digits=1), fill=avdamage)) + 
  geom_tile() + geom_text(color="#FFFFFF") + facet_grid(expt ~ "Fusion rule:"+social) + 
  #scale_fill_viridis() +
  scale_fill_gradientn(
    colours = c("black", "blue", "yellow", "red"),
    values = c(0, 1, 50, 100) / 100
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Targetting rule", y = "Exchange rule", fill = "Proportion\ndamaged\nmtDNA (%)")
fig.4.new

df_sub = df_sub_strip[df_sub_strip$target != "Incomplete",]

my.aov = aov(avdamage ~ target *swap, data=df_sub[df_sub$expt=="BER",])
summary(my.aov)
my.aov = aov(avdamage ~ target *swap, data=df_sub[df_sub$expt=="Template",])
summary(my.aov)
my.aov = aov(avdamage ~ target *swap, data=df_sub)
summary(my.aov)
interaction.plot(df_sub$target, df_sub$swap, df_sub$avdamage)
#my.aov = aov(completes ~ target + swap, data=df_sub)
#TukeyHSD(my.aov)

#my.aov = aov(completes2 ~ target *swap, data=df_sub)
#summary(my.aov)
#interaction.plot(df_sub$target, df_sub$swap, df_sub$completes2)
#my.aov = aov(completes2 ~ target + swap, data=df_sub)
#TukeyHSD(my.aov)

sf = 2

png(paste0("fig-s0-new-v.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.s0)
dev.off()
png(paste0("fig-s1-new-v.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2)
dev.off()
png(paste0("fig-s2-new-v.png", collapse=""), width=600*sf, height=750*sf, res=72*sf)
print(fig.3)
dev.off()
png(paste0("fig-s3-new-v.png", collapse=""), width=600*sf, height=750*sf, res=72*sf)
print(fig.s3)
dev.off()

png(paste0("fig-2-new-v.png", collapse=""), width=600*sf, height=250*sf, res=72*sf)
print(fig.2.new)
dev.off()
png(paste0("fig-3-new-v.png", collapse=""), width=400*sf, height=400*sf, res=72*sf)
print(fig.3.new)
dev.off()
png(paste0("fig-4-new-v.png", collapse=""), width=650*sf, height=450*sf, res=72*sf)
print(fig.4.new)
dev.off()

ggsave("max-fig-1-v.svg", fig.1, width = 7, height=5)
ggsave("max-fig-2-v.svg", fig.2.new.max, width = 6, height=4)
ggsave("max-fig-3-v.svg", fig.3.new, width = 9, height = 4)
ggsave("max-fig-4-v.svg", fig.4.new, width=6.5, height = 5)

ggsave("max-fig-s0-v.svg", fig.s0, width = 10, height=5)
ggsave("max-fig-s1-v.svg", fig.2, width = 7, height=5)
ggsave("max-fig-s2-v.svg", fig.3, width = 7, height=8)
ggsave("max-fig-s3-v.svg", fig.s3, width = 7, height=9)

ggsave("max-fig-1-v.png", fig.1, width = 7, height=5)
ggsave("max-fig-2-v.png", fig.2.new.max, width = 6, height=4)
ggsave("max-fig-3-v.png", fig.3.new, width = 6, height = 6)
ggsave("max-fig-4-v.png", fig.4.new, width=6.5, height = 5)

ggsave("max-fig-s0-v.png", fig.s0, width = 10, height=5)
ggsave("max-fig-s1-v.png", fig.2, width = 7, height=5)
ggsave("max-fig-s2-v.png", fig.3, width = 7, height=8)
ggsave("max-fig-s3-v.png", fig.s3, width = 7, height=9)

