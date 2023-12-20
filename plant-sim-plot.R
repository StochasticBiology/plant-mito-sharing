library(ggplot2)

df = read.csv("k2-out.csv")

ggplot(df, aes(x=n.gene, y=n.mito, fill=log(mean.complex))) + geom_tile() + facet_wrap(~ social)
ggplot(df, aes(x=n.gene, y=n.mito, fill=log(mean.protein))) + geom_tile() + facet_wrap(~ social)
ggplot(df, aes(x=n.gene, y=n.mito, fill=log(mean.rna))) + geom_tile() + facet_wrap(~ social)
ggplot(df, aes(x=n.gene, y=n.mito, fill=prop.zero)) + geom_tile() + facet_wrap(~ social)

