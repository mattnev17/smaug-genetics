3bp<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe3<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
head(oe3)
oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
head(oe5)
oe3.ord<-oe3[order(oe3$ratio),]
oe5.ord<-oe5[order(oe5$ratio),]
oe.ord<-rbind(oe3, oe5)
oe5<-oe5[,-5]
head(oe5)
oe5.ord<-oe5[order(oe5$ratio),]
oe5.ord$rk<-seq(1,631)
oe5.ord$rk<-seq(1,599)
oe3.ord$rk<-seq(1,599)
oe.ord<-rbind(oe3, oe5)
ggplot(oe.ord, aes(x=rk, y=ratio,group=res, colour=res))+geom_point()
require(ggplot2)
ggplot(oe.ord, aes(x=rk, y=ratio,group=res, colour=res))+geom_point()
head(oe.ord)
oe.ord<-rbind(oe3.ord, oe5.ord)
ggplot(oe.ord, aes(x=rk, y=ratio,group=res, colour=res))+geom_point()
ggsave("3bp_vs_5bp.png")
oe.ord$or<-(oe.ord$obs*(100000-oe.ord$exp))/(oe.ord$exp*(100000-oe.ord$obs))
head(oe.ord)
tail(oe.ord)
ggplot(oe.ord, aes(x=rk, y=or,group=res, colour=res))+geom_point()
ggsave("3bp_vs_5bp.png")
oe.ord$logodds<-log(oe.ord$or)
ggplot(oe.ord, aes(x=rk, y=logodds,group=res, colour=res))+geom_point()
ggsave("3bp_vs_5bp.png")
head(oe3.ord)
head(oe5.ord)
oe3.ord[1:15,]
oe5.ord[1:15,]
ls()
head(oe3)
head(oe5)
max(oe3$exp-oe5$exp)
savehistory("3bp_vs_5bp.r")
