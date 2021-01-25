library(tidyverse)
library(data.table)
library(cowplot)
library(wesanderson)


theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
      axis.text.y=element_text(size=15),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      legend.position="none"
    )
}

mex<-fread("MexHigh.allPBE.txt",header=T) %>% mutate(POP="MH",POS=POS/1E6)
gua<-fread("GuaHigh.allPBE.txt",header=T) %>% mutate(POP="GH",POS=POS/1E6)
and<-fread("Andes.allPBE.txt",header=T) %>% mutate(POP="AN",POS=POS/1E6)
swu<-fread("SW_US.allPBE.txt",header=T) %>% mutate(POP="US",POS=POS/1E6)
x<-rbind(swu,mex,gua,and) %>% mutate(POP=factor(POP,levels=c("US","MH","GH","AN")))

gswu<-filter(x,CHROM==3,POS>48,POS<50,POP=="US") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color=wes_palette("Darjeeling1")[3])+
  ylim(-0.5,2)+theme_li()+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=48.76, xmax=48.85, ymin=-Inf, ymax=Inf,alpha=0.2,fill="grey")  

gmex<-filter(x,CHROM==3,POS>48,POS<50,POP=="MH") %>%
  ggplot(aes(y=PBE0,x=POS))+  
  ylim(-0.5,2)+theme_li()+
  geom_point(size=1,color=wes_palette("Darjeeling1")[1])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=48.76, xmax=48.85, ymin=-Inf, ymax=Inf,alpha=0.2,fill="grey")  

ggua<-filter(x,CHROM==3,POS>48,POS<50,POP=="GH") %>%
  ggplot(aes(y=PBE0,x=POS))+
  ylim(-0.5,2)+theme_li()+
  geom_point(size=1,color=wes_palette("Royal2")[3])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=48.76, xmax=48.85, ymin=-Inf, ymax=Inf,alpha=0.2,fill="grey")  


gand<-filter(x,CHROM==3,POS>48,POS<50,POP=="AN") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color="purple")+
  ylim(-0.5,2)+theme_li()+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size=16))+
  xlab("Position (Mb)")+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=48.76, xmax=48.85, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

bob=filter(x,CHROM==3,POS>48,POS<50) %>% 
  ggplot(aes(y=PBE0,x=POS,color=POP))+geom_point(size=2)+ theme(legend.position="bottom")+
  scale_color_manual(values=c(wes_palette("Darjeeling1")[c(3,1)], wes_palette("Royal2")[3], "purple"),name="Populations:")

legend <- get_legend(bob)
pgrid<-plot_grid(gswu,gmex,ggua,gand,ncol=1)
plot_grid(legend,pgrid,  ncol = 1, rel_heights = c(.1, 1))
ggsave("Figure2A.pdf")


#(Chr2: 230276620..230282304)

gswu<-filter(x,CHROM==2,POS>230,POS<230.5,POP=="US") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color=wes_palette("Darjeeling1")[3])+
  ylim(-0.5,1.5)+theme_li()+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=230.276, xmax=230.282, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

gmex<-filter(x,CHROM==2,POS>230,POS<230.5,POP=="MH") %>%
  ggplot(aes(y=PBE0,x=POS))+  
  ylim(-0.5,1.5)+theme_li()+
  geom_point(size=1,color=wes_palette("Darjeeling1")[1])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=230.276, xmax=230.282, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

ggua<-filter(x,CHROM==2,POS>230,POS<230.5,POP=="GH") %>%
  ggplot(aes(y=PBE0,x=POS))+
  ylim(-0.5,1.5)+theme_li()+
  geom_point(size=1,color=wes_palette("Royal2")[3])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=230.276, xmax=230.282, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  


gand<-filter(x,CHROM==2,POS>230,POS<230.5,POP=="AN") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color="purple")+
  ylim(-0.5,1.5)+theme_li()+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size=16))+
  xlab("Position (Mb)")+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=230.276, xmax=230.282, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

bob=filter(x,CHROM==2,POS>230,POS<230.5) %>% 
  ggplot(aes(y=PBE0,x=POS,color=POP))+geom_point(size=2)+ theme(legend.position="bottom")+
  scale_color_manual(values=c(wes_palette("Darjeeling1")[c(3,1)], wes_palette("Royal2")[3], "purple"),name="Populations:")

legend <- get_legend(bob)
pgrid<-plot_grid(gswu,gmex,ggua,gand,ncol=1)
plot_grid(legend,pgrid,  ncol = 1, rel_heights = c(.1, 1))
ggsave("Figure2B.pdf")

#(Chr1: 94099455..94101276)
gswu<-filter(x,CHROM==1,POS>93.5,POS<94.5,POP=="US") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color=wes_palette("Darjeeling1")[3])+
  ylim(-0.5,1.5)+theme_li()+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=94.09, xmax=94.10, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

gmex<-filter(x,CHROM==1,POS>93.5,POS<94.5,POP=="MH") %>%
  ggplot(aes(y=PBE0,x=POS))+  
  ylim(-0.5,1.5)+theme_li()+
  geom_point(size=1,color=wes_palette("Darjeeling1")[1])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=94.09, xmax=94.10, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

ggua<-filter(x,CHROM==1,POS>93.5,POS<94.5,POP=="GH") %>%
  ggplot(aes(y=PBE0,x=POS))+
  ylim(-0.5,1.5)+theme_li()+
  geom_point(size=1,color=wes_palette("Royal2")[3])+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=94.09, xmax=94.10, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  


gand<-filter(x,CHROM==1,POS>93.5,POS<94.5,POP=="AN") %>%
  ggplot(aes(y=PBE0,x=POS))+
  geom_point(size=1,color="purple")+
  ylim(-0.5,1.5)+theme_li()+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size=16))+
  xlab("Position (Mb)")+
  geom_hline(yintercept=0.4, color="red", linetype="dashed")+
  annotate("rect",xmin=94.09, xmax=94.10, ymin=-Inf, ymax=Inf,alpha=0.5,fill="grey")  

bob=filter(x,CHROM==1,POS>93.5,POS<94.5) %>% 
  ggplot(aes(y=PBE0,x=POS,color=POP))+geom_point(size=2)+ theme(legend.position="bottom")+
  scale_color_manual(values=c(wes_palette("Darjeeling1")[c(3,1)], wes_palette("Royal2")[3], "purple"),name="Populations:")

legend <- get_legend(bob)
pgrid<-plot_grid(gswu,gmex,ggua,gand,ncol=1)
plot_grid(legend,pgrid,  ncol = 1, rel_heights = c(.1, 1))
ggsave("Figure2B_1.pdf")