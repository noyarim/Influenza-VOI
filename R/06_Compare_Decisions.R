#Purpose:
#Script to estimate epidemiological outcomes of alternative vaccine decisions 
#Run serially & take arguments from command line as function inputs


library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
#--------------------------------------------------------------------------

setwd("/Users/kyueunlee/Library/CloudStorage/GoogleDrive-leex5499@umn.edu/My Drive/Z Drive/Project_all/2022_FluVOI/Influenza-VOI")
source("R/05_RunSeasonalFluModel.R")

# Emulate ggplot color scheme
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Trend in predominant strain
prev_dt <- read.xlsx("data/Prev_trend.xlsx")
prev_dt$Antigen_test <- factor(prev_dt$Antigen_test, labels=c("Matching","Not matching"))

plt_prev <- function(dt, subtype){
  # subset to one subtype
  this_dt <- dt %>% filter(Subtype == subtype)
  this_dt_t <- melt(this_dt, id.vars = c("Season","Subtype","Strain","Antigen_test"))
  this_dt_t$variable <- factor(this_dt_t$variable, labels=c("When vax decision\n is made", "Between vax decision\n and season","During season"))
  if (subtype == "H1N1"){
    ylab = "Frequency of season-predominant strain(%)"
  }else{
    ylab = ""
  }
  # generate a bar plot
  plt <- ggplot(data=this_dt_t, aes(x=variable,y=value, fill=Antigen_test))+
    ggtitle(subtype)+
    geom_bar(position = "dodge",stat = "identity", width = 0.4)+
    ylab(ylab)+
    facet_wrap(.~Season, ncol=1)+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5))+
    scale_fill_manual(values=c("dark grey","tomato"))
    
  
  return(plt)
}

plt_prev_h1n1 <- plt_prev(prev_dt,"H1N1")
plt_prev_h3n2 <- plt_prev(prev_dt,"H3N2")
plt_prev_byam <- plt_prev(prev_dt,"B/Yam")
plt_prev_bvic <- plt_prev(prev_dt,"B/Vic")

ggarrange(plt_prev_h1n1, plt_prev_h3n2,plt_prev_byam,plt_prev_bvic, 
          #labels=c("A","B","C","D"),
          nrow=1,
          common.legend=TRUE,
          legend = "bottom")

ggsave("output/prev_trend.pdf", width = 18, height = 12)

## 1. Difference in flu dynamic (Prevalence of Infection) ##
## 1a. Overall ##
I_hist = read.csv("output/tot_I(old_VE).csv") # Historical seasons
I_new1_2014 = read.csv("output/tot_I(new_VE1_2014).csv") # New vax decision #1
I_new1_2019 = read.csv("output/tot_I(new_VE1_2019).csv") # New vax decision #1
I_new2_2014 = read.csv("output/tot_I(new_VE2_2014).csv") # New vax decision #2
I_new2_2019 = read.csv("output/tot_I(new_VE2_2019).csv") # New vax decision #2
I_new3_2014 = read.csv("output/tot_I(new_VE3_2014).csv") # New vax decision #3
I_new3_2019 = read.csv("output/tot_I(new_VE3_2019).csv") # New vax decision #3

f_combine_I <- function(I.2014, I.2019, rVE_lb, group){
  
  # Season when vaccine component is updated
  I.2014$season <- '2014/15'
  I.2019$season <- '2019/20'
  I.2014 <- I.2014 %>% filter(X>=5*366 & X<6*366) %>% mutate(days = seq(1,366),value=x) %>% dplyr::select(-X,-x)
  I.2019 <- I.2019 %>% filter(X>=10*366 & X<11*366) %>% mutate(days = seq(1,366),value=x) %>% dplyr::select(-X,-x)
  
  I.comb <- rbind(I.2014, I.2019)
  I.comb$scenario = rVE_lb
  
  if(group == 'historical'){
    I.comb$group <- 'old'
  }else{
    I.comb$group <- 'new'
  }
  return(I.comb)
}

I_old1 <- f_combine_I(I_hist, I_hist, "rVE=1.28",'historical')
I_old2 <- f_combine_I(I_hist, I_hist, "rVE=2",'historical')
I_old3 <- f_combine_I(I_hist, I_hist, "rVE=10",'historical')
I_new1 <- f_combine_I(I_new1_2014, I_new1_2019, "rVE=1.28",'new')
I_new2 <- f_combine_I(I_new2_2014, I_new2_2019, "rVE=2",'new')
I_new3 <- f_combine_I(I_new3_2014, I_new3_2019, "rVE=10",'new')

# generate an epidemic curve graph
I_new <- rbind(I_new1,I_new2,I_new3)
I_old <- rbind(I_old1,I_old2,I_old3)
I_all <- rbind(I_old, I_new)
I_all$season = factor(I_all$season, levels=c("2014/15","2019/20"), labels=c("2014/15 season \n with mis-matching A(H3N2) vaccine component", "2019/20 season \n with mis-matching B/Victoria vaccine component"))
I_all$scenario = factor(I_all$scenario, levels = c("rVE=1.28","rVE=2","rVE=10"))
I_all$group = factor(I_all$group, levels = c("old","new"), labels = c("Historical vax decision","Updated vax decision"))

plt_I <- 
  ggplot(I_all,aes(x=days,y=value*331*10,group=group, color=group))+
  geom_line(size=0.8)+
  facet_grid(scenario~season)+
  #scale_x_continuous(breaks=seq(11*365+1,nrow(I_all),by=365), labels = c("2020/Sep", "2021/Sep", "2022/Sep", "2023/Sep", "2024/Sep"))+
  scale_x_continuous(breaks=seq(1,365,by=30.5), labels = c("Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))+
  #xlab("Year of September")+
  ylab("Number of infected individuals (per 100,000)")+
  scale_color_manual(values=c("dark grey", "tomato"))+
  #guides(color=guide_legend(title="Scenario"))+
  theme_bw()+
  theme(
    # backgrounds
    panel.background = element_blank(), axis.line=element_line(colour = "black"), 
    # axis 
    axis.text.x = element_text(angle = 40,vjust=0.7,size=10),axis.title.x=element_blank(),
    # legend
    legend.key = element_blank(), legend.title = element_blank(),legend.position='bottom',
    legend.text = element_text(size=10),
    # facet labels
    strip.text.x = element_text(size=11), strip.text.y=element_text(size=11)
  )+
  scale_y_continuous(expand = c(0,0), limits = c(0,250)) 
  #scale_color_manual(values=c("#F8766D","#A3A500","#00BF7D","#00B0F6","grey"))
ggsave("output/N_Infected.pdf", width=8, height=6)
ggsave("output/N_Infected.jpeg", width=10, height=6)


## 1b. Prevalence of infection by subtype ##
Is_hist = read.csv("output/tot_I_strain(old_VE).csv") # Historical seasons
Is_new1_2014 = read.csv("output/tot_I_strain(new_VE1_2014).csv") # New vax decision #1
Is_new1_2019 = read.csv("output/tot_I_strain(new_VE1_2019).csv") # New vax decision #1
Is_new2_2014 = read.csv("output/tot_I_strain(new_VE2_2014).csv") # New vax decision #2
Is_new2_2019 = read.csv("output/tot_I_strain(new_VE2_2019).csv") # New vax decision #2
Is_new3_2014 = read.csv("output/tot_I_strain(new_VE3_2014).csv") # New vax decision #3
Is_new3_2019 = read.csv("output/tot_I_strain(new_VE3_2019).csv") # New vax decision #3

# combine 2014 and 2019 and take difference from historical season
f_combine_Is <- function(Is.hist, Is.2014, Is.2019, rVE_lb){
 
  # Subset time periods of interest and take difference from historical season
  dIs.2014 <- Is.2014[Is.2014$X>=5*366 & Is.2014$X<6*366,c(-1)]-Is.hist[Is.hist$X>=5*366 & Is.hist$X<6*366,c(-1)]
  dIs.2019 <- Is.2019[Is.2019$X>=10*366 & Is.2019$X<11*366,c(-1)]-Is.hist[Is.hist$X>=10*366 & Is.hist$X<11*366,c(-1)]
  
  # season label
  dIs.2014$season <- "2014/15"
  dIs.2019$season <- "2019/20"
  
  dI.all <- rbind(cbind(day=seq(1,nrow(dIs.2014)),dIs.2014), cbind(day=seq(1,nrow(dIs.2019)),dIs.2019))
  colnames(dI.all) <- c("days","A(H1N1)","A(H3N2)","B(Yam)","B(Vic)","season")
  dI.all$rVE_lb <- rVE_lb
  
  return(dI.all)
}

dIs_new1 <- f_combine_Is(Is_hist, Is_new1_2014, Is_new1_2019,'rVE=1.28')
dIs_new2 <- f_combine_Is(Is_hist, Is_new2_2014, Is_new2_2019,'rVE=2')
dIs_new3 <- f_combine_Is(Is_hist, Is_new3_2014, Is_new3_2019,'rVE=10')

dIs_all <- rbind(dIs_new1, dIs_new2, dIs_new3)
dIs_all_t <- melt(dIs_all, id.vars=c("days","season","rVE_lb"))
dIs_all_t$season = factor(dIs_all_t$season, levels=c("2014/15","2019/20"), labels=c("2014/15 season \n with mis-matching A(H3N2) vaccine component", "2019/20 season \n with mis-matching B/Victoria vaccine component"))
dIs_all_t$rVE_lb = factor(dIs_all_t$rVE_lb, levels=c("rVE=1.28","rVE=2","rVE=10"))
# plot
plt_I_bystrain <-
  ggplot(data = dIs_all_t, aes(x=days, y=value*331*10, fill=variable))+
  facet_grid(rVE_lb~season)+
  geom_ribbon(aes(ymax=value*331*10,ymin=0),alpha=0.5)+
  #  geom_bar(position = "dodge",stat = "identity", width = 10)+
#  ggtitle("Number of averted hospitalization with alternative vaccine selection")+
  #geom_errorbar(data=rall_lb_df, aes(x=variable, ymin=value*328.2*10^6*0.005,ymax=value_ub*328.2*10^6*0.005), width=0.2)+
  #xlab(xlab)+
  ylab("Difference in the number of infected individuals  \n with updated vs. original vax decision (per 100,000)")+
  theme_bw()+
  theme(
    # backgrounds
    panel.background = element_blank(), axis.line=element_line(colour = "black"), 
    # axis 
    axis.text.x = element_text(angle = 40,vjust=0.7,size=10),axis.title.x=element_blank(),
    # legend
    legend.key = element_blank(), legend.title = element_blank(),legend.position='bottom',
    legend.text = element_text(size=10),
    # facet labels
    strip.text.x = element_text(size=11), strip.text.y=element_text(size=11)
  )+
  scale_x_continuous(breaks=seq(1,365,by=30.5), labels = c("Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))+
  scale_y_continuous(expand=c(0,0), limits=c(min(dIs_all_t$value*331*10)-10,max(dIs_all_t$value*331*10)+10))
ggsave(paste0("output/N_Infected_bystrain.pdf"),width=10,height=6)
# combine the two I plots
ggarrange(plt_I, plt_I_bystrain, labels=c("A","B"),ncol=1)
ggsave(paste0("output/N_Infected_all.pdf"),width=9,height=12)

ggplot()+
  geom_bar()

## 2. Number of hospitalizations ##
## 2a. Total number of hospitalizations ##
season_lb <-c("2012/13", "2013/14","2014/15","2015/16","2016/17","2017/18","2018/19","2019/20")
r_hist <- read.csv("output/ratebystrain(old_VE).csv")
r_new1_2014 <- read.csv("output/ratebystrain(new_VE1_2014).csv")
r_new1_2019 <- read.csv("output/ratebystrain(new_VE1_2019).csv")
r_new2_2014 <- read.csv("output/ratebystrain(new_VE2_2014).csv")
r_new2_2019 <- read.csv("output/ratebystrain(new_VE2_2019).csv")
r_new3_2014 <- read.csv("output/ratebystrain(new_VE3_2014).csv")
r_new3_2019 <- read.csv("output/ratebystrain(new_VE3_2019).csv")

# Combine 2014 and 2019 seasons
f_combine_r <- function(r.2014, r.2019){
  
  # Season when vaccine component is updated
  r.all <- r.2014
  r.all[8,] <- r.2019[8,]
  r.all$All <- rowSums(r.all[,-c(1)])
  colnames(r.all) <- c("season","A(H1N1)","A(H3N2)","B(Vic)","B(Yam)","All")
  r.all$season = season_lb
  
  return(r.all)
}

# Difference in the number of hospitalizations 
r_hist2 <- f_combine_r(r_hist, r_hist)
r_newVE1 <- f_combine_r(r_new1_2014, r_new1_2019)
r_newVE2 <- f_combine_r(r_new2_2014, r_new2_2019)
r_newVE3 <- f_combine_r(r_new3_2014, r_new3_2019)
# total hospitalizationi 
dr_newVE1 <- rowSums(r_hist - r_newVE1)
dr_newVE2 <- rowSums(r_hist - r_newVE2)
dr_newVE3 <- rowSums(r_hist - r_newVE3)
dr_all <- data.frame(season = season_lb, newVE1=dr_newVE1, newVE2=dr_newVE2, newVE3=dr_newVE3)
dr_all_t <- melt(dr_all, id.vars = "season")
ggplot(data = dr_all_t, aes(x=season, y=value/100000*328.2*10^6, fill=variable))+
  geom_bar(position = "dodge",stat = "identity", width = 0.6)+
  #ggtitle("Number of averted hospitalization with alternative vaccine selection")+
  #geom_errorbar(data=rall_lb_df, aes(x=variable, ymin=value*328.2*10^6*0.005,ymax=value_ub*328.2*10^6*0.005), width=0.2)+
  ylab("Number of averted hospitalization with alternative vaccine selection")+
  #xlab(xlab)+
  labs("Relative vaccine effectiveness with matching vaccine strain")+
  theme_bw()+
  theme(
    panel.spacing.x = unit(0,'mm'),
    panel.background = element_blank(), axis.line=element_line(colour = "black"),
    legend.title = element_text(size=12),legend.text = element_text(size=11),
    legend.position = c(0.9,0.8),
    axis.text.x = element_text(size=11, angle = 30, vjust=0.7),
    axis.text.y = element_text(size=11),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14))+
    scale_fill_manual(values = c("darkgrey","seagreen","skyblue","skyblue4","tomato","tomato4"),
                      labels=c("rVE=1.28", "rVE=2", "rVE=10"))+
    scale_x_discrete(limits=c("2014/15","2019/20"))+
    scale_y_continuous(expand=c(0,0), limits=c(-10000, max(dr_all_t$value/100000*328.2*10^6)+5000))
ggsave(paste0("output/num_hospitalizations2.pdf"),width=10,height=6)

##2b. Number of hospitalizations by strain ##
r_hist_t <- melt(r_hist2, id.vars="season")
r_newVE1_t <- melt(r_newVE1,id.vars="season")
r_newVE1_t$rVE <- "rVE=1.28"
r_newVE2_t <- melt(r_newVE2,id.vars="season")
r_newVE2_t$rVE <- "rVE=2"
r_newVE3_t <- melt(r_newVE3,id.vars="season")
r_newVE3_t$rVE <- "rVE=10"
#combine hospitalization outcomes 
r_all_t <- rbind(r_newVE1_t, r_newVE2_t, r_newVE3_t)%>% 
              mutate(hist = rep(r_hist_t$value,3),
                     dr = hist-value,
                     variable = factor(variable, levels=c("All","A(H1N1)","A(H3N2)","B(Vic)","B(Yam)")),
                     rVE = factor(rVE, levels=c("rVE=1.28","rVE=2","rVE=10"))) %>%
              filter(season %in% c("2014/15","2019/20"))


ggplot(data = r_all_t, aes(x=variable, y=dr*328.2*10, fill=variable))+
  geom_bar(position = "dodge",stat = "identity", width = 0.5, alpha=0.9)+
  facet_grid(rVE~season)+
  #ggtitle("Number of averted hospitalization with alternative vaccine selection")+
  #geom_errorbar(data=rall_lb_df, aes(x=variable, ymin=value*328.2*10^6*0.005,ymax=value_ub*328.2*10^6*0.005), width=0.2)+
  ylab("Number of averted hospitalizations \n with alternative vaccine selection")+
  #xlab(xlab)+
  theme_bw()+
  theme(
    panel.spacing.x = unit(0,'mm'),
    panel.background = element_blank(), axis.line=element_line(colour = "black"),
    legend.title = element_text(size=12),legend.text = element_text(size=11),
    legend.position = 'none',#c(0.9,0.8),
    axis.text.x = element_text(size=11, angle = 30, vjust=0.7),
    axis.text.y = element_text(size=11),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14))+
  scale_fill_manual(values = c("grey34",gg_color_hue(4)))
  #scale_x_discrete(limits=c("2014/15","2019/20"))+
  #scale_y_continuous(expand=c(0,0), limits=c(-10, max(r_all_t$dr*328.2*10)))
ggsave(paste0("output/num_hospitalizations2.pdf"),width=8,height=6)
