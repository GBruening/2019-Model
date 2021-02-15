library(ggplot2)
#library(ggforce)
library(ggthemes)
library(lmtest)
library(tidyr)
library(dplyr)
library(multcomp)
library(lme4)
library(knitr)
library(grid)
library(gridExtra)
library("RColorBrewer")
library(wesanderson)
library(extrafont)
library(cowplot)
library(latex2exp)
library(caret)
library(ggfortify)
library(ggResidpanel)
# library(readr)

# loadfonts(device='win')
windowsFonts(Times=windowsFont("TT Times New Roman"))

theme_cowplot <- theme_cowplot(font_size=12,font_family = "Times")
# theme_set(theme_gray(base_family = "Times"))
theme_classic_old <- theme_classic()
theme_classic <- theme_classic(base_family = "Times")+
  theme(text = element_text(color='black'),
        axis.text.x = element_text(color='black'),
        axis.text.y = element_text(color='black'),
        axis.ticks  = element_line(color='black'),
        axis.line = element_line(color='black',linetype='solid'),
        plot.title = element_text(hjust = 0.5))

theme_set(theme_classic(base_family = "Times")+
            theme(text = element_text(color='black'),
                  axis.text.x = element_text(color='black'),
                  axis.text.y = element_text(color='black'),
                  axis.ticks  = element_line(color='black'),
                  axis.line = element_line(color='black',linetype='solid'),
                  plot.title = element_text(hjust = 0.5)))

save_plots = 1
filter_data = 0

# Do Plots? ------------------
do_linfit_plots = 1
do_param_errbar_plots = 1
do_groupparam_plots = 1
do_diag_plots = 0

#Load Data------------------
# Using NeuromechVR
if (dir.exists('D:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Data')){
  # If Neuromech VR folder exists
  parent_fold = ('D:/Users/Gary/Google Drive/Muscle modeling/Min_jerk_files')
  eff_masses_meta <- read.csv('D:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Data/eff_masses.csv',header=0)
} else if (dir.exists('C:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Data/')){
  # If laptop folder exists
  parent_fold = ('C:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Data')
  eff_masses_meta <- read.csv('C:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Lappy Data/eff_masses.csv',header=0)
} else if (dir.exists('D:/Google Drive/Muscle modeling/Min_jerk_files')){
  # If Desky folder exists
  parent_fold = ('D:/Google Drive/Muscle modeling/Min_jerk_files')
  eff_masses_meta <- read.csv('D:/Google Drive/Muscle modeling/Metabolics/Lappy Data/eff_masses.csv',header=0)
} else {
  print('Didn\'t find either folder')
}

mpdata <-read.csv(paste(parent_fold,'/../metabolics/data/met_power_data.csv',sep=''))
sumdata_buff <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_buff.csv',sep=''))
sumdata_nobuff <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_nobuff.csv',sep=''))
sumdata_mj <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_mj.csv',sep=''))

sumdata_buff$buff = 'buff'
sumdata_nobuff$buff = 'nobuff'
sumdata_mj$buff = 'mj'

sumdata = rbind(sumdata_buff,sumdata_nobuff,sumdata_mj)


min_labs = c(TeX('Stress $(N/m^2)$'),
             TeX('$Stress^2$ $(N/m^2)^2$'),
             TeX('Force $(N)$'),
             TeX('$Force^2$ $(N^2)$'),
             TeX('Active State$'),
             TeX('Active $State^2$'),
             TeX('Neural Drive'),
             TeX('Neural $Drive^2$'),
             TeX('Umberger Model $(W)$'))
# min_labs2 = c('bquote(\'Minimizing Stress (N /\'~m^2~\')\')',
#               'bquote(\'Minimizing \'~Stress^2~\'(N /\'~m^2\'~)^2\')',
#               'bquote(\'Minimizing Force~(N)\')',
#               'bquote(\'Minimizing Force\'^2 (N\'^2)\')',
#               'bquote(\'Minimizing Active State\')',
#               'bquote(\'Minimizing Active State\'^2)',
#               'bquote(\'Minimizing Neural Drive\')',
#               'bquote(\'Minimizing Neural Drive\'^2\')',
#               'bquote(\'Minimizing Umberger Energy Model (W)\')')
min_labs2 = c('bquote(\'Minimizing Stress\')',
              'bquote(\'Minimizing \'~Stress^2)',
              'bquote(\'Minimizing Force\')',
              'bquote(\'Minimizing \'~Force^2)',
              'bquote(\'Minimizing Active State\')',
              'bquote(\'Minimizing Active \'~State^2)',
              'bquote(\'Minimizing Neural Drive\')',
              'bquote(\'Minimizing Neural \'~Drive^2)',
              'bquote(\'Minimizing Umberger Energy Model\')')

fitting_labs_split = c(TeX('Torque (Nm)'),
                       TeX('Output Force $(N)$'),
                       TeX('Muscle Force $(N)$'),
                       TeX('Stress $(N/m^2)$'),
                       TeX('Active State$'),
                       TeX('Neural Drive'),
                       TeX('$Torque^2$ $(Nm)^2$'),
                       TeX('Output $Force^2$ $(N)$'),
                       TeX('Muscle $Force^2$ $(N^2)$'),
                       TeX('$Stress^2$ $(N/m^2)^2$'),
                       TeX('Active $State^2$'),
                       TeX('Neural $Drive^2$'),
                       TeX('Umberger Model $(W)$'),
                       TeX('Bhar Model $(W)$'),
                       TeX('Uch Model $(W)$'),
                       TeX('Lich Model $(W)$'),
                       TeX('Marg Model $(W)$'))#,
# TeX('Mine Model $(W)$'),
# TeX('Houd Model $(W)$'))

fitting_labs = c(TeX('Torque (Nm)'),
                 TeX('$Torque^2$ $(Nm)^2$'),
                 TeX('Output Force $(N)$'),
                 TeX('Output $Force^2$ $(N)$'),
                 TeX('Muscle Force $(N)$'),
                 TeX('Muscle $Force^2$ $(N^2)$'),
                 TeX('Stress $(N/m^2)$'),
                 TeX('$Stress^2$ $(N/m^2)^2$'),
                 TeX('Active State$'),
                 TeX('Active $State^2$'),
                 TeX('Neural Drive'),
                 TeX('Neural $Drive^2$'),
                 TeX('Umberger Model $(W)$'),
                 TeX('Bhar Model $(W)$'),
                 TeX('Uch Model $(W)$'),
                 TeX('Lich Model $(W)$'),
                 TeX('Marg Model $(W)$'))#),
# TeX('Mine Model $(W)$'),
# TeX('Houd Model $(W)$'))

rsq_frame = read.csv('rsq_frame_sensitivity.csv')
rsq_frame$buff = ordered(rsq_frame$buff,levels=c('buff','nobuff','mj'))
# rsq_frame[rsq_frame$buff == 'buff','buff'] = 'Buffer'
# rsq_frame[rsq_frame$buff == 'nobuff','buff'] = 'No Buffer'
# rsq_frame[rsq_frame$buff == 'mj','buff'] = 'Minimum Jerk'
# rsq_frame$buff = revalue(rsq_frame$buff,c('buff'='Buffer','nobuff'='No Buffer','mj'='Minimum Jerk'))

choose_vars = c(1,2,5,6,9,10,11,12,13,14,15,16,17)#,18,19)
filt_rsq = filter(rsq_frame, variable %in% choose_vars)
filt_rsq$var_num = rep(1:length(choose_vars),max(filt_rsq$minfunc_num))

RSQ_bar <- ggplot(filt_rsq,
                  aes(fill=factor(minfunc_num),
                      x=var_num,
                      y=rsquared))+
  geom_bar(position='dodge',stat='identity')+
  geom_text(data=filt_rsq,
            aes(x=var_num,
                y=rsquared,
                label=rsquared),
            position=position_dodge(width=0.91), 
            angle=90,
            hjust=-0.1,
            size = 1)+
  # theme_classic(base_family='Times')+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # legend.position = c(.85,.5),
        legend.background = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
  scale_x_discrete(limits = 1:length(choose_vars), 
                   labels=parse(text=fitting_labs[choose_vars]))+
  labs(x='Fitted Variable',
       y='R Squared Value',
       fill=minparams)+
  guides(fill=guide_legend(title="Minimization\nFunction"))+
  scale_fill_manual(labels = parse(text = min_labs), 
                    values = brewer.pal(9,"Spectral"))+
  # coord_cartesian(ylim = c(min(rsq_frame$rsquared)*.98,max(rsq_frame$rsquared)*1.02))+
  facet_grid(rows = 'buff')


choose_vars = c(1,2,5,6,9,10,11,12,13,14,15,16,17)#,18,19)
filt_rsq = filter(rsq_frame, variable %in% choose_vars)
filt_rsq$var_num = rep(1:length(choose_vars),max(filt_rsq$minfunc_num))

filt_rsq[filt_rsq$buff == 'nobuff','rsquared'] = filt_rsq[filt_rsq$buff == 'buff','rsquared'] - filt_rsq[filt_rsq$buff == 'nobuff','rsquared']
filt_rsq[filt_rsq$buff == 'mj','rsquared'] = filt_rsq[filt_rsq$buff == 'buff','rsquared'] - filt_rsq[filt_rsq$buff == 'mj','rsquared']

RSQ_bar_diff <- ggplot(filt_rsq,
                  aes(fill=factor(minfunc_num),
                      x=var_num,
                      y=rsquared))+
  geom_bar(position='dodge',stat='identity')+
  geom_text(data=filt_rsq,
            aes(x=var_num,
                y=rsquared,
                label=rsquared),
            position=position_dodge(width=0.91), 
            angle=90,
            hjust=-0.1,
            size = 1)+
  # theme_classic(base_family='Times')+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # legend.position = c(.85,.5),
        legend.background = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
  scale_x_discrete(limits = 1:length(choose_vars), 
                   labels=parse(text=fitting_labs[choose_vars]))+
  labs(x='Fitted Variable',
       y='R Squared Value',
       fill=minparams)+
  guides(fill=guide_legend(title="Minimization\nFunction"))+
  scale_fill_manual(labels = parse(text = min_labs), 
                    values = brewer.pal(9,"Spectral"))+
  # coord_cartesian(ylim = c(min(rsq_frame$rsquared)*.98,max(rsq_frame$rsquared)*1.02))+
  facet_grid(rows = 'buff')
