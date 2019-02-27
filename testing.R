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
# library(readr)

loadfonts(device='win')
windowsFonts(Times=windowsFont("TT Times New Roman"))

theme_set(theme_cowplot(font_size=12,font_family = "Times"))
theme_set(theme_gray(base_family = "Times"))
theme_set(theme_classic(base_family = "Times"))
#=================== Load Data  ==============================
# Using NeuromechVR
if (dir.exists('D:/Users/Gary/Google Drive/2019 Model')){
  # If Neuromech VR folder exists
  parent_fold = ('D:/Users/Gary/Google Drive/2019 Model')
  eff_masses_meta <- read.csv('D:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Data/eff_masses.csv',header=0)
} else if (dir.exists('C:/Users/Gary/Google Drive/2019 Model')){
  # If laptop folder exists
  parent_fold = ('C:/Users/Gary/Google Drive/2019 Model')
  eff_masses_meta <- read.csv('C:/Users/Gary/Google Drive/Muscle modeling/Metabolics/Lappy Data/eff_masses.csv',header=0)
} else {
  print('Didn\'t find either folder')
}
mpdata <-read.csv(paste(parent_fold,'/met_power_data.csv',sep=''))
sumdata <-read.csv(paste(parent_fold,'/Data_2-5-2019.csv',sep=''))
sumdata_rng <-read.csv(paste(parent_fold,'/Data_2-5-2019_rng.csv',sep=''))


#=================== Do Plots?  ==============================
# Do Plots?
do_linfit_plots = 1
do_param_errbar_plots = 1
do_groupparam_plots = 1

#=================== Add squared vars  ==============================
values <- c('stress','stress2','force','force2','actstate','actstate2','drive','drive2','umberger')
index <- c(1,2,3,4,5,6,7,8,9)
sumdata$minfunc <- values[match(sumdata$minfunc,index)]

sumdata$sumtorque2 = sumdata$sumtorque^2
sumdata$sumforceout2 = sumdata$sumforceout^2
sumdata$sumforcemus2 = sumdata$sumforcemus^2
sumdata$sumstress2 = sumdata$sumstress^2
sumdata$sumactstate2 = sumdata$sumactstate^2
sumdata$sumdrive2 = sumdata$sumdrive^2

sumdata_rng$sumtorque2 = sumdata_rng$sumtorque^2
sumdata_rng$sumforceout2 = sumdata_rng$sumforceout^2
sumdata_rng$sumforcemus2 = sumdata_rng$sumforcemus^2
sumdata_rng$sumstress2 = sumdata_rng$sumstress^2
sumdata_rng$sumactstate2 = sumdata_rng$sumactstate^2
sumdata_rng$sumdrive2 = sumdata_rng$sumdrive^2

#=================== Add effmass2  ==============================
# Adding Column for estimated eff_mass to met data experiment.
index <- c(1,2,3,4)
values <- c(2.47,4.730,6.990,11.50)
subjects = c(1,2,3,4,5,6,7,8)

eff_mass = numeric(length(sumdata$subj))
for (i in 1:length(sumdata$subj)){
  eff_mass[i] = eff_masses_meta[match(sumdata$c[i],index),match(sumdata$subj[i],subjects)]
}

sumdata$effmass <- eff_mass
sumdata_rng$effmass <- eff_mass
#=================== Add Colors  ==============================
# Adding Colors

color1 = c(255,255,0)/255
color2 = c(255,196,0)/255
color3 = c(216,117,1)/255
color4 = c(139,0,0)/255

mass_colors <- c(rgb(color1[1],color1[2],color1[3]),
                 rgb(color2[1],color2[2],color2[3]),
                 rgb(color3[1],color3[2],color3[3]),
                 rgb(color4[1],color4[2],color4[3]))
#=================== Begin loopin  ==============================
minparams = c('stress',
              'stress2',
              'force',
              'force2',
              'actstate',
              'actstate2',
              'drive',
              'drive2',
              'umberger')

fitting_vars = c('sumtorque',
                 'sumtorque2',
                 'sumforceout',
                 'sumforceout2',
                 'sumforcemus',
                 'sumforcemus2',
                 'sumstress',
                 'sumstress2',
                 'sumactstate',
                 'sumactstate2',
                 'sumdrive',
                 'sumdrive2',
                 'sumumber')
fitting_labs = c('Sum of Torque\n(Nm)',
                 'Sum of Torque^2\n(Nm)^2',
                 'Sum of Output Force\n(N)',
                 'Sum of Output Force^2\n(N^2)',
                 'Sum of Muscle Force\n(N)',
                 'Sum of Muscle Force^2\n(N^2)',
                 'Sum of Stress\n(N/m^2)',
                 'Sum of Stress^2\n(N/m^2)^2',
                 'Sum of Active State',
                 'Sum of Active State^2',
                 'Sum of Neural Drive',
                 'Sum of Neural Drive^2',
                 'Sum of Energetics\n(W)')

minfunc1 = values[1]
rsq = c()
minfuncs_rsq = c()
vars_rsq = c()
vars_rsq_labs = c()
expo_rsq = c()

minfuncs_abcd = c()
minfuncnum_abcd = c()
vars_abcd = c()
vars_abcd_labs = c()

a_val = c()
a_ste = c()
b_val = c()
b_ste = c()
c_val = c()
c_ste = c()
d_val = c()
d_ste = c()

fit_plots <- list()
fun.fit1 <- list()
fun.fit2 <- list()
fun.fit3 <- list()
fun.fit4 <- list()

# names(rsq_matrix) = c('min_func','variable','expo','rsq')
rsq_matrix = matrix(,nrow=length(minparams),ncol=length(fitting_labs))
expo=1
minfunc_count = 0

# fitting_vars = fitting_vars[3]
for (minfunc1 in minparams){
  minfunc_count = minfunc_count + 1
  varfit_count=0
  for (var in fitting_vars){
    varfit_count = varfit_count+1
    plotdata=filter(sumdata,minfunc==paste(minfunc1))
    eval(parse(text = paste('fit_',minfunc1,'_',var,'=summary(lm(plotdata$mpowernet ~ plotdata$',var,'))',sep='')))
    
    # print(paste('plot_',minfunc1,'_',var,sep=''))
    int = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[1]',sep='')))
    slope = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))
    eval(parse(text = paste('tempdata = filter(sumdata,minfunc==\'',minfunc1,'\')',sep='')))
    tempdata['temp_fitted'] = int + slope*tempdata[var]
    tempdata['lintransform'] = slope*tempdata[var]
    tempdata['unfitted'] = tempdata[var]
    # tempdata['metabolic'] = 1
    # tempdata['fitted'] = 0
    # tempdata['unfitted'] = 2
    tempdata = rbind(cbind(tempdata$c,
                           tempdata$effmass,
                           tempdata$movedur,
                           tempdata$mpowernet,
                           tempdata$temp_fitted,
                           tempdata$unfitted,
                           tempdata$lintransform))
    tempdata=data.frame(tempdata)
    colnames(tempdata)=c('c','effmass','movedur','mpowernet','fitted','unfitted','lintransform')
    
    rsq_val = eval(parse(text=paste('round(fit_',minfunc1,'_',var,'$r.squared,digits = 3)',sep='')))
    if (do_linfit_plots){
      titlestr = paste('Min ',minfunc1,'\nR^2 = ',rsq_val,sep='')
      string = paste('lmplot_',minfunc1,'_',var,'=
        ggplot(data=plotdata,aes(x=',var,',y=mpowernet))+
        geom_point()+geom_smooth(method=\'lm\')+
        labs(title = \'',titlestr,'\',
          x =\'',fitting_labs[varfit_count],'\',
          y=\'Net Metabolic Power (W)\')+
        theme_classic(base_family=\'Times\')+theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(color=\'black\',size = 1,linetype=\'solid\'))',sep='')
      eval(parse(text = string))
    }
    
    rsq=c(rsq,rsq_val)
    minfuncs_rsq=c(minfuncs_rsq,minparams[minfunc_count])
    vars_rsq = c(vars_rsq,varfit_count)
    vars_rsq_labs = c(vars_rsq_labs,var)
    expo_rsq = c(expo_rsq,expo)
    rsq_matrix[minfunc_count,varfit_count] = rsq_val

#====================== Creating the fitted plots==========================================
    MP_model=nls(mpowernet ~ a1 + 100*a2*(effmass^a3)/(movedur^a4), 
                 data=tempdata,
                 start=list(a1=1,a2=.1,a3=1,a4=1))
    modelsum_met = summary(MP_model)
    a=coef(MP_model)[1]
    b=coef(MP_model)[2]*100
    c=coef(MP_model)[3]
    d=coef(MP_model)[4]
    
    fun.1 <- function(t) a+b*(2.47^c)/(t^d)
    fun.2 <- function(t) a+b*(4.73^c)/(t^d)
    fun.3 <- function(t) a+b*(6.99^c)/(t^d)
    fun.4 <- function(t) a+b*(11.50^c)/(t^d)
    
    Fitted_model=nls(fitted ~ a1 + 100*a2*(effmass^a3)/(movedur^a4), 
                     data=tempdata,
                     start=list(a1=1,a2=.1,a3=1,a4=1))
    modelsum_prox = summary(Fitted_model)
    afit=coef(Fitted_model)[1]
    bfit=coef(Fitted_model)[2]*100
    cfit=coef(Fitted_model)[3]
    dfit=coef(Fitted_model)[4]
    
    eval(parse(text = paste('fun.fit1[[var]] <- function(t) ',afit,'+',bfit,'*(2.47^',cfit,')/(t^',dfit,')',sep='')))
    eval(parse(text = paste('fun.fit2[[var]] <- function(t) ',afit,'+',bfit,'*(4.73^',cfit,')/(t^',dfit,')',sep='')))
    eval(parse(text = paste('fun.fit3[[var]] <- function(t) ',afit,'+',bfit,'*(6.99^',cfit,')/(t^',dfit,')',sep='')))
    eval(parse(text = paste('fun.fit4[[var]] <- function(t) ',afit,'+',bfit,'*(11.50^',cfit,')/(t^',dfit,')',sep='')))
    # if (tempdata$metabolic[1]!="Metabolic"){
    #   index <- c(0,1)
    #   eval(parse(text=paste('values <- c(\'',var,'\',\'Metabolic\')',sep='')))
    #   tempdata$metabolic = values[match(tempdata$metabolic,index)]
    # }
    
    if (do_linfit_plots){
      fit_plots[[var]] <- ggplot(tempdata,aes(x=movedur,y=mpowernet,color=factor(c)))+#,shape=factor(metabolic)))+
        geom_point(size=3)+
        geom_point(data = tempdata,
                   aes(x=movedur,y=mpowernet),
                   shape=21,
                   fill=factor(c),
                   color='black',
                   size=3)#+
        # geom_point(data=filter(tempdata,metabolic==var),
        #            aes(x=movedur,y=mpowernet),
        #            shape=24,
        #            fill=factor(c),
        #            color='black',
        #            size=3)
      
      fit_plots[[var]] <- fit_plots[[var]] + 
        # stat_function(fun=fun.1,size=3,color=mass_colors[1])+
        # stat_function(fun=fun.2,size=3,color=mass_colors[2])+ 
        # stat_function(fun=fun.3,size=3,color=mass_colors[3])+ 
        # stat_function(fun=fun.4,size=3,color=mass_colors[4])+ 
        stat_function(fun=fun.fit1[[var]],size=1.5,color=mass_colors[1])+#, linetype="dashed")+
        stat_function(fun=fun.fit2[[var]],size=1.5,color=mass_colors[2])+#, linetype="dashed")+ 
        stat_function(fun=fun.fit3[[var]],size=1.5,color=mass_colors[3])+#, linetype="dashed")+ 
        stat_function(fun=fun.fit4[[var]],size=1.5,color=mass_colors[4])+#, linetype="dashed")+
        scale_color_manual(values = mass_colors)
      fit_plots[[var]] <- fit_plots[[var]] + labs(y='Metabolic Power (W)',x='Movement Duration (s)',shape='Metabolic',color='Effective\nMass (kg)')
      fit_plots[[var]] <- fit_plots[[var]] + theme_classic(base_family='Times')+theme(plot.title = element_text(hjust = 0.5),
                                     # axis.line = element_line(color='black',size = 1,linetype='solid'),
                                     #text=element_text(family="Arial"),
                                     legend.position='none')
      eval(parse(text = paste('fit_plots[[var]]<-fit_plots[[var]]+labs(title=\'Fitting ',var,'\nMin ',minfunc1,', R^2=',rsq_val,'\')')))
      # dev.off()
      # eval(parse(text = paste('fitplot_',var,'<-g',sep='')))
    }
    #================== Create parameter fits from utility model
    if (var=='sumforcemus2' | var=='sumforcemus' | var=='sumforceout' | var=='sumforceout'){
      scale = 100
    } else if (var=='sumstress2' | var=='sumstress'){
      scale = 100
    } else {
      scale = 10
    }
    Unfitted_model=nls(lintransform ~ a1 + scale*a2*(effmass^a3)/(movedur^a4), 
                     data=tempdata,
                     start=list(a1=1,a2=.1,a3=1,a4=1))
    modelsum_unfitted = summary(Unfitted_model)
    
    minfuncs_abcd=c(minfuncs_abcd,minparams[minfunc_count])
    minfuncnum_abcd=c(minfuncnum_abcd,minfunc_count)
    vars_abcd = c(vars_abcd,varfit_count)
    vars_abcd_labs = c(vars_abcd_labs,var)
    
    a_val = c(a_val,modelsum_unfitted$coefficients[1])
    a_ste = c(a_ste,modelsum_unfitted$coefficients[5])
    b_val = c(b_val,modelsum_unfitted$coefficients[2]*scale)
    b_ste = c(b_ste,modelsum_unfitted$coefficients[6]*scale)
    c_val = c(c_val,modelsum_unfitted$coefficients[3])
    c_ste = c(c_ste,modelsum_unfitted$coefficients[7])
    d_val = c(d_val,modelsum_unfitted$coefficients[4])
    d_ste = c(d_ste,modelsum_unfitted$coefficients[8])
    
  }
  
  minfuncs_abcd=c(minfuncs_abcd,minparams[minfunc_count])
  minfuncnum_abcd=c(minfuncnum_abcd,minfunc_count)
  vars_abcd = c(vars_abcd,varfit_count+1)
  vars_abcd_labs = c(vars_abcd_labs,'Metabolics')
  
  a_val = c(a_val,modelsum_met$coefficients[1])
  a_ste = c(a_ste,modelsum_met$coefficients[5])
  b_val = c(b_val,modelsum_met$coefficients[2]*100)
  b_ste = c(b_ste,modelsum_met$coefficients[6]*100)
  c_val = c(c_val,modelsum_met$coefficients[3])
  c_ste = c(c_ste,modelsum_met$coefficients[7])
  d_val = c(d_val,modelsum_met$coefficients[4])
  d_ste = c(d_ste,modelsum_met$coefficients[8])
  
  # Variable raised to the first
  if (do_linfit_plots){
    string = paste('lmgroupvar1_',minfunc1,
                   '=plot_grid(lmplot_',minfunc1,'_','sumtorque,
                   lmplot_',minfunc1,'_','sumforceout,
                   lmplot_',minfunc1,'_','sumforcemus,
                   lmplot_',minfunc1,'_','sumstress,
                   lmplot_',minfunc1,'_','sumactstate,
                   lmplot_',minfunc1,'_','sumdrive,
                   lmplot_',minfunc1,'_','sumumber,
                   align = \'vh\',
                   labels = c(\'A\',\'B\',\'C\',\'D\',\'E\',\'F\',\'G\'),
                   hjust=-1,
                   nrow=2)',sep='')
    eval(parse(text=string))
    
    setwd(paste(parent_fold,'/2019 Graphs/lm_group_var1',sep=''))
    string = paste('lmgroupvar1_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    string = paste('ggsave(filename,plot=lmgroupvar1_',minfunc1,',useDingbats = FALSE,width=9,height=6,units=\'in\')',sep='')
    eval(parse(text=string))
    
    # Variable Raised to the 2nd
    string = paste('lmgroupvar2_',minfunc1,
                   '=plot_grid(lmplot_',minfunc1,'_','sumtorque2,
                   lmplot_',minfunc1,'_','sumforceout2,
                   lmplot_',minfunc1,'_','sumforcemus2,
                   lmplot_',minfunc1,'_','sumstress2,
                   lmplot_',minfunc1,'_','sumactstate2,
                   lmplot_',minfunc1,'_','sumdrive2,
                   align = \'vh\',
                   labels = c(\'A\',\'B\',\'C\',\'D\',\'E\',\'F\',\'G\'),
                   hjust=-1,
                   nrow=2)',sep='')
    eval(parse(text=string))
    
    setwd(paste(parent_fold,'/2019 Graphs/lm_group_var2',sep=''))
    string = paste('lmgroupvar2_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    string = paste('ggsave(filename,plot=lmgroupvar2_',minfunc1,',useDingbats = FALSE,width=9,height=6,units=\'in\')',sep='')
    eval(parse(text=string))
    
    
    # Fitted Plot Raised to the 1st
    add_legend <- get_legend(fit_plots[[1]]+theme(legend.position='right'))
    string = paste('fittedvar1_',minfunc1,
                   '=plot_grid(fit_plots[[1]],
                   fit_plots[[2]],
                   fit_plots[[3]],
                   fit_plots[[4]],
                   fit_plots[[5]],
                   fit_plots[[6]],
                   fit_plots[[7]],
                   fit_plots[[8]],
                   fit_plots[[9]],
                   fit_plots[[10]],
                   fit_plots[[11]],
                   fit_plots[[12]],
                   fit_plots[[13]],
                   add_legend,
                   align = \'vh\',
                   labels = c(\'A\',\'B\',\'C\',\'D\',\'E\',\'F\',\'G\',\'H\',\'I\',\'J\',\'K\',\'L\',\'M\'),
                   hjust=-1,
                   nrow=4)',sep='')
    eval(parse(text=string))
    # eval(parse(text = paste('fittedvar1_',minfunc1,'<-plot_grid(fittedvar1_',minfunc1,',add_legend,ncol=2,rel_widths=c(.85,.15))',sep='')))
    
    setwd(paste(parent_fold,'/2019 Graphs/fitted_var1',sep=''))
    string = paste('fittedvar1_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    string = paste('ggsave(filename,plot=fittedvar1_',minfunc1,',useDingbats = FALSE,width=16,height=12,units=\'in\')',sep='')
    eval(parse(text=string))
   
    # # Fitted Plot Raised to the 2nd
    # string = paste('fittedvar2_',minfunc1,
    #                '=plot_grid(fitplot_sumtorque2,
    #                fitplot_sumstress2,
    #                fitplot_sumforcemus2,
    #                fitplot_sumactstate2,
    #                fitplot_sumdrive2,
    #                align = \'vh\',
    #                labels = c(\'A\',\'B\',\'C\',\'D\',\'E\'),
    #                hjust=-1,
    #                nrow=2)',sep='')
    # eval(parse(text=string))
    # add_legend <- get_legend(fitplot_sumstress+theme(legend.position='right'))
    # eval(parse(text = paste('fittedvar2_',minfunc1,'<-plot_grid(fittedvar2_',minfunc1,',add_legend,ncol=2,rel_widths=c(.85,.15))',sep='')))
    # 
    # setwd('C:/Users/Gary/Google Drive/2019 Model/2019 Graphs/fitted_var2')
    # string = paste('fittedvar2_',minfunc1,sep='')
    # filename = paste(string,'.pdf',sep='')
    # string = paste('ggsave(filename,plot=fittedvar2_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
    # eval(parse(text=string))
  }
  
}

rsq_frame = data.frame(minfunc = minfuncs_rsq,variable=vars_rsq,expo = expo_rsq,rsquared = rsq)
# rsq_frame$variable = as.numeric(rsq_frame$variable)
# rsq_frame <- rsq_frame[order(rsq_frame$variable),]
# rownames(rsq_frame) = NULL

min_labs = c(TeX('Stress $(N/m^2)$'),
             TeX('$Stress^2$ $(N/m^2)^2$'),
             TeX('Force $(N)$'),
             TeX('$Force^2$ $(N^2)$'),
             TeX('Active State$'),
             TeX('Active $State^2$'),
             TeX('Neural Drive'),
             TeX('Neural $Drive^2$'),
             TeX('Energy Model $(W)$'))
fitting_labs_split = c(TeX('Torque (Nm)'),
                 TeX('Output Force $(N)$'),
                 TeX('Muscle Force $(N)$'),
                 TeX('Stress $(N/m^2)$'),
                 TeX('Active State$'),
                 TeX('Neural Drive'),
                 TeX('$Torque^2$ $(Nm)^2$'),
                 TeX('$Output Force^2$ $(N)$'),
                 TeX('$Muscle Force^2$ $(N^2)$'),
                 TeX('$Stress^2$ $(N/m^2)^2$'),
                 TeX('Active $State^2$'),
                 TeX('Neural $Drive^2$'),
                 TeX('Energy Model $(W)$'))
fitting_labs = c(TeX('Torque (Nm)'),
                 TeX('$Torque^2$ $(Nm)^2$'),
                 TeX('Output Force $(N)$'),
                 TeX('$Output Force^2$ $(N)$'),
                 TeX('Muscle Force $(N)$'),
                 TeX('$Muscle Force^2$ $(N^2)$'),
                 TeX('Stress $(N/m^2)$'),
                 TeX('$Stress^2$ $(N/m^2)^2$'),
                 TeX('Active State$'),
                 TeX('Active $State^2$'),
                 TeX('Neural Drive'),
                 TeX('Neural $Drive^2$'),
                 TeX('Energy Model $(W)$'))

if (do_linfit_plots){
  ggplot(rsq_frame,aes(fill=factor(minfunc),x=reorder(variable,variable),y=rsquared))+
    geom_bar(position='dodge',stat='identity')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:13), labels=parse(text=fitting_labs))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))
  setwd(paste(parent_fold,'/2019 Graphs',sep=''))
  ggsave('Rsq_matrix_expo1.pdf', useDingbats = FALSE,width=7,height=4.5,units='in')
  
  
  rsq_frame$variable = factor(rsq_frame$variable,
                              levels = (1:13),
                              labels = c(1,7,2,8,3,9,4,10,5,11,6,12,13))
  ggplot(rsq_frame,aes(fill=factor(minfunc),x=reorder(variable,variable),y=rsquared))+
    geom_bar(position='dodge',stat='identity')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:13), labels=parse(text=fitting_labs_split))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))
  rsq_frame$variable = factor(rsq_frame$variable,
                              labels = (1:13),
                              levels = c(1,7,2,8,3,9,4,10,5,11,6,12,13))
  setwd(paste(parent_fold,'/2019 Graphs',sep=''))
  ggsave('Rsq_matrix_expo1_split.pdf', useDingbats = FALSE,width=7,height=4.5,units='in')
  
  best_min_func = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc
  rsq_frame_filt = filter(rsq_frame,minfunc==filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc)
  ggplot(rsq_frame_filt,aes(x=variable,y=rsquared))+
    geom_bar(stat='identity',fill='blue')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limits = (1:13), labels=parse(text=fitting_labs_split))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))#+
    # coord_cartesian(ylim = c(.4, .7))
  setwd(paste(parent_fold,'/2019 Graphs',sep=''))
  ggsave('Rsq_matrix_expo1_act2.pdf', useDingbats = FALSE,width=7,height=4.5,units='in')
  
  rsq_frame_filt$variable = factor(rsq_frame_filt$variable,
                              levels = (1:13),
                              labels = c(1,7,2,8,3,9,4,10,5,11,6,12,13))
  ggplot(rsq_frame_filt,aes(x=variable,y=rsquared))+
    geom_bar(stat='identity',fill='blue')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limits = (1:13), labels=parse(text=fitting_labs_split))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))
  rsq_frame_filt$variable = factor(rsq_frame_filt$variable,
                              labels = (1:13),
                              levels = c(1,7,2,8,3,9,4,10,5,11,6,12,13))
  setwd(paste(parent_fold,'/2019 Graphs',sep=''))
  ggsave('Rsq_matrix_expo1_act2_split.pdf', useDingbats = FALSE,width=7,height=4.5,units='in')
}

#=================== Param Error Bar Plots ==============================
# Creating the parameter error bar plots
parameter_frame = data.frame(minfunc = minfuncs_abcd,minfuncnum = minfuncnum_abcd,variable=vars_abcd,var_lab=vars_abcd_labs,a_val=a_val,a_ste=a_ste,b_val=b_val,b_ste=b_ste,c_val=c_val,c_ste=c_ste,d_val=d_val,d_ste=d_ste)


if (do_param_errbar_plots){
# if (0){
  minparams = c('stress',
                'stress2',
                'force',
                'force2',
                'actstate',
                'actstate2',
                'drive',
                'drive2',
                'umberger')
  fitting_vars = c('sumtorque',
                   'sumforceout',
                   'sumforcemus',
                   'sumstress',
                   'sumactstate',
                   'sumdrive',
                   'sumtorque2',
                   'sumforceout2',
                   'sumforcemus2',
                   'sumstress2',
                   'sumactstate2',
                   'sumdrive2',
                   'sumumber',
                   'Metabolics')
  fitting_labs = c(TeX('Torque (Nm)'),
                   TeX('Output Force $(N)$'),
                   TeX('Muscle Force $(N)$'),
                   TeX('Stress $(N/m^2)$'),
                   TeX('Active State$'),
                   TeX('Neural Drive'),
                   TeX('$Torque^2$ $(Nm)^2$'),
                   TeX('$Output Force^2$ $(N^2)$'),
                   TeX('$Muscle Force^2$ $(N^2)$'),
                   TeX('$Stress^2$ $(N/m^2)^2$'),
                   TeX('Active $State^2$'),
                   TeX('Neural $Drive^2$'),
                   TeX('Energy Model $(W)$'),
                   TeX('Metabolic (W)'))
  
  # paste('Sum of Torque\n(Nm)'),
  # paste('Sum of \'',TeX('$Torque^2$'),'\'\n\'',TeX('$(Nm)^2$'),'\''),
  # paste('Sum of Force\n(Nm)'),
  # paste('Sum of ',TeX('$Force^2$'),'\n',TeX('$(N)^2$')),
  # paste('Sum of Stress\n(Nm)'),
  # paste('Sum of ',TeX('$Stress^2$'),'\n',TeX('$(N/m^2)^2$')),
  # paste('Sum of Active\n(Nm)'),
  # paste('Sum of ',TeX('$Active^2$')),
  # paste('Sum of Neural\n(Nm)'),
  # paste('Sum of ',TeX('$Neural^2$')),
  # paste('Sum of Energetics\n(Nm)'),
  # paste('Metabolic (W)'))
  
  util_param_select = c(7,8,9,10,11,12,13)
  fitting_vars = fitting_vars[util_param_select]
  fitting_labs = fitting_labs[util_param_select]
  parameter_frame2 = filter(parameter_frame,var_lab %in% fitting_vars)
  
  util_params = c('a','b','c','d')
  minfunc_count = 0
  for (minfunc1 in minparams){
    minfunc_count = minfunc_count + 1
    
    param_count = 0
      for (param in util_params){
        rm(g)
        param_count = param_count+1
        temp_pframe = filter(parameter_frame2,minfunc==eval(minfunc1))
        string = paste('g=ggplot(temp_pframe,aes(x=factor(variable),y=',param,'_val,color=factor(variable)))+
                       geom_point(size=3.5)+
                       geom_errorbar(aes(ymin=',param,'_val-',param,'_ste,ymax=',param,'_val+',param,'_ste),width=.3,size=1)',sep='')
        
        eval(parse(text = paste(string)))
        
        yminval=temp_pframe[length(fitting_vars),param_count*2+3]-temp_pframe[length(fitting_vars),param_count*2+4]
        ymaxval=temp_pframe[length(fitting_vars),param_count*2+3]+temp_pframe[length(fitting_vars),param_count*2+4]
        # g+
        
        g<-g+
          geom_rect(xmin=0,xmax=length(fitting_vars)+1,ymin=yminval,ymax=ymaxval,alpha=.01)+
          scale_x_discrete(labels=parse(text=fitting_labs))+
          labs(x='Effort Variable',title=eval(param),y=eval(param),color="Effort\nRep")+
          scale_color_manual(labels=parse(text=fitting_labs),
                             values = colorRampPalette(brewer.pal(12,'Dark2'))(12))+
          theme_classic(base_family='Times')+
          theme(legend.position='none',axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(plot.title = element_text(hjust = 0.5))
        
        eval(parse(text = paste('paramplot_',minfunc1,'_',param,'<-g',sep='')))
        
      }
    # 'paramplot_',minfunc1,
    string = paste('g<-plot_grid(paramplot_',minfunc1,'_a,
                   paramplot_',minfunc1,'_b,
                   paramplot_',minfunc1,'_c,
                   paramplot_',minfunc1,'_d,
                   align = \'vh\',
                   labels = c(\'A\',\'B\',\'C\',\'D\'),
                   hjust=-1,
                   nrow=2)',sep='')
    eval(parse(text=string))
    string = paste('add_legend <- get_legend(paramplot_',minfunc1,'_a+theme(legend.position=\'right\'))',sep='')
    eval(parse(text=string))
    
    g<-plot_grid(g,add_legend,nrow=1,rel_widths = c(.8,.2))
    
    
    title <- ggdraw() + 
      draw_label("Parameter Fitting, min stress",
                 fontface = 'bold')
    g<-plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1))
    eval(parse(text = paste('paramplot_',minfunc1,'<-g',sep='')))
    setwd(paste(parent_fold,'/2019 Graphs/param_plots',sep=''))
    filename = paste('paramplot_',minfunc1,'.pdf',sep='')
    string = paste('ggsave(filename,plot=paramplot_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
    eval(parse(text=string))
  }
}
#=================== Every param ==============================
# Create the plots where every min func is on every plot
if (do_groupparam_plots){
  minparams = c('stress',
                'stress2',
                'force',
                'force2',
                'actstate',
                'actstate2',
                'drive',
                'drive2',
                'umberger')
  min_labs = c(TeX('Stress $(N/m^2)$'),
              TeX('$Stress^2$ $(N/m^2)^2$'),
              TeX('Force $(N)$'),
              TeX('$Force^2$ $(N^2)$'),
              TeX('Active State$'),
              TeX('Active $State^2$'),
              TeX('Neural Drive'),
              TeX('Neural $Drive^2$'),
              TeX('Energy Model $(W)$'))
  fitting_vars = c('sumtorque',
                   'sumtorque2',
                   'sumforceout',
                   'sumforceout2',
                   'sumforcemus',
                   'sumforcemus2',
                   'sumstress',
                   'sumstress2',
                   'sumactstate',
                   'sumactstate2',
                   'sumdrive',
                   'sumdrive2',
                   'sumumber',
                   'Metabolics')
  fitting_labs = c(TeX('Torque (Nm)'),
                   TeX('$Torque^2$ $(Nm)^2$'),
                   TeX('Output Force $(N)$'),
                   TeX('$Output Force^2$ $(N^2)$'),
                   TeX('Muscle Force $(N)$'),
                   TeX('$Muscle Force^2$ $(N^2)$'),
                   TeX('Stress $(N/m^2)$'),
                   TeX('$Stress^2$ $(N/m^2)^2$'),
                   TeX('Active State$'),
                   TeX('Active $State^2$'),
                   TeX('Neural Drive'),
                   TeX('Neural $Drive^2$'),
                   TeX('Energy Model $(W)$'),
                   TeX('Metabolic (W)'))
  util_param_select = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  fitting_vars = fitting_vars[util_param_select]
  fitting_labs = fitting_labs[util_param_select]
  temp_pframe = filter(parameter_frame,var_lab %in% fitting_vars)
  
  util_params = c('a','b','c','d')
  
  param_count = 0
  for (param in util_params){
    param_count = param_count+1
    # temp_pframe = filter(parameter_frame,minfunc==eval(minfunc1))
    string = paste('g<-ggplot(temp_pframe,aes(x=factor(variable),y=',param,'_val,color=factor(minfuncnum),group=minfuncnum))+geom_point(size=3, position=position_dodge(width=.8))+geom_errorbar(aes(x=factor(variable),ymin=',param,'_val-',param,'_ste,ymax=',param,'_val+',param,'_ste),width=.3,size=1,position=position_dodge(width=.8))',sep='')
    eval(parse(text = string))
    yminval=temp_pframe[length(temp_pframe$minfunc),param_count*2+3]-temp_pframe[length(temp_pframe$minfunc),param_count*2+4]
    ymaxval=temp_pframe[length(temp_pframe$minfunc),param_count*2+3]+temp_pframe[length(temp_pframe$minfunc),param_count*2+4]
    g<-g+
      geom_rect(xmin=0,xmax=length(fitting_vars)+1,ymin=yminval,ymax=ymaxval,alpha=.005)+
      scale_x_discrete(labels=parse(text=fitting_labs))+
      labs(x='Effort Variable',title=eval(param),y=eval(param),color="Minimization\nFunction")+
      scale_color_discrete(labels=parse(text=min_labs))+
      theme_classic(base_family='Times')+
      theme(legend.position='none',axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(plot.title = element_text(hjust = 0.5))
      # ,values = colorRampPalette(brewer.pal(12,'Dark2'))(12))
    
    eval(parse(text = paste('groupparamplot_',param,'<-g',sep='')))
  }
  string = paste('g<-plot_grid(groupparamplot_a,
                 groupparamplot_b,
                 groupparamplot_c,
                 groupparamplot_d,
                 align = \'vh\',
                 labels = c(\'A\',\'B\',\'C\',\'D\'),
                 hjust=-1,
                 nrow=2)',sep='')
  eval(parse(text=string))
  string = paste('add_legend <- get_legend(groupparamplot_a+theme(legend.position=\'right\'))',sep='')
  eval(parse(text=string))
  
  # make_equation = ggplot()+geom_point(aes(x=0,y=0),size=0)+annotate('text',x=0,y=0,label=TeX('$\\dot{e} = a+\\frac{bm^c}{d}$'),size=40)#+theme_void()
  # add_legend <- plot_grid(make_equation,add_legend,ncol=1,rel_heights = .2,.8)
  g<-plot_grid(g,
               add_legend,
               nrow=1,
               rel_widths = c(.8,.2))
  g<-g+annotate('text',x=.9,y=.85,label=TeX('$\\dot{e} = a+\\frac{bm^c}{d}$'),size=9)
  title <- ggdraw() + 
    draw_label("Parameter Fitting",
               fontface = 'bold')
  g<-plot_grid(title, 
               g, 
               ncol = 1, 
               rel_heights = c(0.1, 1))
  eval(parse(text = paste('groupparamplot<-g',sep='')))
  setwd(paste(parent_fold,'/2019 Graphs/param_plots',sep=''))
  filename = paste('groupparamplot.pdf',sep='')
  string = paste('ggsave(filename,plot=groupparamplot,useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
  eval(parse(text=string))
}

#================================================================
# Extra stuff For making the plots
# g<-plot_grid(fitplot_sumactstate,fitplot_sumactstate2,align = 'vh',labels = c('a','b'),hjust = -1,nrow=1)
# legend_g <- get_legend(fitplot_sumactstate+theme(legend.position='right'))
# p<-plot_grid(g,legend_g,ncol=2,rel_widths=c(.8,.2))


#=================== Compare RNG ==============================
# Comparing Between RNG and Non RNG
sumdata$rng = 0
sumdata_rng$rng = 1
rng_data = rbind(sumdata,sumdata_rng)
