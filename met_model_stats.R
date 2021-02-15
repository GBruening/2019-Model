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
library(rmcorr)

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

save_plots = 0
filter_data = 0

# Do Plots? ------------------
do_linfit_plots = 1
do_param_errbar_plots = 1
do_groupparam_plots = 0
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
sumdata <-read.csv(paste(parent_fold,'/Data/Data_11-11-2020.csv',sep=''))

sumdata_buff <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_buff.csv',sep=''))
sumdata_nobuff <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_nobuff.csv',sep=''))
sumdata_mj <-read.csv(paste(parent_fold,'/2020 Data/Data_11-11-2020_mj.csv',sep=''))

if (filter_data){
  graph_folder = '/graphs_filt/'
  sumdata <- subset(sumdata, mpowernet<max(sumdata$mpowernet)*.99)
} else {
  graph_folder = '/Graphs/'
}
sumdata_rng <-read.csv(paste(parent_fold,'/Data/Data_5-9-2019_rng.csv',sep=''))


# Add squared vars ##############
values <- c('stress','stress2','force','force2','actstate','actstate2','drive','drive2','umberger','uch','bhar','lich','marg')
index <- c(1,2,3,4,5,6,7,8,9)
sumdata$minfunc <- values[match(sumdata$minfunc,index)]

# sumdata$sumtorque2 = sumdata$sumtorque^2
# sumdata$sumforceout2 = sumdata$sumforceout^2
# sumdata$sumforcemus2 = sumdata$sumforcemus^2
# sumdata$sumstress2 = sumdata$sumstress^2
# sumdata$sumactstate2 = sumdata$sumactstate^2
# sumdata$sumdrive2 = sumdata$sumdrive^2

sumdata_rng$sumtorque2 = sumdata_rng$sumtorque^2
sumdata_rng$sumforceout2 = sumdata_rng$sumforceout^2
sumdata_rng$sumforcemus2 = sumdata_rng$sumforcemus^2
sumdata_rng$sumstress2 = sumdata_rng$sumstress^2
sumdata_rng$sumactstate2 = sumdata_rng$sumactstate^2
sumdata_rng$sumdrive2 = sumdata_rng$sumdrive^2

# Add effmass2  ==============================
# Adding Column for estimated eff_mass to met data experiment.
index <- c(1,2,3,4)
values <- c(2.47,4.730,6.990,11.50)
subjects = c(1,2,3,4,5,6,7,8)

eff_mass = numeric(length(sumdata$subj))
for (i in 1:length(sumdata$subj)){
  eff_mass[i] = eff_masses_meta[match(sumdata$c[i],index),match(sumdata$subj[i],subjects)]
}

sumdata$effmass <- eff_mass
# sumdata_rng$effmass <- eff_mass

eff_masses_meta <- read.csv(paste(parent_fold,'/Data/eff_masses_meta.csv',sep=''),header=0)

# Adding Column for estimated eff_mass to preferred experiment.
index <- unique(mpdata$effmass)
eff_mass = numeric(length(mpdata$subj))
for (i in 1:length(mpdata$subject)){
  eff_mass[i] = eff_masses_meta[match(mpdata$effmass[i],index),match(mpdata$subj[i],subjects)]
}

mpdata$effmass2 <- eff_mass
mpdata$effmass <- eff_mass

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

lm_eqn <- function(lin){
  # m <- lm(y ~ x, df);
  # eq <- substitute(italic('MPower') == a~~italic(x) + b*"\n"~~italic(R)^2~"="~r2,
  #                  list(a = format(coef(lin)[2], digits = 3),
  #                       b = format(coef(lin)[1], digits = 3),
  #                       r2 = format(lin$r.squared, digits = 3)))
  # as.character(as.expression(eq));
  coefs = list(a = format(coef(lin)[2], digits = 3),
               b = format(coef(lin)[1], digits = 3),
               r2 = format(lin$r.squared, digits = 3))
  # return(paste('MPower = ',coefs$a,'*x+',coefs$b,' \n ','R^2 = ',coefs$r2,sep=''))
  return(paste('atop(italic(MPower)==~',coefs$a,'~x+',coefs$b,',italic(R)^2==',coefs$r2,')',sep=''))
}

lm_eqn_rmcorr <- function(lin){
  coefs = list(b = lin$model$coefficients[9],
               r2 = format(my.rmc$r, digits = 3))
  return(paste('atop(italic(MPower)==~',coefs$b,'~x+int,italic(R)^2==',coefs$r2,')',sep=''))
}

fitted_eq <- function(lin,a,b,c,d){
  # eq <- substitute(~~italic(R)^2~"="~r2, 
  #                  list(r2 = format(lin$r.squared, digits = 3)))
  # as.character(as.expression(eq));
  coefs = list(a = format(coef(lin)[2], digits = 3),
               b = format(coef(lin)[1], digits = 3),
               r2 = format(lin$r.squared, digits = 3))
  return(paste('atop(italic(R)^2==',coefs$r2,',Model==',
               format(a,digits=3),'+',
               format(b,digits=3),'~m^',
               format(c,digits=3),'/T^',
               format(d,digits=3),')',sep=''))
}

fitted_eq2 <- function(a,b,c,d){
  # eq <- substitute(~~italic(R)^2~"="~r2, 
  #                  list(r2 = format(lin$r.squared, digits = 3)))
  # as.character(as.expression(eq));
  return(paste('Model==',
               format(a,digits=3),'+',
               format(b,digits=3),'~m^',
               format(c,digits=3),'/T^',
               format(d,digits=3),sep=''))
}
# Begin loopin  ------------------
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
                 'sumumber',
                 'sumbhar',
                 'sumuch',
                 'sumlich',
                 'summarg')#,
                 # 'summine',
                 # 'sumhoud')

fitting_labs = c('Sum of Torque (Nm)',
                 'Sum of Torque^2 (Nm)^2',
                 'Sum of Output Force (N)',
                 'Sum of Output Force^2 (N^2)',
                 'Sum of Muscle Force (N)',
                 'Sum of Muscle Force^2 (N^2)',
                 'Sum of Stress (N/m^2)',
                 'Sum of Stress^2 (N/m^2)^2',
                 'Sum of Active State',
                 'Sum of Active State^2',
                 'Sum of Neural Drive',
                 'Sum of Neural Drive^2',
                 'Sum of Umberger (W)',
                 'Sum of Bhargava (w)',
                 'Sum of Uchida (W)',
                 'Sum of Lichtwark (W)',
                 'Sum of Margaria (W)')#,
                 # 'Sum of Minetti (W)',
                 # 'Sum of Houdjick (W)')

min_labs2 = c(('Stress (N/m^2)'),
              ('Stress^2 (N/m^2)^2'),
              ('Force (N)'),
              ('Force^2 (N^2)'),
              ('Active State'),
              ('Active State^2'),
              ('Neural Drive'),
              ('Neural Drive^2'),
              ('Umberger Model (W)'))

fitting_labs2 = c(('Torque (Nm)'),
                  ('Torque^2 (Nm)^2'),
                  ('Output Force (N)'),
                  ('Output Force^2 (N)'),
                  ('Muscle Force (N)'),
                  ('Muscle Force^2 (N^2)'),
                  ('Stress (N/m^2)'),
                  ('Stress^2 (N/m^2)^2'),
                  ('Active State'),
                  ('Active State^2'),
                  ('Neural Drive'),
                  ('Neural Drive^2'),
                  ('Umberger Model (W)'),
                  ('Bhargava Model (W)'),
                  ('Uchida Model (W)'),
                  ('Lichtwark Model (W)'),
                  ('Margaria Model (W)'))#,
                  # ('Minetti Model (W)'),
                  # ('Houdjick Model (W)'))

minfunc1 = values[1]
rsq = c()
minfuncs_rsq = c()
minfuncs_num_rsq = c()
vars_rsq = c()
vars_rsq_labs = c()
expo_rsq = c()

lm_all_rsq = c()
lm_all_minfunc = c()

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

lmplot      <- list()
resid_plots <- list()
fit_plots   <- list()
lm_all_plot <- list()
lm_all_pred <- list()


fun.fit1 <- list()
fun.fit2 <- list()
fun.fit3 <- list()
fun.fit4 <- list()

rmcorrelation = list()

resid_norm_df <- data.frame(minfunc = character(),
                            minfunc_count = double(),
                            var = character(),
                            varfit_count = double(),
                            norm_p = double(),
                            log_p = double(),
                            sqrt_p = double(),
                            stringsAsFactors = FALSE)

# names(rsq_matrix) = c('min_func','variable','expo','rsq')
rsq_matrix = matrix(,nrow=length(minparams),ncol=length(fitting_labs))
expo=1
minfunc_count = 0

# fitting_vars = fitting_vars[3]
for (minfunc1 in minparams){#c('actstate2')){
# for (minfunc1 in c('actstate2')){
  
  minfunc_count = minfunc_count + 1
  varfit_count=0
  
  for (var in fitting_vars){#c('sumbhar')){
  # for (var in c('sumforcemus2')){
    varfit_count = varfit_count+1
    plotdata=filter(sumdata,minfunc==paste(minfunc1))
    eval(parse(text = paste('plotdata$var = plotdata$',var,sep='')))
    eval(parse(text = paste('fit_',minfunc1,'_',var,'=
                            summary(lm(plotdata$mpowernet ~ plotdata$',var,'))',sep='')))
    
    rmcorrelation[[minfunc1]][[var]] = eval(parse(text = paste('rmcorr(subj,',var,',mpowernet,plotdata)',sep='')))
    eval(parse(text = paste('rmcorr_',minfunc1,'_',var,'<-rmcorrelation[[minfunc1]][[var]]',sep='')))
    
    # print(paste('plot_',minfunc1,'_',var,sep=''))
    int = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[1]',sep='')))
    slope = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))
    eval(parse(text = paste('tempdata = filter(sumdata,minfunc==\'',minfunc1,'\')',sep='')))
    tempdata['temp_fitted'] = int + slope*tempdata[var]
    tempdata['lintransform'] = slope*tempdata[var]
    tempdata['unfitted'] = tempdata[var]
    tempdata['error'] = tempdata['mpowernet']-tempdata['temp_fitted']
    # tempdata['metabolic'] = 1
    # tempdata['fitted'] = 0
    # tempdata['unfitted'] = 2
    tempdata = rbind(cbind(tempdata$c,
                           tempdata$speed,
                           tempdata$effmass,
                           tempdata$movedur,
                           tempdata$mpowernet,
                           tempdata$temp_fitted,
                           tempdata$unfitted,
                           tempdata$lintransform,
                           tempdata$error))
    tempdata=data.frame(tempdata)
    colnames(tempdata)=c('c','speed','effmass','movedur','mpowernet','fitted','unfitted','lintransform','error')
    
    # lin model plots ==========================================
    rsq_val = eval(parse(text=paste('round(fit_',minfunc1,'_',var,'$r.squared,digits = 3)',sep='')))
    # rsq_val = round(rmcorrelation[[minfunc1]][[var]]$r,digits = 3)
    if (do_linfit_plots){
      titlestr = paste('Min ',min_labs2[minfunc_count],'\nR^2 = ',rsq_val,sep='')
      
      # # This is Using RMCORR
      # g <- eval(parse(text = paste('ggplot(plotdata,aes(x = var,
      #                                                   y = mpowernet,
      #                                                   group = factor(subj)))+
      #                                   geom_point(aes(color = factor(c)))+
      #                                   geom_line(aes(y = rmcorr_',minfunc1,'_',var,'$model$fitted.values),
      #                                             linetype = 1)',sep='')))
      # g <- g+
      #         labs(title = titlestr,
      #              x = fitting_labs[varfit_count],
      #              y= 'Net Metabolic Power (W)')+
      #         scale_color_manual(values  = mass_colors, name = "Added\nMass (kg)", labels = c("0 lbs", "5 lbs", "10 lbs","20 lbs"))+
      #         theme_classic+
      #         annotate(geom ='text',
      #                  label = eval(parse(text = paste('lm_eqn_rmcorr(rmcorrelation[[minfunc1]][[var]])',sep=''))),
      #                  x=-Inf,
      #                  y=150,
      #                  hjust=-.1,
      #                  vjust=1,
      #                  parse = TRUE,
      #                  size=4,family= theme_get()$text[["family"]])+
      #         ylim(low=0,high=200)+
      #         geom_abline(slope=1,intercept=0,show.legend = NA,linetype = 'dashed')

      # This is colored by massa
      g <- ggplot(data=plotdata,aes(x=var,y=mpowernet))+
            geom_point(data = plotdata,
                       aes(x=var,y=mpowernet,fill=factor(c)),
                       shape=21,
                       color='black',
                       size=3)+
            geom_smooth(data=plotdata,aes(x=var,y=mpowernet), color = 'black',method= 'lm')+
                     labs(title = titlestr,
                          x = fitting_labs[varfit_count],
                          y= 'Net Metabolic Power (W)')+
            scale_fill_manual(values  = mass_colors, name = "Added\nMass (kg)", labels = c("0 lbs", "5 lbs", "10 lbs","20 lbs"))+
            theme_classic+
            annotate(geom='text',
                     label=eval(parse(text = paste('lm_eqn(fit_',minfunc1,'_',var,')',sep=''))),
                     x=-Inf,
                     y=150,
                     hjust=-.1,
                     vjust=1,
                     parse = TRUE,
                     size=4,family= theme_get()$text[["family"]])+
            ylim(low=0,high=200)+
            geom_abline(slope=1,intercept=0,show.legend = NA,linetype = 'dashed')
                                 
      string = paste('lmplot[[\'',minfunc1,'\']][[\'',var,'\']] <- g',sep='')
      eval(parse(text = string))
      eval(parse(text = paste('lmplot[[\'',minfunc1,'\']][[\'',var,'\']]',sep='')))
      
      ggplot(data=tempdata)+
        geom_point(aes(x = movedur,
                       y = error,
                       color = factor(c),
                       shape = factor(speed)),
                   size=5)+
        labs(title = titlestr,
             x = paste('Movement Duration, ',fitting_labs[varfit_count],sep=''),
             y = 'Residual (W)',
             shape = 'Speed Condition')+
        scale_color_manual(values  = mass_colors,
                           name = "Added Mass (kg)",
                           labels = c("0 lbs", "3 lbs", "5 lbs","8 lbs"))+
        theme_classic
    }
    
    eval(parse(text = paste('lin_mod <- lm(plotdata$mpowernet ~ plotdata$',var,')')))
    eval(parse(text = paste('lin_mod_log <- lm(plotdata$mpowernet ~ log(plotdata$',var,'))')))
    eval(parse(text = paste('lin_mod_sqrt <- lm(plotdata$mpowernet ~ sqrt(plotdata$',var,'))')))
    if (do_diag_plots){
      # Residuals plots
      # str = paste('g<-autoplot(lm(plotdata$mpowernet ~ plotdata$',var,'))',sep='')
      # eval(parse(text = str))
      
      title <- ggdraw() +
        draw_label(
          paste(fitting_labs2[varfit_count],', Norm test (p = ',
                signif(shapiro.test(resid(lin_mod))$p.value,3),')',
                sep=''),
          fontface = 'bold', x = 0, hjust = -.25) + theme(plot.margin = margin(0, 0, 0, 7))
      
      # g <- resid_panel(lin_mod, 
      #                  plots = c('resid','qq','hist'),
      #                  nrow = 1,
      #                  smoother = TRUE,
      #                  qqbands=TRUE,
      #                  theme = 'classic')
      g <- plot_grid(resid_compare(models = list(lin_mod),
                              plots= c('resid'),
                              smoother = TRUE,
                              qqbands = TRUE,
                              theme = 'classic'),
                     resid_compare(models = list(lin_mod),
                                   plots= c('qq'),
                                   smoother = TRUE,
                                   qqbands = TRUE,
                                   theme = 'classic'),
                ggplot(data.frame(data = resid(lin_mod)), aes(x=data))+
                  geom_histogram(aes(y=..density..),
                                 color='black',fill='lightgrey',bins=30)+
                  labs(x='Residuals',y='Density')+
                  stat_density(geom = "line", color = 'blue',size = 1)+
                  theme_classic_old+
                  theme(plot.title = element_text(face='bold'))+
                  labs(title = 'Histogram'),
                nrow=1)
      
      resid_plots[[minfunc1]][[var]] <- plot_grid(title, 
                                                  g, 
                                                  nrow=2,
                                                  align = 'v',
                                                  rel_heights = c(.1,.9))
      
      # Resid compare graphs ================
      g_n <- plot_grid(resid_compare(models = list(lin_mod),
                                     plots= c('resid','qq'),
                                     smoother = TRUE,
                                     qqbands = TRUE,
                                     theme = 'classic'),
                       ggplot(data.frame(data = resid(lin_mod)), aes(x=data))+
                         geom_histogram(aes(y=..density..),
                                        color='black',fill='lightgrey',bins=30)+
                         labs(x='Residuals',y='Density')+
                         stat_density(geom = "line", color = 'blue',size = 1)+
                         theme_classic_old+
                         theme(plot.title = element_text(face='bold'))+
                         labs(title = 'Histogram'),
                       nrow=2,
                       rel_heights = c(.66,.33))
      g_log <- plot_grid(resid_compare(models = list(lin_mod_log),
                                       plots= c('resid','qq'),
                                       smoother = TRUE,
                                       qqbands = TRUE,
                                       theme = 'classic'),
                         ggplot(data.frame(data = resid(lin_mod_log)), aes(x=data))+
                           geom_histogram(aes(y=..density..),
                                          color='black',fill='lightgrey',bins=30)+
                           labs(x='Residuals',y='Density')+
                           stat_density(geom = "line", color = 'blue',size = 1)+
                           theme_classic_old+
                           theme(plot.title = element_text(face='bold'))+
                           labs(title = 'Histogram'),
                         nrow=2,
                         rel_heights = c(.66,.33))
      g_sqrt <- plot_grid(resid_compare(models = list(lin_mod_sqrt),
                                        plots= c('resid','qq'),
                                        smoother = TRUE,
                                        qqbands = TRUE,
                                        theme = 'classic'),
                          ggplot(data.frame(data = resid(lin_mod_sqrt)), aes(x=data))+
                            geom_histogram(aes(y=..density..),
                                           color='black',fill='lightgrey',bins=30)+
                            labs(x='Residuals',y='Density')+
                            stat_density(geom = "line", color = 'blue',size = 1)+
                            theme_classic_old+
                            theme(plot.title = element_text(face='bold'))+
                            labs(title = 'Histogram'),
                          nrow=2,
                          rel_heights = c(.66,.33))
      title_main <- ggdraw() +
        draw_label(
          paste(fitting_labs2[varfit_count],
                sep=''),
          fontface = 'bold',x = 0,hjust = -.25,size = 20) + theme(plot.margin = margin(0, 0, 0, 7))
      title_n <- ggdraw() +
        draw_label(
          paste('Norm test (p = ',
                signif(shapiro.test(resid(lin_mod))$p.value,3),')',
                sep=''),
          fontface = 'bold', x = 0,hjust = -.25) +theme(plot.margin = margin(0, 0, 0, 7))
      title_log <- ggdraw() +
        draw_label(paste('LOG Norm test (p = ',
                signif(shapiro.test(resid(lin_mod_log))$p.value,3),')',
                sep=''),
          fontface = 'bold', x = 0,hjust = -.25) + theme(plot.margin = margin(0, 0, 0, 7))
      title_sqrt <- ggdraw() +
        draw_label(
          paste('SQRT Norm test (p = ',
                signif(shapiro.test(resid(lin_mod_sqrt))$p.value,3),')',
                sep=''),
          fontface = 'bold', x = 0, hjust = -.25) + theme(plot.margin = margin(0, 0, 0, 7))
      g <- plot_grid(title_main,
                     plot_grid(title_n,title_log,title_sqrt, align = 'h', nrow = 1),
                     plot_grid(g_n,g_log,g_sqrt,nrow = 1, align = 'h'),
                     rel_heights = c(.05,.05,.9),
                     align = 'v',
                     nrow=3)
      
        # Save histogram plots
        setwd(paste(parent_fold,graph_folder,'resid_plots/',minfunc1,sep=''))
        string = paste('daig_min',minfunc1,'_',var,sep='')
        filename = paste(string,'.pdf',sep='')
        if (do_diag_plots){
          string = paste('ggsave(filename,plot = g ,useDingbats = FALSE,width=12,height=10,units=\'in\')',sep='')
          eval(parse(text=string))
        }
        
        # plot_grid(resid_compare(models = list(lin_mod),
        #                         plots= c('resid','qq','hist'),
        #                         smoother = TRUE,
        #                         qqbands = TRUE,
        #                         theme = 'classic'),
        #           ggplot()+
        #             geom_histogram(data=data.frame(data = resid(lin_mod)),
        #                            aes(x=data),
        #                            color='black',fill='lightgrey')+
        #             labs(x='Residuals',y='Density')+
        #             geom_density(data=data.frame(data = resid(lin_mod)),
        #                          aes(x=data))+ scale_y_continuous(sec.axis = sec_axis(~ . + 10)))
    }
    
    # Adding diagnostic p values to df
    tempdf = data.frame(minfunc1,
                        minfunc_count,
                        var,
                        varfit_count,
                        shapiro.test(resid(lin_mod))$p.value,
                        shapiro.test(resid(lin_mod_log))$p.value,
                        shapiro.test(resid(lin_mod_sqrt))$p.value)
    names(tempdf) = c('minfunc','minfunc_count','var','varfit_count','norm_p','log_p','sqrt_p')
    resid_norm_df = rbind(resid_norm_df,tempdf)
    
    # R squared stuff
    rsq=c(rsq,rsq_val)
    minfuncs_rsq=c(minfuncs_rsq,minparams[minfunc_count])
    minfuncs_num_rsq=c(minfuncs_num_rsq,minfunc_count)
    vars_rsq = c(vars_rsq,varfit_count)
    vars_rsq_labs = c(vars_rsq_labs,var)
    expo_rsq = c(expo_rsq,expo)
    rsq_matrix[minfunc_count,varfit_count] = rsq_val
    
    # nls ==========================================
    MP_model=nls(mpowernet ~ a1 + a2*(effmass^a3)/(movedur^a4),
                 data=tempdata,
                 start=list(a1=1,a2=1,a3=1,a4=1))
    modelsum_met = summary(MP_model)

    if (var=='sumforcemus2' | var=='sumforcemus' | var=='sumforceout' | var=='sumforceout'){
      scale = 100
      param_start = list(a1=1,a2=.1,a3=5,a4=3)
    } else if (var=='sumstress2' | var=='sumstress'){
      scale = 1000
      param_start = list(a1=1,a2=1000,a3=5,a4=3)
    } else {
      scale = 10
      param_start = list(a1=1,a2=.1,a3=1,a4=1)
    }
 
    if (var=='sumforcemus2' | var=='sumforcemus' | var=='sumforceout' | var=='sumforceout'){
      scale = 1000
      param_start = list(a1=1, a2=.1,a3=5,a4=3)
    } else if (var=='sumstress2' | var=='sumstress'){
      scale = 10000
      param_start = list(a1=1, a2=1000,a3=3,a4=4)
    } else {
      scale = 10
      param_start = list(a1=1, a2=20,a3=3,a4=2)
    }

    unfitted_model=nls(unfitted ~ a1 + scale*a2*(effmass^a3)/(movedur^a4),
                 data=tempdata,
                 start=param_start)#,
                 # control = list(maxiter=500))
    modelsum_unfitted = summary(unfitted_model)

    a=coef(unfitted_model)[1]
    # a=0
    b=coef(unfitted_model)[2]*scale
    c=coef(unfitted_model)[3]
    d=coef(unfitted_model)[4]
    
    afit=coef(modelsum_unfitted)[1]
    bfit=coef(modelsum_unfitted)[2]*scale
    cfit=coef(modelsum_unfitted)[3]
    dfit=coef(modelsum_unfitted)[4]
    
    a_val = c(a_val,modelsum_unfitted$coefficients[1])
    a_ste = c(a_ste,modelsum_unfitted$coefficients[5])
    b_val = c(b_val,modelsum_unfitted$coefficients[2]*scale)
    b_ste = c(b_ste,modelsum_unfitted$coefficients[6]*scale)
    c_val = c(c_val,modelsum_unfitted$coefficients[3])
    c_ste = c(c_ste,modelsum_unfitted$coefficients[7])
    d_val = c(d_val,modelsum_unfitted$coefficients[4])
    d_ste = c(d_ste,modelsum_unfitted$coefficients[8])
    
    # if (varfit_count < 13){
      eval(parse(text = paste('fun.fit1[[var]] <- function(t) ',int,'+',slope,'*(',afit,'+',bfit,'*(2.47^',cfit,')/(t^',dfit,'))',sep='')))
      eval(parse(text = paste('fun.fit2[[var]] <- function(t) ',int,'+',slope,'*(',afit,'+',bfit,'*(4.73^',cfit,')/(t^',dfit,'))',sep='')))
      eval(parse(text = paste('fun.fit3[[var]] <- function(t) ',int,'+',slope,'*(',afit,'+',bfit,'*(6.99^',cfit,')/(t^',dfit,'))',sep='')))
      eval(parse(text = paste('fun.fit4[[var]] <- function(t) ',int,'+',slope,'*(',afit,'+',bfit,'*(11.50^',cfit,')/(t^',dfit,'))',sep='')))
    # } else {
    #   eval(parse(text = paste('fun.fit1[[var]] <- function(t) ',afit,'+',bfit,'*(2.47^',cfit,')/(t^',dfit,')',sep='')))
    #   eval(parse(text = paste('fun.fit2[[var]] <- function(t) ',afit,'+',bfit,'*(4.73^',cfit,')/(t^',dfit,')',sep='')))
    #   eval(parse(text = paste('fun.fit3[[var]] <- function(t) ',afit,'+',bfit,'*(6.99^',cfit,')/(t^',dfit,')',sep='')))
    #   eval(parse(text = paste('fun.fit4[[var]] <- function(t) ',afit,'+',bfit,'*(11.50^',cfit,')/(t^',dfit,')',sep='')))
    # }
    
    fun.1 <- function(t) a+b*(2.47^c)/(t^d)
    fun.2 <- function(t) a+b*(4.73^c)/(t^d)
    fun.3 <- function(t) a+b*(6.99^c)/(t^d)
    fun.4 <- function(t) a+b*(11.50^c)/(t^d)
    
    Fitted_model=nls(fitted ~ 10*a2*(effmass^a3)/(movedur^a4), 
                     data=tempdata,
                     start=list(a2=.1,a3=1,a4=1))
    modelsum_prox = summary(Fitted_model)
    
    #================== Create parameter fits from utility model
    print(paste('Min Function: ',minfunc1))
    print(paste('Fitting Variable: ',var))
    
    if (var=='sumforcemus2' | var=='sumforcemus' | var=='sumforceout' | var=='sumforceout2'){
      scale = 100
      param_start = list(a1=1,a2=1,a3=5,a4=3)
    } else if (var=='sumstress2' | var=='sumstress'){
      scale = 1000
      param_start = list(a1=1,a2=1,a3=5,a4=3)
    } else {
      scale = 10
      param_start = list(a1=1,a2=1,a3=3,a4=2)
    }
    
    if (var == 'sumlich'){
      param_start$a3=3
    }
    
    lintransform_model=nls(lintransform ~ a1 + scale*a2*(effmass^a3)/(movedur^a4), 
                 data=tempdata,
                 start=param_start)
    modelsum_lintransform = summary(lintransform_model)
    
    minfuncs_abcd=c(minfuncs_abcd,minparams[minfunc_count])
    minfuncnum_abcd=c(minfuncnum_abcd,minfunc_count)
    vars_abcd = c(vars_abcd,varfit_count)
    vars_abcd_labs = c(vars_abcd_labs,var)

    # a_val = c(a_val,modelsum_lintransform$coefficients[1])
    # a_ste = c(a_ste,modelsum_lintransform$coefficients[5])
    # b_val = c(b_val,modelsum_lintransform$coefficients[2])
    # b_ste = c(b_ste,modelsum_lintransform$coefficients[6])
    # c_val = c(c_val,modelsum_lintransform$coefficients[3])
    # c_ste = c(c_ste,modelsum_lintransform$coefficients[7])
    # d_val = c(d_val,modelsum_lintransform$coefficients[4])
    # d_ste = c(d_ste,modelsum_lintransform$coefficients[8])
    # Fitted plots==========================================
    
    if (do_linfit_plots){
      eval(parse(text = paste('fitted_eq_str=fitted_eq(fit_',minfunc1,'_',var,',',
                              modelsum_unfitted$coefficients[1],',',
                              modelsum_unfitted$coefficients[2],',',
                              modelsum_unfitted$coefficients[3],',',
                              modelsum_unfitted$coefficients[4],')',sep='')))
      # eval(parse(text = paste('fitted_eq_str=fitted_eq(fit_',minfunc1,'_',var,',',
      #                         a_val,',',
      #                         b_val,',',
      #                         c_val,',',
      #                         d_val,')',sep='')))
      
      fit_plots[[var]] <- ggplot(tempdata,
                                 aes(x=movedur,
                                     y=lintransform,
                                     color=factor(c)))+
        geom_point(data = tempdata,
                   aes(x=movedur,
                       y=mpowernet,
                       fill=factor(c)),
                   shape=21,
                   color='black',
                   size=3)+
        # geom_point(data = tempdata,
        #            aes(x=movedur,
        #                y=lintransform,
        #                fill=factor(c)),
        #            shape=21,
        #            color='black',
        #            size=3)+
        scale_fill_manual(values = mass_colors)+
        stat_function(fun=fun.fit1[[var]],size=1.5,color=mass_colors[1])+#, linetype="dashed")+
        stat_function(fun=fun.fit2[[var]],size=1.5,color=mass_colors[2])+#, linetype="dashed")+ 
        stat_function(fun=fun.fit3[[var]],size=1.5,color=mass_colors[3])+#, linetype="dashed")+ 
        stat_function(fun=fun.fit4[[var]],size=1.5,color=mass_colors[4])+#, linetype="dashed")+
        theme_classic+
        theme(legend.position='none')+
        labs(y='Metabolic Power (W)',
             x='Movement Duration (s)',
             fill='Effective\nMass (kg)',
             title = paste('Fitting',fitting_labs2[varfit_count],'\n Min',min_labs2[minfunc_count],', R^2=',rsq_val))+
        ylim(low=0,high=200)+
        annotate('text',
                 label = fitted_eq_str,
                 x=1,
                 y=150,
                 vjust=0,
                 parse=TRUE,
                 size=4,
                 family= theme_get()$text[["family"]])
      
    }
  }
  
  #=================== Doing the caret thing ====================
  # best_min_func = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc
  # best_min_func_num = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc_num
  
  sumdata2 = filter(sumdata, minfunc == minfunc1)
  
  testdata <- sumdata2[,!names(sumdata2) %in% c('c','subj','speed','movedur','minfunc','effmass')]
  testdata2 <- testdata[,1:7]
  split <- createDataPartition(y=testdata$mpowernet,p=0.8,list=FALSE)
  traindata <- testdata2[split,]
  testdata   <- testdata2[-split,]
  
  # testlm <- train(mpowernet~.,data=traindata,method='lm')
  eval(parse(text = paste( 'lm_all_eq_',
                           minfunc1,
                           '<- lm(mpowernet ~ sumtorque + sumforceout + sumforcemus + sumstress + sumactstate + sumdrive,data=sumdata2)',
                           sep = '')))
  testlm2 <- lm(mpowernet ~ sumtorque + sumforceout + sumforcemus + sumstress + sumactstate + sumdrive,data=sumdata2)
  
  pred <- predict(testlm2, sumdata2)
  pred <- data.frame(pred = pred, mpowernet = sumdata2$mpowernet)
  
  lm_all_rsq = c(lm_all_rsq,summary(testlm2)$r.squared)
  lm_all_minfunc = c(lm_all_minfunc,minfunc_count)
  summary(testlm2)$r.squared
  
  lm_all_pred[[minfunc1]] <- ggplot(data = pred)+
    geom_smooth(aes(x=pred,y=mpowernet),method='lm')+
    geom_point(data=pred,aes(x=pred,y=mpowernet))+
    labs(x='Predicted Value', y='Metabolic Power (W)',title='Multiple linear regression')+
    annotate(geom='text',
             label = paste('atop(italic(R)^2==~',format(summary(testlm2)$r.squared, digits = 3),')',sep=''),
             x = -Inf,
             y = 150,
             hjust = -1,
             parse = TRUE)
  
  string = paste('pred_multlin_',minfunc1,sep='')
  filename = paste(string,'.pdf',sep='')
  if (do_diag_plots){
    string = paste('ggsave(filename,plot = lm_all_pred[[minfunc1]],useDingbats = FALSE,width=4,height=4,units=\'in\')',sep='')
    eval(parse(text=string))
  }
  
  all_lm_diag <- plot_grid(resid_compare(models = list(testlm2),
                                         plots= c('resid'),
                                         smoother = TRUE,
                                         qqbands = TRUE,
                                         theme = 'classic'),
                           resid_compare(models = list(testlm2),
                                         plots= c('qq'),
                                         smoother = TRUE,
                                         qqbands = TRUE,
                                         theme = 'classic'),
                           ggplot(data.frame(data = resid(testlm2)), aes(x=data))+
                             geom_histogram(aes(y=..density..),
                                            color='black',fill='lightgrey',bins=30)+
                             labs(x='Residuals',y='Density')+
                             stat_density(geom = "line", color = 'blue',size = 1)+
                             theme_classic_old+
                             theme(plot.title = element_text(face='bold'))+
                             labs(title = 'Histogram'),
                           nrow=1)
  title_n <- ggdraw() +
    draw_label(
      paste('Norm test (p = ',
            signif(shapiro.test(resid(testlm2))$p.value,3),')',
            sep=''),
      fontface = 'bold', x = 0,hjust = -.25) +theme(plot.margin = margin(0, 0, 0, 7))
  lm_all_plot[[minfunc1]] <- plot_grid(title_n,all_lm_diag,rel_heights = c(.1,.9),nrow=2)
  
  setwd(paste(parent_fold,graph_folder,'resid_plots/',minfunc1,sep=''))
  string = paste('daig_multlin_',minfunc1,sep='')
  filename = paste(string,'.pdf',sep='')
  if (do_diag_plots){
    string = paste('ggsave(filename,plot = lm_all_plot[[minfunc1]],useDingbats = FALSE,width=12,height=10,units=\'in\')',sep='')
    eval(parse(text=string))
  }
  
  #========================== ISB 2019 Plots
  plot_vars = c('sumtorque','sumforceout','sumforcemus','sumstress','sumactstate','sumdrive')
  # ISB_fitted <- plot_grid(fit_plots[[plot_vars[1]]],
  #           fit_plots[[plot_vars[2]]],
  #           fit_plots[[plot_vars[3]]],
  #           fit_plots[[plot_vars[4]]],
  #           fit_plots[[plot_vars[5]]],
  #           fit_plots[[plot_vars[6]]],
  #           nrow=2,
  #           align='vh')
  # ISB_linfit <- plot_grid(lmplot[[minfunc1]][[plot_vars[1]]]+
  #                           labs(title = 'Joint Torque (Nm)')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         lmplot[[minfunc1]][[plot_vars[2]]]+
  #                           labs(title = 'Output Force (N)')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         lmplot[[minfunc1]][[plot_vars[3]]]+
  #                           labs(title = 'Muscle Force (N)')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         lmplot[[minfunc1]][[plot_vars[4]]]+
  #                           labs(title = 'Stress (N/m^2')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         lmplot[[minfunc1]][[plot_vars[5]]]+
  #                           labs(title = 'Active State')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         lmplot[[minfunc1]][[plot_vars[6]]]+
  #                           labs(title = 'Neural Drive')+
  #                           theme(plot.title = element_text(face = "bold"),
  #                                 text = element_text(family='Arial',
  #                                                     size = 16)),
  #                         nrow=2,
  #                         align='vh')
  # setwd(parent_fold)
  # setwd('../../School Notes/Grad School/Conferences/ISB 2019')
  # ggsave('ISB_fitted.pdf',plot=ISB_fitted,useDingbats=FALSE,width=13.75,height=8.36,units='in')
  # ggsave('ISB_linfit.pdf',plot=ISB_linfit,useDingbats=FALSE,width=13.75,height=8.36,units='in')
  # setwd(parent_fold)
  
  string = paste('ggsave(filename,plot=lmgroupvar1_',minfunc1,
                 ',useDingbats = FALSE,width=9,height=6,units=\'in\')',sep='')
  
  scale = 100
  # MAKE MPOWERGROSS DUMBO
  MP_model=nls(mpowernet ~ a1 + scale*a2*(effmass^a3)/(movedur^a4), 
               data=tempdata,
               start=list(a1=1,a2=1,a3=1,a4=1))
  modelsum_met = summary(MP_model)
  
  minfuncs_abcd=c(minfuncs_abcd,minparams[minfunc_count])
  minfuncnum_abcd=c(minfuncnum_abcd,minfunc_count)
  vars_abcd = c(vars_abcd,varfit_count+1)
  vars_abcd_labs = c(vars_abcd_labs,'Metabolics')
  
  a_val = c(a_val,modelsum_met$coefficients[1])
  a_ste = c(a_ste,modelsum_met$coefficients[5])
  b_val = c(b_val,modelsum_met$coefficients[2]*scale)
  b_ste = c(b_ste,modelsum_met$coefficients[6]*scale)
  c_val = c(c_val,modelsum_met$coefficients[3])
  c_ste = c(c_ste,modelsum_met$coefficients[7])
  d_val = c(d_val,modelsum_met$coefficients[4])
  d_ste = c(d_ste,modelsum_met$coefficients[8])
  
  a=coef(MP_model)[1]
  # a=0
  b=coef(MP_model)[2]*scale
  c=coef(MP_model)[3]
  d=coef(MP_model)[4]
  
  eval(parse(text = paste('fun.fit1[[\'metabolics\']] <- function(t) ',a,'+',b,'*(2.47^',c,')/(t^',d,')',sep='')))
  eval(parse(text = paste('fun.fit2[[\'metabolics\']] <- function(t) ',a,'+',b,'*(4.73^',c,')/(t^',d,')',sep='')))
  eval(parse(text = paste('fun.fit3[[\'metabolics\']] <- function(t) ',a,'+',b,'*(6.99^',c,')/(t^',d,')',sep='')))
  eval(parse(text = paste('fun.fit4[[\'metabolics\']] <- function(t) ',a,'+',b,'*(11.50^',c,')/(t^',d,')',sep='')))
  
  fitted_eq_str=fitted_eq2(a,b,c,d)
  fit_plots[['metabolics']] <- ggplot(tempdata,aes(x=movedur,y=mpowernet,color=factor(c)))+
    geom_point(data = tempdata,
               aes(x=movedur,y=mpowernet,fill=factor(c)),
               shape=21,
               color='black',
               size=3)+
    scale_fill_manual(values = mass_colors)+
    stat_function(fun=fun.fit1[['metabolics']],size=1.5,color=mass_colors[1])+#, linetype="dashed")+
    stat_function(fun=fun.fit2[['metabolics']],size=1.5,color=mass_colors[2])+#, linetype="dashed")+ 
    stat_function(fun=fun.fit3[['metabolics']],size=1.5,color=mass_colors[3])+#, linetype="dashed")+ 
    stat_function(fun=fun.fit4[['metabolics']],size=1.5,color=mass_colors[4])+#, linetype="dashed")+
    theme_classic+theme(legend.position='none')+
    labs(y='Metabolic Power (W)',
         x='Movement Duration (s)',
         fill='Effective\nMass (kg)',
         title = 'Fitting Metabolics (W)')+
    ylim(low=0,high=200)+
    annotate('text',
             label = fitted_eq_str,
             x=1,
             y=150,
             vjust=0,
             parse=TRUE,
             size=4,
             family= theme_get()$text[["family"]])
  
  # Grouping Vars ============================
  # Variable raised to the first     'lmplot[[',minfunc1,']][[',var,']]
  if (do_linfit_plots){
    get_legend(fit_plots[[1]]+theme(legend.position='right'))
    eval(parse(text = paste('legend_mass = get_legend(lmplot[[\'',minfunc1,'\']][[\'sumtorque\']]+theme(legend.position=\'right\'))',sep='')))
    string = paste('lmgroupvar1_',minfunc1,
                   '=plot_grid(plot_grid(
                                         lmplot[[\'',minfunc1,'\']][[\'sumtorque\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumtorque2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumforceout\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumforceout2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumforcemus\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumforcemus2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumstress\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumstress2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumactstate\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumactstate2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumdrive\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumdrive2\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumumber\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumbhar\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumuch\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumlich\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'summarg\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'summine\']]+theme(legend.position=\'none\'),
                                         lmplot[[\'',minfunc1,'\']][[\'sumhoud\']]+theme(legend.position=\'none\'),
                                         align = \'vh\',
                                         labels = \'AUTO\',
                                         hjust=-1,
                                         ncol=4
                                         ),
                              legend_mass,
                              ncol = 2, rel_widths=c(.85,.15))',sep='')
    eval(parse(text=string))
    
    setwd(paste(parent_fold,graph_folder,'lm_group_var1',sep=''))
    string = paste('lmgroupvar1_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    if (save_plots){
      string = paste('ggsave(filename,plot=lmgroupvar1_',minfunc1,',useDingbats = FALSE,width=12,height=15,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    # Variable Raised to the 2nd
    string = paste('lmgroupvar2_',minfunc1,
                   '=plot_grid(lmplot[[\'',minfunc1,'\']][[\'sumtorque2\']],
                   lmplot[[\'',minfunc1,'\']][[\'sumforceout2\']],
                   lmplot[[\'',minfunc1,'\']][[\'sumforcemus2\']],
                   lmplot[[\'',minfunc1,'\']][[\'sumstress2\']],
                   lmplot[[\'',minfunc1,'\']][[\'sumactstate2\']],
                   lmplot[[\'',minfunc1,'\']][[\'sumdrive2\']],
                   align = \'vh\',
                   labels = "AUTO",
                   hjust=-1,
                   nrow=2)',sep='')
    eval(parse(text=string))
    
    setwd(paste(parent_fold,graph_folder,'lm_group_var2',sep=''))
    string = paste('lmgroupvar2_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    if (save_plots){
      string = paste('ggsave(filename,plot=lmgroupvar2_',minfunc1,',useDingbats = FALSE,width=9,height=6,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    # Residual group plots ===
    
    string = paste('residplots_',minfunc1,
                   '=plot_grid(
                   resid_plots[[\'',minfunc1,'\']][[\'sumtorque\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumtorque2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumforceout\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumforceout2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumforcemus\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumforcemus2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumstress\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumstress2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumactstate\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumactstate2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumdrive\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumdrive2\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumumber\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumbhar\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumuch\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumlich\']],
                   resid_plots[[\'',minfunc1,'\']][[\'summarg\']],
                   resid_plots[[\'',minfunc1,'\']][[\'summine\']],
                   resid_plots[[\'',minfunc1,'\']][[\'sumhoud\']],
                   align = \'vh\',
                   # labels = c(\'Tex{Torque}\',\'Tex{Torque^2}\',\'Force\',\'Force^2\',
                   #            \'Musc Force\',\'Musc Force^2\',\'Stress\'\'Stress^2\',\'ActState\',\'ActState^2\',
                   #            \'Drive\',\'Drive^2\',\'Umber\',\'Bhar\',\'Uch\',\'Lich\',\'Marg\'),
                   hjust=-1,
                   ncol=3)',sep='')
    eval(parse(text=string))
    
    setwd(paste(parent_fold,graph_folder,'resid_plots',sep=''))
    string = paste('resid_plots_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    if (save_plots){
      string = paste('ggsave(filename,plot=residplots_',minfunc1,',useDingbats = FALSE,width=18,height=15,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    
    # Save fitted to new variable name for later
    eval(parse(text = paste('fit_plots_',minfunc1,'<-fit_plots',sep='')))
    
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
                   fit_plots[[14]],
                   fit_plots[[15]],
                   fit_plots[[16]],
                   fit_plots[[17]],
                   fit_plots[[18]],
                   add_legend,
                   align = \'vh\',
                   labels = \'AUTO\',
                   hjust=-1,
                   ncol=4)',sep='')
    eval(parse(text=string))
    # eval(parse(text = paste('fittedvar1_',minfunc1,'<-plot_grid(fittedvar1_',minfunc1,',add_legend,ncol=2,rel_widths=c(.85,.15))',sep='')))
    
    setwd(paste(parent_fold,graph_folder,'fitted_var1',sep=''))
    string = paste('fittedvar1_',minfunc1,sep='')
    filename = paste(string,'.pdf',sep='')
    if (save_plots){
      string = paste('ggsave(filename,plot=fittedvar1_',minfunc1,',useDingbats = FALSE,width=16,height=12,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
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
    # setwd('C:/Users/Gary/Google Drive/2019 Model/Graphs/fitted_var2')
    # string = paste('fittedvar2_',minfunc1,sep='')
    # filename = paste(string,'.pdf',sep='')
    # string = paste('ggsave(filename,plot=fittedvar2_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
    # eval(parse(text=string))
  }
  
}


#=================== Bar Charts ==============================
print(paste('Making bar charts.'))
rsq_frame = data.frame(minfunc = minfuncs_rsq,
                       minfunc_num = minfuncs_num_rsq,
                       variable=vars_rsq,
                       expo = expo_rsq,
                       rsquared = rsq)
setwd(parent_fold)
write.csv(rsq_frame,'rsq_frame.csv')

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

if (do_linfit_plots){
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
          legend.position = c(.15,.8),
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
    coord_cartesian(ylim = c(min(rsq_frame$rsquared)*.98,max(rsq_frame$rsquared)*1.02))
  
  setwd(paste(parent_fold,graph_folder,sep=''))
  if (save_plots){
    ggsave('Rsq_matrix_expo1.pdf',plot = RSQ_bar,useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  rsq_frame$variable = factor(rsq_frame$variable,
                              levels = (1:17),
                              labels = c(1,7,2,8,3,9,4,10,5,11,6,12,13,14,15,16,17))
  
  RSQ_1to2 <- ggplot(rsq_frame,aes(fill=factor(minfunc_num),x=reorder(variable,variable),y=rsquared))+
    geom_bar(position='dodge',stat='identity')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:17), labels=parse(text=fitting_labs_split))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))
  rsq_frame$variable = factor(rsq_frame$variable,
                              labels = (1:17),
                              levels = c(1,7,2,8,3,9,4,10,5,11,6,12,13,14,15,16,17))
  setwd(paste(parent_fold,graph_folder,sep=''))
  if (save_plots){
    ggsave('Rsq_matrix_expo1_split.pdf',plot = RSQ_1to2, useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  best_min_func = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc
  best_min_func_num = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc_num
  
  rsq_frame_filt = filter(filt_rsq,minfunc==filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc)
  
  RSQ_bestmin <- ggplot(rsq_frame_filt,
                        aes(x=var_num,
                            y=rsquared))+
    geom_bar(stat='identity',
             fill=brewer.pal(9,"Spectral")[best_min_func_num])+
    geom_text(data=rsq_frame_filt,
              aes(x=var_num,
                  y=rsquared,
                  label=rsquared),
              position=position_dodge(width=0.91), 
              angle=0,
              vjust= -1,
              size = 3)+
    theme(axis.text.x = element_text(angle = 30, 
                                     hjust = 1),
          plot.title = element_text(hjust = 0.5))+
    scale_x_discrete(limits = (1:length(choose_vars)), 
                     labels=parse(text=fitting_labs[choose_vars]))+
    labs(x='Fitted Variable',
         y='R Squared Value',
         # title = parse(text = min_labs[best_min_func_num]))+
         title = eval(parse(text = min_labs2[best_min_func_num])))+
    scale_fill_manual(values = brewer.pal(9,"Spectral")[best_min_func_num])+
    coord_cartesian(ylim = c(min(rsq_frame_filt$rsquared)*.98,max(rsq_frame_filt$rsquared)*1.02))
  
  setwd(paste(parent_fold,graph_folder,sep=''))
  
  if (save_plots){
    ggsave('Rsq_matrix_expo1_act2.pdf',plot = RSQ_bestmin, useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  rsq_frame_filt$variable = factor(rsq_frame_filt$variable,
                                   levels = (1:17),
                                   labels = c(1,7,2,8,3,9,4,10,5,11,6,12,13,14,15,16,17))
  ggplot(rsq_frame_filt,aes(x=variable,y=rsquared))+
    geom_bar(stat='identity',fill='blue')+
    theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limits = (1:17), labels=parse(text=fitting_labs_split))+
    labs(x='Fitted Variable',y='R Squared Value',fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), values = colorRampPalette(brewer.pal(9,"Spectral"))(9))
  rsq_frame_filt$variable = factor(rsq_frame_filt$variable,
                                   labels = (1:17),
                                   levels = c(1,7,2,8,3,9,4,10,5,11,6,12,13,14,15,16,17))
  setwd(paste(parent_fold,graph_folder,sep=''))
  
  if (save_plots){
    ggsave('Rsq_matrix_expo1_act2_split.pdf', useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  # Residual Diagnostic Plots ===============================
  print('Making Diagnostic Plots')
  resid_normp <- ggplot(resid_norm_df,
         aes(fill=factor(minfunc),
             x=reorder(varfit_count,varfit_count),
             y=norm_p))+
    geom_bar(position='dodge',stat='identity')+
    geom_text(data=resid_norm_df,
              aes(x=reorder(varfit_count,varfit_count),
                  y=norm_p,
                  label=signif(norm_p,3)),
              position=position_dodge(width=0.91), 
              angle=90,
              hjust=-0.1,
              size = 1)+
    geom_hline(yintercept=0.05, color = 'red', linetype = 'dashed')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.15,.8),
          legend.background = element_blank(),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:17), 
                     labels=parse(text=fitting_labs))+
    labs(x='Fitted Variable',
         y='Residual Normality test p-value',
         fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), 
                      values = brewer.pal(9,"Spectral"))+
    scale_y_continuous(limits = c(0,max(resid_norm_df$norm_p)*1.1), expand = c(0,0))
  
  setwd(paste(parent_fold,graph_folder,'resid_plots',sep=''))
  if (save_plots){
    ggsave('resid_normp.pdf',plot = resid_normp, useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  resid_logp <- ggplot(resid_norm_df,
                        aes(fill=factor(minfunc),
                            x=reorder(varfit_count,varfit_count),
                            y=log_p))+
    geom_bar(position='dodge',stat='identity')+
    geom_text(data=resid_norm_df,
              aes(x=reorder(varfit_count,varfit_count),
                  y=log_p,
                  label=signif(log_p,3)),
              position=position_dodge(width=0.91), 
              angle=90,
              hjust=-0.1,
              size = 1)+
    geom_hline(yintercept=0.05, color = 'red', linetype = 'dashed')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.15,.8),
          legend.background = element_blank(),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:17), 
                     labels=parse(text=fitting_labs))+
    labs(x='Fitted Variable',
         y='Log Residual Normality test p-value',
         fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), 
                      values = brewer.pal(9,"Spectral"))+
    scale_y_continuous(limits = c(0,max(resid_norm_df$log_p)*1.1), expand = c(0,0))
  
  setwd(paste(parent_fold,graph_folder,'resid_plots',sep=''))
  if (save_plots){
    ggsave('resid_logp.pdf',plot = resid_logp, useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  resid_sqrtp <- ggplot(resid_norm_df,
                        aes(fill=factor(minfunc),
                            x=reorder(varfit_count,varfit_count),
                            y=sqrt_p))+
    geom_bar(position='dodge',stat='identity')+
    geom_text(data=resid_norm_df,
              aes(x=reorder(varfit_count,varfit_count),
                  y=sqrt_p,
                  label=signif(sqrt_p,3)),
              position=position_dodge(width=0.91), 
              angle=90,
              hjust=-0.1,
              size = 1)+
    geom_hline(yintercept=0.05, color = 'red', linetype = 'dashed')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.15,.8),
          legend.background = element_blank(),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:17), 
                     labels=parse(text=fitting_labs))+
    labs(x='Fitted Variable',
         y='Residual Normality test p-value',
         fill=minparams)+
    guides(fill=guide_legend(title="Minimization\nFunction"))+
    scale_fill_manual(labels = parse(text = min_labs), 
                      values = brewer.pal(9,"Spectral"))+
    scale_y_continuous(limits = c(0,max(resid_norm_df$sqrt_p)*1.1), expand = c(0,0))
  
  setwd(paste(parent_fold,graph_folder,'resid_plots',sep=''))
  if (save_plots){
    ggsave('resid_sqrtp.pdf',plot = resid_sqrtp, useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
  
  # Caret Bar Plots ===============================
  lmall_rsq_frame = data.frame(minfunc_num = lm_all_minfunc,
                         rsquared = lm_all_rsq)
  
  RSQ_mult_lm <- ggplot(lmall_rsq_frame,
                    aes(fill=factor(minfunc_num),
                        x=reorder(minfunc_num,minfunc_num),
                        y=rsquared))+
    geom_bar(position='dodge',stat='identity')+
    geom_text(data=lmall_rsq_frame,
              aes(x=reorder(minfunc_num,minfunc_num),
                  y=rsquared,
                  label=format(rsquared, digits = 3)),
              position=position_dodge(width=0.91), 
              angle=90,
              hjust=-0.1,
              size = 3)+
    # theme_classic(base_family='Times')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'None',
          legend.background = element_blank(),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          legend.key.size = unit(0.3, "cm"))+#,axis.line = element_line(color='black',size = 1,linetype='solid'))+
    scale_x_discrete(limits = (1:max(lmall_rsq_frame)), 
                     labels=parse(text=min_labs[1:max(lmall_rsq_frame)]))+
    labs(x='Minimization Function',
         y='R Squared Value',
         title = 'Multiple Linear Regression r-squared',
         fill=minparams)+
    # scale_fill_manual(labels = parse(text = min_labs), 
    #                   values = brewer.pal(9,"Spectral"))+
    coord_cartesian(ylim = c(min(lmall_rsq_frame$rsquared)*.98,max(lmall_rsq_frame$rsquared)*1.02))
  setwd(paste(parent_fold,graph_folder,sep=''))
  if (save_plots){
    ggsave('RSQ_mult_lm.pdf',plot = RSQ_mult_lm,useDingbats = FALSE,width=7,height=4.5,units='in')
  }
  
}

#=================== Param Error Bar Plots ==============================
print('Making the parameter fits plots.')
# Creating the parameter error bar plots
parameter_frame = data.frame(minfunc = minfuncs_abcd,
                             minfuncnum = minfuncnum_abcd,
                             variable=vars_abcd,
                             var_lab=vars_abcd_labs,
                             a_val=a_val,
                             a_ste=a_ste,
                             b_val=b_val,
                             b_ste=b_ste,
                             c_val=c_val,
                             c_ste=c_ste,
                             d_val=d_val,
                             d_ste=d_ste)


## =================== Adding p-values ==============================

minfunc_count = 0

p_vals <- function(mu1,se1,mu2,se2,df){
  T = (abs(mu1-mu2))/sqrt(se1+se2)
  return(2*pt(-T,df=184))
}

# fitting_vars = fitting_vars[3]
parameter_frame$a_pval = NA
parameter_frame$b_pval = NA
parameter_frame$c_pval = NA
parameter_frame$d_pval = NA
parameter_frame$a_pval2 = NA
parameter_frame$b_pval2 = NA
parameter_frame$c_pval2 = NA
parameter_frame$d_pval2 = NA

parameter_frame$prox_or_model = NA

f_temp = fitting_vars
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
                 'sumbhar',
                 'sumuch',
                 'sumlich',
                 'summarg',
                 'Metabolics')

for (minfunc1 in minparams){
  minfunc_count = minfunc_count + 1
  varfit_count=0
  
  temp_frame = filter(parameter_frame,minfunc==minfunc1)
  df1 = length(mpdata$speed) - 4
  tc = abs(qt(0.05/2,df1))
  
  for (var in fitting_vars){
    varfit_count = varfit_count+1
    
    m_av = temp_frame$a_val[temp_frame$var_lab=='Metabolics']
    m_as = temp_frame$a_ste[temp_frame$var_lab=='Metabolics']
    m_bv = temp_frame$b_val[temp_frame$var_lab=='Metabolics']
    m_bs = temp_frame$b_ste[temp_frame$var_lab=='Metabolics']
    m_cv = temp_frame$c_val[temp_frame$var_lab=='Metabolics']
    m_cs = temp_frame$c_ste[temp_frame$var_lab=='Metabolics']
    m_dv = temp_frame$d_val[temp_frame$var_lab=='Metabolics']
    m_ds = temp_frame$d_ste[temp_frame$var_lab=='Metabolics']
    
    for (param in c('a','b','c','d')){
      eval(parse(text = paste(param,'v = temp_frame$',param,'_val[temp_frame$var_lab==\'',var,'\']',sep='')))
      eval(parse(text = paste(param,'s = temp_frame$',param,'_ste[temp_frame$var_lab==\'',var,'\']',sep='')))
      
      eval(parse(text = paste('param_p_val = p_vals(',param,'v,',param,'s,','m_',param,'v,m_',param,'s,df1)',sep='')))
      
      parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,paste(param,'_pval',sep='')] = param_p_val
      
      if (varfit_count<12.5){
        parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,'prox_or_model'] = 1
      } else if (varfit_count<17.5){
        parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,'prox_or_model'] = 2
      } else {
        parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,'prox_or_model'] = 3
      }
      
      # Wolfe and Hanley 2002
      # Reward feedback accelerates motor learning
      # abs(mean_a-mean_b)>2*sqrt(SE1^2+SE2^2)
      # param_p_val2 = (2*sqrt(as^2+m_as^2))/abs(av-m_av)
      
      param_p_val2 = eval(parse(text = paste('(2*sqrt(',param,'s^2+m_',param,'s^2))/abs(',param,'v-m_',param,'v)',sep='')))
      if (param_p_val2<1){
        parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,paste(param,'_pval2',sep='')] = 0
      } else {
        parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,paste(param,'_pval2',sep='')] = 1
      }
      if (param == 'a'){
        param_p_val2 = eval(parse(text = paste(param,'v-1.96*',param,'s',sep='')))
        if (param_p_val2<0){
          parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,paste(param,'_pval2',sep='')] = 0
        } else {
          parameter_frame[colSums(t(parameter_frame[,c("minfunc", "var_lab")]) == c(minfunc1,var))==2,paste(param,'_pval2',sep='')] = 1
        }
      }
    }
  }
}
fitting_vars = f_temp
if (do_param_errbar_plots){
  # if (1){
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
  
  min_labs3 = c('$Stress (N/m^2)$',
                '$Stress^2 (N/m^2)^2$',
                '$Muscle Force (N)$',
                '$Muscle Force^2 (N^2)$',
                '$Active State$',
                '$Active State^2$',
                '$Neural Drive$',
                '$Neural Drive^2$',
                '$Umberger Model (W)$')
  
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
                   'sumbhar',
                   'sumuch',
                   'sumlich',
                   'summarg',
                   'Metabolics')
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
                   TeX('Marg Model $(W)$'),
                   TeX('Metabolics ($W$)'))
  
  util_param_select = c(1,2,5,6,9,10,11,12,13,14,15,16,17,18)
  fitting_vars = fitting_vars[util_param_select]
  fitting_labs = fitting_labs[util_param_select]
  parameter_frame2 = filter(parameter_frame,var_lab %in% fitting_vars)
  
  util_params = c('a','b','c','d')
  minfunc_count = 0
  
  paramplot = list()
  paramplot_pval = list()
  paramplot_pval2 = list()
  paramplot_grouped = list()
  
  for (minfunc1 in minparams){
  # for (minfunc1 in c('actstate2')){
    minfunc_count = minfunc_count + 1
    
    param_count = 0
    for (param in util_params){
    # for (param in util_params){
      param_count = param_count+1
      
      rm(g)
      temp_pframe = filter(parameter_frame2,minfunc==eval(minfunc1))
      temp_pframe$var2 = 1:length(temp_pframe$variable)
      
      string = paste('g=ggplot(temp_pframe,
                               aes(x=factor(variable),
                               y=',param,'_val,
                               color=factor(prox_or_model)))+
                     geom_point(size=3.5)+
                     geom_errorbar(aes(ymin=',param,'_val-',param,'_ste,ymax=',param,'_val+',param,'_ste),width=.3,size=1)',sep='')
      
      eval(parse(text = paste(string)))
      
      yminval=temp_pframe[length(fitting_vars),param_count*2+3]-temp_pframe[length(fitting_vars),param_count*2+4]
      ymaxval=temp_pframe[length(fitting_vars),param_count*2+3]+temp_pframe[length(fitting_vars),param_count*2+4]
      
      g<-g+
        geom_rect(xmin=0,
                  xmax=length(fitting_vars)+1,
                  ymin=yminval,
                  ymax=ymaxval,
                  alpha=.01,
                  color = '#56B4E9')+
        scale_x_discrete(labels=parse(text=fitting_labs))+
        labs(x='Effort Variable',title=eval(param),y=eval(param),
             color="Representations")+
        # scale_color_manual(labels=parse(text=fitting_labs),
        #                    values = brewer.pal(12,'Set3'))+
        scale_color_manual(labels = c('Effort Model','Metabolic Model','Metabolic data'),
                           values = c('#E69F00','#009E73','#56B4E9'))+
        geom_vline(xintercept=8.5,size=1)+
        geom_vline(xintercept=13.5,size=1)+
        theme_classic+
        theme(legend.position='none',
              axis.text.x = element_text(angle = 30, 
                                         hjust = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x=element_blank())
      
      for (row in 1:(nrow(temp_pframe)-1)){
        if (temp_pframe[row,paste(param,'_pval2',sep='')]>.5){
          g <- g + geom_text(x=row,y=temp_pframe[row,paste(param,'_val',sep='')]+max(temp_pframe[,paste(param,'_val',sep='')]/8),
                             label = "*",
                             size = 10,
                             color = 'black')
        }
      }
      
      # eval(parse(text = paste('paramplot_',minfunc1,'_',param,'<-g',sep='')))
      paramplot[[minfunc1]][[param]] <- g
      
      # Pvalue Bar chart
      rm(g)
      string = paste('g=ggplot(data = filter(temp_pframe,var_lab != \'Metabolics\'),
                               aes(x=factor(variable),
                               y=',param,'_pval,
                               fill=factor(prox_or_model)))+
                        geom_bar(position=\'dodge\',stat=\'identity\')',sep='')
      
      eval(parse(text = paste(string)))
      ylimit = eval(parse(text = paste('max(filter(temp_pframe,var_lab != \'Metabolics\')$',param,'_pval)',sep='')))+.05
      g<-g+
        geom_hline(yintercept = 0.05,color='red',linetype='dashed')+
        scale_x_discrete(labels=parse(text=fitting_labs))+
        labs(x='Effort Variable',
             title=eval(param),
             y=paste('P-Value compared to Met: ',param,sep=''),
             # fill="Effort\nProxy")+
             fill="Representation")+
        # scale_fill_manual(labels=parse(text=fitting_labs),
        # values = brewer.pal(12,'Set3'))+
        scale_fill_manual(labels = c('Effort Model','Metabolic Model'),
                          values = c('#E69F00','#009E73'))+
        geom_vline(xintercept=8.5,size=1)+
        # theme_classic+
        theme(legend.position='none',
              axis.text.x = element_text(angle = 30, 
                                         hjust = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x=element_blank())+ 
        scale_y_continuous(expand = c(0, 0), limits = c(0, ylimit))
      
      # eval(parse(text = paste('paramplot_pval_',minfunc1,'_',param,'<-g',sep='')))
      paramplot_pval[[minfunc1]][[param]] <- g
      
      # Pvalue 0 1 Bar chart
      rm(g)
      string = paste('g=ggplot(data = filter(temp_pframe,var_lab != \'Metabolics\'),
                               aes(x=factor(variable),
                               y=',param,'_pval2,
                               fill=factor(prox_or_model)))+
                        geom_bar(position=\'dodge\',stat=\'identity\')',sep='')
      
      eval(parse(text = paste(string)))
      ylimit = eval(parse(text = paste('max(filter(temp_pframe,var_lab != \'Metabolics\')$',param,'_pval2)',sep='')))+.05
      g<-g+
        geom_hline(yintercept = 0.05,color='red',linetype='dashed')+
        scale_x_discrete(labels=parse(text=fitting_labs))+
        labs(x='Effort Variable',
             title=eval(param),
             y=paste('P-Value compared to Met: ',param,sep=''),
             # fill="Effort\nProxy")+
             fill="Representation")+
        # scale_fill_manual(labels=parse(text=fitting_labs),
        # values = brewer.pal(12,'Set3'))+
        scale_fill_manual(labels = c('Effort Model','Metabolic Model'),
                          values = c('#E69F00','#009E73'))+
        geom_vline(xintercept=8.5,size=1)+
        # theme_classic+
        theme(legend.position='none',
              axis.text.x = element_text(angle = 30, 
                                         hjust = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x=element_blank())+ 
        scale_y_continuous(expand = c(0, 0), limits = c(0, ylimit))
      
      # eval(parse(text = paste('paramplot_pval2_',minfunc1,'_',param,'<-g',sep='')))
      paramplot_pval2[[minfunc1]][[param]] <- g
      
      # Combined Bar and error bar chart
      rm(g)
      pval_max = eval(parse(text = paste('sort(temp_pframe$',param,'_pval)[length(temp_pframe$a_pval)-1]',sep='')))
      param_val_max = eval(parse(text = paste('max(temp_pframe$',param,'_val)',sep='')))
      met_param_val = eval(parse(text = paste('temp_pframe$',param,'_val[length(temp_pframe$',param,'_val)]',sep='')))
      scale = met_param_val
      
      string = paste('g=ggplot()+
               geom_bar(data = filter(temp_pframe,var_lab != \'Metabolics\'),
                        aes(x=var2,
                            y=',param,'_pval*',scale,',
                            fill=factor(prox_or_model)),
                        width = 0.5,
                        position=\'dodge\',
                        stat=\'identity\')+
               geom_point(data=temp_pframe,
                          aes(x=var2,
                              y=',param,'_val,
                              color=factor(prox_or_model)),
                              size=2.5)+
               geom_errorbar(data=temp_pframe,
                             aes(x=var2,
                                 ymin=',param,'_val-',param,'_ste,
                                 ymax=',param,'_val+',param,'_ste,
                                 color=factor(prox_or_model)),
                             width=.5,size=1)',sep='')
      
      eval(parse(text = paste(string)))
      
      g<-g+
        geom_rect(data=temp_pframe,
                  xmin=0,
                  xmax=length(fitting_vars)+1,
                  ymin=yminval,
                  ymax=ymaxval,
                  alpha=.01,
                  color = '#56B4E9')+
        scale_x_continuous(breaks=c(1:length(temp_pframe$var2)),
                           labels=parse(text=fitting_labs))+
        labs(x='Effort Variable',
             title=param,
             y=paste('Parameter: ',param,sep=''),
             color="Representations",
             fill="Representations")+
        geom_hline(yintercept = 0.05*scale,linetype='dashed',color='red')+
        scale_color_manual(labels = c('Effort Model','Metabolic Model','Metabolic data'),
                           values = c('#E69F00','#009E73','#56B4E9'))+
        scale_fill_manual(labels = c('Effort Model','Metabolic Model'),
                          values = c('darkgrey','darkgrey'),
                          guide = 'none')+
        geom_vline(xintercept=8.5,size=1)+
        geom_vline(xintercept=13.5,size=1)+
        theme_classic+
        theme(legend.position='none',
              axis.text.x = element_text(angle = 30, 
                                         hjust = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x=element_blank())
      
      eval(parse(text = paste('g<-g+scale_y_continuous(name = paste(\'Parameter: \',param,sep=\'\'),
                                                      sec.axis = sec_axis(trans=~./',scale,',
                                                                          name=\'P-Value Compared to Metabolic Data\',
                                                                          breaks = c(0,.25,.5,.75,1),
                                                                          labels = round(seq(0,
                                                                                             1,#pval_max
                                                                                             length.out = 5),
                                                                                         digits = 3)))',sep='')))
      
      for (row in 1:(nrow(temp_pframe)-1)){
        if (temp_pframe[row,paste(param,'_pval2',sep='')]>.5){
          y_val = temp_pframe[row,paste(param,'_val',sep='')]+max(temp_pframe[,paste(param,'_val',sep='')]/8)
          eval(parse(text = paste('g <- g + geom_text(aes(x=',row,',y=',y_val,'),
                             label = \'*\',
                             size = 10,
                             color = \'black\')',sep='')))
        }
      }
      # eval(parse(text = paste('paramplot_grouped_',minfunc1,'_',param,'<-g',sep='')))
      paramplot_grouped[[minfunc1]][[param]] <- g
      
    }
    # Parameter Estimated Plot
    # paramplot[[minfunc1]]$a,
    g<-plot_grid(paramplot[[minfunc1]]$c,
                 paramplot[[minfunc1]]$d,
                 align = 'vh',
                 labels = c('A','C','D'),
                 hjust=-1,
                 nrow=3)
    add_legend <- get_legend(paramplot[[minfunc1]]$a+theme(legend.position='right'))
    
    g<-plot_grid(g,add_legend,nrow=1,rel_widths = c(.9,.1))
    
    title <- ggdraw() + 
      draw_label(parse(text = TeX(paste("Parameter Estimates\nMinimizing: ",min_labs3[minfunc_count],sep = ''))),
                 fontface = 'bold',
                 fontfamily = 'Times')
    g<-plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1))+
      annotate('text',x=.92,y=.85,label=TeX('$\\dot{e} = a+\\frac{bm^c}{T^d}$'),size=9,family="Times")
    eval(parse(text = paste('paramplot_',minfunc1,'<-g',sep='')))
    
    if (save_plots){
      setwd(paste(parent_fold,graph_folder,'param_plots',sep=''))
      filename = paste('paramplot_',minfunc1,'.pdf',sep='')
      string = paste('ggsave(filename,plot=paramplot_',minfunc1,',useDingbats = FALSE,width=7,height=9,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    # P_val Bar Chart
    g<-plot_grid(paramplot_pval[[minfunc1]]$a,
                 paramplot_pval[[minfunc1]]$b,
                 paramplot_pval[[minfunc1]]$c,
                 paramplot_pval[[minfunc1]]$d,
                 align = 'vh',
                 labels = c('A','B','C','D'),
                 hjust=-1,
                 nrow=2)
    add_legend <- get_legend(paramplot_pval[[minfunc1]]$a+theme(legend.position='right'))
    
    g<-plot_grid(g,add_legend,nrow=1,rel_widths = c(.9,.1))
    
    title <-ggdraw() + 
      draw_label(parse(text = TeX(paste("Parameter P-Values\nMinimizing: ",min_labs3[minfunc_count],sep = ''))),
                 fontface = 'bold',
                 fontfamily = 'Times')
    g<-plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1))+
      annotate('text',x=.92,y=.9,label=TeX('$\\dot{e} = a+\\frac{bm^c}{T^d}$'),size=9,family="Times")
    eval(parse(text = paste('paramplot_pval_',minfunc1,'<-g',sep='')))
    
    if (save_plots){
      setwd(paste(parent_fold,graph_folder,'param_plots/pvals',sep=''))
      filename = paste('paramplot_pval_',minfunc1,'.pdf',sep='')
      string = paste('ggsave(filename,plot=paramplot_pval_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    # P_val2 Bar Chart
    g<-plot_grid(paramplot_pval2[[minfunc1]]$a,
                 paramplot_pval2[[minfunc1]]$b,
                 paramplot_pval2[[minfunc1]]$c,
                 paramplot_pval2[[minfunc1]]$d,
                 align = 'vh',
                 labels = c('A','B','C','D'),
                 hjust=-1,
                 nrow=2)
    add_legend <- get_legend(paramplot_pval2[[minfunc1]]$a+theme(legend.position='right'))
    
    g<-plot_grid(g,add_legend,nrow=1,rel_widths = c(.9,.1))
    
    title <- ggdraw() + 
      draw_label(parse(text = TeX(paste("Parameter Significance\nMinimizing: ",min_labs3[minfunc_count],sep = ''))),
                 fontface = 'bold',
                 fontfamily = 'Times')
    g<-plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1))+
      annotate('text',x=.92,y=.9,label=TeX('$\\dot{e} = a+\\frac{bm^c}{T^d}$'),size=9,family="Times")
    eval(parse(text = paste('paramplot_pval2_',minfunc1,'<-g',sep='')))
    
    if (save_plots){
      setwd(paste(parent_fold,graph_folder,'param_plots/pvals01',sep=''))
      filename = paste('paramplot_pval2_',minfunc1,'.pdf',sep='')
      string = paste('ggsave(filename,plot=paramplot_pval2_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
      eval(parse(text=string))
    }
    
    # Grouped Bar and error bar
    g<-plot_grid(paramplot_grouped[[minfunc1]]$a,
                 paramplot_grouped[[minfunc1]]$b,
                 paramplot_grouped[[minfunc1]]$c,
                 paramplot_grouped[[minfunc1]]$d,
                 align = 'vh',
                 labels = c('A','B','C','D'),
                 hjust=-1,
                 nrow=2)
    add_legend <- get_legend(paramplot_grouped[[minfunc1]]$a+theme(legend.position='right'))
    
    g<-plot_grid(g,add_legend,nrow=1,rel_widths = c(.9,.1))
    
    title <- ggdraw() + 
      draw_label(parse(text = TeX(paste("Parameter Est and P-Val\nMinimizing: ",min_labs3[minfunc_count],sep = ''))),
                 fontface = 'bold',
                 fontfamily = 'Times')
    g<-plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1))+
      annotate('text',x=.92,y=.90,label=TeX('$\\dot{e} = a+\\frac{bm^c}{T^d}$'),size=9,family="Times")
    eval(parse(text = paste('paramplot_grouped_',minfunc1,'<-g',sep='')))
    
    if (save_plots){
      setwd(paste(parent_fold,graph_folder,'param_plots/combined',sep=''))
      filename = paste('paramplot_grouped_',minfunc1,'.pdf',sep='')
      string = paste('ggsave(filename,plot=paramplot_grouped_',minfunc1,',useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
      eval(parse(text=string))
    }
  }
}

#=================== Every param ==============================

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
                 'sumbhar',
                 'sumuch',
                 'sumlich',
                 'summarg',
                 'Metabolics')
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
                 TeX('Marg Model $(W)$'),
                 TeX('Metabolics ($W$)'))
# Create the plots where every min func is on every plot
if (do_groupparam_plots){
  util_param_select = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
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
      theme_classic+
      theme(legend.position='none',axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x=element_blank())
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
               rel_widths = c(.9,.1))
  g<-g+annotate('text',x=.92,y=.85,label=TeX('$\\dot{e} = a+\\frac{bm^c}{d}$'),size=9,family="Times")
  title <- ggdraw() + 
    draw_label("Parameter Fitting",
               fontface = 'bold',fontfamily = 'Times')
  g<-plot_grid(title, 
               g, 
               ncol = 1, 
               rel_heights = c(0.1, 1))
  eval(parse(text = paste('groupparamplot<-g',sep='')))
  setwd(paste(parent_fold,graph_folder,'param_plots',sep=''))
  filename = paste('groupparamplot.pdf',sep='')
  if (save_plots){
    string = paste('ggsave(filename,plot=groupparamplot,useDingbats = FALSE,width=14,height=9,units=\'in\')',sep='')
    eval(parse(text=string))
  }
}

# ================== Extra stuff For making the plots ===========
# 
# g<-plot_grid(fit_sumactstate,
#              fit_sumactstate2,
#              align = 'vh',
                  
#              labels = c('a','b'),
#              hjust = -1,nrow=1)
# legend_g <- get_legend(fitplot_sumactstate+theme(legend.position='right'))
# p<-plot_grid(g,legend_g,ncol=2,rel_widths=c(.8,.2))

#=================== Met power Density  =========================
print('Making metabolic power density.')
d = density(sumdata$mpowernet)

metpower_kdf <-ggplot(sumdata,aes(x=mpowernet))+
  geom_density(color = 'darkblue',
               fill='lightblue')+
  labs(x = 'Net Metabolic Power (W)', 
       y = 'Kernal Density Function')+
  geom_vline(aes(xintercept=mean(sumdata$mpowernet))
             ,color = 'red',
             linetype = 'dashed')+
  annotate('text',
           x = .5*max(d$x),
           y = .75*max(d$y),
           label = paste('Mean Value = ',mean(d$x), sep =''))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'',sep=''))
  ggsave('metpower_kdf.pdf',plot = metpower_kdf,useDingbats = FALSE,width=4,height=3,units='in')
  setwd(paste(parent_fold,graph_folder,'resid_plots',sep=''))
  ggsave('metpower_kdf.pdf',plot = metpower_kdf,useDingbats = FALSE,width=4,height=3,units='in')
}

#=================== Grouping for paper =========================
print('Making paper plots.')
## RSQaured
rsq_group <- plot_grid(RSQ_bestmin,
                       RSQ_bar,
                       labels='AUTO',
                       nrow=2,
                       align = 'vh')

if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('rsq_group.pdf',plot = rsq_group,useDingbats = FALSE,width=7,height=9,units='in')
}

## Umberger Fit
# best_min_func = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc
# best_min_func_num = filter(rsq_frame,rsquared == (max(rsq_frame$rsquared)))$minfunc_num

for (mod in c('umber','bhar','uch','lich','marg','torque2','actstate2','drive2')){
  eval(parse(text = paste(mod,'_bestmin_lm  <- lmplot[[\'',best_min_func,'\']][[\'sum',mod,'\']]+
                                               theme(legend.position = \'none\')',sep = '')))
  eval(parse(text = paste(mod,'_bestmin_fit <- fit_plots_',best_min_func,'$sum',mod,sep='')))
}

comb_umb <- plot_grid(umber_bestmin_lm,umber_bestmin_fit,
                      labels='AUTO',
                      nrow=1,
                      align = 'vh')
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('comb_umb.pdf',plot = comb_umb,useDingbats = FALSE,width=9,height=4.5,units='in')
}

legend_mass = get_legend(fit_plots[['sumtorque']]+theme(legend.position = 'right'))
## lm with fit all met models
comb_met <- plot_grid(plot_grid(umber_bestmin_lm,
                                umber_bestmin_fit,
                                bhar_bestmin_lm,
                                bhar_bestmin_fit,
                                uch_bestmin_lm,
                                uch_bestmin_fit,
                                lich_bestmin_lm,
                                lich_bestmin_fit,
                                marg_bestmin_lm,
                                marg_bestmin_fit,
                                fit_plots[['metabolics']],
                                nrow = 3,
                                align = 'vh',
                                labels = c('A','','B','','C','','D','','E','','F')
                                ),
                      legend_mass,
                      ncol = 2, rel_widths = c(.9,.1))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('comb_met.pdf',plot = comb_met,useDingbats = FALSE,width=12,height=9,units='in')
}

## lm with fit all met models
lm_met <- plot_grid(plot_grid(umber_bestmin_lm,
                                bhar_bestmin_lm,
                                uch_bestmin_lm,
                                lich_bestmin_lm,
                                marg_bestmin_lm,
                                nrow = 2,
                                align = 'vh',
                                labels = 'AUTO'),
                      legend_mass,
                      ncol = 2, rel_widths = c(.9,.1))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('lm_met.pdf',plot = lm_met,useDingbats = FALSE,width=12,height=8,units='in')
}

legend_mass = get_legend(fit_plots[['sumtorque']]+theme(legend.position = 'right'))
## Fit plots with met models
fits_met <- plot_grid(
              plot_grid(
                plot_grid(plot.new(),
                          fit_plots[['metabolics']],
                          plot.new(),
                          rel_widths = c(.2,.6,.2),
                          nrow=1,
                          labels = c('','A','')),
                plot_grid(umber_bestmin_fit,
                          bhar_bestmin_fit,
                          uch_bestmin_fit,
                          lich_bestmin_fit,
                          marg_bestmin_fit,
                          nrow = 2,
                          align = 'vh',
                          labels = c('B','C','D','E','F')),
                nrow = 2,
                align = 'h',
                labels = c(),
                rel_heights = c(.4,.8)),
              legend_mass,
              ncol = 2, rel_widths = c(.9,.1))

if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('fits_met.pdf',plot = fits_met,useDingbats = FALSE,width=12,height=12,units='in')
}

# fits_all <- plot_grid(
#     plot_grid(plot.new(),
#               fit_plots[['metabolics']],
#               plot.new(),
#               plot.new(),plot.new(),plot.new(),
#               umber_bestmin_fit,
#               bhar_bestmin_fit,
#               uch_bestmin_fit,
#               lich_bestmin_fit,
#               marg_bestmin_fit,
#               plot.new(),
#               plot.new(),plot.new(),plot.new(),
#               torque2_bestmin_fit,
#               actstate2_bestmin_fit,
#               drive2_bestmin_fit,
#               nrow = 6,
#               rel_heights = c(1,.2,1,1,.2,1),
#               align = 'vh',
#               labels = c('','A','','','','','B','C','D','E','F','','','','','G','H','I')),
#       legend_mass,
#       align = 'vh',
#       ncol = 2, rel_widths = c(.9,.1))

fits_all <- plot_grid(
  plot_grid(fit_plots[['metabolics']],
            umber_bestmin_fit,
            bhar_bestmin_fit,
            uch_bestmin_fit,
            lich_bestmin_fit,
            marg_bestmin_fit,
            torque2_bestmin_fit,
            actstate2_bestmin_fit,
            drive2_bestmin_fit,
            nrow = 3,
            align = 'vh',
            labels = 'AUTO'),
  legend_mass,
  align = 'vh',
  ncol = 2, rel_widths = c(.9,.1))

if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('fits_all.pdf',plot = fits_all,useDingbats = FALSE,width=12,height=12,units='in')
}

## Lm plots
eval(parse(text = paste('plots = lmplot[[\'',best_min_func,'\']]',sep='')))
legend_mass = get_legend(plots[['sumtorque']]+theme(legend.position = 'right'))
var1_lm_paper = plot_grid(plot_grid(
                          plots[['sumtorque']]+theme(legend.position = 'none'),
                          plots[['sumforceout']]+theme(legend.position = 'none'),
                          plots[['sumforcemus']]+theme(legend.position = 'none'),
                          plots[['sumstress']]+theme(legend.position = 'none'),
                          plots[['sumactstate']]+theme(legend.position = 'none'),
                          plots[['sumdrive']]+theme(legend.position = 'none'),
                          nrow=2,
                          labels='AUTO',
                          align = 'vh'),
                        legend_mass,
                        rel_widths = c(.9,.1))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('var1_lm.pdf',plot = var1_lm_paper,useDingbats = FALSE,width=3*3,height=3*2,units='in')
}

eval(parse(text = paste('plots = lmplot[[\'',best_min_func,'\']]',sep='')))
var2_lm_paper = plot_grid(plot_grid(plots[['sumtorque2']]+theme(legend.position = 'none'),
                                    plots[['sumforceout2']]+theme(legend.position = 'none'),
                                    plots[['sumforcemus2']]+theme(legend.position = 'none'),
                                    plots[['sumstress2']]+theme(legend.position = 'none'),
                                    plots[['sumactstate2']]+theme(legend.position = 'none'),
                                    plots[['sumdrive2']]+theme(legend.position = 'none'),
                                    nrow=2,
                                    labels='AUTO',
                                    align = 'vh'),
                        legend_mass,
                        rel_widths = c(.9,.1))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('var2_lm.pdf',plot = var2_lm_paper,useDingbats = FALSE,width=3*3,height=3*2,units='in')
}

eval(parse(text = paste('plots = lmplot[[\'',best_min_func,'\']]',sep='')))
var_lm_paper = plot_grid(plot_grid(
                                   plots[['sumtorque']]+theme(legend.position = 'none'),
                                   plots[['sumtorque2']]+theme(legend.position = 'none'),
                                   plots[['sumforcemus']]+theme(legend.position = 'none'),
                                   plots[['sumforcemus2']]+theme(legend.position = 'none'),
                                   plots[['sumactstate']]+theme(legend.position = 'none'),
                                   plots[['sumactstate2']]+theme(legend.position = 'none'),
                                   plots[['sumdrive']]+theme(legend.position = 'none'),
                                   plots[['sumdrive2']]+theme(legend.position = 'none'),
                                   nrow=4,
                                   labels='AUTO',
                                   align = 'vh'),
                          legend_mass,
                          rel_widths = c(.9,.1))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('var_lm.pdf',plot = var_lm_paper,useDingbats = FALSE,width=2*4,height=4*4,units='in')
}

## Other fits
eval(parse(text = paste('plots = fit_plots_',best_min_func,sep='')))
add_legend <- get_legend(fit_plots[[1]]+theme(legend.position='right'))
# title <- ggdraw() + draw_label("Model Fits raised to the fi", fontface='bold')
var1_paper = plot_grid(plot_grid(plots[['sumtorque']],
                                plots[['sumforceout']],
                                plots[['sumforcemus']],
                                plots[['sumstress']],
                                plots[['sumactstate']],
                                plots[['sumdrive']],
                                nrow=2,
                                labels='AUTO',
                                align = 'vh'),
                      add_legend,
                      rel_widths = c(.9,.1),
                      nrow=1)
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('var1_fits.pdf',plot = var1_paper,useDingbats = FALSE,width=(3*3)/.9,height=3*2,units='in')
}

add_legend <- get_legend(fit_plots[[1]]+theme(legend.position='right'))
var2_paper = plot_grid(plot_grid(plots[['sumtorque2']],
                                 plots[['sumforceout2']],
                                 plots[['sumforcemus2']],
                                 plots[['sumstress2']],
                                 plots[['sumactstate2']],
                                 plots[['sumdrive2']],
                                 nrow=2,
                                 labels='AUTO',
                                 align = 'vh'),
                       add_legend,
                       rel_widths = c(.9,.1),
                       nrow=1)
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('var2_fits.pdf',plot = var2_paper,useDingbats = FALSE,width=(3*3)/.9,height=3*2,units='in')
}
# abcd param plot
eval(parse(text = paste('save_param_plot = paramplot_',best_min_func,sep='')))
if (save_plots){
  setwd(paste(parent_fold,graph_folder,'Paper_graphs',sep=''))
  ggsave('paramplot.pdf',plot = save_param_plot,useDingbats = FALSE,height=(4.5*3)/.9,width=9,units='in')
}

# Make tables of correlation coefficients -----------------

library(psychometric)
expo=1
minfunc_count = 0
slope_table = list()

for (minfunc1 in minparams){
  
  minfunc_count = minfunc_count + 1
  varfit_count=0
  
  slope_table[[minfunc1]] = matrix(,nrow = length(fitting_vars)-1,ncol = 12)
  
  for (var in fitting_vars[1:length(fitting_vars)-1]){
    varfit_count = varfit_count + 1
    
    slope_table[[minfunc1]][varfit_count,1] = rmcorrelation[[minfunc1]][[var]]$r
    slope_table[[minfunc1]][varfit_count,2] = rmcorrelation[[minfunc1]][[var]]$CI[1]
    slope_table[[minfunc1]][varfit_count,3] = rmcorrelation[[minfunc1]][[var]]$CI[2]
    
    slope_table[[minfunc1]][varfit_count,4] = rmcorrelation[[minfunc1]][[var]]$model$coefficients[['Measure1']]
    
    slope_table[[minfunc1]][varfit_count,5] = eval(parse(text = paste('fit_',minfunc1,'_',var,'$r.squared',sep='')))
    slope_table[[minfunc1]][varfit_count,6] = CI.Rsq(eval(parse(text = paste('fit_',minfunc1,'_',var,'$r.squared',sep=''))),188,1)$LCL
    slope_table[[minfunc1]][varfit_count,7] = CI.Rsq(eval(parse(text = paste('fit_',minfunc1,'_',var,'$r.squared',sep=''))),188,1)$UCL
    
    slope_table[[minfunc1]][varfit_count,8]  = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))
    slope_table[[minfunc1]][varfit_count,9]  = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))-eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[4]',sep='')))
    slope_table[[minfunc1]][varfit_count,10] = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))+eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[4]',sep='')))
    
    slope_table[[minfunc1]][varfit_count,11]  = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[2]',sep='')))
    slope_table[[minfunc1]][varfit_count,12]  = eval(parse(text = paste('fit_',minfunc1,'_',var,'$coefficients[4]',sep='')))
    
  }
  rownames(slope_table[[minfunc1]]) = fitting_vars[1:17]
  colnames(slope_table[[minfunc1]]) = c('RMCorr-Rsq','RMCorr-LB','RMCorr-uB','RMCorr-Slope','RSQ','RSQ-LB','RSQ-UB','Slope','Slope-LB','Slope-UB','Slope','Slope-SE')
}


# #=================== Compare RNG ==============================
# # Comparing Between RNG and Non RNG
# sumdata$rng = 0
# sumdata_rng$rng = 1
# rng_data = rbind(sumdata,sumdata_rng)



# ctrl<-trainControl(method = 'cv',number = 10)
# lmCVFit <- train(mpowernet~.,data=testdata,method='lm',trControl=ctrl,metric="Rsquared")
# 
# library(brnn)
# testNNlm <- train(mpowernet~.,data=testdata,method='brnn',neurons = 20)