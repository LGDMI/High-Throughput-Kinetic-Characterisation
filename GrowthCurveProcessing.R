####Preamble####
rm(list = ls())
library(ggplot2)
library(xlsx)
library(minpack.lm)


####User Function declaration####
# This section contains all functions used in the high-throughput 
# data processing function (GrowthCurveFit).

#Growth curve functions

# Logistic growth model:
# Most basic growth curve model, growth slows down until carrying capacity (A) is reached
#  log(OD/OD0) = A / (1 + exp(4*mu/A * (lambda - t) + 2))  (Zwietering1990)
# with:
#  A = Carrying capacity
#  lambda = lag time
#  mu = growth rate

Logistic <- function(t,parameters){
  mu <- parameters$mu
  A <- parameters$A
  lambda <- parameters$lambda
  y <- A / (1 + exp(4*mu/A*(lambda-t)+2))
  return(y) 
}

# Richards growth curve model:
# Expanded logisitic model with shape factor, usually resulting in better fits (Zwietering, 1990)
# log(OD/OD0) = A*(1+v*exp(1+v)*exp(mu/A*(1+v)^(1+1/v)*(lambda-t)))^(-1/v)
# with:
#  A, lambda and mu the same as for the logistic model
#  v = a shape factor allowing modification of inflection point

Richards <- function(t,parameters){
  mu <- parameters$mu
  A <- parameters$A
  lambda <- parameters$lambda
  v <- parameters$v
  y <-A*(1+v*exp(1+v)*exp(mu/A*(1+v)^(1+1/v)*(lambda-t)))^(-1/v)
  return(y) 
}

# Gompertz model:
# Third commonly used type of growth curve function, containing only 3 parameters
# instead of 4. The reduction in parameter number increases identifiability of each
# parameter, and also results in all parameters having a mechanistic meaning.
# log(OD/OD0) = A * exp(-exp(mu*exp(1)/A*(lambda-t)+1))
#  with A, mu and lambda having the same meaning as for logistic and Richard's model

Gompertz <- function(t,parameters){
  mu <- parameters$mu
  A <- parameters$A
  lambda <- parameters$lambda
  y <- A * exp(-exp(mu*exp(1)/A*(lambda-t)+1))
  return(y) 
}

#Fitting functions

#nls.lm fitting
fit_mdl_lm <- function(d,modelfunction,par_init){
  # we need this temporary function for nls.lm:
  residual_fun <- function(parameters,observed){
    model_data <- modelfunction(observed$t,parameters)
    residuals <- model_data - observed$y
    # add a penalty if one of the parameters is negative
    x <- unlist(parameters)
    residuals <- residuals + sum(x[x<0])^2
    return(residuals)
  }
  fit <- nls.lm(par = par_init,
                fn = residual_fun,
                observed = d,
                control = nls.lm.control())
  print(fit)
  predicted <- data.frame('t'=d$t,
                          'y_fit'=modelfunction(d$t,fit$par))
  return(predicted)
}

#Plotting function to show fit for 1 well
single_plot <- function(fit_data,title){
  max_yf<-ceiling(max(fit_data$y_fit))
  lim<-max(max_yf,4)
  ggplot(data = fit_data, aes(x=t, y=y), environment=environment()) +
    geom_point() +
    geom_line(aes(x=t,y=y_fit), col="black",lty=2) +
    ylab(expression(ln~(OD/OD[0]))) + 
    #ylim(-1,lim) +
    xlab('t (h)') +
    expand_limits(y=0) +
    scale_y_continuous(expand=c(0,0),limits=c(-1,lim))+
    expand_limits(x=0) +
    scale_x_continuous(expand=c(0,0),limits=c(0,100))+
    ggtitle(title) +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_blank())
}

# Functions to allow the use of different growth models in GrowthCurveFit

## EmptyParDF feeds an empty data frame to CurveFitting, adjusted for each model type.
## This dataframe is then filled with the parameters obtained from fitting
## the selected model to the data from each well.

EmptyParDF <- function(model){
  {if(model=="Logistic"){
      df_par<-data.frame(mu=double(),
                    A=double(),
                    lambda=double(),
                    model=double())
  }
    else if(model=="Gompertz"){
      df_par<-data.frame(mu=double(),
                      A=double(),
                      lambda=double(),
                      model=double())
    }
    else if(model=="Richards"){
      df_par<-data.frame(mu=double(),
                      A=double(),
                      lambda=double(),
                      v=double(),
                      model=double())
    }
  else{
    stop("Error: Unknown model supplied")
  }}
return(df_par)
}

## InitialPar assembles the different parameters in a list, to be used in fit_mdl_lm
InitialPar <- function(model,mu,A,lambda){
  {if(model=="Logistic"){
    initPar<-list('mu'=mu,'A'=A,'lambda'=lambda)
  }
  else if(model=="Gompertz"){
    initPar<-list('mu'=mu,'A'=A,'lambda'=lambda)     
  }
  else if(model=="Richards"){
    initPar<-list('mu'=mu,'A'=A,'lambda'=lambda,'v'=1) 
  }
    else{
      stop("Error: Unknown model supplied")
    }}
  return(initPar)
}

## OutputPar extracts the parameters after fitting from the fitted object
## and returns them in a single vector
OutputPar <- function(model,fit){
  {if(model=="Logistic"){
    outputPar<-c(fit$par_fit$mu, fit$par_fit$A, fit$par_fit$lambda)
  }
  else if(model=="Gompertz"){
    outputPar<-c(fit$par_fit$mu, fit$par_fit$A, fit$par_fit$lambda)
  }
  else if(model=="Richards"){
    outputPar<-c(fit$par_fit$mu, fit$par_fit$A, fit$par_fit$lambda, fit$par_fit$v)
  }
    else{
      stop("Error: Unknown model supplied")
    }}
  return(outputPar)
}

####Automated growth curve fitting####

# This function automatically reads an Excel file with OD-data (corrected with blank),
# processes the data and returns the parameters of the optimised model fit for each well.
# An example dataset can be found on github.
# 
# Template-structure:
#   - Column 1: Time
#   - Column 2-end: OD-data corrected with blank; Count-data can also be used here.
#   - Row 1: ID of each experiment (e.g. Well or experiment-specific ID)
# 
# This function log-transforms the corrected OD-data, 
# fits the selected model (Richard's as default) to the log-transformed data, 
# and returns the fitted parameters for all growth curves in a single dataframe.
#
# Each log-transformed growth curve, along with the fitted curve, is printed to a PDF file,
# with one experiment per figure, allowing to control the quality of the fits.
# If the fit quality is not good enough, selecting a different model type can often
# result in a better fit (best fits obtained using Richards and Gompertz models).
# While the function makes an initial estimation of the parameters to be fitted,
# an initial estimation of the parameters can also be supplied to the function.
# This can further improve the quality of the fit.
#
# For transitioning the processed output to the kinetic modelling framework:
#  - Make an aggregate table for all experiments with the conditions applied in each well 
#     (e.g. substrate concentrations)
#  - Convert it to .csv
# Template is also available on github
# 

GrowthcurveFit<-function(filename,sheet=1,grMod="Richards",par_est=list('mu'=NULL,'A'=NULL,'lambda'=NULL)){
  #Preliminary declarations to run function
  model<-eval(parse(text=grMod))
  par<-EmptyParDF(grMod)
  title<-unlist(strsplit(filename,split='.',fixed=TRUE))[1]
  
  #Reading and pre-processing data
  d<-read.xlsx(filename,1,header=T)
  time<-d[,1]
  ID<-colnames(d)[2:ncol(d)]
  
  #Fitting of all curves to selected model type
  pdf(file=paste("Output ",title," ",grMod,".pdf",sep=""))
  for(i in 1:(ncol(d)-1)){
    OD<-d[,i+1]
    OD<-replace(OD,OD<0.0001,0.0001)
    lOD<-log(OD/OD[1])    #Log-transformation of data
    d_well<-data.frame(t=time,
                  y=lOD)
    
    #Pre-estimation of parameters to get better fits of the growth curve models
    par_est_i<-par_est
    ## ?
    {if(is.null(par_est$mu)){
      par_est_i$mu<-max(predict(smooth.spline(time,lOD),time,deriv=1)$y)
      }
    else{
      par_est_i$mu<-par_est$mu
      }
    }
    ## lambda
    ## Here, either the time at which lOD becomes greater than 0.5 (start of growth), or,
    ## half of the total length of the experiment is used as estimation for lambda
    {if(is.null(par_est$lambda)){
      par_est_i$lambda<-time[min(which(lOD>0.5))]
      if (is.na(par_est_i$lambda)){
        par_est_i$lambda<-max(time)/2
        }
      }
    else{
      par_est_i$lambda<-par_est$lambda
      }}
    
    ## A
    ## A is estimated with either the maximum lOD, or, 
    ## the lOD of a fully grown experiment with 10% incoulum (=log(10))
    {if(is.null(par_est$A)){
      par_est_i$A<-max(lOD)
      if (is.nan(par_est_i$A)){
        par_est_i$A<-log(10)
      }
      if (par_est_i$A==0){
        par_est_i$A<-log(10)
      }
    }
    else{
      par_est_i$A<-par_est$A
      }
    }
    
    #Fitting of model to data
    fit<- fit_mdl_lm(d_well,model,par_init=InitialPar(grMod,mu=par_est_i$mu,A=par_est_i$A,lambda=par_est_i$lambda))
    par[i,]<-c(OutputPar(grMod,fit),grMod)
    d_fit<- cbind(d_well, fit$y_fit)
    colnames(d_fit) <- c("t","y","y_fit")
    
    #Plot of fit to well i
    fig<-single_plot(d_fit,paste(grMod," Model Fit to Well ",ID[i],sep=""))
    print(fig)
  }
  dev.off()
  rownames(par)<-ID
  return(par)
}

####Using the functions####

#Example Dataset
Fit_Example<-GrowthcurveFit("ExampleData.xlsx")
write.xlsx(Fit_Example,"Output.xlsx",sheetName="Growth parameters_Rich",col.names=T,row.names=T)

Fit_Gompertz<-GrowthcurveFit("ExampleData.xlsx","Gompertz")
write.xlsx(Fit_Gompertz,"Output.xlsx",sheetName="Growth parameters_Gomp",col.names=T,row.names=T,append=T)

Fit_Log<-GrowthcurveFit("ExampleData.xlsx","Logistic")
write.xlsx(Fit_Log,"Output.xlsx",sheetName="Growth parameters_Log",col.names=T,row.names=T,append=T)

#### Bibliography ####

# Begot C, Desnier I, Daudin JD, Labadie JC, Lebert A. 1996. Recommendations for calculating
# growth parameters by optical density measurements. J. Microbiol. Methods 25(3):225–232.

# Zwietering MH, Jongenburger I, Rombouts FM, Van ’t Riet K. 1990. Modeling of the Bacterial 
# Growth Curve 56:1875–1881.

