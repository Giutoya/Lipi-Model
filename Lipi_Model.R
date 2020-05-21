library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#List of Parameters of the Model:
  alpha=0.5 #Kaldor-Verdoorn parameter
  beta=0.2 #Institutional Innertia
  gamma=2 #growth elasticity to RER
  #mux=0.06 #mum=0.02  #gamma=1-mux+mum
  pi=3 # Income elasticity of demand for imports
  u=vartheta=1 #Speed of adjustment of the RER
  v=1 #Speed of ajustment of the quadractic RER
  qd=qt=3 #Desired RER for developmentalist state
  yf=0.02 #Foreign growth rate
  tau=0.2 #Share of domestic goods
  
#Initial values  
  epsilon0=2 #Initial income elasticity of exports
  q0=2 #Initial RER

# Parameters for conflicting claims  
  qg=2.1 #Desired government RER
  qw=1.3 #Desired workers RER
  h=0.3
  j=0.4
 
  
  
  ####################################################################
  #################### MAIN STRUCTURE OF THE MODEL ###################
  ####################################################################
  
  #y=(epsilon/pi)*yf+(gamma/pi)*qdot
  #epsilondot=(alpha/(1+beta))*(gamma/pi)*qdot

  #P=z*a*W
  #Phat=zhat+ahat+what
  #q=ln*(Pf*E/P)
  #qdot=Ehat+(zfhat-zhat)+(afhat-ahat)+(wfhat-what)
  #y=(e/pi)*yf+(gamma/pi)*(Ehat+(zfhat-zhat)+(afhat-ahat)+(wfhat-what))
  #epsilondot=(alpha/(1+beta))*(gamma/pi)*(Ehat+(zfhat-zhat)+(afhat-ahat)+(wfhat-what))
  
  #a=L/Y
  #sigma=(z-1)/z
  #######################################################################
  #######################################################################
  
  ######################## CASES CLOSURE ################################
  
  ### Which case? 
  ### 1: developmentalist state. 2: Conflicting claims. 3: financialization
  case=1
if (case==1)  
{

  ################ Developmentalist State #####################
  
  ############### CHARACTERISTICS  OF MODEL ################### 
  # zfhat=zhat
  # what=-ahat
  # wfhat=-afhat
  #     Then
  # qdot=Ehat
  # Ehat=vartheta*(qd-q)
  
  #y=(e/pi)*yf+(gamma/pi)*vartheta*(qd-q)
  #epsilondot=(alpha/(1+beta))*(gamma/pi)*vartheta*(qd-q)
 
  
  developm <- function (t, y, parms) {
    with(as.list(y), {
      #Differential Equations
      dq=vartheta*(qt-q)
      list(c(dq)) })
  }
  #Initial Conditions
  yini <- c(q=q0)
  #Computing the dynamic equation (ODE)
  times <- seq(from = 0, to = 99, by = 1)
  out <- ode(y = yini, times = times, func = developm,
             parms = NULL)
  #Correctly graph the equations
  q=out[,2]
 
  #Integrating epsilon in relationship to q
  #epsilond=epsilon0+(alfa/(1+beta))*(gamma/pi)*vartheta*(qt*qvar-(((qvar)^2)/2)) #Check and compare those two
  epsilond=epsilon0+(alpha/(1+beta))*(gamma/pi)*vartheta*((q-q0)^2)/2
  #ed=e0*(alfa/(1+beta))*(gamma/pi)*v*(qt*qvar-(((qvar)^2)/2))
  #y=(epsilon/pi)*yf+(gamma/pi)*vartheta*(qd-q)
  #epsilondot=(alpha/(1+beta))*(gamma/pi)*vartheta*(qd-q)
  
  #Make the graphs:
  y0=(epsilon0/pi)*yf
  y1=(1/pi)*((epsilon0)*yf+(gamma/pi)*vartheta*(qt-q0))
  yd=(yf/pi)*(epsilon0+(alpha/(1+beta))*(gamma/pi)*vartheta*((qd-q0)^2))
  
  yfvar=(0:10)/100
  y0=(epsilon0/pi)*yfvar
  y1=(1/pi)*((epsilon0)*yfvar)+(gamma/pi)*vartheta*(qd-q0)
  yd=(yfvar/pi)*(epsilon0+(alpha/(1+beta))*(gamma/pi)*vartheta*((qd-q0)^2))
  
  #plot(y0, yfvar, lwd = 2, type='l', main="", xlab="q", ylab="epsilon")
  #plot(y1, yfvar, lwd = 2, type='l', main="", xlab="q", ylab="epsilon")
  
  plot(y0,type="l",col="red")
  lines(yd,type="l", col="green")
  lines(y1,type="l", col="black")

  #Graphs
  #par(mfrow=c(2,2))
  #par(mar=c(2, 4, 1, 1))
  par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  plot(q, epsilond, lwd = 2, type='l', main="", xlab="q", ylab="epsilon")
  plot(q, lwd = 2, type='l', main="", xlab="time", ylab="q")
}  
  
if (case==1.2)  
{
  ###############################################################################
 ###########Case with a quadractic equation#####################################
  #Colocar equação dinâmica nova
  #Ehat=vartheta*(u*(qd-q)-v*(qd-q)^2)
  u=0.2
  v=0.1
  
  developm <- function (t, y, parms) {
    with(as.list(y), {
      #Differential Equations
      dq=vartheta*(u*(qd-q)-v*(qd-q)^2)
      list(c(dq)) })
  }
  #Initial Conditions
  yini <- c(q=q0)
  #Computing the dynamic equation (ODE)
  times <- seq(from = 0, to = 99, by = 1)
  out <- ode(y = yini, times = times, func = developm,
             parms = NULL)
  #Correctly graph the equations
  q=out[,2]
  
  epsilond=epsilon0+(alpha/(1+beta))*(gamma/pi)*vartheta*((u*(q-q0)^2)/2-v*((q-q0)^3)/3)
  depsilon= epsilond-epsilon0
  
  # Growth changes
  y0=(epsilon0/pi)*yfvar
  y1=(1/pi)*((epsilon0)*yfvar)+(gamma/pi)*vartheta*(u*(qd-q0)-v*(qd-q0)^2)
  yd=(yfvar/pi)*(epsilon0+(alpha/(1+beta))*(gamma/pi)*vartheta*((u*(qd-q0)^2)/2-v*((qd-q0)^3)/3))
  plot(y0,type="l",col="red")
  lines(yd,type="l", col="green")
 
  plot(q, lwd = 2, type='l', main="", xlab="time", ylab="q")
  plot(epsilond, lwd = 2, type='l', main="", xlab="time", ylab="epsilondot")

  #qt-d0<3u/2v for a positive impact on income elasticity of exports
}
  
  ################# CONFLICTING CLAIMS ########################

  ############### CHARACTERISTICS  OF MODEL ################### 
  #P=(P^tau)*(Pf*E)^(1-tau)
  #omega=1/z*a
  #what=varsigma*(Ln(1/(z*(e^(qw*(1-tau)))))-Ln(1/(z*(e^(q*(1-tau))))))
  #what=h*(q-qw) #h=varsigma*(1-tau)
  #Ehat=j*(qg-q)
  #qdot=j*(qg-q)-h*(q-qw)
  #being h=j-1
  #qdot=j*(qg-q)-(j-1)*(q-qw)
  
  j=0.6
  qg=3
  qw=2
  
    if (case==2)
{  
 
   conflic <- function (t, y, parms) {
    with(as.list(y), {
      #Differential Equations
      dq=j*(qg-q)-(j-1)*(q-qw)
      list(c(dq)) })
  }
  #Initial Conditions
  yini <- c(q=q0)
  #Computing the dynamic equation (ODE)
  times <- seq(from = 0, to = 99, by = 1)
  out <- ode(y = yini, times = times, func = conflic,
             parms = NULL)
  #Correctly graph the equations
  q=out[,2]
    
#Conflicting Claims and the RER
  
  y=(epsilon/pi)*yf+(gamma/pi)*(j*(qg-q)-(j-1)*(q-qw))
  epsilondot=(alpha/(1+beta))*(gamma/pi)*(j*(qg-q)-(j-1)*(q-qw))
  epsilone=epsilon0+(alpha/(1+beta))*(gamma/pi)*(-(qw+qg)*((1/2)*(((qw*(1-j)+qg*j)^2)-q0^2))+(qw*(1-j)+qg*j-q0)*(qw*(1-j)+qg*j)) #test
  #epsilone=epsilon0+(alpha/(1+beta))*(gamma/pi)*((qe^2)/2+(j-1)*qw*qe-j*qg*qe-(q0^2)/2-(j-1)*qw*q0+j*qg*q0) #My version - test
  #test with j=1 and qg=qd
}   
  
  ############## Financialization and neoliberal coalition #######################
if (case==3)
  {
  #rdot=rho0*(what-theta)-rho1*r
  #rdot=rho0*(h*(q-qw)-theta)-rho1*r
  
  
    #rdot=-g+rho*q-rho1*r
    #qdot=phi*(rf-r)
  
    rho0=0.02 #adjustment speed to inflation target
    h=0.8 # adjustment speed between real exchange rate
    rho=rho0*h
    theta=0.01 # Inflation target
      g=rho*qw-theta
      
    rho1=0.6 # Sensitivity of real interest rate changes to itself
    rho1=0 # Case with cycles
    
    rf=0.02  # International real interest rate
    r0=0.05
    phi=0.01 # Adjustment parameter
      
    financ <- function (t, y, parms) {
      with(as.list(y), {
        #Differential Equations
        dr=-g+rho*q-rho1*r
        dq=phi*(rf-r)
        list(c(dq, dr)) })
    }
    #Initial Conditions
    yini <- c(q=q0, r=r0)
    #Computing the dynamic equation (ODE)
    times <- seq(from = 0, to = 999, by = 1)
    out <- ode(y = yini, times = times, func = financ,
               parms = NULL)
    #Correctly graph the equations
    q=out[,2]    
    r=out[,3]
    
    plot(q,type="l",col="red")
    plot(r,type="l", col="green")
    
      
    #Steady State
  #re=fr
  #qe=(g+rho1*rf)/rho
  rf  
  (g+rho1*rf)/rho
  #case rho1=0
  
  #Case path dependence
  edot=(alpha/(1+beta1)*(gamma)/pi) #if qdot>0
  edot=(alpha/(1+beta2)*(gamma)/pi) #if qdot<0
  #beta1>beta2
  }
  
  
  
  
######################### Garbage ###############################
  

#################### Conflicting Claims #########################
# 
# confli <- function (t, y, parms) {
#   with(as.list(y), {
# 
#     #Differential Equations
#     dq=j*(qg-q)-h*(q-qw)
#     de=(alfa/(1+beta))*(gamma/pi)*(j*(qg-q)-h*(q-qw))
#     list(c(dq, de)) })
# }
# 
# #Initial Conditions
# yini <- c(q=q1, e=e1)
# 
# 
# #Computing the dynamic equation (ODE)
# times <- seq(from = 0, to = 30, by = 0.1)
# out2 <- ode(y = yini, times = times, func = confli,
#            parms = NULL)
# 
# 
# #Correctly graph the equations
# qvar2=out2[,2]
# evar2=out2[,3]
# 
# plot(qvar2, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(qvar2, qvar2, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(ycc, lwd = 2, type='l', main="", xlab="time", ylab="y")
# plot(evar2, lwd = 2, type='l', main="", xlab="time", ylab="e")
# 
# 
# 
# ycc1=(e1/pi)*yc
# ycc1
# ycc2=(1/pi)*((evar2*yc)+(gamma)*(j*(qg-qvar2)-h*(qvar2-qw)))
# ycc2
# ycc3=((tail(evar2, n=1))/pi)*yc
# ycc3
# 
# ycc=c(ycc1,ycc2,ycc3)
# 
# ######################## Create PDF with Graphs ###########
# 
# pdf("rplot.pdf") 
# par(mfrow=c(2,2))
# par(mar=c(2, 4, 1, 1))
# plot(qvar, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(qvar, qvar, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(y, lwd = 2, type='l', main="", xlab="time", ylab="y")
# plot(evar, lwd = 2, type='l', main="", xlab="time", ylab="e")
# 
# plot(qvar2, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(qvar2, qvar2, lwd = 2, type='l', main="", xlab="time", ylab="q")
# plot(ycc, lwd = 2, type='l', main="", xlab="time", ylab="y")
# plot(evar2, lwd = 2, type='l', main="", xlab="time", ylab="e")
# dev.off()
# 
# y1
# y3
# ycc1
# ycc3
