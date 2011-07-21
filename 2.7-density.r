#X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
#Derivative-solving library
#All parameters must be positive. Remember, X, Y and mu all refer to numbers; rho ≤ 1 as it is a probability.
library(deSolve)
sircarrierdens<-function(time, state, parameters)     #Solve XYZ equations
{ with(as.list(c(state, parameters)),
   {
       #       N <- S + E + I + R
       dS <- mu - beta*S*I - epsilon*beta*S*C - mu*S
       dI <- beta*S*I + epsilon*beta*S*C - gamma*I - mu*I
       dC <- gamma*q*I - GAMMA*C - mu*C
       dR <- gamma*I - gamma*q*I + GAMMA*C - mu*R

       return(list(c(dS, dI, dC, dR)))
   })
}
parms<- c(
          mu = 1/(70*365),
          beta = (.0043),    # transmission rate Ro=beta/gamma
          gamma = (1/500),  # recovery rate
          epsilon =  0.1, # proportion reduction in trans from carriers
                       # compared to standard infected individuals
          q = 0.84, #proportion of infected indivs that become carriers
                #rather than fully recover
          GAMMA =  0.0001) #recovery rate associated w/ carriers (1/GAMMA is
                  #the avg time in the carrier class

N = 1.0           # total pop
S = .92           # initial proportion susceptible
I = .08           # initial proportion infected
C = .0          # initial proportion carriers
R = 0           # initial proportion recovered

yini<- c(S=S,I=I,C=C,R=R)

ND<- seq(0,60*365,by=1)             #time sequence to use

out <- as.data.frame(ode(func = sircarrierdens, y=yini, parms = parms, times=ND))
out

                                        #Plot values and add legend
graph.colors = c("green", "red","orange","blue")
matplot(out[,1], out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "SIS Model: Dens", bty = "l", lwd=1, col = graph.colors)
legend(40, 0.7,c("Susceptible", "Infected","Carriers","Recovered"), pch=1, col=graph.colors, xjust=0, yjust=2)