                                       #X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
                                        #Derivative-solving library
                                        #All parameters must be positive. Remember, X, Y and mu all refer to numbers; rho ≤ 1 as it is a probability.
library(deSolve)
seirdens<-function(time, state, parameters)     #Solve XYZ equations
{ with(as.list(c(state, parameters)),
   {
                                        #       N <- S + E + I + R
       dS <- mu - beta*I*S - mu*S
       dE <- beta*S*I - sigma*E - mu*E
       dI <- sigma*E -  gamma*I - mu*I
       return(list(c(dS, dE, dI)))
   })
}
parms<- c(
          mu = 1/(70*365),
          beta = (.03),    # transmission rate Ro=beta/gamma
          sigma = (1/500),
          gamma= (1/500))  # recovery rate


N = 1.0           # total pop
S = .92           # initial no./density susceptible
E = .08
I = .0          # initial no./density infected

yini<- c(S=S,E=E,I=I)

ND<- seq(0,60*365,by=1)             #time sequence to use

out <- as.data.frame(ode(func = seirdens, y=yini, parms = parms, times=ND))
out

                                        #Plot values and add legend
graph.colors = c("green", "orange","red")
matplot(out[,1], out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "SIS Model: Dens", bty = "l", lwd=1, col = graph.colors)
legend(40, 0.7,c("Susceptible", "Exposed","Infected"), pch=1, col=graph.colors, xjust=0, yjust=2)