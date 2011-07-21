                                        #X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
                                        #Derivative-solving library
                                        #All parameters must be positive. Remember, X, Y and mu all refer to numbers; rho ¡Ü 1 as it is a probability.
library(deSolve)
seirfreq<-function(time, state, parameters)     #Solve XYZ equations
{ with(as.list(c(state, parameters)),
   {
       N <- S + E + I + R
       dS <- mu - beta*I*S/N - mu*S
       dE <- beta*S*I/N - sigma*E - mu*E
       dI <- sigma*E -  gamma*I - mu*I
       dR <- gamma*I - mu*R
       return(list(c(dS, dE, dI, dR)))
   })
}
parms<- c(
          mu = 1/(70*365),
          beta = (.03),    # transmission rate Ro=beta/gamma
          sigma = (1/500),
          gamma= (1/500))  # recovery rate


N = 10           # total pop
S = 9           # initial no./density susceptible
E = 1
I = 0          # initial no./density infected
R = 0

yini<- c(S=S,E=E,I=I,R=R)

ND<- seq(0,60*365,by=1)             #time sequence to use

out <- as.data.frame(ode(func = seirfreq, y=yini, parms = parms, times=ND))
out

                                        #Plot values and add legend
graph.colors = c("green", "orange","red","blue")
matplot(out[,1], out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "SIS Model: Dens", bty = "l", lwd=1, col = graph.colors)
legend(40, 0.7,c("Susceptible", "Exposed","Infected", "Recovered"), pch=1, col=graph.colors, xjust=0, yjust=2)