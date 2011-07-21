   #X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
                                        #Derivative-solving library
                                        #All parameters must be positive. Remember, X, Y and mu all refer to numbers; rho ≤ 1 as it is a probability.
library(deSolve)
sisdens<-function(time, state, parameters)     #Solve XYZ equations
{ with(as.list(c(state, parameters)),
   {
       N <- S + I
       dS <- gamma*I - (beta*S*I)
       dI <- (beta*S*I) - (gamma*I)
       return(list(c(dS, dI)))
   })
}
parms<- c(
          beta= .64,    # transmission rate Ro=beta/gamma
          gamma= 0.5)  # recovery rate


N = 1.0           # total pop
S = .999           # initial no./density susceptible
I = .111          # initial no./density infected

yini<- c(S=S,I=I)

ND<- seq(0,100,by=1)             #time sequence to use

out <- as.data.frame(ode(func = sisdens, y=yini, parms = parms, times=ND))
out

                                        #Plot values and add legend
graph.colors = c(rgb(0,1,0), rgb(1,0,0))
matplot(out[,1], out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "SIS Model: Dens", bty = "l", lwd=1, col = graph.colors)
legend(40, 0.7,c("Susceptible", "Infected"), pch=1, col=graph.colors, xjust=0, yjust=2)