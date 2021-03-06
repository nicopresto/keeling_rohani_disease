                                        #X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
                                        #Derivative-solving library
                                        #All parameters must be positive. Remember, X, Y and mu all refer to numbers; rho ¡Ü 1 as it is a probability.
library(deSolve)
sisfreq<-function(time, state, parameters)     #Solve XYZ equations
{ with(as.list(c(state, parameters)),
   {
       N <- X + Y
       dX <- gamma*Y - (beta*X*Y)/N
       dY <- (beta*X*Y)/N - (gamma*Y)
       return(list(c(dX, dY)))
   })
}
parms<- c(
          beta= .64,    # transmission rate Ro=beta/gamma
          gamma= 0.5)  # recovery rate


N = 10           # total pop
X = 9           # initial no./density susceptible
Y = 1          # initial no./density infected

yini<- c(X=X,Y=Y)
ND<- seq(0,100,by=1)             #time sequence to use

out <- as.data.frame(ode(func = sisfreq, y=yini, parms = parms, times=ND))
out

                                        #Plot values and add legend
graph.colors = c(rgb(0,1,0), rgb(1,0,0))
matplot(out[,1], out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "SIS Model: Freq", bty = "l", lwd=1, col = graph.colors)
legend(40, 0.7,c("Susceptible", "Infected"), pch=1, col=graph.colors, xjust=0, yjust=2)