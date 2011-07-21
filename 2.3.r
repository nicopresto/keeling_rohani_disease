#X, Y, and Z are the initial *totals* of susceptible, infected, and recovered number/density of individuals
#Derivative-solving library
#All parameters must be positive. Remember, X, Y and ν all refer to numbers; ρ ≤ 1 as it is a probability.
library(deSolve) 

xyz_density<-function(time, state, parameters)  	#Solve XYZ equations
{
with(as.list(c(state, parameters)), {
    	dX <- mu - (beta*X*Y) - (mu*X)
    	dY <- (beta*X*Y) - ((gamma+mu)/(1-rho))*Y
    	dZ <- (gamma*Y) - (mu*Z)
    	return(list(c(dX, dY, dZ)))
  		})}

parms<-	c(beta= 520/365,        # transmission rate  (Ro=beta/gamma)
          gamma= .14286,			# recovery rate
	  mu= 1/(70*365),		# birth/mortality rate
          rho=0.5)            # mortality probability

N = 20           # total pop
X = 18           # initial no./density susceptible
Y =2          # initial no./density infected
Z =0             # initial no./density recovered

yini<- c(X=X,Y=Y,Z=0)
	
times<- seq(0,20,by=1)   	    #time sequence to use		

out <- as.data.frame(ode(func = xyz_density, y=yini, parms = parms, times=times))
out

#Plot values and add legend
  		graph.colors = c(rgb(0,1,0), rgb(1,0,0), rgb(0,0,1))
  		matplot(times, out[,-1], lty=1, type = "l", xlab = "Time", ylab = "no./density of Individuals", main = "XYZ Model: Density", bty = "l", lwd=1, col = graph.colors)
  		legend(par("usr")[1],par("usr")[3],c("Susceptible", "Infected", "Recovered"), pch=1, col=graph.colors, xjust=0, yjust=2)