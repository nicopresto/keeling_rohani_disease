#S, I, and R are the initial *totals* proportions of susceptible, infected, and recovered individuals
#Derivative-solving library
library(deSolve) 

sir_density<-function(time, state, parameters)  #Solve SIR equations
	{
  	with(as.list(c(state, parameters)), 
	{
    	dS <- mu - (beta*S*I) - (mu*S)
    	dI <- (beta*S*I) - (gamma*I) - (mu*I)
    	dR <- (gamma*I) - (mu*R)
    	return(list(c(dS, dI, dR)))
  		})}

parms<-	c(beta=2, 				# transmission rate  (Ro=beta/gamma)
		gamma= 0.14,			# recovery rate
		mu = 0.1)			# birth/mortality rate
		
yini<- 	c(S=1-.000001, 				# initial proportion susceptible
		I=.000001, 			# initial proportion infected
		R=0) 				# initial proportion recovered
	
times<- seq(0,100,by=1)   	   		#time sequence to use		

out <- as.data.frame(ode(func = sir_density, y=yini, parms = parms, times=times))
out

graph.colors = c(rgb(0,1,0), rgb(1,0,0), rgb(0,0,1)) #Plot values and add legend
matplot(times, out[,-1], lty=1, type = "l", xlab = "Time", ylab = "% of Individuals", main = "SIR Model: Density", bty = "l", lwd=1, col = graph.colors)
legend(par("usr")[1],par("usr")[3],c("Susceptible", "Infected", "Recovered"), pch=1, col=graph.colors, xjust=0, yjust=2)