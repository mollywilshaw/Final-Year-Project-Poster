# install.packages('tidyverse')
# install.packages('ggplot')

# Load required packages (use above code to install packages if not done before)
library(tidyverse)
library(ggplot2)

# Define SIRS model as follows:
sirs <- function(N, beta, gamma, xi, t.end, dt){
  # Arguments:
  # - N = population size
  # - beta = rate of infection, gamma = rate of recovery, xi = rate of loss of immunity
  # - t.end = end of time period
  # - dt = time increment size
  
  # Initialise sizes of compartments:
  St <- N    
  It <- 1
  Rt <- 0 
  
  # Initialise vectors to store the values at each iteration:
  Susceptible <- c(N)
  Infected <- c(1)
  Recovered <- c(0)
  
  # Create vector of time iterations:
  t <- seq(0, t.end, dt)
  
  # Loop over the time iterations vector:
  for (i in 2:length(t)){
    # These are the solutions to the SIRS differential equations
    St <- St + (xi*Rt - beta*St*It/N)*dt
    It <- It + (beta*St*It/N - gamma*It)*dt
    Rt <- Rt + (gamma*It - xi*Rt)*dt
    
    # Add the new values to the results vector:
    Susceptible <- append(Susceptible, St)
    Infected <- append(Infected, It)
    Recovered <- append(Recovered, Rt)
  }
  
  return(data.frame(t, Susceptible, Infected, Recovered)) 
}

# Here we plot some example data with N = 67100000 (population of the UK):
data <- sirs(67100000, 0.3, 1/10, 1/152, 365, 1)
View(data)
data1 <- pivot_longer(data, cols=c(2,3,4), names_to='cat', values_to = 'n')
ggplot(data1, aes(x=t, y=n, color=cat)) + geom_line(size=1.2) + xlab('Time t') + ylab('Number of individuals in compartment') + labs(color='Compartment') + theme(legend.position='bottom')
