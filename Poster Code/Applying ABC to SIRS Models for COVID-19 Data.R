# install.packages('ggplot2')
# install.packages('tidyverse')
# install.packages('grid')
# install.packages('janitor')

# Load required packages (use above code to install packages if not done before)
library(ggplot2)
library(tidyverse)
library(grid)
library(janitor)

# Load combined data set
data <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data%20Cleaning/Combined%20Data.csv')
data$date <- as.Date(data$date) # Make sure R sees 'date' variable as date rather than text

# Create new data sets for each of the 3 areas considered:

# South West - Tier 1
swdata <- data|>
  filter(area_name=='Cornwall'|area_name=='Devon'|area_name=='Dorset'|area_name=='Somerset'|area_name=='Gloucestershire'|area_name=='Wiltshire'|area_name=='West of England')|>
  group_by(date)|>
  summarise(area='South West', population=sum(population), new_cases=sum(new_cases), cumulative_cases=sum(cumulative_cases))

# North East - Tier 2
nedata <- data|>
  filter(area_name=='County Durham'|area_name=='Tyne and Wear'|area_name=='Tees Valley'|area_name=='Northumberland')|>
  group_by(date)|>
  summarise(area='North East', population=sum(population), new_cases=sum(new_cases), cumulative_cases=sum(cumulative_cases))

# Merseyside - Tier 3
msdata <- data|>
  filter(area_name=='Merseyside')|>
  group_by(date)|>
  summarise(area='Merseyside', population=sum(population), new_cases=sum(new_cases), cumulative_cases=sum(cumulative_cases))

# Plot cases for all dates in data set, coloured based on the area
ggplot(rbind(swdata, nedata, msdata), aes(x=date, y=new_cases, color=area)) + geom_line() + xlab('Date') + ylab('New Daily Cases') + scale_color_discrete(name='Area', breaks=c('South West', 'North East', 'Merseyside'), type=c('red', 'yellow', 'deepskyblue'))

# Define SIRS function, slightly adapted from previous code:
sirs <- function(beta, gamma, xi, d.start, d.end, I0, R0, data){
  # We now additionally specify the start date d.start and the initial conditions I0 and R0
  # We extract the population N from the input data 'data', and dt is now 1 day, so we do not include it as an argument
  
  N <- data$population[1]

  St <- N-R0-I0 # Initial condition for S, calculated using I0 and R0 
  It <- I0
  Rt <- R0 
  
  Susceptible <- c(St) 
  Infected <- c(It)
  Recovered <- c(Rt)
  
  # Time iteration vector is now a list of dates:
  date <- seq(as.Date(d.start,'%d/%m/%Y')-1, as.Date(d.end,'%d/%m/%Y'), by='days')  
  
  for (i in 2:length(date)){
    St <- St + (xi*Rt - beta*St*It/N)
    It <- It + (beta*St*It/N - gamma*It)
    Rt <- Rt + (gamma*It - xi*Rt)
    
    Susceptible <- append(Susceptible, St)
    Infected <- append(Infected, It)
    Recovered <- append(Recovered,Rt)
  }
  
  # To compare directly with the real COVID-19 data, which is given in terms of new daily cases, we calculate the daily change in the size of the infective compartment:
  new_cases <- c()
  for(j in 2:length(Infected)){
    new_cases[j-1] <- Infected[j] - Infected[j-1]
  }
  
  # Remove the 1st date - this was only here to calculate the new cases on our desired start date:
  date <- date[-1]
  Susceptible <- Susceptible[-1]
  Infected <- Infected[-1]
  Recovered <- Recovered[-1]
  
  return(data.frame(date, Susceptible, Infected, Recovered, new_cases)) 
}

# Define ABC function as below:
abc <- function(N, epsilon, prior.beta, prior.I0, prior.R0, dates, inpt.data){ 
  # Arguments are:
  # - Number of iterations, N
  # - Tolerance, epsilon
  # - c(lower bound/upper bound) of the uniform priors for beta/I0/R0, prior.beta/prior.I0/prior.R0
  # - c(start date, end date), dates
  # - input data, inpt.data
  
  # Extract start/end date from 'dates' argument:
  d.start <- as.Date(dates[1],'%d/%m/%Y')
  d.end <- as.Date(dates[2],'%d/%m/%Y')
  
  # Create data frame of randomly generated prior values for beta, I0 and R0:
  prior <- data.frame(beta = runif(N, prior.beta[1], prior.beta[2]),
                      I0 = runif(N, prior.I0[1], prior.I0[2]),
                      R0 = runif(N, prior.R0[1], prior.R0[2]))
  # Create an empty data frame to hold the posterior values, as well as the average distance between observed/simulated data points, 'diff' 
  post <- data.frame('beta'=character(), 'I0'=character(), 'R0'=character(), 'diff'=character())
  
  # Loop over N iterations:
  for(i in 1:N){
    # Select the ith prior values as proposed parameters:
    beta <- prior[i,1] 
    I0 <- prior[i,2]
    R0 <- prior[i,3]
    
    # Run the SIRS model with these parameters:
    data <- sirs(beta, 1/10, 1/152, d.start, d.end, I0, R0, inpt.data)
    # Note: gamma and xi chosen according to the following:
    # 10 days = recovery period (via https://www.gov.uk/guidance/people-with-symptoms-of-a-respiratory-infection-including-covid-19)
    # 152 days = immunity period (via https://www.gov.uk/government/news/past-covid-19-infection-provides-some-immunity-but-people-may-still-carry-and-transmit-virus)
    
    # Calculate the absolute difference between the simulated and observed data:
    diff <- abs(as.numeric(data$new_cases) - filter(inpt.data, between(inpt.data$date, d.start, d.end))$new_cases)
    
    j <- 0 # Set counter
    
    # Count the number of simulated points that are within the chosen tolerance of the 'observed' data: 
    for(x in diff){
      if(x < epsilon){
        j <- j + 1
      } 
    }
    # If all points are within the tolerance, we accept these parameters into the posterior distribution: 
    if(j == length(diff)){
      post <- rbind(post, data.frame(beta=beta, I0=I0, R0=R0, diff=mean(diff)))
    }
  }
  return(post)
}

# First we need to choose the prior bounds:

# I0, R0:
# We will start the model on 14/10/2020 for each area, although we actually look at the data from the day before due to the 'day zero' counting convention used for the infective period:
sw.start <- swdata$cumulative_cases[which(swdata$date==as.Date('13/10/2020','%d/%m/%Y'))]
ne.start <- nedata$cumulative_cases[which(nedata$date==as.Date('13/10/2020','%d/%m/%Y'))]
ms.start <- msdata$cumulative_cases[which(msdata$date==as.Date('13/10/2020','%d/%m/%Y'))]

# People who first became infected 10 days before this (3/10/2020) are the active cases on 14/10/2020 - we use the cumulative case numbers to calculate this number:
sw.inf <- sw.start - swdata$cumulative_cases[which(swdata$date==as.Date('3/10/2020','%d/%m/%Y'))]
ne.inf <- ne.start - nedata$cumulative_cases[which(nedata$date==as.Date('3/10/2020','%d/%m/%Y'))]
ms.inf <- ms.start - msdata$cumulative_cases[which(msdata$date==as.Date('3/10/2020','%d/%m/%Y'))]
sw.inf
ne.inf
ms.inf
# To allow the same prior for each area, we set the lower bound of the I0 prior to be 6000. 

# We now attempt to estimate the upper bound for I0 using data on the proportion of positive results in samples of test swabs:

# Load and clean data:
testresults <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/20230127covid19infectionsurveydatasetsengland%20%5B9%5D.csv')[c(5,18:19),]|>
  row_to_names(1)|>
  pivot_longer(2:21, names_to = 'area', values_to = 'test_numbers')
testresults$area <- gsub('Nu', '/Nu', testresults$area)
testresults$area <- gsub('To', '/To', testresults$area)
testresults <- separate(testresults, 2, sep='/', into=c('Area','2'))
testresults$Area <- gsub('England ', 'England', testresults$Area)
testresults <- pivot_wider(testresults, names_from = 3, values_from = 4)|>
  filter(Area=='South West\n'|Area=='North East\n'|Area=='North West\n')
testresults$`Total number of tests in sample ` <- as.numeric(gsub(',', '', testresults$`Total number of tests in sample `))
testresults <- testresults|>  
  group_by(Area)|>
  summarise(n_pos = sum(as.numeric(`Number of tests positive for COVID-19`)), n_tests = sum(`Total number of tests in sample `))

# Add population data:
populations <- data.frame('Area' = c('South West\n', 'North East\n', 'North West\n'), population = c(swdata$population[1], nedata$population[1], msdata$population[1]))

# Add variables for the proportions of positive tests and estimated number of cases:
testresults <- full_join(testresults, populations)|>
  mutate(pos_prop = n_pos/n_tests)|>
  mutate(estimated_pos = population*pos_prop)|>
  summarise(Area,population,pos_prop,estimated_pos)
testresults$Area[2] <- 'Merseyside'
View(testresults)
# Highest here is over 30000 - since this is still only from a sample, we still have some uncertainty. We express this by extending the upper bound to 50000.

# R0 = number of individuals who have had COVID, are no longer active cases and had it less than 5 months ago (i.e. those who had COVID between 13/5/2020 and 3/10/2020)
sw.rec <- swdata$cumulative_cases[which(swdata$date==as.Date('3/10/2020','%d/%m/%Y'))] - swdata$cumulative_cases[which(swdata$date==as.Date('13/5/2020','%d/%m/%Y'))]
ne.rec <- nedata$cumulative_cases[which(nedata$date==as.Date('3/10/2020','%d/%m/%Y'))] - nedata$cumulative_cases[which(nedata$date==as.Date('13/5/2020','%d/%m/%Y'))]
ms.rec <- msdata$cumulative_cases[which(msdata$date==as.Date('3/10/2020','%d/%m/%Y'))] - msdata$cumulative_cases[which(msdata$date==as.Date('13/5/2020','%d/%m/%Y'))]
sw.rec
ne.rec
ms.rec
# Set lower bound for R0 to be 9000
# R0 can't be any larger than I0, so we also set the upper bound for R0 to be 50000, as there is very little other data on recovery numbers to draw from.

# We now need to decide on the prior for beta:

# First we define a slightly different ABC algorithm, which will tell us which parameter values would be rejected and accepted:

abc.priors <- function(N, epsilon, prior.beta, prior.I0, prior.R0, dates, inpt.data){ 
  
  d.start <- as.Date(dates[1],'%d/%m/%Y')
  d.end <- as.Date(dates[2],'%d/%m/%Y')
  
  # Add indicator variable 'post' to show whether prior values are in the posterior distribution or not:
  prior <- data.frame(beta = runif(N, prior.beta[1], prior.beta[2]), I0 = runif(N, prior.I0[1], prior.I0[2]), R0 = runif(N, prior.R0[1], prior.R0[2]), post = rep(NA,N))
  
  for(i in 1:N){
    beta <- prior[i,1] 
    I0 <- prior[i,2]
    R0 <- prior[i,3]
    
    data <- sirs(beta, 1/10, 1/152, d.start, d.end, I0, R0, inpt.data)
    diff <- abs(as.numeric(data$new_cases) - filter(inpt.data, between(inpt.data$date, d.start, d.end))$new_cases)
    j <- 0
    
    for(x in diff){
      if(x < epsilon){
        j <- j + 1
      } 
    }
    
    # Indicator variable is 1 if the values are accepted, 0 if not:
    ifelse(j == length(diff), prior$post[i] <- 1, prior$post[i] <- 0)
  }
  return(prior)
}

# Run ABC algorithm with twice the number of iterations than we will do for the results later:

set.seed(1) # Run to reproduce poster plots exactly

abc.sw <- abc.priors(20000, 1000, c(0,1), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), swdata)|>
  mutate(area = 'South West')
abc.ne <- abc.priors(20000, 1000, c(0,1), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), nedata)|>
  mutate(area = 'North East')
abc.ms <- abc.priors(20000, 1000, c(0,1), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), msdata)|>
  mutate(area = 'Merseyside')
abc.all <- rbind(abc.sw, abc.ne, abc.ms)

# Plot histograms of attempted beta values
prior.plot <- ggplot(abc.all, aes(x=beta, color=factor(post), fill=factor(post))) + geom_histogram(binwidth=0.01) + xlab('Attempted Beta Value') + ylab('Frequency') + scale_fill_discrete(name='', breaks=c('0','1'), labels=c('Not accepted','Accepted'), type = c('darkgrey','red')) + scale_color_discrete(type=c('darkgrey','black')) + guides(color='none') + facet_wrap(~factor(area, levels=c('South West','North East','Merseyside')), nrow=1, ncol=3)
prior.plot <- prior.plot + geom_vline(xintercept=0.1, linetype='dotted', color='black')
prior.plot <- prior.plot + geom_vline(xintercept=0.2, linetype='dotted', color='black')
prior.plot <- prior.plot + annotate('text', x=0.08, y=-10, size=2, label='0.1', color='black')
prior.plot <- prior.plot + annotate('text', x=0.18, y=-10, size=2, label='0.2', color='black')
prior.plot <- prior.plot + theme(legend.position='bottom')

# The following code is used to change the colours of the facet titles based on the region (lines 242-251 from https://github.com/tidyverse/ggplot2/issues/2096)
plot.col <- ggplot_gtable(ggplot_build(prior.plot))
striprt <- which( grepl('strip-1', plot.col$layout$name) | grepl('strip-t', plot.col$layout$name))
fills <- c('deepskyblue','yellow','red')
k <- 1
for (i in striprt){
  j <- which(grepl('rect', plot.col$grobs[[i]]$grobs[[1]]$childrenOrder))
  plot.col$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(plot.col)
# We can see that all accepted beta values are within [0.1,0.2], hence we narrow the prior between these bounds. 

# With these priors, run the ABC algorithm for each of the 3 areas:
set.seed(1) # Run to reproduce poster plots exactly
abc.sw <- abc(10000, 1000, c(0.1,0.2), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), swdata)
abc.ne <- abc(10000, 1000, c(0.1,0.2), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), nedata)
abc.ms <- abc(10000, 1000, c(0.1,0.2), c(6000,50000), c(9000,50000), c('14/10/2020','5/11/2020'), msdata)

# Create simulated data sets using the parameter estimates with the lowest average distance from the observed data:
simdata.sw <- sirs(abc.sw[which.min(abc.sw$diff),1], 1/10, 1/152, '14/10/2020', '5/11/2020', abc.sw[which.min(abc.sw$diff),2], abc.sw[which.min(abc.sw$diff),3], swdata)|>
  mutate(cat='simulated', area='South West')
simdata.ne <- sirs(abc.ne[which.min(abc.ne$diff),1], 1/10, 1/152, '14/10/2020', '5/11/2020', abc.ne[which.min(abc.ne$diff),2], abc.ne[which.min(abc.ne$diff),3], nedata)|>
  mutate(cat='simulated', area='North East')
simdata.ms <- sirs(abc.ms[which.min(abc.ms$diff),1], 1/10, 1/152, '14/10/2020', '5/11/2020', abc.ms[which.min(abc.ms$diff),2], abc.ms[which.min(abc.ms$diff),3], msdata)|>
  mutate(cat='simulated', area='Merseyside')
# Join for all areas
simdata <- rbind(simdata.sw, simdata.ne, simdata.ms)

# Filter the observed data between the dates of the first tier system (ready for plotting)
swdata.1 <- filter(swdata, between(date, as.Date('14/10/2020','%d/%m/%Y'), as.Date('5/11/2020','%d/%m/%Y')))|>
  mutate(cat='real', area='South West')
nedata.1 <- filter(nedata, between(date, as.Date('14/10/2020','%d/%m/%Y'), as.Date('5/11/2020','%d/%m/%Y')))|>
  mutate(cat='real', area='North East')
msdata.1 <- filter(msdata, between(date, as.Date('14/10/2020','%d/%m/%Y'), as.Date('5/11/2020','%d/%m/%Y')))|>
  mutate(cat='real', area='Merseyside')
# Join for all areas
obsdata <- rbind(swdata.1, nedata.1, msdata.1)

# Combine the simulated and observed data ready for plotting:
plot.data <- rbind(simdata[,c(1,5:7)], obsdata[,c(1,2,4,6)])

# Create plots of simulated and observed data for each area:
p <- ggplot(plot.data, aes(x=date, y=new_cases, color=cat, linetype=cat)) + geom_line() + xlab('Date') + ylab('New Daily Cases') + scale_linetype_manual(name='Data type', labels=c('Observed','Simulated'), values=c('solid','longdash')) + scale_color_discrete(name='Data type', type=c('black','red'), labels=c('Observed','Simulated')) + facet_wrap(~factor(area, levels=c('South West','North East','Merseyside')), nrow=3, ncol=1, scales='free_y') + theme(legend.position='bottom')

# Again we use the code from https://github.com/tidyverse/ggplot2/issues/2096 to colour the facet titles according to the area:
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-1', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- c('deepskyblue','yellow','red')
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

# Plot histograms of the marginal posterior distributions for beta:
par(mfrow=c(3,1))
hist(abc.sw$beta, col='deepskyblue', main='South West (Tier 1)', xlab='Beta')
hist(abc.ne$beta, col='yellow', main='North East (Tier 2)', xlab='Beta')
hist(abc.ms$beta, col='red', main='Merseyside (Tier 3)', xlab='Beta')

# Compare directly the best parameter estimates for each area:
best.sw <- abc.sw[which.min(abc.sw$diff),-4]|>
  mutate(area='sw')
best.ne <- abc.ne[which.min(abc.ne$diff),-4]|>
  mutate(area='ne')
best.ms <- abc.ms[which.min(abc.ms$diff),-4]|>
  mutate(area='ms')
best <- rbind(best.sw, best.ne, best.ms)
View(best)

# What parameter values would we get if we considered the entire population of England together?
engdata <- data|>
  group_by(date)|>
  summarise(area = 'England', population = sum(population), new_cases = sum(new_cases), cumulative_cases = sum(cumulative_cases))

abc.eng <- abc(10000, 10000, c(0.1,0.2), c(100000,500000), c(300000,600000), c('14/10/2020','5/11/2020'), engdata)
best.eng <- abc.eng[which.min(abc.eng$diff),-4]|>
  mutate(area='all')

# Based on initial conditions from ABC of whole of England and beta from individual tiers, what would the data look like?

# Create simulated data under these conditions:
sim.eng.sw <- sirs(best.sw[1,1], 1/10, 1/152, '14/10/2020', '6/1/2021', best.eng[1,2], best.eng[1,3], engdata)|>
  mutate(area='sw')
sim.eng.ne <- sirs(best.ne[1,1], 1/10, 1/152, '14/10/2020', '6/1/2021', best.eng[1,2], best.eng[1,3], engdata)|>
  mutate(area='ne')
sim.eng.ms <- sirs(best.ms[1,1], 1/10, 1/152, '14/10/2020', '6/1/2021', best.eng[1,2], best.eng[1,3], engdata)|>
  mutate(area='ms')
sim.eng.all <- sirs(best.eng[1,1], 1/10, 1/152, '14/10/2020', '6/1/2021', best.eng[1,2], best.eng[1,3], engdata)|>
  mutate(area='all')

simulations <- rbind(sim.eng.sw, sim.eng.ne, sim.eng.ms, sim.eng.all)

# Plot
sims <- ggplot(simulations, aes(x=date, y=new_cases, colour=area)) + geom_line(size=1.1) + xlab('Date') + ylab('New Daily Cases') + scale_color_discrete(name='Area on which beta is based', breaks=c('sw','ne','ms','all'), labels=c('South West','North East','Merseyside','Whole of England'), type=c('black','red','yellow','deepskyblue')) #+theme(legend.position = 'none') # Legend not included on poster for space
sims <- sims + geom_hline(yintercept=max(sim.eng.sw$new_cases), linetype='dashed', color='deepskyblue', size=2)
sims <- sims + geom_hline(yintercept=max(sim.eng.ne$new_cases), linetype='dashed', color='yellow', size=2)
sims <- sims + geom_hline(yintercept=max(sim.eng.ms$new_cases), linetype='dashed', color='red', size=2)
sims <- sims + geom_hline(yintercept=max(sim.eng.all$new_cases),linetype='dashed',color='black',size=2)
sims
