# install.packages('ggplot2')
# install.packages('tidyverse')
# install.packages('rgeos')
# install.packages('rgdal')
# install.packages('maptools')

# Load required packages (use above code to install packages if not done before)
library(ggplot2)
library(tidyverse)
library(rgeos)
library(rgdal)
library(maptools)

# Load the required data
data <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data%20Cleaning/Combined%20Data.csv') # Large combined data set - contains council codes as well as names, which we will need to join with the boundary data later
changes <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/1st%20system%20tier%20changes%20%5B1%5D%5B2%5D%5B3%5D.csv') # Data set containing the tier classification of each council on the dates when the regulations were changed

# Extract the column from the data containing all council names
la <- changes$la_name

# Here we re-use the sorting algorithm used for larger areas while making the combined data set, but for individual Local Authorities.

# Create vectors to store the sorted council names:
tier1 <- c() 
tier2 <- c() 
tier3 <- c() 
change12 <- c() 
change13 <- c() 
change23 <- c() 
change123 <- c() 

for(x in la){ # Loop over each council
  perarea <- filter(changes, changes$la_name==x) # Select row corresponding to considered council
  t <- unique(as.numeric(perarea[1,2:8])) # List unique values in that row, i.e. all tiers which this council area was in
  # If there is only one unique value, the council stayed in 1 tier over the whole period - we sort based on which tier this is
  if(length(t)==1){
    if(t==1){
      tier1 <- append(tier1,x)
    }
    else if(t==2){
      tier2 <- append(tier2,x)
    }
    else if(t==3){
      tier3 <- append(tier3,x)
    } 
  }
  # Similarly if there are 2 unique values, the council could have moved from tier 1->2, 2->3 or 1->3 (note in this system no council area moved down a tier)
  if(length(t)==2) {
    if( (t[1]==1 & t[2]==2)| 
        (t[1]==2 & t[2]==1)){
      change12 <- append(change12,x)
    }
    else if((t[1]==1 & t[2]==3)|
            (t[1]==3 & t[2]==1)){
      change13 <- append(change13,x)
    }
    else if((t[1]==2 & t[2]==3)|
            (t[1]==3 & t[2]==2)){
      change23 <- append(change23,x)
    }
  }
  # If there are 3 unique values, the council must have moved from tier 1->2->3
  if(length(t)==3){
    change123 <- append(change123,x)
  }
}
# We now have several lists of council names, sorted by their tier status over the first tiered lockdown period

# Add a variable for this status to the data set:
# Note that we subset the combined data set, as we only need the council names and codes.
data1 <- unique(data[,c(1,2,3)])|>
  mutate(status=ifelse(la_name %in% tier1,'tier 1',
                       ifelse(la_name %in% tier2,'tier 2',
                       ifelse(la_name %in% tier3,'tier 3',
                       ifelse(la_name %in% change12,'changed 1 to 2',
                       ifelse(la_name %in% change23,'changed 2 to 3',NA))))))

# Load boundary data:
bdys <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/LAD_DEC_2020_UK_BUC%20%5B8%5D.csv')

# Join the boundary data set with our sorting data and prepare the new data set for plotting:
datawmap <- left_join(data1, bdys, by=c('la_code' ='id'))
datawmap <- arrange(datawmap)

# Plot a map of the councils in England, filled based on their status:
ggplot(data=datawmap, aes(x=long, y=lat, group=group, fill=factor(status))) + geom_polygon() + coord_equal() + theme_void() + scale_fill_discrete(name='Status over first system', breaks=c('tier 1', 'tier 2', 'tier 3', 'changed 1 to 2', 'changed 2 to 3'), labels=c('Medium', 'High', 'Very High', 'Changed from Medium to High', 'Changed from High to Very High'), type=c('chartreuse', 'darkorange', 'deepskyblue', 'yellow', 'red')) + theme(legend.position='bottom', legend.direction='vertical')
