# install.packages('janitor')
# install.packages('tidyverse')

# Load required packages (use above code to install packages if not done before)
library(janitor)
library(tidyverse)

# First we upload and clean the data sets we will be combining:

# Data on tier allocations of Local Authorities during the first tier system
tiers1 <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/1st%20system%20tier%20data%20%5B1%5D%5B2%5D%5B3%5D.csv')
tiers1 <- tiers1[1:312,]
tiers1 <- pivot_longer(tiers1, -1, names_to='date', values_to='tier') 
tiers1$date <- as.Date(tiers1$date, format='X%d.%m.%Y')

# Data on tier allocations of Local Authorities during the second tier system
tiers2 <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/2nd%20sytem%20tier%20data%20%5B4%5D.csv')
tiers2 <- tiers2[1:312,1:37]
tiers2 <- pivot_longer(tiers2, -1, names_to='date', values_to='tier')
tiers2$date <- as.Date(tiers2$date, format='X%d.%m.%Y')

# Here we add the period between these two systems (when there was a general national lockdown):
dates <- data.frame(date = seq(as.Date('5/11/2020','%d/%m/%Y'), as.Date('1/12/2020','%d/%m/%Y'), by='days'))
names <- data.frame(la_name = unique(tiers1$la_name))
tiers <- data.frame(tier='National Lockdown')
natlock <- full_join(dates, names, by=character())
natlock <- full_join(natlock, tiers, by=character())

# Here we add 'tier' values for data before/after our tier list data - this is to allow better analysis later (due to lag in infection data)

# Note - the early days of the case data do not record cases for all councils - we want to find the first date when every Local Authority's cases were recorded to start our data set from, which we do as follows:
# Download case data
cases <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/ltla_2023-01-12%20%5B6%5D.csv')|>
  filter(grepl('E', areaCode) == TRUE) # Filter out data for Scotland/Wales/NI
cases$date <- as.Date(cases$date, format='%Y-%m-%d') # Make sure that R recognises the 'date' variable of the cases data as a date

names <- unique(cases$areaName) # Create a vector of all area names

cases1 <- cases[ ,c(2,4)] # Select only date and area names - these are the only variables we are interested in for this purpose
dates <- seq(as.Date('30/1/2020','%d/%m/%y'), as.Date('13/10/2021','%d/%m/%Y'), by='days') # Create a vector of all dates before the first tier system was implemented
for(i in 1:length(dates)){
  areas <- unique(cases1[cases1$date==dates[i],1]) # List the areas for which we have data on each date
  if(length(areas)==length(names)) return(print(dates[i])) # Return the first value for which we have data for all of the areas
}
# This date turns out to be 2020-03-31, so this is where we will start our data set.

names <- data.frame(la_name = unique(tiers1$la_name))
predates <- data.frame(date=seq(as.Date('31/3/2020','%d/%m/%Y'), as.Date('13/10/2020','%d/%m/%Y'), by='days'))
postdates <- data.frame(date=seq(as.Date('7/1/2021','%d/%m/%Y'), as.Date('6/2/2021','%d/%m/%Y'), by='days'))
pretiers <- data.frame(tier='Before tier system')
posttiers <- data.frame(tier=4)
pre <- full_join(names, predates, by=character())
pre <- full_join(pre, pretiers, by=character())
post <- full_join(names, postdates, by=character())
post <- full_join(post, posttiers, by=character())
pp <- rbind(pre, post)

# Combine all of the tier data into one data set:
tiersall <- rbind(tiers1, tiers2, natlock, pp)

# We now filter the COVID case data for England over the time period we want to consider:
cases <- cases[ ,-c(3,6)]|>
  filter(between(date, as.Date('31/3/2020','%d/%m/%Y'), as.Date('6/2/2021','%d/%m/%Y')) == TRUE)

# There are discrepancies between data sets for areas in Buckinghamshire - we fix this issue as follows:
cases2 <- cases|>
  filter(areaName=='South Bucks'|areaName=='Wycombe'|areaName=='Chiltern'|areaName=='Aylesbury Vale')|>
  select(c(3,4,5))|>
  group_by(date)|>
  summarise('cumCasesBySpecimenDate' = sum(cumCasesBySpecimenDate), 'newCasesBySpecimenDate' = sum(newCasesBySpecimenDate))

cases3 <- cbind(data.frame('areaCode'=rep('E06000060',313), 'areaName'=rep('Buckinghamshire',313), cases2))

cases <- rbind(cases, cases3)|>
  filter(areaName!='South Bucks' & areaName!='Wycombe' & areaName!='Chiltern' & areaName!='Aylesbury Vale')

# Download the population data for each area:
population <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/ukpopestimatesmid2020on2021geography%20%5B7%5D.csv') 
population <- population[-c(1:6), c(1:3,34)]|>
  row_to_names(row_number=1)|>
  filter(grepl('E', Code)==TRUE & Geography!='Country' & Geography!='Region')|>
  select(c(1,4))

# Download data on larger areas to which our Local Authorities belong (this data set was created using 'Local Authority to Area Matchups' code)
lawareas <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data%20Cleaning/la_with_areas.csv')

# We now join all of our data sets:
casesntiers <- full_join(tiersall, cases, by=c('date'='date','la_name'='areaName'))
caseswtiersnpop <- left_join(casesntiers, population, by=c('areaCode'='Code'))
caseswtiersnpop <- left_join(caseswtiersnpop, lawareas)

# Create a variable for the cases as a proportion of the population of each Local Authority (in order to compare them more directly):
caseswtiersnpop$`2020 Mid-year estimate` <- as.numeric(gsub(',', '', caseswtiersnpop$`2020 Mid-year estimate`))
caseswtiersnpop <- caseswtiersnpop|>
  mutate(new_cases_proportion = (newCasesBySpecimenDate/`2020 Mid-year estimate`)*10000)|>
  rename('population'=`2020 Mid-year estimate`)|>   # Here we also change some column names to be more intuitive. 
  rename('la_code'=areaCode)|>
  rename('new_cases'=newCasesBySpecimenDate)|>
  rename('cumulative_cases'=cumCasesBySpecimenDate)|>
  rename('area_name'= CTY21NM)

# Rearrange the columns to a 'nicer' order:
caseswtiersnpop <- caseswtiersnpop[ ,c(4,1,8,7,2,6,9,5,3)]

# We now want to add a variable which tells us how the tier status of each larger area changed over the entirety of each tiered lockdown period.
# To do so we use the sorting algorithms below:

# Load tier data for the first system (only on days where legislation changed) and join with larger area data:
tier1changes <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/1st%20system%20tier%20changes%20%5B1%5D%5B2%5D%5B3%5D.csv')
tier1wcounties <- left_join(tier1changes, lawareas)

# Create a list of each unique larger area:
areas <- unique(tier1wcounties$CTY21NM)

# First we create vectors to store the sorted area names:

tier1 <- c() # Areas that stayed in tier 1 for the whole period
tier2 <- c() # Areas that stayed in tier 2 for the whole period
tier3 <- c() # Areas that stayed in tier 3 for the whole period
change12 <- c() # Areas that moved from tier 1 to tier 2
change13 <- c() # Areas that moved from tier 1 to tier 3
change23 <- c() # Areas that moved from tier 2 to tier 3
change123 <- c() # Areas that moved between all tiers
notwholecounty <- c() # Areas in which different Local Authorities moved tiers on different days

# Now we sort the areas:
for(x in areas){ # Loop over all unique areas
  perarea <- filter(tier1wcounties, tier1wcounties$CTY21NM==x) # Filter the data to the area being sorted
  j <- 0 # Set counter
  # The following loop checks each column and increases the counter by 1 if there is only one unique value across all rows for that column:
  for(i in 2:8){
    if(length(unique(perarea[,i]))==1){
      j <- j+1 
    }
    else{
      j <- j
    }
  }
  # 'j' is a counter of how many times the entire area moved tier together - the maximum number j could be is 7
  if(j!=7){
    notwholecounty <- append(notwholecounty,x)
    # j!=7 means that there was at least one day on which only some Local Authorities within this area changed tier, while the others didn't
    # We therefore sort this area into 'notwholecounty' if this is true
  }
  else if(j==7){
    # j==7 means that the entire area moved tiers as a whole - we now find the tier(s) that the area moved to/stayed in
    # Note that each row is identical, so we only need to consider one row (arbitrarily we choose the first)
    t <- unique(as.numeric(perarea[1,2:8])) # We find the unique values in each row
    # If there is only one unique value, the area must have stayed in one tier for the whole period - we sort it into tier1, tier2 or tier3 based on what the unique number is
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
    # If there are 2 unique values, the tier must have changed at some point - we sort the area into tier12 or tier23 based on the values
    if(length(t)==2) {
      if((t[1]==1&t[2]==2)|(t[1]==2&t[2]==1)){
        change12 <- append(change12,x)
      }
      else if((t[1]==1&t[2]==3)|(t[1]==3&t[2]==1)){
        change13 <- append(change13,x)
      }
      else if((t[1]==2&t[2]==3)|(t[1]==3&t[2]==2)){
        change23 <- append(change23,x)
      }
    }
    if(length(t)==3){
      change123 <- append(change123,x)
    }
  }
}
# Note that change13 and change123 both return NULL values

# We can now add a variable to our data set to tell us how the area behaved over the entire first tier system:
caseswtiersnpop <- caseswtiersnpop|>
  mutate(first_system=ifelse(area_name %in% tier1, 'tier 1',
                      ifelse(area_name %in% tier2, 'tier 2',
                      ifelse(area_name %in% tier3, 'tier 3',
                      ifelse(area_name %in% change12, 'changed 1 to 2',
                      ifelse(area_name %in% change23, 'changed 2 to 3',
                      ifelse(area_name %in% notwholecounty, 'different within councils', NA)))))))

# The sorting algorithm for the second tier system proceeds similarly:
tier2changes <- read.csv('https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/2nd%20sytem%20tier%20changes%20%5B4%5D.csv')
tier2changes <- tier2changes[1:312,]
tier2wcounties <- left_join(tier2changes, lawareas)
areas <- unique(tier2wcounties$CTY21NM)

# We now include additional vectors to account for the addition of tier 4:
tier1 <- c()
tier2 <- c()
tier3 <- c()
tier4 <- c()
change12 <- c()
change13 <- c()
change14 <- c()
change23 <- c()
change24 <- c()
change34 <- c()
change123 <- c()
change124 <- c()
change134 <- c()
change234 <- c()
change1234 <- c()
notwholecounty <- c()

# We similarly include additional if() statements to sort areas which involved tier 4.
# Since the entire country went into tier 4 from 6th January, we exclude the final column in our sorting.
for(x in areas){
  perarea <- filter(tier2wcounties, tier2wcounties$CTY21NM==x)
  j <- 0
  for(i in 2:7){
    if(length(unique(perarea[,i]))==1){
      j <- j+1
    }
    else{
      j <- j
    }
  }
  if(j!=6){
    notwholecounty <- append(notwholecounty,x)
  }
  else if(j==6){
    t <- unique(as.numeric(perarea[1,2:7]))
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
      else if(t==4){
        tier4 <- append(tier4,x)
      } 
    }
    if(length(t)==2) {
      if((t[1]==1&t[2]==2)|(t[1]==2&t[2]==1)){
        change12 <- append(change12,x)
      }
      else if((t[1]==1&t[2]==3)|(t[1]==3&t[2]==1)){
        change13 <- append(change13,x)
      }
      else if((t[1]==1&t[2]==4)|(t[1]==4&t[2]==1)){
        change14 <- append(change14,x)
      }
      else if((t[1]==2&t[2]==3)|(t[1]==3&t[2]==2)){
        change23 <- append(change23,x)
      }
      else if((t[1]==2&t[2]==4)|(t[1]==4&t[2]==2)){
        change24 <- append(change24,x)
      }
      else if((t[1]==3&t[2]==4)|(t[1]==4&t[2]==3)){
        change34 <- append(change34,x)
      }
    }
    if(length(t)==3){
      if(1 %in% t & 2 %in% t & 3 %in% t){
        change123 <- append(change123,x)
      }
      else if(1 %in% t & 2 %in% t & 4 %in% t){
        change124 <- append(change124,x)
      }
      else if(1 %in% t & 3 %in% t & 4 %in% t){
        change134 <- append(change134,x)
      }
      else if(2 %in% t & 3 %in% t & 4 %in% t){
        change234 <- append(change234,x)
      }
    }
    if(length(t)==4){
      change1234 <- append(change1234,x)
    }
  }
}

# Note that tier1, tier2, tier4, change12, change13, change14, change124, change 1234 return NULL values

# Again we add a new variable to our data set to show this information:
caseswtiersnpop <- caseswtiersnpop|>
  mutate(second_system=ifelse(area_name %in% tier3, 'tier 3',
                       ifelse(area_name %in% change23, 'changed between tiers 2 and 3',
                       ifelse(area_name %in% change24, 'changed between tiers 2 and 4',
                       ifelse(area_name %in% change34, 'changed between tiers 3 and 4',        
                       ifelse(area_name %in% change123, 'changed between tiers 1, 2 and 3',
                       ifelse(area_name %in% change134, 'changed between tiers 1, 3 and 4',
                       ifelse(area_name %in% change234, 'changed between tiers 2, 3 and 4',
                       ifelse(area_name %in% notwholecounty, 'different within councils', NA)))))))))

# Save combined data set
#write.csv(caseswtiersnpop, "Combined Data.csv", row.names = FALSE)
