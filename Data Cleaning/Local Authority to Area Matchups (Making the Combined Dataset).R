# install.packages("tidyverse")

# Load required package (use above code to install package if not done before)
library(tidyverse)

# Load our tier data - this gives us the Local Authorities whose general area we wish to know.
data <- read.csv("https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/1st%20system%20tier%20changes%20%5B1%5D%5B2%5D%5B3%5D.csv")

# Load area data:
counties <- read.csv("https://github.com/mollywilshaw/Final-Year-Project-Poster/raw/main/Data/Local_Authority_District_to_County_(April_2021)_Lookup_in_England%20%5B5%5D.csv")

# Join these data sets and select the relevant columns
lawcounties <- left_join(data, counties, by=c("la_name"="LAD21NM"))
lawcounties <- lawcounties|>
  select(c(1,12))

# Note that there are several missing values here. We fix these manually as below:

lawcounties[lawcounties$la_name=="Hartlepool"|lawcounties$la_name=="Middlesbrough"|lawcounties$la_name=="Darlington"|lawcounties$la_name=="Redcar and Cleveland"|lawcounties$la_name=="Stockton-on-Tees",2] <- "Tees Valley"
lawcounties[lawcounties$la_name=="York",2] <- "North Yorkshire"
lawcounties[lawcounties$la_name=="Halton"|lawcounties$la_name=="Warrington"|lawcounties$la_name=="Cheshire East"|lawcounties$la_name=="Cheshire West and Chester",2] <- "Cheshire"
lawcounties[lawcounties$la_name=="Blackpool"|lawcounties$la_name=="Blackburn with Darwen",2] <- "Lancashire"
lawcounties[lawcounties$la_name=="East Riding of Yorkshire"|lawcounties$la_name=="Kingston upon Hull, City of",2] <- "East Riding of Yorkshire"
lawcounties[lawcounties$la_name=="North Lincolnshire"|lawcounties$la_name=="North East Lincolnshire",2] <-"Lincolnshire"
lawcounties[lawcounties$la_name=="Derby"|lawcounties$la_name=="Barms"|lawcounties$la_name=="Blackbrook"|lawcounties$la_name=="Burbage"|lawcounties$la_name=="Buxton Central"|lawcounties$la_name=="Chapel East"|lawcounties$la_name=="Chapel West"|lawcounties$la_name=="Corbar"|lawcounties$la_name=="Cote Heath"|lawcounties$la_name=="Dinting"|lawcounties$la_name=="Gamesley"|lawcounties$la_name=="Hadfield North"|lawcounties$la_name=="Hadfield South"|lawcounties$la_name=="Hayfield"|lawcounties$la_name=="Hope Valley"|lawcounties$la_name=="Howard Town"|lawcounties$la_name=="Limestone Peak"|lawcounties$la_name=="New Mills East"|lawcounties$la_name=="New Mills West"|lawcounties$la_name=="Old Glossop"|lawcounties$la_name=="Padfield"|lawcounties$la_name=="St John's"|lawcounties$la_name=="Sett"|lawcounties$la_name=="Simmondley"|lawcounties$la_name=="Stone Bench"|lawcounties$la_name=="Temple"|lawcounties$la_name=="Tintwistle"|lawcounties$la_name=="Whaley Bridge"|lawcounties$la_name=="Whitfield",2] <- "Derbyshire"
lawcounties[lawcounties$la_name=="Leicester",2] <- "Leicestershire"
lawcounties[lawcounties$la_name=="Rutland",2] <- "Rutland"
lawcounties[lawcounties$la_name=="Nottingham",2] <- "Nottinghamshire"
lawcounties[lawcounties$la_name=="Herefordshire, County of",2] <- "Herefordshire"
lawcounties[lawcounties$la_name=="Shropshire"|lawcounties$la_name=="Telford and Wrekin",2] <- "Shropshire"
lawcounties[lawcounties$la_name=="Stoke-on-Trent",2] <- "Staffordshire"
lawcounties[lawcounties$la_name=="Bristol, City of"|lawcounties$la_name=="South Gloucestershire"|lawcounties$la_name=="Bath and North East Somerset",2] <- "West of England"
lawcounties[lawcounties$la_name=="North Somerset",2] <- "Somerset"
lawcounties[lawcounties$la_name=="Plymouth"|lawcounties$la_name=="Torbay",2] <- "Devon"
lawcounties[lawcounties$la_name=="Swindon"|lawcounties$la_name=="Wiltshire",2] <- "Wiltshire"
lawcounties[lawcounties$la_name=="Peterborough",2] <- "Cambridgeshire"
lawcounties[lawcounties$la_name=="Bedford"|lawcounties$la_name=="Central Bedfordshire"|lawcounties$la_name=="Luton",2] <- "Bedfordshire"
lawcounties[lawcounties$la_name=="Southend-on-Sea"|lawcounties$la_name=="Thurrock",2] <- "Essex"
lawcounties[lawcounties$la_name=="Medway",2] <- "Kent"
lawcounties[lawcounties$la_name=="Slough"|lawcounties$la_name=="Bracknell Forest"|lawcounties$la_name=="West Berkshire"|lawcounties$la_name=="Reading"|lawcounties$la_name=="Windsor and Maidenhead"|lawcounties$la_name=="Wokingham",2] <- "Berkshire"
lawcounties[lawcounties$la_name=="Milton Keynes"|lawcounties$la_name=="Buckinghamshire",2] <- "Buckinghamshire"
lawcounties[lawcounties$la_name=="Brighton and Hove",2] <- "East Sussex"
lawcounties[lawcounties$la_name=="Portsmouth"|lawcounties$la_name=="Southampton",2] <- "Hampshire"
lawcounties[lawcounties$la_name=="Isle of Wight",2] <- "Isle of Wight"
lawcounties[lawcounties$la_name=="County Durham",2] <- "County Durham"
lawcounties[lawcounties$la_name=="Cornwall and Isles of Scilly",2] <- "Cornwall"
lawcounties[lawcounties$la_name=="Northumberland",2] <- "Northumberland"
lawcounties[lawcounties$la_name=="Dorset"|lawcounties$la_name=="Bournemouth, Christchurch and Poole",2] <- "Dorset"
lawcounties[lawcounties$la_name=="Corby"|lawcounties$la_name=="Kettering"|lawcounties$la_name=="Wellingborough"|lawcounties$la_name=="North Northamptonshire"|lawcounties$la_name=="East Northamptonshire",2] <- "North Northamptonshire"
lawcounties[lawcounties$la_name=="Northampton"|lawcounties$la_name=="Daventry"|lawcounties$la_name=="South Northamptonshire"|lawcounties$la_name=="West Northamptonshire",2] <- "West Northamptonshire"
lawcounties[lawcounties$la_name=="Hackney and City of London",2] <- "Inner London"

# Save the new data set
#write.csv(lawcounties,"la_with_areas.csv", row.names = FALSE)
