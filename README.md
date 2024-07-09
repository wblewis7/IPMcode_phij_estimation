# IPMcode_phij_estimation
Data and R/NIMBLE code for running an integrated population model inferring juvenile survival from breeding productivity and post-breeding population abundance and age/sex ratio data. William B. Lewis, Chloe R. Nater, Justin A. Rectenwald, D. Clay Sisson, and James A. Martin.
2024. Use of breeding and post-breeding data to estimate juvenile survival with integrated population models. International Statistical Ecology Conference.
Contact: william.lewis27@uga.edu


---

# Metadata

# S.GA.NOBO.IPM.data.simple.gzip
Data for running the integrated population model are stored in this gzip file. Demographic data were collected from a population of northern bobwhites in southern Georgia, USA, from 1998 - 2022. There are four main data sources used in the model: 
1) radiotelemetry data for estimating survival
2) nest monitoring data for estimating breeding productivity
3) post-breeding population surveys
4) age/sex ratio of post-breeding harvested birds.
   
## nyears
The number of years of data (25).
## years
Actual years of data.
## n.b.month
Number of months during the breeding season (6, April - September).
## n.breeding.months
Number of months during which juveniles enter population (4, June - September).
## Harv.sexage
The number of adult females (AF), first-year females (JF), and males (M) harvested in November in each year of the study. Birds were aged and sexed based on plumage characteristics.
## Harv.N
The total number of harvested birds in November of each year. Equal to row sums of Harv.sexage.
## chicks.nest.mean and chicks.nest.sd
Nests were monitored each month of the breeding season (June - September) to determine the nest success and number of chicks hatched. A small number of nests were labeled as successful but did not have data recorded on the number of chicks produced. We estimated the number
of chicks hatched from these nests based on the clutch size and the average hatch rate for a given month/year. The mean and standard deviation of these posterior distribution of total number of chicks produced are given in chicks.nest.mean and chicks.nest.sd, respectively. 
Rows correspond to years while columns correspond to breeding months.
## N.tracked
The number of adult females from which the productivity data (chicks.nest.mean and chicks.nest.sd) was collected. Data is structured by year (rows) and breeding month June - September (columns).
## count
Total number of calling bobwhite coveys (aggregations of bobwhites in the fall) detected in each year of the study. Covey counts were performed using a 4-person quadrat sampling method.
## effort
Covey count surveys occurred at a variable portion of the total study area in each year.
## csize.mean, csize.sd, csize.min, csize.max
Covey counts were related to population abundance through estimating average covey size. In each year, a subset of detected coveys were flushed to count size. The mean size (csize.mean) and sd (csize.sd) of flushed coveys is given for each year, while csize.min and 
csize.max give the 2.5% and 97.5% quantiles of observed covey sizes across the entire study duration.
## avail.mean, avail.sd, avail.min, avail.max, pdet_mean, pdet_sd, pdet_min
The detection process was decomposed into calling availability (probability that covey was within survey area and calling during survey) and conditional detection (probability of at least one of four observers detecting a covey conditional on it being available). No 
covariates were collected to model these processes so we used informative priors. Minimum, maximum, and yearly mean and sd of availability (logit scale) are given by avail.min, avail.max, avail.mean, and avail.sd, respectively. Values were generated based on observed
neighbor density in each year (covey count - 1) and coefficients in Wellendorf, S. D., W. E. Palmer, and P. T. Bromley (2004). Estimating calling rates of northern bobwhite coveys and measuring abundance. Mean, sd, and min values of conditional detection probability (logit
scale) are given by pdet_mean, pdet_sd, and pdet_min, respectively. Values were estimated based on data from Howell, P. E., T. M. Terhune, and J. A. Martin (2021). Edge density affects demography of an exploited grassland bird.
## CH.state
Birds were captured using baited traps each year in the spring and fall and fitted with radiotransmitters. Tagged birds were tracked every few days until mortality or tag failure. CH.state gives the formatted capture history data from radiotelemetry. Rows correspond to
birds while columns correspond to days of the study. Values of 1 respresent the bird was located alive while values of 0 represent the transmitter was recovered and assumed dead.
## CH.age
Capture history of CH.state structured by age and season. The breeding season extended from April - September while the non-breeding season extended from October - March. Subadults (young birds) transition to adults at the start of April. Values of 0 represent transmitter
recovered (assumed dead), 1 represent detected alive nonbreeding subadult, 2 represent detected alive nonbreeding adult, and 3 represent alive breeding adult. 
## CH.year
Vector giving the year of the study for every tracking occasion (columns of CH.state and CH.age).
## nCH
Number of birds with telemetry data (number of rows in CH.state and CH.age).
## trackdates
Matrix giving, for each bird (rows), the columns from CH.state which have active tracking data. Excludes the first tracking column from CH.state, as this is the release date. A small number birds were tracked until radio failure, caught in a subsequent year, and were 
redployed with another transmitter.
## trackperiods
The number of active tracking periods for each bird. Corresponds to the length of non-NA values in each row of trackdates.


<br />
<br />

# S.GA.NOBO.IPM.Simple.Code.R
R code for running the integrated population model using demographic data from northern bobwhites (Colinus virginianus). The model is female-specific and structured by age (juvenile, subadult, or adult), season (breeding (June - September) or non-breeding (October -
March)), and year. Adult survived each month based on the breeding-season survival rate. Juveniles enter the population in each month June - September based on the number of surviving adults and the monthly per-capita productivity rates. Juveniles survive until the 
non-breeding season based on the juvenile survival rate, after which they transition to subadults. Subadults and adults survive through the non-breeding season based on the age-specific non-breeding survival rate, with subadults transitioning to adults at the start of
the breeding season in April. Vital rates vary yearly based on a Normal distribution with a global mean and standard deviation.
Survival of subadults and adults were assessed throughout the year via radiotelemetry. Monthly per-capita productivity rates were assessed each month June - September from nest monitoring data. The model incorporates uncertainty in the total number of chicks produced
from monitorred nests from nest with missing data. Post-breeding data was collected in November. Covey counts were related to total abundance through accounting for survey effort, average covey size, calling availability, and conditional detection. Harvest data was used
to assess age and sex ratios in November. This model is used to infer juvenile survival in the absence of direct data. Productivity data was used to estimate the number of chicks produced in the breeding season while the covey counts and age/sex ratio were used to estimate the number of young birds surviving in November, with the difference reflective of juvenile survival rates.
