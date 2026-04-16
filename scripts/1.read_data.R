rm(list = ls())

set.seed(123)

library(readxl)

df <- read_excel("data (raw).xlsx")

names(df) <- make.names(c("Date",                                                                                                                 
                          "Excess.mortality",                                                                                                     
                          "Schools",                                                                                                              
                          "Face.masks",                                                                                                           
                          "Lockdown.severity",                                                                                                    
                          "Majority.COVID.19.variant",
                          "Flights.7.day.moving.average",                                                                                       
                          "OpenTable.restaurant.bookings.London.index",                                                                         
                          "Google.homeworking.Greater.London.mobility.index",                                                                   
                          "Google.workplace.Greater.London.mobility.index",                                                                     
                          "Apple.walking.London.mobility.index",                                                                                
                          "Google.parks.Greater.London.mobility.index",                                                                         
                          "Google.retail.recreation.Greater.London.mobility.index",                                                           
                          "Google.grocery.pharmacy.Greater.London.mobility.index",                                                            
                          "Google.transit.stations.mobility.index",                                                                               
                          "TfL.Tube.mobility.index",                                                                                              
                          "TfL.Bus.mobility.index",                                                                                               
                          "Citymapper.journeys.mobility.index",                                                                                   
                          "Season",                                                                                                               
                          "PCR.tests",                                                                                                            
                          "PCR.tests.capacity",                                                                                                   
                          "Antibody.tests",                                                                                                       
                          "Antibody.tests.capacity",                                                                                              
                          "Pillar.1.NHS.and.UKHSA.capacity",                                                                                    
                          "Pillar.2.UK.Government.capacity",                                                                                    
                          "Pillar.3.Antibody.capacity",                                                                                         
                          "Pillar.4.Surveillance.capacity",                                                                                     
                          "Pillar.1.NHS.and.UKHSA.tests",                                                                                       
                          "Pillar.2.UK.Government.tests",                                                                                       
                          "Pillar.3.Antibody.tests",                                                                                            
                          "Pillar.4.Surveillance.tests",                                                                                        
                          "Tests.across.all.4.Pillars",                                                                                           
                          "New.cases",                                                                                            
                          "New.infections",                                                                                                       
                          "Reinfections",                                                                                                         
                          "New.cases.specimen.date.7.day.change.",                                                                               
                          "New.cases.specimen.date.7.day...change.",                                                                             
                          "Hospital.admissions",                                                                                                  
                          "Patients.in.hospital",                                                                                                 
                          "Patients.in.MVBs",                                                                                            
                          "Vaccinations.total",                                                                                                 
                          "Vaccinations.1st.dose",                                                                                              
                          "Vaccinations.2nd.dose",                                                                                              
                          "Vaccinations.3rd.dose",                                                                                              
                          "First.dose.uptake",                                                                                                      
                          "Second.dose.uptake",                                                                                                      
                          "Third.dose.uptake",                                                                                                      
                          "Deaths.within.28.days.of.COVID.19",                                                                                    
                          "COVID-19.deaths.on.certificate"))

df <- df[ , !(names(df) %in% c("New.cases.specimen.date.7.day.change.",                                                                               
                               "New.cases.specimen.date.7.day...change.", 
                               "Deaths.within.28.days.of.COVID.19"))]

df <- df[order(df$Date, decreasing = TRUE),]
char_cols <- sapply(df, is.character)
df[char_cols] <- lapply(df[char_cols], factor)

df <- droplevels(df)

raw <- as.character(df$Apple.walking.London.mobility.index)
raw[ raw == "NA" ] <- NA
num <- as.numeric(raw)
cat("NAs after coercion:", sum(is.na(num)), "\n")
df$Apple.walking.London.mobility.index <- num

date_column <- df$Date
df$Date <- NULL

