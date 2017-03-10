#Create function that extracts daily mean,min,max, DTR from 
first2 <- function (df) {
 #Convert to C
  df[,2] <- (df[,2] - 32)/1.8
 #Calculate daily values
  Temp.av <- mean(df[,2])
  Temp.min <- min(df[,2])
  Temp.max <- max(df[,2])
  Temp.DTR <- max(df[,2]) - min(df[,2])
  RH.av <- mean(df[,1])
  RH.min <- min(df[,1])
  RH.max <- max(df[,1])
  RH.DTR <- max(df[,1]) - min(df[,1])
  data <- c(i, Temp.av, Temp.min, Temp.max, Temp.DTR,RH.av,
                     RH.min,RH.max,RH.DTR)
  return(data)  
} 

#create vector of file name numbers
vec <- c(150615:150630,150701:150731,150801:150831,150901:150930,
         151001:151013)

#Create empty matrix
weather <- matrix(nrow=length(vec), ncol=9)
#Run a for loop for the file names of the above that puts the data into 
#the weather matrix
for (i in 1:length(vec)){
  file <- paste("http://weather.ggy.uga.edu/data/csv/",vec[i],".csv", sep="")
  df <- read.csv(file, header=F, colClasses = c(rep("NULL", 3),rep("numeric", 2), 
                                                rep("NULL", 7)), skip=6)
  
  weather[i,] <- first2(df)
}
weather<-as.data.frame(weather)

#Name columns
colnames(weather) <- c("Date", "Temp.av", "Temp.min", "Temp.max", "Temp.DTR"
                        ,"RH.av", "RH.min","RH.max","RH.DTR")
#Change Date to actual dates
weather$Date <- as.Date(as.character(vec),"%y%m%d" )
#save as CSV
write.csv(weather, file="UGAStation.csv")

