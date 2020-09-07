appsPackages = c(
  "rstudioapi",
  "DT",
  'gnm',
  'tidyverse',
  'demography',
  'plot3D',
  'nloptr',
  'rgl',
  'car',
  'plotly',
  'scales',
  'pracma',
  'data.table',
  'dfoptim',
  'forecast',
  'StMoMo',
  'rmarkdown',
  'webshot'
)

for (package in appsPackages) {
  if (!require(package, character.only = T, quietly = T)) {
    if (Sys.getenv(c("R_PLATFORM")) == "x86_64-redhat-linux-gnu") {
      install.packages(package, contriburl = "file:///apps/rpp/")
    } else{
      install.packages(package)
    }
  }
  library(package, character.only = T)
}




#########  Data Preperation ##############

  current_path <- getActiveDocumentContext()$path
  setwd(dirname(current_path))
  
  #read the text files
  
  DEATH <- 'C:\\Users\\Silas.Amarasinghe\\Desktop\\mpt_repo\\www\\data\\US Death Rates (1933-2016)2.txt'
    # read.table(InputDeath$datapath,header = input$csvHeader1, skip = input$csvSkip1 - 1)
    # fread('D:\\Folder\\Filename.csv', stringsAsFactors = F, sep=",", na.strings="NA")
  EXPO <- 'C:\\Users\\Silas.Amarasinghe\\Desktop\\mpt_repo\\www\\data\\US Pop Exposure (1933-2016)2.txt'
    #InputExposure$datapath
    # read.table(InputExposure$datapath,header = input$csvHeader2, skip = input$csvSkip2 - 1)
    # fread('D:\\Folder\\Filename.csv', stringsAsFactors = F, sep=",", na.strings="NA")
  
  #choosing Gender
  Gender <- "Male" #choices = c("Male", "Female", "Total") 
  
  mortality_data <- read.demogdata(DEATH, EXPO, type = "mortality", label = "Mx")
  
     
      M_data <- as.data.frame(mortality_data[["rate"]][[tolower(isolate(Gender))]])
      M_data$x <- rownames(M_data)
      last_column <- ncol(M_data)
      M_data <- M_data[, c(last_column , (1:(last_column - 1)))]
      M_data <- M_data %>% gather(y, "M", 2:last_column)
      
      # Get E (Exposure) data from demogdata based on gender chosen
      
      ex_data <- as.data.frame(mortality_data[["pop"]][[tolower(isolate(Gender))]])
      ex_data$x <- rownames(ex_data)
      last_column <- ncol(ex_data)
      ex_data <- ex_data[, c(last_column , (1:(last_column - 1)))]
      ex_data <- ex_data %>% gather(y, "E", 2:last_column)
      
      # Join M and Ex
      
      data <- left_join(M_data, ex_data, by = c("x", "y"))
      
      # Compute D (No of Deaths) - colnames = x, y, M, E, D
      
      data$D <- data$M * data$E
      
      # Rename M to mx for version purpose
      
      colnames(data)[3] <- "mx"
      
      df <- data
      
      #Remove "110+" Age band (We should remember to review deaths in the age band - not done)
      df = df[which(df$x != "110+"), ]
      #Convert ages and years to numeric
      df$x <- as.numeric(as.character(df$x))
      df$y <- as.numeric(as.character(df$y))
      df <- df[order(df$x), ]
      df <- as.data.table(df)
      
      # adding year of birth
      df[, yob := as.numeric(y - x)]
      
      #reorder column sequence
      setcolorder(df, c("x", "y", "yob", "E", "D", "mx"))
      
      # CDR is for Crude Death rates (Death over Exposure) & log CDR for the analysis
      # df[,mx := D/E]
      df[, logmx := log(mx)]
      
      # adding the q(x,t)
      df[, qx := 1 - exp(-mx)]
      
      # adding the mortality improvement
      df[, MI := 1 - qx / c(0, qx[-.N]) , by = x]
      # Number of years for rolling average
      smoothMI <- 9
      
      df[, smoothed_MI := movavg(MI, smoothMI, "s")]
      
      # # adding the mean of mortality improvement by ages
      # df[, meanMI := mean(MI, na.rm = TRUE),by = x]
      
      # adding the Mortality Reduction Factor
      df[, Rx := mx / c(0, mx[-.N]) , by = x]
      
      # adding Death times logarithmic of Exosure
      df[, DlogE := D * log(E)]
      
      # adding logarithmic of Death factorial
      df[, logDfac := 0.5 * log(2 * pi * D) + D * log(D) - D]
      
      df <- df[order(df$x), ]
      
      # write.csv(df, "df.csv", row.names = FALSE)
      # Select Ages and Years (sliderInput)
      Ages <- c(30, 100)
      Years <- c(1933, 2016)
      
      df <- df %>% filter(x %in% seq(isolate(Ages[1]), isolate(Ages[2]), 1)) %>% filter(y %in% seq(isolate(Years[1]), isolate(Years[2]), 1))
      
      df <- as.data.frame(df)
    
    # Data Summary - DataBase 
      df <- as.data.table(df) 
      df <- df[, c(1:11)]
      names(df)[names(df) == "x"] = "Age"
      names(df)[names(df) == "y"] = "Calendar Year"
      names(df)[names(df) == "yob"] = "Year of Birth"
      names(df)[names(df) == "D"] = "Deaths"
      names(df)[names(df) == "E"] = "Exposure"
      names(df)[names(df) == "mx"] = "Death Rates"
      names(df)[names(df) == "logmx"] = "Log Death Rates"
      names(df)[names(df) == "qx"] = "Mortality Rate"
      names(df)[names(df) == "MI"] = "MI Rate"
      names(df)[names(df) == "smoothed_MI"] = "Smoothed MI Rate"
      names(df)[names(df) == "Rx"] = "Mortality Reduction Factor"
    

    # M1 - Fitting Lee Carter Model using Maximum Likelihood Estimation
    # Maximum likelihood estimation function with StMoMo
      
      data <- StMoMoData(mortality_data, series = tolower(Gender))
      ### M1
      model <-
        suppressWarnings(
          fit(
            lc(link = "log"),
            data = data,
            verbose = FALSE,
            options(warn= -1)
          )
        )
      M1_mle <- model
      AIC(M1_mle)
      plot(M1_mle)
      
    

      
      wxt = genWeightMat(
        ages = c(30,35,40,45,50),
        years =  seq(Years[1], Years[2], 1),
        clip = 5
      )
      
      eta1 <- 1 
      
      ###M2
      
      model <-
        suppressWarnings(
          fit(
            rh(
              link = "log",
              cohortAgeFun = if (eta1 == "Estimated") {"NP"} else {"1"},
              approxConst = F
            ),
            wxt = wxt,
            data = data,
            start.ax = M1_mle$ax,
            start.bx = M1_mle$bx,
            start.kt = M1_mle$kt,
            verbose = FALSE,
            options(warn= -1)
          )
        )
    
      M2_mle <- model
      AIC(M2_mle)
      plot(M2_mle)
      
      ###M3
      
      model <-
        suppressWarnings(
          fit(
            apc(link = "log"),
            data = data,
            wxt=wxt,
            start.ax = M1_mle$ax,
            start.kt = M1_mle$kt,
            verbose = FALSE,
            options(warn= -1) 
          )
        )
      
      M3_mle <- model
      AIC(M3_mle)
      plot(M3_mle)
      
      
      ###M5
      model <-
        suppressWarnings(
          fit(
            cbd(link = "log"),
            data = data,
            verbose = FALSE,
            options(warn= -1)
          )
        )
      M5_mle <- model
      plot(M5_mle)
      AIC(M5_mle)    
      
      ###M6
      model <-
        suppressWarnings(
          fit(
            m6(link = "log"),
            data = data,
            verbose = FALSE,
            options(warn= -1)
          )
        )
      M6_mle <- model
      
      AIC(M6_mle)
      
      
      model <-
        suppressWarnings(
          fit(
            m7(link = "log"),
            data = data,
            verbose = FALSE,
            options(warn= -1)
          )
        )
      M7_mle <- model
      AIC(M7_mle)
     
      
      model <-
        suppressWarnings(
          fit(
            m8(link = "log", xc = Ages[2]),
            data = data,
            ages.fit = seq(Ages[1], Ages[2], 1),
            years.fit =  seq(Years[1], Years[2], 1),
            verbose = FALSE,
            options(warn= -1)
          )
        )
      M8_mle <- model
      AIC(M8_mle)
   
    
     
      
      
    