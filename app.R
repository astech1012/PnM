#install.packages("shiny")
#install.packages("shinythemes")
library(shiny)
library(shinythemes)

#LINEAR REGRESSION MODELING WITH THE USE OF MATRICES
ols=function(data){
  
  indeps=data[,-1]
  dep=data[,1]
  n=nrow(indeps)
  if (is.null(ncol(indeps))){p=1}
  else{p=ncol(indeps)}
  
  x<-cbind(constant = 1, as.matrix(indeps))
  y<-as.matrix(dep)
  
  #beta estimation
  betas=solve(crossprod(x), crossprod(x,y))
  
  # Computation of standard errors
  s2 <- sum((y - x%*%betas)^2)/(nrow(x) - ncol(x))
  VCV <- as.numeric(s2)*solve(t(x)%*%x)
  SE <- matrix(sqrt(diag(VCV)))
  
  # Computation of t-statistics
  t <- betas/SE
  # Computation of p-values
  np <- 2*pt(abs(t),nrow(x) - ncol(x), lower.tail = FALSE)#PROBABILITY THAT Z IS LESS THAN THE Z-SCORE
  betas<-format(round(betas, 4), nsmall = 4)
  SE<-format(round(SE, 4), nsmall = 4)
  t<-format(round(t, 3), nsmall = 3)
  np<-format(round(np, 4), nsmall = 4)
  if (is.null(ncol(indeps))){mat<-matrix(c("CONSTANT",colnames(data)[2],betas,SE,t,np), byrow=F, nrow=p+1, ncol=5)}
  else{mat<-matrix(c("CONSTANT",colnames(data[1,-1]),betas,SE,t,np), byrow=F, nrow=p+1, ncol=5)}
  df<-as.data.frame(mat)
  colnames(df)<-c("VARIABLES","COEFFICIENTS","STANDARD ERRORS","t","p-value")
  df
}
#ESTIMATION OF R2 AND MEAN SQUARE PREDICTION ERROR
fit=function(data){
  indeps=data[,-1]
  dep=data[,1]
  n=nrow(indeps)
  if (is.null(ncol(indeps))){p=1}
  else{p=ncol(indeps)}
  x<-cbind(constant = 1, as.matrix(indeps))
  y<-as.matrix(dep)
  #beta estimation
  betas=solve(crossprod(x), crossprod(x,y))
  #prediction of the dependent variables based on computed value of regression coefficients
  Y_hat=x%*%betas
  #sum of squares of residuals
  SSr=sum((y - Y_hat)^2)
  mspe=SSr/nrow(data)
  SSt <- sum((y - mean(y))^2) #among sum of squares
  R2 <- 1 - (SSr/SSt)
  adj.R2 <- 1 - ((1 - R2)*(nrow(x) - 1))/(nrow(x) - ncol(x[,-1]) - 1)
  if (p==1){adj.R2=R2}
  cat("MEAN SQUARE PREDICTION ERROR(MSPE):\t\t",mspe)
  cat("\nCOEFFICIENT OF DETERMINATION(R2):\t\t",R2)
  cat("\nADJUSTED COEFFICIENT OF DETERMINATION(Adj-R2):\t",adj.R2)
}

#checking of the normality of the data
normality = function(data,alpha=0.05){
  #It is theoretically invalid to evaluate normality with less than three observations; thus this should be checked first.	
  if(length(data)<3){
    cat("Cannot evaluate normality of data with less than three observations.\n")
  }
  
  #For data with more than 50 observations, Anderson-Darling Test shall be used.
  #Anderson-Darling Test is an improved version of Kolmogorov-Smirnov test.
  #This tests the null hypothesis that the data are normally distributed and rejects the null hypothesis when the test statistic exceeds the critical value.
  if(length(data)>50){
    cat("DESCRIPTIVE STATISTICS\n")
    cat(" Mean: ",mean(data),"\t","Standard error: ",sd(data),"\n")
    cat(" Minimum: ",min(data),"\t","Maximum: ",max(data),"\n\n")	
    
    testing = ad.test(data)
    pval = testing$p.value		
    
    cat("ANDERSON-DARLING TEST FOR NORMALITY\n")
    cat(" Test statistic, A: ",testing$statistic,"\n")
    cat(" p-value: ",pval,"\n")
    if(pval<alpha){cat(" Conclusion: At ", alpha*100,"% level of significance, data are not normally distributed.\n")}
    else{cat(" Conclusion: At ", alpha*100,"% level of significance, data are normally distributed.\n")}			
  }
  
  #Shapiro-Wilk Test for Normality is used when the number of observsations is between 3 and 50, inclusive.
  #Similar to Anderson-Darling Test, this tests the null hypothesis that the data are normally distributed. This test rejects the null hypothesis when the test statistic is less than the critical value.
  if((length(data)<=50)  && (length(data)>2)){
    cat("DESCRIPTIVE STATISTICS\n")
    cat(" Mean: ",mean(data),"\t","Standard error: ",sd(data),"\n")
    cat(" Maximum: ",max(data),"\t","Minimum: ",min(data),"\n\n")
    
    testing = shapiro.test(data)	
    pval = testing$p.value	
    
    cat("SHAPIRO-WILK TEST FOR NORMALITY\n")
    cat(" Test statistic, W: ",testing$statistic,"\n")
    cat(" p-value: ",pval,"\n")
    if(pval<alpha){cat(" Conclusion: At ", alpha*100,"% level of significance, data are  not normally distributed.\n")}
    else{cat(" Conclusion: At ", alpha*100,"% level of significance, data are normally distributed.\n")}	
  }
}

#ftest to provide evidence on the lack of fit with regards to the generated model
ftest = function(data,alpha){
  indeps=data[,-1]
  dep=data[,1]
  n=nrow(indeps)
  if (is.null(ncol(indeps))){p=1}
  else{p=ncol(indeps)}
  
  x<-cbind(constant = 1, as.matrix(indeps))
  y<-as.matrix(dep)
  
  #beta estimation
  #betas=solve(crossprod(x), crossprod(x,y))
  betas=solve(crossprod(x), crossprod(x,y))
  #prediction of the dependent variables based on computed value of regression coefficients
  Y_hat=x%*%betas
  c<-ncol(data)-1
  n<-nrow(data)
  r<-1
  #SS Regression
  SSR <- sum((y - mean(y))^2)
  
  #SS Lack of Fit
  SSLF<-0
  for(i in 1:n){
    SSLF<-SSLF+(mean(y)-Y_hat[i])^2
  }
  
  #SS Pure Error
  SSPE<- 0
  for(i in 1:n){
    SSPE<-SSLF+(y[i]-mean(y))^2
  }
  
  SSE<-SSLF+SSPE
  SSTO<-SSR+SSE
  
  MSR<-SSR/1
  MSE<-SSE/(n-2)
  MSLF<-SSLF/(c-2)
  MSPE<-SSPE/(n-c)
  FReg<-MSR/MSE
  FLack<-MSLF/MSPE
  
  cat("LACK OF FIT TEST\n\n")
  
  cat("OUTPUT\n")
  cat("Source\t \t df\t SS\t \t MS\t \t Fvalue\t\n")
  cat("Regression\t",r,"\t",SSR,"\t", MSR,"\t",FReg,"\t\n")
  cat("Residual Error\t",(n-2),"\t",SSE,"\t", MSR,"\t\t\n")
  cat("  Lack of Fit\t",(c-2),"\t",SSLF," \t",MSLF,"\t",FLack,"\t\n")
  cat("  Pure Error\t",(n-c),"\t",SSPE,"\t",MSPE,"\t\t\n")
  cat("  Total\t \t",(n-1),"\t",SSTO,"\t",MSPE,"\t\t\n\n")
  
  Ftab<-qf(1-alpha,c-2,n-c)
  cat("  Ftab\t \t",Ftab,"\n")
  if(FLack>Ftab){
    cat("Interpretation: There is lack of fit in the model.\n")
  }else{
    cat("Interpretation: There is no lack of fit in the model.\n")
  }  
  
}
#This function assumes that there is NO grouping variable in the dataset "data" because the program has its own grouping mechanism.
#The data here is the dataset which contains ALL the variables to be used in evaluating homogeneity of variances. 

hmgnt = function(data,alpha=0.05){
  
  #The program should not allow dataset which contain only one column (variable), so it is checked first.
  if(ncol(data)<2){
    cat("Dataset should have two or more variables for comparison of variances.\n")
  }
  else{
    group_vector = numeric() #A vector that indicates the group of the observations
    data_vector = numeric() #Vector of all the data
    #The following are for the descriptive statistics of the variables.
    
    for(i in 1:length(data)){
      #This identifies the group of a particular observation
      group_vector = c(group_vector,rep(i,length(data[,i])))
      #This appends all the data into just one vector
      data_vector = c(data_vector,data[,i])
      #These are done for the data to be appropriate to be sujected to Bartlett's test
    }
    dataset = data.frame(data_vector,group_vector)
    testing = bartlett.test(data_vector ~ group_vector, data = dataset)	
    
    pval = testing$p.value
    cat("BARTLETT'S TEST FOR HOMOGENEITY OF VARIANCES\n")
    cat(" Test statistic, K-squared: ",testing$statistic,"~ df =",testing$parameter,"\n")
    cat(" p-value: ",pval,"\n")
    if(pval<alpha){cat(" Conclusion: At ", alpha*100,"% level of significance, variances of the data are not homogenous.\n")}
    else{cat(" Conclusion: At ", alpha*100,"% level of significance, variances of the data are homogenous.\n")}
  }	
}
#test for the bivariate normality prior to the usage of Pearson correlation analysis
biv<-function(data,rho){
  X=data[,2]
  Y=data[,1]
  i<-mean(X)
  id<-var(X)
  j<-mean(Y)
  jd<-var(Y)	
  n<-c(X, Y)
  n1<-mean(n)
  n1
  n2<-var(n)
  n2
  a=j+(jd/id)*rho*(X-i)
  a
  b=sqrt((i-rho^2)*jd^2)
  b
  if(n1==a && n2==b){
    #ggplot(data, aes(X, Y)) + geom_point()
    return(1)
  }
  else{
    #ggplot(data, aes(X, Y)) + geom_point()	
    return(0)
  }
}

#computation of correlation
correlate=function(data){
  dependent=data[,1]
  independent=data[,-1]
  r=cov(dependent,independent)/(sqrt(var(dependent))*sqrt(var(independent)))
}

#testing the significance of the Pearson correlation coefficient
csigtest=function(data,alternative,alpha,r) {
  dependent=data[,1]
  independent=data[,-1]
  #COMPUTING TEST STAT (t DISTRIBUTED)
  k=length(dependent)-2
  tc=r/sqrt((1-r^2)/k)
  
  #TRANSFORMING TEST STAT TO Z-SCORE
  a=((8*k)+1)/((8*k)+3)
  b=sqrt(k*log(1+((tc^2)/k)))
  zc=a*b
  
  #COMPUTING THE PROBABILITY OF COMPUTING Z-SCORE
  const=c(0.319381530,-0.356563782,1.781477937,-1.821255978,1.330274429)
  z=abs(zc)
  f=exp(-(z^2)/2)/sqrt(2*pi)
  t=1/(1+(0.2316419*z))
  sum=0
  for(i in 1:length(const)) {
    sum=sum+(const[i]*(t^i))
  }
  cat(" [Hypotheses]\nNull(Ho):\t  The two variables have no significant correlation.\n")
  if (alternative=="not equal"){
    cat("Alternative(Ha):  They have a significant correlation.\n")
  } else if (alternative=="greater"){
    cat("Alternative(Ha):  They have a significant positive correlation.\n")
  } else {
    cat("Alternative(Ha):  They have a significant negative correlation.\n")
  }
  cat("Decision Rule:    Reject the null hypothesis if p-value <",alpha,"\n")
  cat("                  Otherwise, fail to reject the null hpothesis.")
  prob=1-f*sum	#PROBABILITY THAT Z IS LESS THAN THE Z-SCORE
  
  if(alternative == "greater") {
    if(tc > 0) {
      p_value = 1-prob
    }
    else {
      p_value= prob
    }
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant positive linear relationship between the two variables.")	
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant positive correlation between the two variables.")	
    }
  }
  else if(alternative == "less") {
    if(tc > 0) {
      p_value = prob
    }
    else {
      p_value = 1-prob
    }
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant negative linear relationship between the two variables.")	
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant negative correlation between the two variables.")	
    }
  }
  else {
    p_value=2*min(prob,1-prob)
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant linear relationship between the two variables.")
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant correlation between the two variables.")
    }
  }
  #cat("\n r = ",r,"\n t = ",tc,"\n df = ", k, "\n zc = ",zc,"\n prob = ",prob,"\n\n\n") 
}
#estimation of the spearman correlation coefficient when the bivariate normality is not satisfied
spearman<-function(data){
  #cor(data[1],data[2],method="spearman")
  x<-data[,1]
  y<-data[,2]
  rankX<-rank(x)					#ranks elements in vector x
  rankY<-rank(y)					#ranks elements in vector y
  d<-numeric()
  i<-0
  j<-0
  n<-length(x)
  for(i in 1:n){
    d[i]<-(rankX[i]-rankY[i])^2			#difference in ranks
  }
  D<-sum(d)						#D(included in formula)
  
  tieXdf<-as.data.frame(table(rankX))
  tieX<-numeric()
  tieX<-tieXdf[,"Freq"]				#vector of ties in ranks of X
  
  tieYdf<-as.data.frame(table(rankY))
  tieY<-numeric()
  tieY<-tieYdf[,"Freq"]				#vector of ties in ranks of Y
  
  if((length(tieX)<length(rankX))||(length(tieY)<length(rankY))){	#if any of the variable has tie ranks
    s<-0
    t<-0
    n<-length(x)
    for(i in 1:length(tieX)){
      s<-s+((tieX[i]^3)-tieX[i])
    }
    for(i in 1:length(tieY)){
      t<-t+((tieY[i]^3)-tieY[i])
    }
    C1<-((s+t)/(2*n*((n^2)-1)))
    C2<-sqrt((1-(s/(n*((n^2)-1))))*(1-(t/(n*((n^2)-1)))))
    
    rs<-(1-(6*D/(n*((n^2)-1)))-C1)/C2
    rs
    
  }else{						#if there is no tie in ranks for the both variable
    rs<-1-(6*D/(n*((n^2)-1)))
    rs
  }
}
#testing of the significance of the computed spearman correlation
ssigtest=function(data,alternative,alpha,rs){
  cat(" [Hypotheses]\nNull(Ho):\t  The two variables have no significant association.\n")
  if (alternative=="not equal"){
    cat("Alternative(Ha):  They have a significant association.\n")
  } else if (alternative=="greater"){
    cat("Alternative(Ha):  They have a significant positive association.\n")
  } else {
    cat("Alternative(Ha):  They have a significant negative association.\n")
  }
  cat("Decision Rule:    Reject the null hypothesis if p-value <",alpha,"\n")
  cat("                  Otherwise, fail to reject the null hpothesis.")
  
  zc=rs*sqrt(nrow(data)-1)
  #COMPUTING THE PROBABILITY OF COMPUTING Z-SCORE
  const=c(0.319381530,-0.356563782,1.781477937,-1.821255978,1.330274429)
  z=abs(zc)
  f=exp(-(z^2)/2)/sqrt(2*pi)
  t=1/(1+(0.2316419*z))
  sum=0
  for(i in 1:length(const)) {
    sum=sum+(const[i]*(t^i))
  }
  prob=1-f*sum	#PROBABILITY THAT Z IS LESS THAN THE Z-SCORE
  if(alternative == "greater") {
    if(Zc > 0) {
      p_value = 1-prob
    }
    else {
      p_value= prob
    }
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant positive association between the two variables.")	
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant positive association between the two variables.")	
    }
  }
  else if(alternative == "less") {
    if(tc > 0) {
      p_value = prob
    }
    else {
      p_value = 1-prob
    }
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant negative association between the two variables.")	
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant negative association between the two variables.")	
    }
  }
  else {
    p_value=2*min(prob,1-prob)
    if(p_value < alpha) {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value <", alpha, ", then reject Ho.\n")
      cat("\nThere is a significant association between the two variables.")
    }
    else {
      if(p_value == 0) {
        cat("\np-value:\t  <0.001 \n")
      }
      else {
        cat("\np-value:\t  ", p_value, "\n")
      }
      cat("Decision: \t  Since p-value >", alpha, ", then fail to reject Ho.\n")
      cat("\nThere is no significant association between the two variables.")
    }
  }
}
#SHINY SERVER
server=shinyServer(function(input,output,session){
  #file reading from the data uploaded by the user
  data=reactive({
    file1=input$file
    if(is.null(file1)){return ()}
    read.table(file=file1$datapath,sep=input$sep,header=TRUE)
  })
  #acquiring the desired dependent and independent variable
  observe({
    s_options<-colnames(data())
    updateSelectInput(session, "yvar",choices = s_options)
    updateSelectInput(session, "xvar",choices = s_options)
  })
  #creation of data designated for summary statistics
  selectedSUMData <- reactive({
    data()[, c(input$yvar, input$xvar)]
  })
  #generation of random seed for reproducibility of the sample for cross validation
  .Random.seed
  selectedMData <- reactive({
    x=sample(1:nrow(data()),size=nrow(data())-input$training,replace=FALSE)
    if(length(x)!=0){data()[-x, c(input$yvar, input$xvar)]}
    else{data()[, c(input$yvar, input$xvar)]}
  })
  #these data catches the sample removed from the training data set selectedMData()
  Observed=reactive({
    x=sample(1:nrow(data()),size=nrow(data())-input$training,replace=FALSE)
    data()[x,input$yvar]
  })
  #This outputs the selected dependent and independent variables for the user to see the selected data
  output$table=renderTable({
    if(is.null(input$xvar)){return ()}
    selectedSUMData()
  })
  
  observe({
    s_options<-colnames(data())
    updateSelectInput(session, "svar",choices = s_options)
  })
  #this catches the single variable to be summarised and tested
  selectedSData=reactive({
    data()[, c(input$svar)]
  })
  #outputs the summary statistics
  output$summary=renderPrint({
    if(is.null(data())){return ()}
    cat("\nMean:\t\t",mean(selectedSData()))
    cat("\nMinimum:\t",min(selectedSData()))
    cat("\nFirst Quartile:\t",quantile(selectedSData(),0.25))
    cat("\nMedian:\t\t",quantile(selectedSData(),0.50))
    cat("\nThird Quartile:\t",quantile(selectedSData(),0.75))
    cat("\nMaximum:\t",max(selectedSData()))
    cat("\n\n\n")
  })
  #this visualizes the data from the selected single variable
  output$graph=renderPlot({
    if(is.null(input$xvar)){return ()}
    boxplot(selectedSData(),main=" ",xlab=input$svar,col="darkslategray",horizontal=TRUE)
  })
  #outputs the dependent variable
  output$model_desc_y=renderPrint({
    cat(input$yvar)
  })
  #outputs the independent variable
  output$model_desc_x=renderPrint({
    cat(input$xvar)
  })
  #outputs the model via ordinary least squares
  output$model=renderTable({
    if(is.null(input$xvar)){return ()}
    ols(selectedMData())
  })
  #outputs the adj.R2, R2,MSPE
  output$fit=renderPrint({
    if(is.null(input$xvar)){return ()}
    fit(selectedMData())
  })
  #calls the normality function which test the normality of the errors
  output$normality=renderPrint({
    if(is.null(input$alpha)){return ()}
    normality(selectedSData(),input$alpha)
  })
  #calls the hmgnt function which test the homogeneity of the error variances
  output$hmg=renderPrint({
    if(is.null(input$alpha)){return ()}
    hmgnt(selectedMData(),input$alpha)
  })
  
  observe({
    s_options<-colnames(data())
    updateSelectInput(session, "one",choices = s_options)
    updateSelectInput(session, "two",choices = s_options)
  })
  #creation of the data subjected to correlation analysis
  corData=reactive({
    data()[, c(input$one,input$two)]
  })
  #outputs the interpretation from the bivariate normality tests
  output$bnorm=renderPrint({
    if(is.null(corData())){return ()}
    if(biv(corData(),correlate(corData()))==1){cat("The data follows a bivariate normal distribution. With this, Pearson correlation coefficient is valid, and will be used as a measure of association.")}
    else {cat("The data does not follow bivariate normal distribution. With this, Pearson correlation coefficient is invalid, and Spearman will be used as a measure of association.")}
  })
  #calculates the association coefficients accounting parametric and nonparametric methods
  output$association=renderText({
    if(is.null(corData())){return ()}
    if(biv(corData(),correlate(corData()))==1){correlate(corData())}
    else{spearman(corData())}
  })
  #outputs the interpretation of the association coefficient from the results generated spearman and correlation function
  output$assocint=renderPrint({
    if(is.null(corData())){return ()}
    if(biv(corData(),correlate(corData()))==1){
      if(correlate(corData())==0){
        cat("There is no correlation between the two.",sep=" ")
      } else if(correlate(corData())<0){
        if(abs(correlate(corData()))>0 && abs(correlate(corData()))<0.2){
          cat("There is a very weak negative correlation between the two.",sep='space',fill=TRUE)
        } else if(abs(correlate(corData()))>=0.2 && abs(correlate(corData()))<0.4){
          cat("There is a weak negative correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.4 && abs(correlate(corData()))<0.6){
          cat("There is a moderate negative correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.6 && abs(correlate(corData()))<0.8){
          cat("There is a strong negative correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.8 && abs(correlate(corData()))<1){
          cat("There is a very strong negative correlation between the two.",sep=" ",fill=TRUE)
        } else {
          cat("There is a perfect negative correlation between the two.",sep=" ",fill=TRUE)
        }
      } else {
        if(abs(correlate(corData()))>0 && abs(correlate(corData()))<0.2){
          cat("There is a very weak positive correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.2 && abs(correlate(corData()))<0.4){
          cat("There is a weak positive correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.4 && abs(correlate(corData()))<0.6){
          cat("There is a moderate positive correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.6 && abs(correlate(corData()))<0.8){
          cat("There is a strong positive correlation between the two.",sep=" ",fill=TRUE)
        } else if(abs(correlate(corData()))>=0.8 && abs(correlate(corData()))<1){
          cat("There is a very strong positive correlation between the two.",sep=" ",fill=TRUE)
        } else {
          cat("There is a perfect positive correlation between the two.",sep=" ",fill=TRUE)
        }
      }
    } else {
      if(spearman(corData())==0){
        cat("There is no correlation between the two.",sep=" ")
      } else if(spearman(corData())<0){
        if(abs(spearman(corData()))>0 && abs(spearman(corData()))<0.2){
          cat("There is a very weak negative association between the two.",sep='space',fill=TRUE)
        } else if(abs(spearman(corData()))>=0.2 && abs(spearman(corData()))<0.4){
          cat("There is a weak negative association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.4 && abs(spearman(corData()))<0.6){
          cat("There is a moderate negative association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.6 && abs(spearman(corData()))<0.8){
          cat("There is a strong negative association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.8 && abs(spearman(corData()))<1){
          cat("There is a very strong negative association between the two.",sep=" ",fill=TRUE)
        } else {
          cat("There is a perfect negative association between the two.",sep=" ",fill=TRUE)
        }
      } else {
        if(abs(spearman(corData()))>0 && abs(spearman(corData()))<0.2){
          cat("There is a very weak positive association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.2 && abs(spearman(corData()))<0.4){
          cat("There is a weak positive association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.4 && abs(spearman(corData()))<0.6){
          cat("There is a moderate positive association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.6 && abs(spearman(corData()))<0.8){
          cat("There is a strong positive association between the two.",sep=" ",fill=TRUE)
        } else if(abs(spearman(corData()))>=0.8 && abs(spearman(corData()))<1){
          cat("There is a very strong positive association between the two.",sep=" ",fill=TRUE)
        } else {
          cat("There is a perfect positive association between the two.",sep=" ",fill=TRUE)
        }
      }
    }
  })
  #tests the significance of the association coefficients accounting the parametric and nonparametric methods
  output$sigcor=renderPrint({
    if(is.null(corData())){return ()}
    r=correlate(corData())
    rs=spearman(corData())
    if(biv(corData(),correlate(corData()))==1){csigtest(corData(),input$dir,input$alphacor,r)}
    else {ssigtest(corData(),input$dir,input$alphacor,rs)}
  })
  #graphs the scatter plot of the two variables to visually assess the association between the two
  output$scatter=renderPlot({
    if(is.null(corData())){return ()}
    plot(corData()[,1], corData()[,2],xlab = input$two, ylab = input$one, pch = 16, col = "darkslategray", cex = 1) 
    abline(lm(corData()[,1]~corData()[,2]), col="grey", lwd = 2) 
  })
  #allows the user to modify the traning data set
  output$training=renderUI({
    sliderInput("training","Size of Training Dataset:",step=nrow(data())%%2+1,min=0,max=nrow(data()),value=0.5*nrow(data()))
  })
  observe({
    updateSliderInput(session,"training",step=nrow(data())%%2+1,min=0,max=nrow(data()),value=0.5*nrow(data()))
  })
  #outputs the size of the training data set
  output$ntraining=renderPrint({
    cat(input$training)
  })
  #outputs the size of the validating data set
  output$nvalidating=renderText({
    validating=nrow(data())-input$training
  })
  #outputs the data sampled and classified as one of the validating data set and their predicted values
  output$tablevalidate=renderTable({
    if(is.null(input$xvar)){return ()}
    indeps=selectedMData()[,-1]
    dep=selectedMData()[,1]
    n=nrow(indeps)
    if (is.null(ncol(indeps))){p=1}
    else{p=ncol(indeps)}
    x<-cbind(constant = 1, as.matrix(indeps))
    y<-as.matrix(dep)
    #beta estimation
    betas=solve(crossprod(x), crossprod(x,y))
    Predicted=as.vector(x%*%betas)
    Observed=Observed()
    data.frame(Observed,Predicted)
  })
  #function for computing p-value on detecting difference between the actual and predicted values
  pvalcomp=function(data1,data2){
    indeps=data1[,-1]
    dep=data1[,1]
    n=nrow(indeps)
    if (is.null(ncol(indeps))){p=1}
    else{p=ncol(indeps)}
    x<-cbind(constant = 1, as.matrix(indeps))
    y<-as.matrix(dep)
    #beta estimation
    betas=solve(crossprod(x), crossprod(x,y))
    Predicted=as.vector(x%*%betas)
    diff=c()
    for(i in 1:length(Observed())){ #computes for the difference
      diff[i]=data2[i]-Predicted[i]
    }
    nprime=length(diff[diff!=0])
    cval=length(diff[diff>0]) #determines the critical value for the sign test
    p_value=binom.test(cval,nprime,alternative="two.sided",p=0.5)$p.value
    return(p_value)
  }
  #test for the difference of the predicted and the actual values of the dependent variable using binomial exact test
  output$testvalidate=renderPrint({
    if(is.null(input$xvar)){return ()}
    p_value=pvalcomp(selectedMData(),Observed())
    cat("p-value=",p_value)
  })
  #outputs the conclusion derived from the test of the difference between the actual and predicted values of the dependent variables
  output$testvalidateint=renderPrint({
    if(is.null(input$xvar)){return ()}
    p_value=pvalcomp(selectedMData(),Observed())
    if(p_value<input$alpha){
      cat("There is a significant difference between the actual and predicted values of",input$yvar)
    }else{
      cat("There is no significant difference between the actual and predicted values of",input$yvar)
    }
  })
  #calculates the vif needed for model diagnostics
  output$vif=renderPrint({
    if(is.null(input$xvar)){return ()}
    indeps=selectedMData()[,-1]
    dep=selectedMData()[,1]
    n=nrow(indeps)
    if (is.null(ncol(indeps))){p=1}
    else{p=ncol(indeps)}
    x<-cbind(constant = 1, as.matrix(indeps))
    y<-as.matrix(dep)
    #beta estimation
    betas=solve(crossprod(x), crossprod(x,y))
    #prediction of the dependent variables based on computed value of regression coefficients
    Y_hat=x%*%betas
    SSr=sum((y - Y_hat)^2)
    mspe=SSr/nrow(data)
    SSt <- sum((y - mean(y))^2)
    R2 <- 1 - (SSr/SSt)
    adj.R2 <- 1 - ((1 - R2)*(nrow(x) - 1))/(nrow(x) - ncol(x[,-1]) - 1)
    if (p==1){adj.R2=R2}
    vif=1/(1-adj.R2)
    cat("MULTICOLLINEARITY\n")
    cat("VIF = ",vif,"\n")
    if(vif>5){
      cat("Interpretation: With the variation inflation factor greater than 5, multicollinearity exists among the independent variable.\n")  
    }else{
      cat("Interpretation: With the variation inflation factor less than 5, there exist no multicollinearity among the predictors.\n")  
    }
  })
  #tests the significance of the model using ftest
  output$ftest=renderPrint({
    if(ncol(selectedMData()[,-1])<3){return ()}
    ftest(selectedMData(),input$alphacor) 
  })
  #outputs the normality of errors
  output$normal=renderPrint({
    if(is.null(input$xvar)){return ()}
    indeps=selectedMData()[,-1]
    dep=selectedMData()
    n=nrow(indeps)
    if (is.null(ncol(indeps))){p=1}
    else{p=ncol(indeps)}
    
    x<-cbind(constant = 1, as.matrix(indeps))
    y<-as.matrix(dep)
    #beta estimation
    #betas=solve(crossprod(x), crossprod(x,y))
    betas=solve(crossprod(x), crossprod(x,y))
    #prediction of the dependent variables based on computed value of regression coefficients
    Y_hat=x%*%betas
    errors=y-Y_hat
    normality(errors,input$alpha)
  })
})

ui=shinyUI(navbarPage(
  "Correlation |  Association  |  Regression Analysis",
  theme = shinythemes::shinytheme("flatly"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition="input.panel==1",
                       fileInput("file",label=h4("File Input")),
                       tags$hr(style="border-color: gray;"),
                       radioButtons(inputId='sep',label='Data Separator',choices=c(Comma=',',Semicolon=';',Tab='\t',Space=' ')),
                       tags$hr(style="border-color: gray;"),
                       selectInput("yvar","Dependent variable:",names(data())),
                       selectInput("xvar","Independent variable/s:",multiple=TRUE,names(data()))
      ),
      conditionalPanel(condition="input.panel==2",
                       h1("Summary Statistics",align='center'),
                       selectInput("svar","Variable:",names(data())),
                       verbatimTextOutput("summary")
      ),
      conditionalPanel(condition="input.panel==3",
                       h1("Correlation",align='center'),
                       selectInput("one","First variable:",names(data())),
                       selectInput("two","Second variable:",names(data())),
                       tags$hr(style="border-color: gray;"),
                       h4("Coefficient:",align='center'),tags$i(uiOutput("association",align='center')),uiOutput("assocint",align='center'),
                       tags$hr(style="border-color: gray;"),
                       selectInput("dir","Direction[Alternative Hypothesis]",choices=c("not equal","less","greater")),
                       numericInput("alphacor","Level of significance",value=0.05,min=0,max=1)
      ),
      conditionalPanel(condition="input.panel==4",
                       h2("Via Ordinary Least Squares method",align='center'),
                       tags$hr(style="border-color: gray;"),
                       h4("Modeling: "),
                       verbatimTextOutput("model_desc_y"),
                       h4("Based on:"),
                       verbatimTextOutput("model_desc_x",placeholder=TRUE),
                       uiOutput("training")
      ),
      conditionalPanel(condition="input.panel==5",
                       h3("Model Validation",align='center'),
                       tags$hr(style="border-color: gray;"),
                       h5("Size of Training Dataset:"),
                       verbatimTextOutput("ntraining"),
                       h5("Size of Validating Dataset:"),
                       verbatimTextOutput("nvalidating",placeholder=TRUE),
                       tags$hr(style="border-color: gray;"),
                       h4("Test for the Difference between the Actual and Predicted Values",align='center'),
                       uiOutput("testvalidate",align='center'),
                       tags$hr(style="border-color: gray;"),
                       uiOutput("testvalidateint",align='center')
      ),
      conditionalPanel(condition="input.panel==6",
                       h1("Model Assumptions",align='center'),
                       tags$hr(style="border-color: gray;"),
                       #h5("Absence of Outliers",align='center'),
                       h5("Normality of Errors",align='center'),
                       h5("Independence of Errors",align='center'),
                       h5("Homogeneity of Error Variances",align='center'),
                       h5("Non-multicollinearity",align='center'),
                       tags$hr(style="border-color: gray;"),
                       numericInput("alpha","Level of significance",value=0.05,min=0,max=1)
      )
    ),
    mainPanel(
      tabsetPanel(id="panel",
                  tabPanel("Data",value=1,uiOutput("table")),
                  tabPanel("Summary",value=2,plotOutput("graph")),
                  tabPanel("Correlation",value=3,tags$hr(style="border-color: gray;"),h4("Test on the Bivariate Normality of the Two Variables",align='center'),uiOutput("bnorm",align='center'),tags$hr(style="border-color: gray;"),h3("Test on the Significance of the Association",align='center'),verbatimTextOutput("sigcor"),plotOutput("scatter")),
                  tabPanel("Model Building",value=4,tags$hr(style="border-color: gray;"),verbatimTextOutput("fit"),tags$hr(style="border-color: gray;"),uiOutput("model",align='center')),
                  tabPanel("Cross Validation",value=5,uiOutput("tablevalidate")),
                  tabPanel("Diagnostics",value=6,verbatimTextOutput("normal"),verbatimTextOutput("hmg"),verbatimTextOutput("vif"),verbatimTextOutput("ftest"))
                  #tabPanel("Diagnostics",value=6,verbatimTextOutput("normal"),verbatimTextOutput("hmg"),verbatimTextOutput("vif"))
        )
      )
    )
  )
)

shinyApp(ui,server)