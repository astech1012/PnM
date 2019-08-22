#install.packages("shiny")
#install.packages("shinythemes")

#  some inspo from https://rich.shinyapps.io/regression/
# odbc()
# con <- DBI::dbConnect(odbc::odbc(),
#                       Driver       = "/opt/amazon/redshift/lib/libamazonredshiftodbc.dylib",
#                       servername   = "rtbcluster-west2.cjgay0xd9jwq.us-west-2.redshift.amazonaws.com",
#                       database     = "rtb",
#                       UID          = "rtb_read",
#                       PWD          = "RTB_password1",
#                       Port         = 5439)

#  Spend Installs CPI GOAL (LIST SIZE)




library(shiny)
library(shinythemes)
library(googleVis)
#Date functions
#st1 <- tail(King$date_cet)[1]
#en1 <- head(King$date_cet)[1]
#subset(King, date_cet >= st1 & date_cet <= en1)
suppressPackageStartupMessages(library(googleVis))
King <- read.csv(file="/Users/Aarki/Desktop/KingModel.csv",header=TRUE)
King$date_cet <- as.Date(King$date_cet)
King<-King[-nrow(King),]
King <- as.data.frame(King)


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
} # need F-STAT and p-value

coeffs = function(data){
  indeps=data[,-1]
  dep=data[,1]
  n=nrow(indeps)
  if (is.null(ncol(indeps))){p=1}
  else{p=ncol(indeps)}
  x<-cbind(constant = 1, as.matrix(indeps))
  y<-as.matrix(dep)
  #beta estimation
  betas=solve(crossprod(x), crossprod(x,y))
  cat(colnames(data)[1], " = ", betas[1]) #colnames(data[1:nrow(indeps),-1])
  for(i in 1:ncol(indeps)) {
    #cat("FUCK")}
    if(betas[i+1]<0){
      cat(" ",betas[i+1],"*",colnames(indeps[i]))
    }
    else{
      cat(" + ", betas[i+1],"*",colnames(indeps[i]))
    }
  }
}

#cat(colnames(indeps[i]))}
#dat <- read.csv(file="/Users/Aarki/Documents/Project/sampledata.csv",header=TRUE)
#ols(dat)
#fit(dat)
#model <- lm(CRUDE ~ INTEREST + FOREIGN + DJIA  + GNP + PURCHASE + CONSUMER,data = King)
#summary(model)
#SHINY SERVER
#SHINY SERVER

server=shinyServer(function(input, output,session) {
  # Filter data based on selections
  data <- King
  selectedData<-reactive({
    #Date
    #if (input$DateRange[,2] != en1) {
    #  subset(data, date_cet >= start & date_cet <= end)
    #}
    if (input$app != "All") {
      data <- data[data$app == input$app,]
    }
    if (input$spl != "All") {
      data <- data[data$campaign_split == input$spl,]
    }
    data$date_cet <- as.character(data$date_cet)
    data
  })
  selectedMData <- reactive({
    selectedData()[, c(input$yvar, input$xvar)]
  })
  ## pwede if else
  observe({
    s_options<-colnames(selectedData())
    updateSelectInput(session, "yvar",choices = s_options)
    updateSelectInput(session, "xvar",choices = s_options)
  })
  output$tableset <- renderTable({selectedData()[1:15,]})
  output$summary <- renderPrint({
    if(is.null(selectedData())){return ()}
    summary(cbind(selectedData()[input$yvar], selectedData()[input$xvar]))
  })
  output$model=renderTable({
    #selectedMData()
    #if(is.null(input$xvar)){return ()}
    ols(selectedMData())
  })
  output$fit=renderPrint({
    if(is.null(input$xvar)){return ()}
    fit(selectedMData())
  })
  output$betas=renderPrint({
    if(is.null(input$xvar)){return ()}
    coeffs(selectedMData())
  })

})


ui=shinyUI(navbarPage(
  "KING RT",
  theme = shinythemes::shinytheme("cosmo"),
  tabPanel("Model Building",
           sidebarLayout(
             sidebarPanel(
               conditionalPanel(condition="input.panel==1",
                                #radioButtons(inputId='select',label='Data Selection',choices=c(King='king',"Input Own"='own'),inline=TRUE),
                                #fileInput("file",label=h4("File Input")),
                                #tags$hr(style="border-color: gray;"),
                                #radioButtons(inputId='sep',label='Data Separator',choices=c(Comma=',',Semicolon=';',Tab='\t',Space=' ')),
                                # tags$hr(style="border-color: gray;"),
                                h1("Select Data",align='center'),
                                selectInput("app","App:", c("All",unique(as.character(King$app)))),
                                selectInput("spl","Campaign Split:", c("All",unique(as.character(King$campaign_split)))),
                                selectInput("yvar","Dependent variable:",names(data())),
                                selectInput("xvar","Independent variable/s:",multiple=TRUE,names(data()))
               ),
               conditionalPanel(condition="input.panel==2",
                                h1("Summary Statistics",align='center')
               ),
               conditionalPanel(condition="input.panel==3",
                                h1("Correlation",align='center'),
                                selectInput("one","First variable:",names(data())),
                                selectInput("two","Second variable:",names(data())),
                                tags$hr(style="border-color: gray;"),
                                h4("Coefficient:",align='center'),tags$i(uiOutput("association",align='center')),uiOutput("assocint",align='center')   ),
               conditionalPanel(condition="input.panel==4",
                                h2("Linear Regression",align='center'),
                                tags$hr(style="border-color: gray;")
                                #h4("Model: CONSTANT + XVAR * ")
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
                           tabPanel("Data",value=1,uiOutput("tableset")),
                           tabPanel("Summary",value=2,verbatimTextOutput("summary")),
                           tabPanel("Correlation",value=3,htmlOutput("corr")),
                           tabPanel("Model Building",value=4,tags$hr(style="border-color: gray;"),verbatimTextOutput("fit"),tags$hr(style="border-color: gray;"),uiOutput("model",align='center'),tags$hr(style="border-color: gray;"),uiOutput("betas",align='center')),
                           tabPanel("Cross Validation",value=5,uiOutput("tablevalidate")),
                           tabPanel("Diagnostics",value=6,verbatimTextOutput("normal"),verbatimTextOutput("hmg"),verbatimTextOutput("vif"),verbatimTextOutput("ftest"))
                           #tabPanel("Diagnostics",value=6,verbatimTextOutput("normal"),verbatimTextOutput("hmg"),verbatimTextOutput("vif"))
               )
             )
           )
  ),
  tabPanel("More", img(src='/Users/Aarki/Desktop/Ty.png'))
)
)

shinyApp(ui,server)
