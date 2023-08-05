library(shiny)
library(statmod)
library(EQL);library(esaddle)

# Define UI for random distribution app ----
ui <- fluidPage(
  
  withMathJax(),
  
  # App title ----
  titlePanel(h3("Inverse Gaussian Distribution")),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      

      
      # Input: Select the random distribution type ----
      div("Probability density function:",style="text-indent:20px;font-size:125%;"),

      withMathJax(helpText("$$f(x;\\mu,\\lambda)= \\left(\\frac{\\lambda}{2\\pi x^3}\\right)^\\frac{1}{2} \\exp \\left[ \\frac{-\\lambda(x-\\mu)^2}{2\\mu^2 x} \\right]$$",style="text-indent:20px;font-size:100%;")),
      
      # div(textOutput("fdp"),style="text-indent:20px;font-size:125%;"),
      
      # br() element to introduce extra vertical spacing ----
      br(),
      

      
      sliderInput("sliderMu", 
                  "\\(\\mu\\) parameter value:", 
                  min = 1, 
                  max = 10, 
                  value = 1, 
                  animate = animationOptions(interval = 300)),
      
        sliderInput("sliderLambda", 
                    "\\(\\lambda\\) parameter value:", 
                    min = 1, 
                    max = 30, 
                    value = 2, 
                    animate = animationOptions(interval = 200))
      

    ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Density", plotOutput("densityplot")),
                  tabPanel("Density Approximations", align="center",
                           br(),
                           numericInput("Nn", 
                                        label="Number of random variables in the sum",
                                        value=10,min=1,max=30),
                           
                           plotOutput("edgeworthplot"))
      )
      
    )
  )
)

server <- function(input, output) {
  

  output$densityplot <- renderPlot({
        mu <- input$sliderMu
       lambda<-input$sliderLambda
    
    xrange <- seq(from=0.001 , to=2*1+5/lambda ,by =0.009)
    xprob <- dinvgauss (xrange , mean = mu , shape = lambda )
    mainlabel <- expression (paste ( "Probability Density IG ( " ,mu,", ",lambda," )", sep = " " ) )
    plot(xrange , xprob , type = "l" , main = mainlabel , cex=.9, xlab="x" , ylab = "density",lty=1,lwd=3,col="darkblue")
    legend("topright",cex=1.7, bty="n",
           legend = c(sapply(paste("",mu,"   "), function(x) as.expression(substitute(mu==B,list(B=as.name(x))))),
                      sapply(paste(lambda), function(x) as.expression(substitute(lambda==C,list(C=as.name(x)))))) )
    
  
    })
  
  
  output$edgeworthplot <- renderPlot({
    mu <- input$sliderMu # mu <- 1
    lambda<-input$sliderLambda # lambda <- 3
    n <- input$Nn
    xrange <- seq (from=0.001 , to=2*mu+5/lambda ,by =0.009)

    #Cumulants
    K1 <- mu
    K2 <- mu^3/lambda
    K3 <- 3*mu^5/lambda^2
    K4 <- 15*mu^7/lambda^3

    K3_st <- K3/K2^(3/2) #assim        rho3
    K4_st <- K4/K2^2 #kurt             rho4

    variancia<-mu^3/lambda

    # n<-30
    R=100000 ; simu <- matrix(NA, ncol=n, nrow=R)

    for(i in 1:n){
      simu[,i] <- rinvgauss(R,mean=mu,shape=lambda)
    }

    Srows <- apply(simu,1,sum)
    Sn <- Srows/n
    
    mainLABEL=expression (paste ( "Density Approximations of  Z = ",frac(S[n],n), sep = " " ) )
    plot(density(Sn), lwd=2, main=mainLABEL, xlab=expression(paste(S[n],"=",Y[1],"+...+",Y[n])), ylab=c("Density",paste("n=",n)),cex.lab = 1.1) 
    
    edgew<-edgeworth(xrange, n, rho3=K3_st, rho4=K4_st, mu=mu, sigma2=variancia, type="mean")
    lines(xrange,edgew$approx, col="red",lwd=2)
    
    sddpoint <- dsaddle(y=xrange, X=Sn,decay=0.05,log=T)
    lines(xrange, exp(sddpoint$llk), col="darkgreen",lwd=1.5)

    lines(xrange, n*dnorm(n*xrange,mean=n*mu,sd=sqrt(n*variancia)),col="blue",lwd=2)
    
    
    legend("topright", lty=1,lwd=2,col=c("black","red","blue","darkgreen"), legend=c("Empirical   ","Edgeworth   ","Normal   ","Saddle Point  "),bty="n",pt.bg = "white",cex=1.5,y.intersp = 1.3)
    # 
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)





