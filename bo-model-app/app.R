#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggiraph)

# rs_df <- read_csv('../test.csv')
water_df <- read_csv("data/water.csv") 
# water_df <- read_csv("bo-model-app/data/water.csv") 
abs_water = data.frame(wl = water_df$lambda_nm,
                       abs_water = water_df$absorption_m1)
backs_water = data.frame(wl = water_df$lambda_nm,
                         backs_water = 0.5*(water_df$scattering_m1))

source("R/app-functions.R")

in1 <- sliderInput("absCHL443",
                   "chl absorb at 443:",
                   min = 6.12e-3,
                   max = 9.04e-1,
                   value = 0.143)
in2 <- sliderInput("absNAP443",
                   "NAP absorb at 443:",
                   min = 2.05e-3,
                   max = 1.73,
                   value = 0.139)
in3 <- sliderInput("absCDOM443",
                   "CDOM absorb at 443:",
                   min = 9.27e-5,
                   max = 1,
                   value = 5.01e-1)
in4 <- sliderInput("sNAP",
                   "NAP spectral slope:",
                   min = 1.1e-2,
                   max = 1.5e-2,
                   value = 1.3e-2)
in5 <- sliderInput("sCDOM",
                   "CDOM spectral slope:",
                   min = 1.5e-2,
                   max = 1.9e-2,
                   value = 1.7e-2)
in6 <- sliderInput("sbCHL",
                   "CHL scatter slope:",
                   min = 5.23e-5,
                   max = 2,
                   value = 1.01)
in7 <- sliderInput("sbSED",
                   "sed scatter slope:",
                   min = 9.38e-5,
                   max = 2,
                   value = 1)
in8 <- sliderInput("bbCHL",
                   "chl bb:",
                   min = 1e-8,
                   max = 3.94e-1,
                   value = 1.73e-2)
in9 <- sliderInput("bbSED",
                   "sed bb:",
                   min =  1e-8,
                   max =  6.59e-1,
                   value =  4.64e-2)
in10 <- sliderInput("mu0",
                   "zenith angle:",
                   min =  0,
                   max =  1.4,
                   step = 0.1,
                   value =  1)
in11 <- sliderInput("scatterCHL443",
                   "chl scatter at 443:",
                   min =  1.81e-5,
                   max =  2.08e1,
                   value =  1.85)
in12 <- sliderInput("scatterSED443",
                   "sed scatter at 443:",
                   min =  1.04e-4,
                   max =  2.06e1,
                   value =  1.85)
# outputs
out1 <- ggiraphOutput("distPlot")
out2 <- ggiraphOutput("backsPlot")
out3 <- ggiraphOutput("absPlot")
# out2 <- tableOutput("mydf")


# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Simulate Reflectance Spectra"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            tabsetPanel(
                tabPanel("set1", in1, in2, in3, in4, in5, in6),
                tabPanel("set2", in7, in8, in9, in11, in12, in10))),
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Reflectance", out1),
                tabPanel("Absorbance", out3),
                tabPanel("Backscattering", out2))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    abs_df <- reactive({
        
        ABSnap_spectra <- calc_absorb_spectra(input$absNAP443, aSF = input$sNAP) %>%
            rename(wl = wavelength, abs_nap = abs_m1) %>%
            dplyr::select(wl, abs_nap)
        ABScdom_spectra <- calc_absorb_spectra(input$absCDOM443, aSF = input$sCDOM) %>%
            rename(wl = wavelength, abs_cdom = abs_m1) %>%
            dplyr::select(wl, abs_cdom)
        
        abs_df <- abs_water %>%
            left_join(ABSnap_spectra) %>%
            left_join(ABScdom_spectra) %>%
            mutate(abs_total = abs_water + abs_nap + abs_cdom)
        
        abs_df
    })
    
    BSchl_spectra <- reactive({
        calc_scatter_spectra(scatterREF = input$scatterCHL443, 
                             bSF = input$sbCHL, 
                             bratio = input$bbCHL/input$scatterCHL443) %>%
            rename(wl = wavelength, bb_chl = backscatter_m1) %>%
            dplyr::select(wl, bb_chl)
        
    })    
    bratio_sed <- reactive({
        input$bbSED/input$scatterSED443
    })
    output$bratio1 <- renderText({bratio_sed()})
    BSsed_spectra <- reactive({
        calc_scatter_spectra(scatterREF = input$scatterSED443, 
                             bSF = input$sbSED, 
                             bratio = bratio_sed()) %>%
        rename(wl = wavelength, bb_sed = backscatter_m1) %>%
        dplyr::select(wl, bb_sed)
        })
    
    backs_df <- reactive({
        backs_df <- backs_water %>%
            left_join(BSchl_spectra()) %>%
            left_join(BSsed_spectra()) %>%
            mutate(backs_total = backs_water + bb_chl + bb_sed)
        backs_df
    })
    
    rs_df <- reactive({
        C_mu0 = -0.629 * input$mu0 + 0.975
        rs_df <- abs_df() %>% 
            left_join(backs_df()) %>%
            mutate(rD = 0.544 * C_mu0 * (backs_total/ (abs_total + backs_total)))
        rs_df

    })
    # output$mydf <- renderTable({
    #     head(backs_df())
    # })
    output$distPlot <- renderGirafe({

        rs1 <- rs_df() %>%
            # dplyr::filter(wl %in% water_wls) %>%
            ggplot(aes(x = wl, y = rD)) +
            geom_line() +
            theme_bw() +
            ggtitle("Reflectance")

        gg1 <- rs1 +
            geom_point_interactive(aes(tooltip = wl, data_id = wl), size = 1)
        x <- girafe(code = print(gg1))
        x    })
    
    output$absPlot <- renderGirafe({
        rs3 <- rs_df() %>% 
            ggplot(aes(x = wl, y = abs_total)) +
            geom_line(lwd = 1) +
            geom_line(aes(y = abs_nap), col = "brown") +
            geom_line(aes(y = abs_cdom), col = "orange") +
            geom_line(aes(y = abs_water), col = "blue") +
            theme_bw() +
            ggtitle("Absorbance")
        gg3 <- rs3 +
            geom_point_interactive(aes(tooltip = abs_total, data_id = abs_total), size = 1)
        x3 <- girafe(code = print(gg3))
        x3   
    })
    output$backsPlot <- renderGirafe({
        
        rs2 <-  rs_df() %>% 
            # dplyr::filter(wl %in% water_wls) %>%
            ggplot(aes(x = wl, y = backs_total)) +
            geom_line(lwd = 1) +
            geom_line(aes(y = bb_chl), col = "green") +
            geom_line(aes(y = bb_sed), col = "orange") +
            geom_line(aes(y = backs_water), col = "blue") +
            theme_bw() +
            ggtitle("Backscattering")
        
        gg2 <- rs2 +
            geom_point_interactive(aes(tooltip = backs_total, data_id = backs_total), size = 1)
        x2 <- girafe(code = print(gg2))
        x2    })
}

# Run the application 
shinyApp(ui = ui, server = server)
