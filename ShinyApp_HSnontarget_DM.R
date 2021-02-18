library(tidyverse)
library(shiny)
library(plotly)
# loading in the first list of the result from homol.search() [nontarget package; Martin Loos]
homol_test <- read.csv("homol[[1]]_posmode_3er_8-250mz_05-4RT.csv", header = T)

homo_tibble <- tibble::as_tibble(homol_test) 

homo_tibble <- homo_tibble %>%
  mutate(hsid=homo_tibble$`HS.IDs`) %>%
  mutate(mzsplit=homo_tibble$'m.z.increment') %>%
  tidyr::separate_rows(c('hsid', 'mzsplit'), sep="/") %>%
  mutate(hsid=as.factor(hsid)) 

#summary(homo_tibble)
homo_filter <- homo_tibble %>%
  filter(homo_tibble$"HS.IDs" != "0")

homo_filter$mzsplit <- as.numeric(homo_filter$mzsplit)

# building the ui for the shiny app

ui <-   fluidPage(
  titlePanel("m/z increment"),
  sidebarLayout(sidebarPanel(helpText("Choose the range of m/z increment of interest"),
                             numericInput("min", "Minimum", 0, step = 0.1, value = 16.2),
                             numericInput("max", "Maximum", 400, step = 0.1, value = 18.3),
                             ),
  mainPanel(textOutput("range"), plotlyOutput(outputId = 'homo_plot', height = "600px"))
  ))

# building the server for the shiny app

server <-   function(input, output) {
  output$range <- renderText({ paste("You have chosen a range that goes from", input$min, "to", input$max, "\t", "This example uses the first list of the output from homol.search from the 
  <nontarget: Detecting Isotope, Adduct and Homologue Relations in LC-MS Data> from Martin Loos; 
  settings: data,
  isotopes,	elements = c(C, H, O, N, P, S), use_C = T, 
  minmz=8, 	maxmz=250,
  minrt=0.5,  maxrt=4, 
  ppm=TRUE,
  mztol=4,  rttol=0.5, 
  minlength=3, 
  mzfilter=F,
  spar=.45, 	R2=.98,
  plotit=FALSE")})
  output$homo_plot <- renderPlotly({  
    homo_filter %>% 
      filter(`mzsplit` >= input$min & `mzsplit` <= input$max) %>%
      
      ggplot(aes(x = mz, y = RT, text = paste(
        "m/z increment: ", `m.z.increment`, "\n",
        "HS ID ", `HS.IDs`, "\n",
        sep = ""))) +
      geom_point()+
      geom_line(aes(x=mz, y=RT, group = hsid, alpha = 0.4))
  })
}


# run the shiny app
shinyApp(ui = ui, server = server)
