#| standalone: true
#| packages: [ "shiny", "shinydashboard", "genpwr", "munsell", "scales", "colorspace" ]

library(shiny)
library(genpwr)
library(shinydashboard)

# Interfaz de Usuario (UI)
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(title = "Calculadora GWAS"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculadora", tabName = "calc", icon = icon("calculator")),
      menuItem("Información", tabName = "info", icon = icon("info-circle"))
    ),
    
    hr(),
    h4(" Parámetros", style = "margin-left: 15px;"),
    
    selectInput("tipo_calculo", "Cálculo deseado:",
                choices = c("Tamaño Muestral (N)" = "n", 
                            "Potencia Estadística" = "power")),
    
    selectInput("model", "Modelo del Test:",
                choices = c("Aditivo" = "Additive", 
                            "Dominante" = "Dominant", 
                            "Recesivo" = "Recessive")),
    
    sliderInput("maf", "MAF (Frec. Alelo Menor):", 0.01, 0.5, 0.1, 0.01),
    sliderInput("or", "Odds Ratio (Efecto):", 1.05, 3.0, 1.25, 0.05),
    numericInput("alpha", "Significancia (Alpha):", value = 5e-8),
    
    conditionalPanel(
      condition = "input.tipo_calculo == 'n'",
      sliderInput("target_power", "Potencia deseada:", 0.5, 0.95, 0.8, 0.05)
    ),
    
    conditionalPanel(
      condition = "input.tipo_calculo == 'power'",
      numericInput("target_n", "N total disponible:", value = 1000, min = 1)
    )
  ),
  
  dashboardBody(
    tabItems(
      # Pestaña Principal
      tabItem(tabName = "calc",
              fluidRow(
                # Caja de Resultados - Usamos shinydashboard::box para evitar conflictos
                shinydashboard::box(
                  title = "Resultados del Análisis", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 8,
                  helpText("Asumiendo ratio casos/controles = 1 (k=1)"),
                  uiOutput("resultadoText")
                ),
                
                # Caja de ayuda rápida
                shinydashboard::box(
                  title = "Guía rápida", 
                  status = "warning", 
                  width = 4,
                  "Modifica los parámetros de la izquierda para ver cómo impactan en el tamaño de muestra o la potencia."
                )
              )
      ),
      
      # Pestaña de Información
      tabItem(tabName = "info",
              h2("Sobre esta herramienta"),
              p("Esta calculadora utiliza el paquete 'genpwr' para estimar la potencia y el tamaño de muestra en estudios de asociación genómica (GWAS).")
      )
    )
  )
)

# Servidor
server <- function(input, output) {
  
  output$resultadoText <- renderUI({
    req(input$model, input$maf, input$or, input$alpha)
    
    res <- tryCatch({
      if (input$tipo_calculo == "n") {
        genpwr::ss.calc(
          OR = as.numeric(input$or), k = 1,
          MAF = as.numeric(input$maf), Alpha = as.numeric(input$alpha),
          power = as.numeric(input$target_power), Test.Model = input$model
        )
      } else {
        genpwr::power.calc(
          OR = as.numeric(input$or), k = 1,
          N = as.numeric(input$target_n),
          MAF = as.numeric(input$maf), Alpha = as.numeric(input$alpha),
          Test.Model = input$model
        )
      }
    }, error = function(e) return(NULL))
    
    if(is.null(res)) return(HTML("<b style='color:red'>Error en el cálculo. Revisa los parámetros.</b>"))
    
    get_val <- function(modelo_real) {
      fila <- which(tolower(res$True.Model) == tolower(modelo_real))
      if(length(fila) == 0) return("N/A")
      
      if (input$tipo_calculo == "n") {
        valor <- res$N_total[fila]
        return(paste0("<span class='label label-success' style='font-size:14px'>", format(round(valor), big.mark = ".", decimal.mark = ","), " individuos</span>"))
      } else {
        valor <- res$Power[fila]
        return(paste0("<span class='label label-info' style='font-size:14px'>", round(valor * 100, 2), "% de potencia</span>"))
      }
    }
    
    HTML(paste0(
      "<div style='font-size: 16px; line-height: 2.5;'>",
      "<p>Resultados según la arquitectura genética verdadera:</p>",
      "<b>Modelo Aditivo:</b> ", get_val("Additive"), "<br>",
      "<b>Modelo Dominante:</b> ", get_val("Dominant"), "<br>",
      "<b>Modelo Recesivo:</b> ", get_val("Recessive"),
      "</div>"
    ))
  })
}

shinyApp(ui, server)