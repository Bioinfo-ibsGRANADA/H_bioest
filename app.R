#| standalone: true
#| packages: [ "shiny", "bslib", "genpwr", "ggplot2", "nleqslv", "MASS", "cli", "gtable", "isoband", "lifecycle", "rlang", "S7", "scales", "vctrs", "withr", "glue", "cpp11", "farver", "labeling", "R6", "RColorBrewer", "viridisLite", "munsell" ]

library(shiny)
library(bslib)
library(genpwr)

ui <- page_sidebar(
  title = "Calculadora de Potencia Genómica",
  theme = bs_theme(bootswatch = "flatly"),
  
  sidebar = sidebar(
    h5("Selección de Método"),
    selectInput("metodo", "Tipo de estudio:",
                choices = c("GWAS (Genómica de Asociación)" = "gwas",
                            "RNA-Seq (Expresión Génica)" = "rnaseq")),
    hr(),
    h5("Parámetros"),
    
    # Parámetros GWAS
    conditionalPanel(
      condition = "input.metodo == 'gwas'",
      selectInput("tipo_calculo", "Cálculo deseado:",
                  choices = c("Tamaño Muestral (N)" = "n",
                              "Potencia Estadística" = "power")),
      selectInput("model", "Modelo del Test:",
                  choices = c("Aditivo" = "Additive",
                              "Dominante" = "Dominant",
                              "Recesivo" = "Recessive")),
      sliderInput("maf", "MAF (Frec. Alelo Menor):", 0.01, 0.5, 0.1, 0.01),
      sliderInput("or", "Odds Ratio (Efecto):", 1.05, 3.0, 1.25, 0.05),
      numericInput("alpha_gwas", "Significancia (Alpha):", value = 5e-8),
      conditionalPanel(
        condition = "input.tipo_calculo == 'n'",
        sliderInput("target_power", "Potencia deseada:", 0.5, 0.95, 0.8, 0.05)
      ),
      conditionalPanel(
        condition = "input.tipo_calculo == 'power'",
        numericInput("target_n", "N total disponible:", value = 1000, min = 1)
      )
    ),
    
    # Parámetros RNA-Seq
    conditionalPanel(
      condition = "input.metodo == 'rnaseq'",
      numericInput("depth", "Profundidad de Lectura (Depth):", value = 100, min = 1),
      helpText("Promedio de lecturas por gen."),
      numericInput("cv", "Coeficiente de Variación (CV):", value = 0.7, min = 0.01, max = 1, step = 0.1),
      helpText("Variabilidad biológica (0.4 - 0.7 es común en humanos)."),
      numericInput("effect", "Tamaño del Efecto (Fold Change):", value = 1.25, min = 1.0, step = 0.05),
      helpText("Ej: 1.25 significa un cambio del 25%."),
      numericInput("power_rna", "Potencia deseada:", value = 0.8, min = 0.1, max = 0.99, step = 0.05),
      numericInput("alpha_rna", "Nivel de Significancia (Alpha):", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
      actionButton("btn_calcular", "CALCULAR", class = "btn-primary", style = "width: 100%;")
    )
  ),
  
  navset_tab(
    nav_panel("Calculadora",
      uiOutput("panel_calculadora")
    ),
    nav_panel("Información",
      card(card_body(
        h3("Sobre esta herramienta"),
        p("Esta aplicación combina dos calculadoras de potencia estadística para estudios genómicos:"),
        tags$ul(
          tags$li(strong("GWAS:"), " Usa el paquete 'genpwr' para estudios de asociación genómica."),
          tags$li(strong("RNA-Seq:"), " Cálculo basado en distribución Binomial Negativa, equivalente a RNASeqPower.")
        )
      ))
    )
  )
)

server <- function(input, output) {
  
  # --- GWAS reactivo (se actualiza solo) ---
  resultado_gwas <- reactive({
    req(input$metodo == "gwas", input$model, input$maf, input$or, input$alpha_gwas)
    tryCatch({
      if (input$tipo_calculo == "n") {
        genpwr::ss.calc(
          OR = as.numeric(input$or), k = 1,
          MAF = as.numeric(input$maf), Alpha = as.numeric(input$alpha_gwas),
          power = as.numeric(input$target_power), Test.Model = input$model
        )
      } else {
        genpwr::power.calc(
          OR = as.numeric(input$or), k = 1,
          N = as.numeric(input$target_n),
          MAF = as.numeric(input$maf), Alpha = as.numeric(input$alpha_gwas),
          Test.Model = input$model
        )
      }
    }, error = function(e) return(NULL))
  })
  
  # --- RNA-Seq reactivo (solo al pulsar botón) ---
  resultado_rna <- eventReactive(input$btn_calcular, {
    req(input$depth, input$cv, input$effect, input$alpha_rna, input$power_rna)
    tryCatch({
      z_alpha <- qnorm(1 - input$alpha_rna / 2)
      z_beta  <- qnorm(input$power_rna)
      n <- ((z_alpha + z_beta)^2 * (1/input$depth + input$cv^2)) / (log(input$effect)^2)
      ceiling(n)
    }, error = function(e) return(NULL))
  })
  
  # --- Panel principal que cambia según el método ---
  output$panel_calculadora <- renderUI({
    
    if (input$metodo == "gwas") {
      
      res <- resultado_gwas()
      
      if (is.null(res)) {
        return(card(card_body(HTML("<b style='color:red'>Error en el cálculo. Revisa los parámetros.</b>"))))
      }
      
      get_val <- function(modelo_real) {
        fila <- which(tolower(res$True.Model) == tolower(modelo_real))
        if (length(fila) == 0) return("N/A")
        if (input$tipo_calculo == "n") {
          valor <- res$N_total[fila]
          return(paste0("<span class='badge bg-success' style='font-size:14px'>",
                        format(round(valor), big.mark=".", decimal.mark=","),
                        " individuos</span>"))
        } else {
          valor <- res$Power[fila]
          return(paste0("<span class='badge bg-primary' style='font-size:14px'>",
                        round(valor*100, 2), "% de potencia</span>"))
        }
      }
      
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header("Resultados GWAS"),
          card_body(
            helpText("Asumiendo ratio casos/controles = 1 (k=1)"),
            HTML(paste0(
              "<div style='font-size:16px; line-height:2.5;'>",
              "<p>Resultados según la arquitectura genética verdadera:</p>",
              "<b>Modelo Aditivo:</b> ", get_val("Additive"), "<br>",
              "<b>Modelo Dominante:</b> ", get_val("Dominant"), "<br>",
              "<b>Modelo Recesivo:</b> ", get_val("Recessive"),
              "</div>"
            ))
          )
        ),
        card(
          card_header("Guía rápida"),
          card_body(p("Modifica los parámetros de la izquierda para ver cómo impactan en el tamaño de muestra o la potencia."))
        )
      )
      
    } else {
      
      # Panel RNA-Seq
      val <- if (input$btn_calcular > 0) resultado_rna() else NULL
      
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header("Resultados RNA-Seq"),
          card_body(
            if (is.null(val)) {
              h5("Introduce los datos y pulsa CALCULAR.")
            } else {
              div(style = "background-color:#f7f7f9; padding:20px; border-radius:10px; border:1px solid #e1e1e8;",
                  h4("Tamaño muestral requerido:"),
                  h1(paste(val, "sujetos por grupo"), style = "color:#2c3e50; font-weight:bold;"),
                  h3(paste("Total del estudio:", val * 2, "sujetos"), style = "color:#e74c3c;")
              )
            }
          )
        ),
        card(
          card_header("Guía rápida"),
          card_body(
            tags$ul(
              tags$li(strong("n por grupo:"), " Pacientes necesarios en cada brazo del estudio."),
              tags$li(strong("Total:"), " Suma de los dos grupos."),
              tags$li("Cálculo basado en distribución Binomial Negativa (estándar para RNA-Seq).")
            )
          )
        )
      )
    }
  })
}

shinyApp(ui, server)