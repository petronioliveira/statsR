
####################################################################
# Função sombrear área curva normal                                #
# Calcula a probabilidade da área sombreada em uma curva normal    #          #
# Função criada pelo Professor Petrônio Fagundes de Oliveira Filho #
# Criada em R 4.5.1                                                #
# ##################################################################


# Criar função:

normal_intervalo_sombreado <- function(a, b, mean = 0, sd = 1) {
  library(ggplot2)
  
  # Dados da curva
  x <- seq(mean - 4*sd, mean + 4*sd, length.out = 1000)
  y <- dnorm(x, mean = mean, sd = sd)
  df <- data.frame(x = x, y = y)
  
  # Área sombreada
  x_fill <- seq(a, b, length.out = 500)
  y_fill <- dnorm(x_fill, mean = mean, sd = sd)
  df_fill <- data.frame(x = x_fill, y = y_fill)
  
  # Probabilidade
  prob <- pnorm(b, mean, sd) - pnorm(a, mean, sd)
  
  # Gráfico
  ggplot(df, aes(x = x, y = y)) +
    geom_line(color = "black", size = 1) +
    geom_area(data = df_fill, aes(x = x, y = y), fill = "steelblue", alpha = 0.4) +
    geom_vline(xintercept = c(a, b), linetype = "dashed", color = "red", size = 0.8) +
    scale_y_continuous(expand = expansion(add = c(0, 0.05))) +
    annotate("text", x = (a + b)/2, y = max(y)*0.9,
             label = paste0("P(", round(a, 3), " ≤ X ≤ ", round(b, 3), ") = ", round(prob, 4)),
             size = 6, color = "red", fontface = "plain") +
    labs(
      title = "Probabilidade da Área Sombreada",
      x = "Valor",
      y = "Densidade de Probabilidade"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 13, face = "plain"),
      axis.title.y = element_text(size = 13, face = "plain"),
      axis.text = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10)
    )
}
