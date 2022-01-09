## an example code to plot the variance estimations
## different line types for b/w printing
## facets
## math facet labels

# load("simu83.RData")
library(tidyverse)
library(reshape2)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

toplot <- melt(simu83$errdata1, id.vars = c("X", "Lower", "Upper", "g", "n"))

nlevels <- c(expression("n == 50"),
                expression("n == 100"))

glevels <- c(paste0("g[", 1, "]"),
             paste0("g[", 2, "]"),
             paste0("g[", 3, "]"))

toplot$n <- factor(toplot$n)
levels(toplot$n) <- nlevels


toplot$g <- factor(toplot$g)
levels(toplot$g) <- glevels

(toplot %>% mutate(Curve = variable) %>% 
  ggplot(aes(x = X, y = value, group = Curve, col = Curve)) + 
  geom_line(aes(lty = Curve), lwd = 1) + 
  geom_ribbon(aes(min = Lower, max = Upper), fill = "grey", alpha = 0.35,
              col = "grey") + 
  geom_line(aes(lty = Curve), lwd = 1) +
  geom_line(aes(lty = Curve), lwd = 1) + facet_grid(n ~ g, labeller = label_parsed) + 
  theme_bw() + scale_color_manual(values = cbPalette[2:3]) + 
  xlab("X") + ylab(expression(sigma(epsilon))) + 
  theme(text = element_text(size = 15),
        legend.position = c(0.07,0.87))) %>% ggsave(plot = .,
                                                   filename = "example.pdf",
                                                   width = 10,
                                                   height = 5)
