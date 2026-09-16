
# Tutorium 06 -------------------------------------------------------------

# Datum: 04.07.2023
# Autor: Helen Alber


# Aufgabe 17 ---------------------------------------------------------------

# Lade notwendige Packages:
library(ggplot2)


# Modellerstellung --------------------------------------------------------

mod_air <- lm(Ozone ~ Temp + Wind + Solar.R, data = airquality)
# Das Modell mod_air enthält keine Interaktion oder Kovariablen die beispielsweise 
# mit Splines flexibel aufgenommen werden. Vorherige Tutorien haben Indizien 
# geliefert, dass dies evtl. notwendig wäre. Wir sollten also Verletzungen bei 
# den Modellannahmen feststellen. 


# Datenaufbereitung Modelldiagnostik --------------------------------------

# Der Datensatz dat_air enthält verschiedene Residuen-Arten, die gefitteten 
# Werte sowie die Kovariablen aus dem Modell.

# Datensatzerstellung für Modelldiagnstik-Plots von mod_air:
dat_air <- data.frame(
  # Verschiedene Residuen
  resid = residuals(mod_air),
  resid_stand = rstandard(mod_air),
  resid_stand_abs = abs(rstandard(mod_air)),
  
  # Gefittete Werte
  fitted = mod_air$fitted.values,
  
  # Datengrundlage
  mod_air$model
)


# Plot-Funktion -----------------------------------------------------------

# Die Funktion plot_resid nimmt einen Datensatz wie dat_air als Input und 
# enthält einige weitere Funktionsargumente, mit denen der Output modifiziert 
# werden kann. Der Output ist ein Plot (üblicherweise Residuen vs. Fitted 
# oder Residuen vs. Kovariable).

# Input: Datensatz mit Variablen zur Modelldiagnostik
# Output: Scatterplot
plot_resid <- function(data,          # Datensatz 
                       x = "fitted",  # x-Variable (als Charcater)
                       y = "resid",   # y-Variable (als Character)
                       ylim = NULL,   # y-Achsenlimits (als Vektor: c(unten, oben))
                       line = FALSE,  # Smoothline durch die Daten (als Boolean)
                       shape = 16,    # Punkt-Form (als Integer 1-24)
                       alpha = 1) {   # Punkt-Transparenz (als Numeric 0-1)
  
  
  # Erstelle Basis-Plot:
  p <- ggplot(data = data, mapping = aes_string(x = x, y = y)) +
    theme_bw() + # Theme verändern
    geom_point(shape = shape, # Punkte hinzufügen
               alpha = alpha) +
    geom_hline(yintercept = 0) + # Line bei y = 0 hinzufügen
    ggtitle(paste(y, "vs.", x)) + # Titel hinzufügen
    theme(plot.title = element_text(hjust = 0.5)) # Titel zentrieren
  
  
  # Falls 'line' true ist, dann füge eine Smoothline hinzu:
  if(line == TRUE) {
    p <- p + 
      geom_smooth(size = 0.4, 
                  alpha = 0.1)
  }
  
  
  # Falls y-Limits gesetzt werden, nutze diese:
  if(!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }
  
  
  # Gebe finalen Plot aus:
  return(p)
  
}



# Modelldiagnostik: 'mod_air' ---------------------------------------------

# Diagnostik Teil A: Residuen vs. Fitted/Kovariable
plot_resid(dat_air)

plot_resid(data = dat_air, 
           x = "Temp")

plot_resid(data = dat_air, 
           x = "Wind")

plot_resid(data = dat_air, 
           x = "Solar.R")



# Diagnostik Teil B: abs. stand. Residuen vs. Fitted/Kovariable
plot_resid(dat_air, 
           y = "resid_stand_abs")

plot_resid(data = dat_air, 
           x = "Temp",
           y = "resid_stand_abs")

plot_resid(data = dat_air, 
           x = "Wind",
           y = "resid_stand_abs")

plot_resid(data = dat_air, 
           x = "Solar.R",
           y = "resid_stand_abs")



# Diagnostik Teil C: Residuen vs. Index & Durbin-Watson-Test
library(dplyr)
library(car)

# Erstelle Datensatz mit Index
dat_air_mit_index <- data.frame(
  index = as.numeric(rownames(airquality)), 
  airquality) %>% 
  filter(!is.na(Ozone), !is.na(Wind), !is.na(Temp), 
         !is.na(Solar.R)) %>% 
  mutate(resid = mod_air$residuals,
         fitted = mod_air$fitted.values)

# Plotte Residuals vs. Index
plot_resid(data = dat_air_mit_index,
           x = "index",
           line = T)


# Durbin-Watson-Test
durbinWatsonTest(mod_air)



# Diagnostik Teil D: QQ-Plot
qqnorm(dat_air$resid_stand)
qqline(dat_air$resid_stand)


# Alternative Erstellung via ggplot:
dat_air_qqplot <- data.frame(
  resid_stand_emp = dat_air$resid_stand,
  resid_stand_theo = qnorm(ecdf(dat_air$resid_stand)(dat_air$resid_stand),
                           mean = mean(dat_air$resid_stand),
                           sd = sd(dat_air$resid_stand)))
  
ggplot(dat_air_qqplot, aes(x = resid_stand_theo,
                           y = resid_stand_emp)) +
  theme_bw() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle("Normal-QQ-Plot")


# Diagnostik Teil E: Varianz-Inflations-Faktor
library(car)
vif(mod_air)


