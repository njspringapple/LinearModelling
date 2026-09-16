# Blatt 2

# Aufgabe 4
# Laden von benötigten Paketen:
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(visreg)
library(effects)

# Definition eines ggplot-Themes:
theme <- theme_classic() +
  theme(text = element_text(size = 12), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text.y = element_text(size = 12), legend.text.align = 0,
        strip.placement = "outside", strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))

# a) 
# Einlesen der Daten und Überblick 
load("Uebung26/Daten/wdi.Rdata")
head(wdi)
str(wdi)
summary(wdi)

# Filtern für das Jahr 2013
wdi <- wdi %>% filter(year == 2013)

# Graphische Visualisierung der Verteilungen von Bevölkerung (pro Tausend Einwohner) und CO2-Emissionen:
wdi$Population1000 <- wdi$Population / 1000
wdi_long <- gather(wdi, variable, value, c(CO2emission, Population1000))
gg_wdi <- ggplot(data = wdi_long, mapping = aes(x = value)) +
  geom_density() + facet_wrap(~variable, scale = "free") +
  xlab("") + theme
grid.arrange(gg_wdi + ggtitle("Uneingeschränkter Datenbereich"),
             gg_wdi + ggtitle("Eingeschränkter Datenbereich") +
               xlim(0, 1e5))

# Bevölkerung und CO2 Emission extrem rechtsschief bzw linkssteil (dh. ein Großteil 
# der Länder weist eine vergleichsweise kleine Bevölkerungszahl sowie vergleichsweise 
# niedrige CO2-Emissionen auf.

# Durchführung einer log10 Transformation
wdi <- wdi %>%
  mutate(log10_Population = log10(Population),
         log10_CO2emission = log10(CO2emission))
head(wdi)

# Graphische Visualisierung der Verteilungen von Bevölkerung und CO2-Emissionen:

# Originalvariablen:
gg_streu <- ggplot(data = wdi,
                   mapping = aes(x = Population * 1000, y = CO2emission,
                                 label = country)) +
  geom_point() + xlab("Bevölkerung") + ylab("Co2-Emissionen [kt]") +
  geom_text(mapping = aes(label = ifelse(test = CO2emission >=
                                           sort(CO2emission,
                                                decreasing = TRUE)[3] |
                                           Population >=
                                           sort(Population,
                                                decreasing = TRUE)[3],
                                         yes = as.character(country), no = "")),
            hjust = 1, vjust = 1, size = 3) + theme
# Transformierte Variablen:
gg_streu_log10 <- gg_streu + wdi +
  aes(x = log10_Population, y = log10_CO2emission) +
  xlab("log10(Bevölkerung)") + ylab("log10(CO2-Emissionen) [kt]") +
  geom_text(mapping = aes(label = ifelse(test = log10_CO2emission >=
                                           sort(log10_CO2emission,
                                                decreasing = TRUE)[3] |
                                           log10_Population >=
                                           sort(log10_Population,
                                                decreasing=TRUE)[3],
                                         yes = as.character(country), no = "")),
            hjust = 1, vjust = 1, size = 3)
# Gemeinsame Visualisierung:
grid.arrange(gg_streu + ggtitle("Originaldaten"),
             gg_streu_log10 + ggtitle("Transformierte Daten"), nrow = 1)

# Ein (annähernd) linearer Zusammenhang ist nur auf der Ebene der log-transformierten Daten gegeben.




# b)
# Lineares Regressionsmodell ohne Transformation der Daten 
lm_original <- lm(formula = CO2emission ~ Population, data = wdi)
summary(lm_original)

# Lineares Regressionsmodell mir log10 Transformierten Daten 
lm_log10 <- lm(formula = log10_CO2emission ~ log10_Population, data = wdi)
summary(lm_log10)

# Graphische Visualisierung der beiden Regressionsmodelle
grid.arrange(gg_streu %+% wdi + geom_smooth(method = "lm", se = FALSE) +
               ggtitle("Originaldaten") ,
             gg_streu_log10 %+% wdi +
               geom_smooth(method = "lm", se = FALSE) +
               ggtitle("Transformierte Daten"),
             nrow = 1)

# Beim Modell ohne Transformation fällt auf, dass die Varianz der Residuen nicht konstant
# ist. Dies stellt eine Verletzung der Annahme der Varianzhomogenität dar.



# c) 
# Graphische Visualisierung der Dichtekurven der einzelnen Variablen
ggplot(data = gather(wdi, variable, value, GDP:PopulationRural),
       mapping = aes(x = value)) +
  geom_density(aes()) +
  facet_wrap(facets = ~variable, scales = "free") + theme


# Transformation einzelner Variablen 
# Bemerkung: log(x+1) bei rechtsschiefen Variablen, die ebenfalls den Wert 0 enthalten
wdi <- wdi %>%
  mutate(Area = log10(Area), CO2emission = log10(CO2emission),
         Population = log10(Population),
         PopulationRural = log10(PopulationRural + 1), GDP = log10(GDP)) 

# Graphische Visualisierung der Dichtekurven nach durchgeführten Transformationen:
ggplot(data = gather(wdi, variable, value, GDP:PopulationRural),
       mapping = aes(x = value)) +
  geom_density(aes()) +
  facet_wrap(facets = ~variable, scales = "free") + theme


# Graphische Visualisierung des Zusammenhangs der Variablen Area, Population,
# Livestock und Employees.M.IND mit den CO2-Emissionen der Länder in Streudiagrammen:
ggplot(data = gather(wdi, variable, value, Area, Population, Livestock, GDP),
         mapping = aes(y = CO2emission, x=value)) + geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  facet_wrap(facets = ~variable, scale = "free_x") + theme

# Für die Fläche, die Einwohnerzahl und das GDP besteht ein positiver Zusammenhang
# mit den CO2-Emissionen. Zwischen dem Livestock-Index und den CO2-Emissionen ist
# kein ausgeprägter Zusammenhang erkennbar.


# d) 
# Schätzung eines multiplen linearen Regressionsmodells unter Verwendung der Einussgröÿen
# Area, Population, Livestock und Employees.M.IND:
lm_co2_mult <- lm(formula = CO2emission ~ Area + Population + Livestock +
                      Employees.M.IND, data = wdi, na.action = na.exclude)
summary(lm_co2_mult)

# e) 
# Konditionelle Darstellung: 
par(mfrow = c(2, 2))
visreg(fit = lm_co2_mult, type = "conditional",
       ylim = range(wdi$CO2emission, na.rm = TRUE))


# Kontrastdarstellung 
par(mfrow = c(2, 2))
visreg(fit = lm_co2_mult, type = "contrast", ylim = c(-2, 2))

# Alternative: effects
plot(allEffects(lm_co2_mult), ylim = range(wdi$CO2emission, na.rm = TRUE))


# f) 
# Graphische Visualisierung des univariaten Zusammenhangs zwischen der Landbevölkerung
# und den CO2-Emissionen:
ggplot(data = wdi, mapping = aes(x = PopulationRural, y = CO2emission)) +
  geom_point() + geom_smooth(method = "lm") + theme

# Es ist ein klar positiver Zusammenhang zwischen der Landbevölkerung und den CO2-
# Emissionen erkennbar.


# Schätzung eines multiplen linearen Regressionsmodells mit zusätzlicher Einflussgröße PopulationRural:
lm_co2_rur <- lm(formula = CO2emission ~ Area + Population + Livestock +
                   Employees.M.IND + PopulationRural,
                 data = wdi, na.action = na.exclude)
summary(lm_co2_rur)

# Graphische Visualiserung der individuellen Effekte:
par(mfrow = c(2, 3))
visreg(fit = lm_co2_rur, type = "contrast", ylim = c(-2, 2))


# g) 
# Schätzung eines multiplen linearen Regressionsmodells mit zusätzlicher Einussgröÿe GDP:
lm_co2_gdp <- lm(formula = CO2emission ~ Area + Population + Livestock +
                   Employees.M.IND + PopulationRural + GDP,
                 data = wdi, na.action = na.exclude)
summary(lm_co2_gdp)


# Graphische Visualisierung der individuellen Eekte:
par(mfrow = c(2, 3))
visreg(fit = lm_co2_gdp, type = "contrast", ylim = c(-2, 2))
