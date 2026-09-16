# Blatt 2

# Aufgabe 5

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

# b) 
# Laden der Daten und Überblick über Datensatz:
dat_dssal <- read.csv("Uebung26/Daten/ds_salaries.csv")
head(dat_dssal)
summary(dat_dssal)

# what happened to the outlier in the salary column
# dat_dssal %>% filter(salary ==  30400000) # -> currency CLP

# only interested in Germany, US, GB
dat_dssal <- dat_dssal %>% filter(company_location %in% c("DE", "GB", "US"))

# Graphische Visualisierung der Zielgröße über die Zeit je Land
# aggregierte Daten pro Jahr zur Analyse inflationsbedingter Erhöhungen
dat_dssal_yrly <- dat_dssal %>%
  group_by(work_year, company_location) %>%
  summarise(median_salary_in_usd = median(salary_in_usd, na.rm = TRUE))
# Jahresgehalt grafisch darstellen, um die inflationsbedingte Erhöhung zu analysieren
ggplot(data = dat_dssal_yrly,
       aes(x=work_year, y=median_salary_in_usd, col=company_location)) +
  geom_line(linewidth = 0.75) +
  xlab("") +
  ylab("Gehalt in USD")+
  ylim(c(0, 150000))


# Graphische Visualisierung der Zielgröÿe je Land:
# Jahresgehalt grafisch darstellen, um die inflationsbedingte Erhöhung zu analysieren
ggplot(data = dat_dssal,
       aes(x=company_location, y=salary_in_usd, fill = company_location)) +
  geom_violin() +
  geom_boxplot(width = 0.25, color = "black") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Gehalt in USD")+
  scale_y_continuous(labels=scales::comma)


# Graphische Visualisierung der log-transformierten Zielgröße je Land:
ggplot(data = dat_dssal,
       aes(x=company_location, y=log(salary_in_usd), fill = company_location)) +
  geom_violin() +
  geom_boxplot(width = 0.25, color = "black") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Logarithmiertes Gehalt")+
  scale_y_continuous(labels=scales::comma)



# c)
dat_dssal_adj <- read.csv("Uebung26/Daten/ds_salaries_adj.csv")
model <- lm(formula = salary_in_usd ~ company_location +
              company_size +
              year,
            data = dat_dssal_adj)
summary(model)

# Lineares Regressionsmodells mit log-Transformation der Daten und Inflationsbereinigung
model_adj <- lm(formula = log_salary_adj ~ company_location +
                  company_size +
                  year,
                data = dat_dssal_adj)
summary(model_adj)


# d) 
# Einschränkung des Trainingsdatensatzes auf Beobachtungen vor 2023:
dat_dssal_adj_train <- dat_dssal_adj %>%
  filter(work_year<2023)
dat_dssal_adj_test <- setdiff(dat_dssal_adj, dat_dssal_adj_train)

# Schätzung eines multiplen linearen Regressionsmodells:
model_d <- lm(formula = log_salary_adj ~ company_location +
                company_size +
                experience_level,
              data = dat_dssal_adj_train)
summary(model_d)


# Berechnung der Prognosegüte in R (transformierte Daten):
# Funktion zur Berechnung des RMSE in Abhängigkeit der wahren Werte y und der
# prädiktierten Werte y_hat:
rmse <- function(y, y_hat) {
  sqrt(sum((y - y_hat)^2) / length(y))
}
# Prädiktion der Beobachtungen im Testdatensatz:
prediction <- predict(object = model_d, newdata = dat_dssal_adj_test)
# Berechnung des RMSE:
rmse_test <- rmse(dat_dssal_adj_test$log_salary_adj, prediction)
rmse_test


# Graphischer Vergleich zwischen prädiktierten und tatsächlich beobachteten Werten auf der
# Originalskala:
ggplot()+
  geom_point(aes(x = dat_dssal_adj_test$salary_adj, y = exp(prediction))) +
  geom_line(aes(x = dat_dssal_adj_test$salary_adj,
                y = dat_dssal_adj_test$salary_adj), col = "red",
            linewidth = 1) + xlim(50000, 450000) + ylim(50000, 450000) +
  geom_text(aes(x = 300000,
                y = 50000,
                label = paste("RMSE", round(rmse_test, 2)))) +
  theme

# Im graphischen Vergleich erkennt man, dass sehr hohe Gehälter deutlich unterschätzt werden.




# e) 
# „employment_type“ umfasst mehrere Kategorien, die nicht miteinander vergleichbar sind, 
# z. B. Teilzeit vs. Vollzeit.
table(dat_dssal_adj$employment_type)

# Datensatz filtern und nur ganztags Beschäftigte aufnehmen
dat_dssal_adj_neu <- dat_dssal_adj %>%
  filter(employment_type == "FT")
model_e1 <- lm(formula = log_salary_adj ~ company_location +
                 company_size +
                 experience_level,
               data = dat_dssal_adj_neu)
summary(model_e1)

# alternativ: ins Modell aufnehmen, mit "FT" als Referenz
employment_type_factor <- factor(dat_dssal_adj$employment_type, levels = c("FT", "CT", "FL", "PT"))
model_e2 <- lm(formula = log_salary_adj ~ company_location +
                 company_size +
                 experience_level +
                 employment_type_factor,
               data = dat_dssal_adj)
summary(model_e2)
