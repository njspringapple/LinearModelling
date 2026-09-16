# Tutorium 04 -------------------------------------------------------------

# Datum: 06.06.2020
# Autor: Helen Alber


# Aufgabe 3 a) -----------------------------------------------------------

head(PlantGrowth)
?PlantGrowth
summary(PlantGrowth)
str(PlantGrowth)

# => metrische Zielvariable
# => diskrete Einflussvariable mit 3 Gruppen á 10 Beobachtungen
# => group ist als Factor spezifiziert.

# Aufgabe 3 b) -----------------------------------------------------------

mod_1 <- lm(weight ~ group, data = PlantGrowth)
summary(mod_1)

# Aufgabe 3 b) i.) -------------------------------------------------------

# Die Summary enthält einen Intercept und Group-Treat-1- und
# Group-Treat-2-Effekt.
# Die Kontrollgruppe ist also nicht mit in der Modellgleichung.
# => Es wird die Referenz-Kodierung mit 'ctrl' als Referenz
# verwendet.
# (Als Referenz wird stets das erste Level der Einflussvariable
# verwendet.)
PlantGrowth$group
levels(PlantGrowth$group)

# Die R-Hilfe kann hier etwas verwirrend sein, da die lm-Funktion
# eine sehr grundlegende Funktion ist, die mit extrem vielen
# Gegebenheiten umgehen kann.
?lm
# -> Weiterleitung zu: model.matrix -> contrasts -> options


# Aufgabe 3 b) ii.) ------------------------------------------------------

# Interpretation der Regressionskoeffizienten:
# Intercept:      Eine Pflanze der Kontroll-Gruppe gibt im
#                 Mittel 5.03 (kg?) Ertrag.
# grouptrt1:      Eine Pflanze der Gruppe trt1 gibt im Mittel
#                 0.37 (kg?) weniger Ertrag als eine Pflanze
#                 der Kontrollgruppe.
# grouptrt2:      Eine Pflane der Gruppe trt2 gibt im Mittel
#                 0.49 (kg?) mehr Ertrag als eine Pflanze
#                 der Kontrollgruppe.

# Weitere Kommentare:
# Das Modell ist nicht sonderlich gut. Nur 26% der Streuung von
# Y können durch die Gruppen erklärt werden. Es gibt also
# an sich schon eine sehr hohe Streuung beim Ertrag bei
# selber Behandlung.
# Auch die hohen p-Werte sind ein Indikator für sehr breite
# Konfidenzintervalle. Selbst, wenn man nicht bezüglich
# multiplem Testen korrigiert, sind keine signifikanten
# Unterschiede zwischen den Gruppen festzustellen.
# Das heißt aber nicht, dass keine Ertragsunterschiede vorhanden 
# sind. Wir sind uns aufgrund der hohen Streuung nur nicht 
# sicher, ob es welche gibt, oder nicht.
# Für genauere Analysen empfiehlt sich ein größerer Stichproben-
# Umfang. 
# So kann man aus den Koeffizienten maximal Tendenzen
# ableiten, die genauer geprüft werden müssen.


# Aufgabe 3 b) iii) ------------------------------------------------------

# Prädiktion für den Erwartungswert einer neuen Beobachtung aus:
# - der Kontrollgruppe:
# y_hut = Intercept = 5.03

# - der Treatment-1-Gruppe:
# y_hut = Intercept + grouptrt1 = 5.03 - 0.37 = 4.66



# Zusätzliche Kommentare:
# - Nicht immer sind komplizierte Modelle die beste Wahl. Ein
#   einfacher Boxplot gibt auch schon viel Aufschluss und
#   zeigt die Tendenz zu höherem Ertrag mit 'trt2'.

plot(weight ~ group, data = PlantGrowth)

# - Weitere Deskriptive Analyse zeigt:
library(dplyr)
PlantGrowth %>% 
  group_by(group) %>% 
  summarise(Mean = mean(weight), 
            Variance = var(weight))
#     - Die Schätzung vom Intercept ist tatsächlich der 
#       Mittelwert der 'ctrl'-Gruppe
#     - Auch die Prädiktion für 'trt1' stimmt mit dem
#       Mittelwert von 'trt1' überein (4.66).
#     - Die Varianz-Homogenitätsannahme ist hier evtl. 
#       nicht gegeben, da z.B. 'trt1' eine höhere Varianz
#       aufweist als z.B. 'trt2'. Die Stichprobe ist 
#       allerdings recht klein und die Abweichung noch 
#       nicht überzubewerten.
#       Die Boxplots haben auch tendenziell ähnlich 
#       breite Boxen, also auch eher ähnliche Varianz.



# Aufgabe 3 c) i.) -------------------------------------------------------

# Erstellen einer Variable mit geänderter Level-Reihenfolge
PlantGrowth$group_trt2 <- relevel(PlantGrowth$group, 
                                  ref = "trt2")

# Prüfe, ob es geklappt hat:
PlantGrowth$group
PlantGrowth$group_trt2

# Fitte das Modell mit der geänderten group-Variable:
mod_trt2 <- lm(weight ~ group_trt2, data = PlantGrowth)
summary(mod_trt2)

# Nun ist 'trt2' der neue Intercept vom Modell und
# die anderen Koeffizienten geben die Unterschiede
# zu dieser Referenz-Gruppe an.


# Aufgabe 3 c) ii.) ------------------------------------------------------

# Ausgabe der ersten Summary:
summary(mod_1)

# 1.) Berechnung der Gruppenmittelwerte:
# mu_ctrl = 5.03 (Refernz-Gruppe)
# mu_trt1 = mu_ctrl + grouptrt1 = 5.03 - 0.37 = 4.66
# mu_crtl = mu_ctrl + grouptrt2 = 5.03 + 0.49 = 5.52

# 2.) Mittelwert der Gruppenmittelwerte:
# Intercept = (mu_crtl + mu_trt1 + mu_trt2)/3 = 5.07

# 3.) Effekte berechnen:
# tau_ctrl = mu_crtl - Intercept = 5.03 - 5.07 = -0.04
# tau_trt1 = mu_trt1 - Intercept = 4.66 - 5.07 = -0.41
# tau_trt2 = mu_trt2 - Intercept = 5.52 - 5.07 =  0.45

# (Das geht natürlich auch alles direkt in R, wird aber
# schnell unübersichtlich bei den ganzen Modellparameter-
# Zugriffen)

# Aufgabe 3 d) -----------------------------------------------------------

# Erstelle Datensatz:
MeineDaten <-data.frame(score = 1:9, 
                        gruppe = c(1, 1, 1, 2, 2, 2, 3, 3, 3))

# Zeige Datensatz:
MeineDaten

# Fitte Modell und gibt die Summary aus:
mod_2 <-lm(score ~ gruppe, data = MeineDaten)
summary(mod_2)

# Problem:
# Anstelle der zwei erwarteten Parameter für den Unterschied
# der TUM und LMU zur HM gibt es nur einen Parameter.
# Ursache ist, dass die Codierung der Hochschulen mit Zahlen
# von R als Zahlen interpretiert und somit als metrisch
# behandelt werden werden.

# Lösung:
# Die Variable 'gruppe' muss zu einem Factor umgewandelt werden:

MeineDaten$gruppe_fact <- factor(MeineDaten$gruppe, 
                                 levels = 1:3,
                                 labels = c("TUM", "HM", "LMU"))

# Nun kann das Modell gefittet werden:
mod_2_fact <- lm(score ~ gruppe_fact, data = MeineDaten)
summary(mod_2_fact)

# Jetzt passt alles!
# (Vergleiche auch mit dem Plot:)
plot(score ~ gruppe_fact, data = MeineDaten)
