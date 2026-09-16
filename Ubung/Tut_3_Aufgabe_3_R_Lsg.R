# Tutorium 03 -------------------------------------------------------------

# Datum: 19.05.2026

# Aufgabe 3 ---------------------------------------------------------------

# Aufgabe 3 b) ------------------------------------------------------------

# Überblick über den Datensatz verschaffen:
?airquality

head(airquality)

summary(airquality)

str(airquality)

# Lineares Modell fitten:
mod <- lm(Ozone ~ Solar.R + Temp + Wind, data = airquality)

# Summary aufrufen
summary(mod)

# Aufgabe 3 d) ------------------------------------------------------------

# Erstelle die Matrix A
A <- matrix(c(0, 1, 3, 0,
              0, 0, 0, 1), 
            nrow = 2, 
            byrow = T)

# Erstelle den Vektor c
c <- c(8, 0)

# Teste die Hypothese.
library(car)
linearHypothesis(mod, hypothesis.matrix = A, rhs = c)

# Der Wert der F-Statistik zu der Hypothese ist 14.388,
# der p-Wert 2.9*10^-6. Da dieser kleiner als 0.05 ist, 
# kann die Nullhypothese abgelehnt werden.


# ZUSATZAUFGABEN
# Aufgabe 3 e) ------------------------------------------------------------

#### Overall-Test ####

# Zunächst muss der Datensatz gefiltert werden, sodass
# für das Gesamt-Modell und das restringierte Modell
# der gleiche Datensatz verwendet wird.
# Hintergrund ist, dass die lm-Funktion automatisch
# alle Beobachtungen rauswirft, die ein NA haben.
# Das tut sie allerdings nur für Variablen, die
# auch im Modell vorkommen.
# Sprich im Modell unter H0 werden nur Beobachtungen
# entfernt, die bei Ozone ein NA haben.
# Für das komplette Modell werden aber auch Beobachtungen
# entfernt, die z.B. auch bei Solar.R ein NA haben.
# => Die Modelle werden mit unterschiedlichen Daten gefittet 
# und können daher nicht mehr mit einem F-Overall-Test 
# verglichen werden.

# Die summary() zeigt, dass außer Ozone nur noch die 
# Solar.R-Variable NAs hat.
summary(airquality)

# Diese filtern wir nun auch für die Schätzung des 
# restringierten Modells manuell heraus.
air_filter <- airquality[!is.na(airquality$Solar.R), ]
# (Zugriff auf die Datensatz: MeinDatensatz[Zeilen, Spalten]
# hier: Zugriff auf die Zeilen, die nicht ('!') NA sind und
# alle Spalten (hinter dem Komma keine Einschränkung))


# Modell unter Nullhypothese (nur ein Intercept):
mod_H0_Overall <- lm(Ozone ~ 1, data = air_filter)
# (Hinweis: Für Modelle mit Kovariablen braucht der Intercept
# nicht explizit aufgeführt werden und wird automatisch 
# hinzugefügt.)

# Wald-Test (Overall-Test):
anova(mod, mod_H0_Overall)

# Auch die anova gibt einen F-Wert von 54.83 aus und einen
# zugehörigen p-Wert von 2.2*10^-16. Vergleiche:
summary(mod) # ganz unten

#### t-Test ####

# (Da die Solar.R-Variable in beiden Modellen vorkommt, muss
# hier kein gesondert gefilterter Datensatz verwendet werden.)

# Modell unter Nullhypothese (ohne Temp)
mod_H0_tTest_Temp <- lm(Ozone ~ Solar.R + Wind, data = airquality)

anova(mod, mod_H0_tTest_Temp)

# Der F-Wert der anova ist 42.463 und der zugehörige 
# p-Wert 2.42*10^-9.
# Vergleich mit der Summary:
summary(mod)
# Hier wird ein t-Wert von 6.516 angegeben. Der p-Wert ist
# allerdings identisch. 
# t-Test und F-Test führen zu äquivalenten Ergebnissen.
# Quadrieren der t-Statistik ergibt den F-Statistik-Wert.
# Beide Teststatistiken können in die Wald-Statistik 
# überführt werden. 