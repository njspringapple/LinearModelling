# Tutorium 02 -------------------------------------------------------------

# Datum: 05.05.2026


# Tipp: Mit den grauen Pfeilen links können Abschnitte ein- 
#       und ausgeklappt werden.


# Aufgabe 1 ---------------------------------------------------------------

# Aufgabe 1 c) -----------------------------------------------------------
# Mach dich mit dem Datensatz mtcars vertraut.

?mtcars
# Die Hilfe zum Datensatz gibt eine Beschreibung des Datensatzes
# mit Variablenerklärungen, Quelle, etc.

# Hier kann man auch ablesen, dass die benötigen Variablen
# 'mpg' (Reichweite) und 'wt' (Gewicht) sind.

head(mtcars)
# Liefert die ersten 5 Zeilen des Datensatzes und gibt einen 
# groben Überblick wie der Datensatz aufgebaut ist, ob 
# alle Variablen richtig eingelesen sind, etc.

summary(mtcars)
# Liefert für alle im Datensatz enthaltenen Variablen eine
# kurze Zusammenfassung.
# Die Art der Zusammenfassung hängt vom Variablentyp 
# (nominal, metrisch, ...) ab.
# In diesem Fall sind alle Variablen metrisch und es wird
# je die 5-Punkte-Zusammenfassung und der Mittelwert ausgegeben.
# Die Summary gibt einen weiteren Einblick in den Datensatz und
# kann auch Hinweise auf Ausreißer geben (z.B. wenn das 
# Minimum vom Gewicht -5 wäre). Auch 'NAs' werden hier ggf.
# mit aufgeführt.

str(mtcars)
# Die str-Funktion gibt den Datentyp der Variablen wieder.
# Später im Kurs ist es z.B. wichtig, dass bestimmte Variablen
# als Faktor codiert sind. Dies kann u.a. mit der str-Funktion
# überprüft werden.

# Aufgabe 1 d) -----------------------------------------------------------
# Fitte das Regressionsmodell mit dem die Reichweite anhand 
# des Gewichts vorhergesagt werden kann und speichere das 
# Objekt in einer Variable.

mod <- lm(mpg ~ wt, data = mtcars)

# Die Funktion zum Berechnen eines linearen Modells ist 'lm'.
# Als Argument wird eine formula nach dem Schema:
# y ~ x1 + x2 + x3
# benötigt. Des Weiteren muss der Datensatz spezifiziert werden,
# wenn man nur die Variablen-Bezeichnungen verwendet.

# Das von der Funktion erzeugte Objekt kann einer Variable
# (hier: 'mod') mit dem Pfeil '<-' zugeordnet werden.

# Aufgabe 1 e) -----------------------------------------------------------
# Gib die Summary des Modells aus.

summary(mod)

# Aufgabe 1 e) i.) -------------------------------------------------------
# Erkläre alle in der Ausgabe angegebenen Werte (mit Ausnahme 
# der letzten Zeile).

?summary.lm
# Gibt genaue Informationen zu dem Output der summary-Funktion.
# Siehe z.B. Abschnitt 'Details' in der Hilfe.
# Details zur Schätzung selber findet man jedoch unter ?lm.


# 'Call':
# Hier wird die zur Berechnung verwendete Funktion nochmal
# ausgegeben.

# 'Residuals':
# Hier wird die 5-Punkte-Zusammenfassung der geschätzten 
# Residuen epsilon_dach ausgegeben.
# Bei Normalverteilungsannahme sollten diese um 0
# symmetrisch sein.
# Auch kann man anhand der Quartile sehen, wie stark die
# Residuen streuen und mit Min/Max, ob es extreme
# Beobachtungen (Ausreißer?) gibt.

# 'Coefficients':
# Alle Zahlen in der ersten Zeile beziehen sich auf
# den Intercept (also beta_dach_0).
# Alle Zahlen in der zweiten Zeile beziehen sich auf 
# auf die Kovariable 'wt' (also beta_dach_gewicht).

# 'Estimate' sind die Schätzung der Regressions-Koeffizienten
# beta_dach.

# 'Std. Error' ist der geschätzte Standardfehler der 
# Regressions-Koeffizienten. Also sigma_dach(beta_dach).

# 't value' ist der Wert der t-Statistik zur 
# Nullhypothese beta_hut_j = 0.
# 'Pr(>|t|)' ist der zugehörige p-Wert.
# Die Sterne geben an zu welchem Signifikanz-Niveau
# die jeweilige Nullhypothese abgelehnt werden kann.
# Erklärung der Sterne ist in der Zeile darunter.

# 'Residual standard error' ist der Standardfehler
# der Residuen. Also sigma_dach - die Schätzung der 
# Varianz der epsilon_i.
# (Zusätzlich sind die entsprechenden Freiheitsgrade
# der Stanardfehler-Schätzung angegeben)

# 'Multiple R-suqared' und 'Adjusted R-squared' sind
# das R^2 bzw. adjustierte R^2 des Modells.

# 'F-statistic' und der 'p-value' sind der Wert der F-Statistik
# und der entsprechende p-Wert zur Nullhypothese: 'Ein Modell
# welches nur einen Intercept enthält, ist genauso gut'.
# (Also die Mittelwerts-Schätzung der Zielvariable als 
# Prädiktion für jede Beobachtung gibt keine schlechteren
# Prädiktion als das aufgestellte Modell. Nochmal anders 
# ausgedrückt: Nullhypothese Das Modell bringt keinen Mehrwert.)

# Aufgabe 1 e) ii.) ------------------------------------------------------
# Gib für alle Größen eine Formel zur Berechnung an 
# (mit Ausnahme von t-value und Pr($> |t|$), sowie der 
# letzten Zeile).

# 'Residuals':
# epsilon_dach = Y - Y_dach = Y - X beta_dach
# Von diesem Vektor der geschätzten n Residuen können
# die Quantile bzw. Min/Max bestimmt werden.

# 'Coefficients':
# 'Estimate'
# beta_dach = (X'X)^-1 X'Y
# (Die Einträge des beta_dach-Vektors sind die Estimates.)

# 'Std. Error':
# sigma_dach(beta_dach)^2 = sigma_dach^2 (X'X)^-1
# (Die Wurzel aus den Einträgen auf der Diagonalen)
# Das benötigte sigma_dach ist das der Residuen 
# (siehe 'Residual Standard Error').

# 'Residual Standard Error':
# sigma_dach = Wurzel(Var(epsilon_dach))
# Die Freiheitsgrade ergeben sich via n-(p+1).

# 'Multiple R-Squared':
# R^2 = SSM/SST = 1 - SSE/SST = 1 - Var(epsilon_dach)/Var(Y)

# Aufgabe 1 f) -----------------------------------------------------------
# Mach dich mit dem Modell-Objekt vertraut.

?lm
# mod$ -> Vorschläge von Objekt-Einträgen

# Aufgabe 1 f) i.) -------------------------------------------------------
# Rufe die Residuen im Modell-Objekt auf.
mod$residuals

# Aufgabe 1 f) ii.) ------------------------------------------------------
# Rufe die prädiktierten Werte im Modell-Objekt auf.
mod$fitted.values

# Aufgabe 1 f) iii.) -----------------------------------------------------
# Rufe die für das Modell verwendeten Daten über das 
# Modell-Objekt auf.
mod$model

# ZUSATZAUFGABEN
# Aufgabe 1 j) -----------------------------------------------------------
# Unter den Annahmen 2.1)-2.5) (siehe VL-Skript) ist
# beta_hut ~ N(beta, (X'X)^-1 * sigma^2)

# Das wahre beta und sigma^2 kennen wir nicht, können beide
# aber schätzen und somit auch die Verteilung von beta_hut
# schätzen.

# Aufgabe 1 j) i.) -------------------------------------------------------

# Die Schätzungen sind die Estimates:
mod$coefficients

# Aufgabe 1 j) ii.) ------------------------------------------------------

# Die Schätzungen für die Varianzen sind die QUADRIERTEN 
# Std. Error aus der Summary.

(summary(mod)$coefficients[ , "Std. Error"])^2

# Aufgabe 1 j) iii.) -----------------------------------------------------

# 95%-Konfidenz-Intervalle für normalverteilte Zufallsvariablen
# berechen sich nach:
# Erwartungswert +- 1.96 * sqrt(Varianz)
# +- 1.96 ist der Wert des Normalverteilungs-Quantils (beidseitig, 5%)

KI_unten <- mod$coefficients - 1.96 * summary(mod)$coefficients[ , "Std. Error"]
KI_oben <- mod$coefficients + 1.96* summary(mod)$coefficients[ , "Std. Error"]

round(KI_unten, 2)
round(KI_oben, 2)

# Interpretation salopp:
# Der wahre Regressionskoeffizient für das Gewicht in 1000 Pfund
# liegt zu 95% zwischen -6.44 und -4.25.

# (Eigentlich richtig: Das KI [-6.44; -4.25] überdeckt mit 
# Wahrscheinlichkeit 95% den wahren Parameter beta_Gewicht 
# - sprich in dieser Formulierung ist das KI zufällig und 
# der Koeffizient fix, was eigentlich richtig ist.)
