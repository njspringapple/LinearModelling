############################# (a) ##############################################

simulation <- function(n = 100, n.x = 10, n.p = 3, n.extra = 0, beta0 = 100,
                       beta = c(3, -0.5, 0.05), sigma.eps = 1, sigma.extra = 2) {
  # Erzeugung der Designmatrix:
  X1 <- replicate(n = n.p, expr = rnorm(n, 0, 1))
  X2 <- replicate(n = n.x - (n.p + n.extra), expr = rnorm(n, 0, 1))
  X3 <- if (n.extra != 0) {
    replicate(n = n.extra, expr = rnorm(n, 0, sigma.extra))
  } else {
    NULL
  }
  X <- cbind(X1, X2, X3)
  colnames(X) <- paste0("x", seq_len(n.x))
  # Erzeugung der Zielvariablen über lineares Regressionsmodell:
  beta <- c(beta, rep(0, n.x - length(beta)))
  y <- beta0 + X %*% beta + rnorm(n, 0, sigma.eps)
  df <- cbind.data.frame(y = y, X)
  return(df)
}

Nsim <- 1000
set.seed(1)
load("../Daten/simulation_data.Rdata")
head(df_list[[1]])


############################# (b) ##############################################

# Bestimmung des AIC für eine Formel
calculate_aic <- function(formula, data) {
  model <- lm(formula = formula, data = data)
  aic <- AIC(model)
  attr(aic, "form") <- formula
  return(aic)
}

# Auswahl des Modells mit dem besten AIC
getbest <- function(aic_list) {
  ind1 <- sapply(X = aic_list, FUN = which.min)
  best_sub <- lapply(
    X = seq_along(ind1),
    FUN = function(z) {
      aic_list[[z]][[ind1[z]]]
    }
  )
  best_ind <- which.min(unlist(best_sub))
  return(best_sub[[best_ind]])
}

# Durchführung von Best Subset-Variablenselektion
best_simple <- function(data) {
  x_vars <- grep(pattern = "^x.*", x = colnames(data), value = TRUE)
  all_forms <- lapply(X = seq_along(x_vars), FUN = function(z) {
    rhs <- combn(x = x_vars, m = z, FUN = paste0, collapse = "+")
    lapply(X = rhs, FUN = function(i) {
      as.formula(paste0("y~", i))
    })
  })
  res_all <- lapply(X = all_forms, FUN = function(z) {
    lapply(X = z, FUN = calculate_aic, data = data)
  })
  model <- lm(
    formula = attr(x = get_best(res_all), which = "form"),
    data = data
  )
  return(model)
}

# Anwendung der Variablenselektionsverfahren
selection_stats <- function(data) {
  lwr <- y ~ 1
  upr <- paste0("y ~", paste0("x", seq_len(ncol(data) - 1), collapse = "+"))
  m_lower <- lm(formula = lwr, data = data)
  m_upper <- lm(formula = upr, data = data)
  # Durchführung der Variablenselektionsverfahren:
  model_list <- list(
    forward = MASS::stepAIC(
      object = m_lower,
      scope = list(upper = upr, lower = lwr),
      direction = "forward"
    ),
    backward = MASS::stepAIC(
      object = m_upper,
      scope = list(upper = upr, lower = lwr),
      direction = "backward"
    ),
    step_forward = MASS::stepAIC(
      object = m_lower,
      scope = list(upper = upr, lower = lwr),
      direction = "both"
    ),
    step_backward = MASS::stepAIC(
      object = m_upper,
      scope = list(upper = upr, lower = lwr),
      direction = "both"
    ),
    best_simple = best_simple(data)
  )
  # Zusammenfassung aller Ergebnisse in einem data.frame:
  result <- reshape2::melt(data = lapply(X = model_list, FUN = function(z) {
    do.call(cbind.data.frame, as.list(coef(z)))
  }))
  return(result)
}

load("../Daten/simulation_stats.Rdata")
head(stats_list[[1]])

results <- do.call(rbind, stats_list) %>%
  rename(Verfahren = L1) %>%
  mutate(variable = factor(
    x = variable,
    levels = c("(Intercept)", paste0("x", 1:10))
  ))
head(results)

ggplot(
  data = filter(results, variable != "(Intercept)"),
  mapping = aes(x = variable, y = value)
) +
  geom_boxplot() +
  facet_wrap(~Verfahren) +
  theme_minimal()

xtabs(formula = ~ Verfahren + variable, data = results)

correct_models <- function(results) {
  results <- lapply(X = split(x = results, f = results$ID), FUN = function(z) {
    lapply(X = split(x = z, f = z$Verfahren), FUN = function(i) {
      !any(is.na(match(paste0("x", 1:3), i$variable))) &
        all(is.na(match(paste0("x", 4:10), i$variable)))
    })
  })
  return(results)
}

overview <- reshape2::melt(correct_models(results))
xtabs(formula = ~ value + L2, data = overview)


############################# (c) ##############################################

rmse <- function(observed, predicted) {
  return(sqrt(sum((observed - predicted)^2) / length(observed)))
}

get_folds <- function(n, K = 10) {
  cv_folds <- cvTools::cvFolds(n = n, K = K)
  cv_folds <- split(x = cv_folds$subsets[, 1], f = cv_folds$which)
  return(cv_folds)
}

do_cv <- function(formula, data, folds) {
  # Berechnung des Mittelwertes der RMSEs aller Teildatensätze:
  cv_mod <- mean(vapply(X = folds, FUN = function(i) {
    model <- lm(formula = formula, data = data[-i, ])
    rmse(
      observed = data[i, "y"],
      predicted = predict(object = model, newdata = data[i, ])
    )
  }, FUN.VALUE = numeric(1)))
  # Ausgabe des Kreuzvalidierungs-RMSEs:
  attr(x = cv_mod, which = "formula") <- formula
  return(cv_mod)
}

do_cvs <- function(data, folds) {
  # Generierung aller möglichen Modellspezifikationen:
  x_vars <- grep(pattern = "^x.*", x = colnames(data), value = TRUE)
  all_forms <- lapply(X = seq_along(x_vars), FUN = function(z) {
    rhs <- combn(x = x_vars, m = z, FUN = paste0, collapse = "+")
    lapply(X = rhs, FUN = function(i) {
      as.formula(paste0("y~", i))
    })
  })
  # Durchführung der Kreuzvalidierung für die einzelnen Modellspezifikationen:
  res_all <- lapply(X = all_forms, FUN = function(z) {
    lapply(X = z, FUN = do_cv, data = data, folds = folds)
  })
  # Extraktion des gemäß des RMSEs besten Modells:
  model_best <- lm(
    formula = attr(x = get_best(res_all), which = "formula"),
    data = data
  )
  return(coef(model_best))
}

set.seed(1)
folds <- get_folds(n = 100)
set.seed(1)

load("../Daten/simulation_cv.RData")

results_cv <- reshape2::melt(cv, id.vars = "ID") %>%
  mutate(variable = factor(x = variable, levels = paste0("x", 1:10)))
head(results_cv)

ggplot(
  data = filter(results_cv, variable != "(Intercept)"),
  mapping = aes(x = variable, y = value)
) +
  geom_boxplot() +
  theme_minimal()

xtabs(formula = ~variable, data = results_cv)

correct_cv <- sapply(
  X = split(x = results_cv, f = results_cv$ID),
  FUN = function(z) {
    !any(is.na(match(x = paste0("x", 1:3), table = z$variable))) &
      all(is.na(match(x = paste0("x", 4:10), table = z$variable)))
  }
)
table(correct_cv)
