# List of required packages
required_packages <- c("shiny", "ggplot2", "MASS", "ggpubr", "dplyr", "shinyBS")

# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load packages
lapply(required_packages, library, character.only = TRUE)

library(shiny)
library(ggplot2)
library(MASS)
library(ggpubr)
library(dplyr)
library(shinyBS)

ui <- fluidPage(
  titlePanel("Simulation Study: Linear Model Assumptions"),
  
  withMathJax(),
  
  h3("Purpose and Instructions"),
  
  helpText(
    strong("This tool illustrates how violations of classical linear model assumptions affect the behavior of OLS estimators."),
    "It uses simulation to show how the sampling distribution, confidence intervals, and coverage rates behave under different model scenarios.",
    br(), br(),
    
    strong("The linear regression model assumes:"),
    tags$ol(
      tags$li("Linearity: \\( Y_i = \\beta_0 + \\beta_1 x_{1i} + \\beta_2 x_{2i} + \\varepsilon_i \\)"),
      tags$li("Zero-mean errors: \\( \\mathbb{E}[\\varepsilon_i] = 0 \\)"),
      tags$li("Homoscedasticity: \\( \\operatorname{Var}(\\varepsilon_i) = \\sigma^2 \\)"),
      tags$li("Independence: \\( \\varepsilon_1, \\ldots, \\varepsilon_n \\text{ are independent} \\)"),
      tags$li("Normality: \\( \\varepsilon_i \\sim \\mathrm{N}(0, \\sigma^2) \\)")
    ),
    
    # br(),
    # strong("How to use it:"),
    # tags$ul(
    #   tags$li("Choose a scenario representing a specific assumption violation."),
    #   tags$li("Set the sample size, number of repetitions, and confidence level."),
    #   tags$li("Set the true values for the regression coefficients and correlation between covariates."),
    #   tags$li("In Scenario F, toggle inclusion of a hidden confounder to test misspecification."),
    #   tags$li("Interpret the results in the table and plots to evaluate bias, spread, and coverage."),
    #   tags$li("Read the Scenario-Specific Interpretation section for guidance.")
    # ),
    br(),
    "This tool helps develop intuition about how regression inference depends on these assumptions — particularly in finite samples."
  ),
  
  
  
  sidebarLayout(
    sidebarPanel(
      actionButton("runSim", "Run Simulation",
                   style = "color: white; background-color: #007BFF; border-color: #007BFF; font-weight: bold; font-size: 16px; padding: 10px 20px;"),
      bsTooltip("runSim", "Click to start the simulation", "right"),
      
      h4("Model Parameters"),
      
      numericInput("beta0", "Intercept (β₀):", value = -1),
      bsTooltip("beta0", "The true intercept of the linear model", "right"),
      
      numericInput("beta1", "Slope for x₁ (β₁):", value = 2),
      bsTooltip("beta1", "Coefficient for the first covariate x₁", "right"),
      
      numericInput("beta2", "Slope for x₂ (β₂):", value = 0.5),
      bsTooltip("beta2", "Coefficient for the second covariate x₂", "right"),
      
      numericInput("rho", "Correlation between x₁ and x₂:", value = 0.3, min = -0.99, max = 0.99),
      bsTooltip("rho", "Controls how strongly x₁ and x₂ are linearly related", "right"),
      
      h4("Simulation Parameters"),
      
      numericInput("conf_level", "Confidence level:", value = 0.95, min = 0.5, max = 0.99),
      bsTooltip("conf_level", "Confidence level for confidence intervals", "right"),
      
      numericInput("n", "Sample Size (n):", value = 200, min = 10),
      bsTooltip("n", "Number of observations in each simulation", "right"),
      
      numericInput("reps", "Number of Repetitions:", value = 1000, min = 100),
      bsTooltip("reps", "How many simulation runs to perform", "right"),
      
      selectInput("scenario", "Choose Scenario:",
                  choices = c(
                    "A: All assumptions are met" = "A",
                    "B: Gamma errors (non-normal)" = "B",
                    "C: Autocorrelated errors\n(dependent errors)" = "C",
                    "D: Heteroscedastic errors\n(depends on x₁)" = "D",
                    "E: Systematic errors depending on x₁ (E[ε|x₁] ≠ 0)" = "E",
                    "F: Model misspecification (true relationship is not Y = β₀ + β₁·x₁ + β₂·x₂)" = "F"
                  )),
      bsTooltip("scenario", "Select the type of assumption violation", "right"),
      conditionalPanel(
        condition = "input.scenario == 'F'",
        checkboxInput("include_confounder", "Include confounder in true model (Scenario F only)", value = TRUE)
      ),
      bsTooltip("include_confounder", "If checked, the true model includes an unobserved confounder term.", "right")
      
    ),
    
    mainPanel(
      tableOutput("summary_table"),
      h4("Mean Estimate Explanation"),
      helpText("The 'Mean_Estimate' column shows the average estimated coefficient across all simulation runs. 
      If the model is unbiased, this average should be close to the true set beta value. 
      Strong deviations from the true value suggest estimation bias."),
      
      h4("Coverage Rate Explanation"),
      helpText(
        "The 'Coverage_rate' column shows the proportion of simulations in which the computed confidence interval for a beta coefficient 
   contains the true beta value. It is calculated as the number of times the true value lies within the interval 
   divided by the total number of simulation repetitions.",
        
        "If model assumptions hold and the estimator is valid, the empirical coverage rate should approximately match the 
   chosen confidence level (e.g., 0.95 for a 95% interval).",
        
        "A coverage rate noticeably below the confidence level suggests that the intervals are too narrow or systematically 
   miss the true value."
      ),
      
      h4("Sampling Distributions of Coefficients"),
      plotOutput("coef_plot"),
      h4("Sampling Distributions of Coefficients"),
      helpText(
        "Each plot shows the sampling distribution of an estimated regression coefficient (β₀, β₁, β₂) across all simulation repetitions.\n",
        
        "The histogram represents the frequency of estimated values. The blue density curve shows the empirical distribution 
   of estimates, while the dashed red line shows the ideal normal distribution based on the empirical mean and standard deviation.\n",
        
        "The vertical dotted line marks the true beta value used to generate the data. If the estimator is unbiased, 
        the histogram should be centered near the true value (dotted line).",
        
        "The sampling distribution of the estimator is often approximately normal, especially when the sample size is large, even if the error distribution is not normal. 
      Deviations from normality tend to appear only when the normality assumption is violated **and** the Central Limit Theorem doesn't apply reliably (e.g. too small sample size or non randomness).

      If the confidence intervals are too narrow (i.e., the empirical coverage is lower than the target confidence level), this is often reflected in a sampling distribution 
      that is **too concentrated** (i.e., has too small a standard deviation). This means that estimated coefficients vary less than they should given the uncertainty, 
      causing the intervals to miss the true value more often than expected."
        
      ),
      
      
      bsCollapse(id = "resid_section", open = NULL,
                 bsCollapsePanel("Error Analysis (click to expand/collapse)", style = "primary",
                                 tabsetPanel(
                                   tabPanel("ε vs. Index", 
                                            plotOutput("eps_index_plot")),
                                   tabPanel("ε vs. x₁", 
                                            plotOutput("eps_x1_plot")),
                                   tabPanel("Estimated Error − True Error", 
                                            plotOutput("resid_minus_eps_plot"))
                                 )
                 )
      ),
      
      h4("Scenario-Specific Interpretation"),
      uiOutput("scenario_notes")
      
    )
  )
)

server <- function(input, output) {
  simulate <- eventReactive(input$runSim, {
    beta <- c(input$beta0, input$beta1, input$beta2)
    p <- length(beta)
    rho <- input$rho
    n <- input$n
    reps <- input$reps
    scenario <- input$scenario
    conf_level <- input$conf_level
    
    Sigma <- matrix(c(1, input$rho, input$rho, 1), 2, 2)
    results <- matrix(NA, reps, p)
    coverage <- numeric(p)
    eps_first <- NULL
    residuals_first <- NULL
    x1_first <- NULL
    
    for (r in 1:reps) {
      X_cov <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
      x1 <- X_cov[,1]
      x2 <- X_cov[,2]
      X <- cbind(1, x1, x2)
      
      eps <- switch(scenario,
                    "A" = rnorm(n, mean = 0, sd = 1),
                    "B" = rgamma(n, rate = 1, shape = 1) - 1,
                    "C" = as.numeric(arima.sim(n = n, list(ar = 0.9))),
                    "D" = rnorm(n, mean = 0, sd = abs(x1)),
                    "E" = rnorm(n, mean = 0.5 * x1, sd = 1),
                    "F" = rnorm(n, mean = 0, sd = 1))
      
      if (scenario == "F") {
        interaction <- x1 * x2
        x1_squared <- x1^2
        
        if (input$include_confounder) {
          confounder <- rnorm(n, mean = 0, sd = 1)
          x2 <- x2 + 0.5 * confounder
          X_cov[,2] <- x2
          Y <- X %*% beta + 1.5 * interaction + 3 * confounder + x1_squared + eps
        } else {
          Y <- X %*% beta + 1.5 * interaction + x1_squared + eps
        }
      } else {
        Y <- X %*% beta + eps
      }
      
      fit <- lm(Y ~ x1 + x2)
      b_hat <- coef(fit)
      results[r, ] <- b_hat
      
      ci <- confint(fit, level = conf_level)
      for (j in 1:p) {
        if (ci[j, 1] < beta[j] && ci[j, 2] > beta[j]) {
          coverage[j] <- coverage[j] + 1
        }
      }
      
      if (r == 1) {
        eps_first <- eps
        residuals_first <- resid(fit)
        x1_first <- x1
      }
    }
    
    df <- as.data.frame(results)
    names(df) <- paste0("beta_", 0:(p - 1))
    
    summary_df <- data.frame(
      Beta = c("β₀", "β₁", "β₂"),
      Mean_Estimate = round(colMeans(df), 3),
      #SD = round(apply(df, 2, sd), 3),
      Coverage_rate = round(coverage / reps, 3)
    )
    
    list(
      df = df,
      summary_df = summary_df,
      beta = beta,
      eps = eps_first,
      residuals = residuals_first,
      x1 = x1_first
    )
  })
  
  output$summary_table <- renderTable({
    simulate()$summary_df
  })
  
  output$coef_plot <- renderPlot({
    sim <- simulate()
    df <- sim$df
    beta <- sim$beta
    
    plots <- lapply(1:3, function(j) {
      ggplot(df, aes(x = .data[[names(df)[j]]])) +
        geom_histogram(bins = 30, fill = c("orange", "skyblue", "lightgreen")[j],
                       color = "black", aes(y = after_stat(density)), alpha = 0.5) +
        geom_density(color = "darkblue") +
        geom_function(fun = dnorm,
                      args = list(mean = mean(df[[j]]), sd = sd(df[[j]])),
                      color = "red", linetype = "dashed") +
        geom_vline(xintercept = beta[j], color = "black", linetype = "dotted") +
        labs(title = paste("Sampling Distribution of", names(df)[j]),
             x = names(df)[j], y = "Density") +
        theme_bw()
    })
    
    ggarrange(plotlist = plots, ncol = 3, common.legend = TRUE, legend = "bottom")
  })
  
  output$eps_index_plot <- renderPlot({
    sim <- simulate()
    df_eps <- data.frame(Index = 1:length(sim$eps), Eps = sim$eps)
    
    ggplot(df_eps, aes(x = Index, y = Eps)) +
      geom_point(color = "gray40") +
      labs(title = "True Error (ε) vs. Index i", y = "ε") +
      theme_bw() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })
  
  output$eps_x1_plot <- renderPlot({
    sim <- simulate()
    df_eps <- data.frame(x1 = sim$x1, Eps = sim$eps)
    
    ggplot(df_eps, aes(x = x1, y = Eps)) +
      geom_point(color = "gray40") +
      labs(title = "True Error (ε) vs. x₁", x = "x₁", y = "ε") +
      theme_bw() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })
  
  output$resid_minus_eps_plot <- renderPlot({
    sim <- simulate()
    df_res <- data.frame(Index = 1:length(sim$eps),
                         Estimated_Error_Minus_True = sim$residuals - sim$eps)
    
    ggplot(df_res, aes(x = Index, y = Estimated_Error_Minus_True)) +
      geom_point(color = "gray40") +
      labs(title = expression(hat(epsilon) - epsilon ~ "vs. Index"),
           y = expression(hat(epsilon) - epsilon),
           x = "Index") +
      theme_bw() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })
  
  output$scenario_notes <- renderUI({
    notes <- switch(input$scenario,
                    "A" = HTML("<b>Scenario A:</b><br>
      • All assumptions are met.<br>
      • We expect unbiased estimators for the beta coefficients and that the coverage rate is equal to the set confidence level.<br>
      • We also see that the empirical distribution matches the normal distribution. This even holds true for small sample size n.<br><br><br><br><br><br><br>"),
                    
                    "B" = HTML("<b>Scenario B:</b><br>
      • Error distribution is skewed (Gamma).<br>
      • Try small vs. large sample sizes n and compare the empirical distribution to the normal distribution.<br>
      • We see that the beta coefficients are estimated correctly and the coverage rate matches the set confidence level for small and large n.<br>
                               but the empirical distribution does not match the normal distribution for small n. For large sample size n it matches approximatly, due to the central limit theorem.<br>
                               Can you guess (from small n), what distribution the estimates actually follow?<br><br><br><br><br><br><br>"),
                    
                    "C" = HTML("<b>Scenario C:</b><br>
      • Variance depends on x₁ (heteroscedastic).<br>
      • We expect to see that the estimates for the beta coefficients are still unbiased, but the coverage rate (for beta_1) doesn't match the set confidence level.<br>
      • Why is that?.<br>
      • Try increasing the sample size n to see, if this solves the issue.<br>
      • Take away: We still estimate the coefficients unbiased, but conducting inference can be problematic.<br>
      • In what plot can you see the heteroscedasticity?<br><br><br><br><br><br><br>"),
                    
                    "D" = HTML("<b>Scenario D:</b><br>
      • Errors are autocorrelated (next error dependents on the previous).<br>
      • We expect to see that the estimates for the beta coefficients are still unbiased, but the coverage rate (for beta_0) doesn't match the set confidence level.<br>
      • In what plot can you see the autocorrelation?<br><br><br><br><br><br><br>"),
                    
                    "E" = HTML("<b>Scenario E:</b><br>
      • E[ε|x₁] ≠ 0, i.e. systematic error based on the value of x₁.<br>
      • We observe biased estimates of the beta coefficients (and therefore the coverage rate doesn't match the set confidence level as well).<br><br><br><br><br><br><br>"),
                    
                    "F" = HTML("<b>Scenario F:</b><br>
      • Model misspecification is somewhat tricky, because the results highly depend on the kind of misspecification.<br>
      • In this example, the actual Y depends on x₁ and x₂ but also on other terms (quadratic and interaction which you will learn later). Especially there is the option to include a confounder variable, i.e. a variable that has an influence on Y and on (in this case) x₂.<br>
      • You got the option to include a confounder variable or leave it out. How does this influence the model estimates? .<br>
      • Without confounder, we observe that the beta coefficients of x₁ and x₂ are still estimated correclty, but the intercept is not. Can you explain why? The coverage rate does not match the set confidence level.<br>
      • With confounder, we observe that the beta coefficient of x₂ is now biased. The same holds true for the beta coefficient of x₁. Can you explain why the beta coefficient of x₁ is biased as well, although the confounder is constructed to influence x₂.<br>
      • In the case of included confounder, try to set the correlation of x₁ and x₂ close to 0 (or 0) and see what happens to the coefficient estimator of x₁. Why does this happen?.<br>
      • Take away: Correct model specification is generally the most crucial assumption. If the true underlying relationship is not what we think it is, this will almost always result in wrong inference statements and (in case of missed confounders) in very wrong estimated effects of interest.<br><br><br><br><br><br><br>"),
                    
                    "Unknown scenario."
    )
    return(notes)
  })
  
}

shinyApp(ui = ui, server = server)