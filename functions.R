
n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}

# CORRELATIONS =====
plot_corrplot <- function(data, grd, lg) {
  # Function to create the grade- and language-specific correlation plots of
  # all relevant Multitudes and WJ/WM measures
  # df: always the main df
  # grd: the grade for which to plot the corrplot (K, G1, G2)
  # lang: the language for which to plot the corrplot (English, Spanish)
  
  df.plot <- data |> 
    filter(grade == grd
           # model_lang == lg,
           # task_lang == lg
    ) |> 
    select(c(student_id, grade, task, score)) |> 
    mutate(task = names[task]) |>
    pivot_wider(names_from = task, values_from = score, values_fn = mean) |>
    filter(!duplicated(student_id))
  
  df.plot <- df.plot |> 
    # Create correlation matrix and format for plotting
    select(where(is.numeric)) |>
    correlate() |>
    stretch() |>
    group_by(x, y) |>
    # add n
    mutate(
      n = sum(!is.na(df.plot[[x]]) & !is.na(df.plot[[y]]))
    ) |>
    ungroup() |>
    filter(!is.na(r)) |>
    # Add p-values
    mutate(
      t_stat = (r * sqrt(n-2))/sqrt(1-r^2),
      p_value = 2 * pt(abs(t_stat), n-2, lower.tail = FALSE)
    )
  
  # Create the heatmap
  fig.heatmap <- df.plot |> 
    filter(as.numeric(factor(x)) < as.numeric(factor(y))) |>
    ggplot(aes(x = x, y = y, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "white",
      mid = "white",
      high = "darkgreen",
      midpoint = 0,
      limits = c(-0.1, 1)
    ) +
    geom_text(aes(
      label = round(r, 2),
      color = p_value < 0.05  # color based on statistical significance
    ),
    size = 3
    ) +
    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "lightgray"),
      guide = "none"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = c(.8, .3)
    ) +
    labs(
      fill = "  r"
    ) +
    scale_x_discrete(position = "top")
  
  # Create the n-chart
  fig.nchart <- df.plot |> 
    filter(as.numeric(factor(x)) < as.numeric(factor(y))) |>
    mutate(n100 = if_else(n >= 100, "100+", "< 100")) |> 
    ggplot(aes(x = x, y = y, fill = n100)) +
    geom_tile() +
    scale_fill_manual(values = nchart_colours) +
    geom_text(aes(
      label = n,
      color = p_value < 0.05  # color based on statistical significance
    ),
    size = 3
    ) +
    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "lightgray"),
      guide = "none"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = c(.8, .3)
    ) +
    labs(
      fill = "  n"
    ) +
    scale_x_discrete(breaks = NULL)
  
  plot <- ggarrange(plotlist = list(fig.heatmap, fig.nchart),
                    nrow = 2,
                    heights = c(.62, .38),
                    labels = "AUTO"
  )
  
  
  return(plot)
}

# BORUTA FEATURE SELECTION =======

boruta_task_selection <- function(data, grd, lg) {
  
  if (lg == "All") {
    df.boruta <- data |> 
      group_by(student_id, task) %>%
      slice_head(n = 1) %>%  # remove 2nd+ administrations of tasks
      ungroup() %>%
      filter(grade == grd,
             # model_lang == lg,
             # task_lang == lg,
             !grepl("wc", task),
             !is.na(risk)
      ) |> 
      select(c(student_id, task, score, risk)) |> 
      pivot_wider(names_from = task, values_from = score) |> 
      # recode outcomes for boruta algorithm
      mutate(risk = case_when(risk == 0 ~ -1,
                              risk == 1 ~ 1
      )
      ) |> 
      select(-student_id)
  } else {
    df.boruta <- data |> 
      group_by(student_id, task) %>%
      slice_head(n = 1) %>%  # remove 2nd+ administrations of tasks
      ungroup() %>%
      filter(grade == grd,
             model_lang == lg,
             task_lang == lg,
             !grepl("wc", task),
             !is.na(risk)
      ) |> 
      select(c(student_id, task, score, risk)) |> 
      pivot_wider(names_from = task, values_from = score) |> 
      # recode outcomes for boruta algorithm
      mutate(risk = case_when(risk == 0 ~ -1,
                              risk == 1 ~ 1
      )
      ) |> 
      select(-student_id)
  }
  
  
  # 1. impute missing data (multivariate imputation by chained equations)
  imp <- mice(df.boruta, print = FALSE, seed = 1)
  
  # 2. run boruta on each of the 5 imputed datasets
  boruta_results <- vector("list", 5)
  for (i in 1:5) {
    imputed_df <- complete(imp, i)
    boruta_results[[i]] <- Boruta(risk ~ .,
                                  data = imputed_df,
                                  doTrace = 0,
                                  # maxRuns = 500
                                  # ntrees = 100000
    )
  }
  
  # 3. aggregate results (both mean importance and selection rate)
  process_boruta_results <- function(boruta_list) {
    # Extract full ImpHistory
    full_imp_history <- boruta_list %>%
      map(~ .x$ImpHistory) %>%
      set_names(paste0("Imputation_", seq_along(.)))
    
    # Process importance data
    importance_data <- boruta_list %>%
      imap_dfr(~ {
        imp_history <- as.data.frame(.x$ImpHistory)
        imp_history %>%
          pivot_longer(everything(), names_to = "Feature", values_to = "Importance") %>%
          mutate(imputation = .y, run = rep(1:nrow(.x$ImpHistory), ncol(.x$ImpHistory)))
      })
    
    # Process decision data
    decision_data <- boruta_list %>%
      map_dfr(~ tibble(
        Feature = names(.x$finalDecision),
        Decision = as.character(.x$finalDecision)
      ))
    
    # Calculate summary statistics
    importance_summary <- importance_data %>%
      group_by(Feature) %>%
      summarise(
        MeanImportance = mean(Importance, na.rm = TRUE),
        SE = sd(Importance, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        CI_lower = MeanImportance - qt(0.975, df = n - 1) * SE,
        CI_upper = MeanImportance + qt(0.975, df = n - 1) * SE
      )
    
    # Aggregate results
    aggregated_results <- importance_summary %>%
      left_join(
        decision_data %>%
          group_by(Feature) %>%
          summarise(
            SelectionRate = mean(Decision == "Confirmed", na.rm = TRUE),
            .groups = 'drop'
          ),
        by = "Feature"
      ) %>%
      arrange(desc(MeanImportance))
    
    list(
      full_imp_history = full_imp_history,
      importance_data = importance_data,
      aggregated_results = aggregated_results
    )
  }
  
  # aggregate and return the Boruta results
  results <- process_boruta_results(boruta_results)
}

plot_boruta <- function(borutadata, disaggregated = TRUE, type = "A") {
  # ADD DESCRIPTION
  # borutadata: aggregated_results output from the process_boruta_results function, either single occurrence or multiple rbound ones with group variable
  # disaggregated: boolean indicating whether faceting by group (ELPD) is desired
  # type: 'A' returns a faceted (by group) horizontal bar chart, 'b' returns a faceted (by task) vertical bar chart
  
  # use proper task names
  borutadata <- borutadata |> 
    mutate(Feature = factor(Feature, levels = names(names), labels = names),
           # group = case_when(group == "EL" ~ "English Learners",
           #                   group == "EO" ~ "English-only Students",
           #                   group == "All" ~ "All",
           #                   TRUE ~ group
           #                   )
           ) #|> 
    # user full ELPD labels
    # mutate(group = factor(group, levels = c("All", "English Learners", "English-only Students")))
  
  # default to type a for non-disaggregated view
  # if (!'group' %in% colnames(borutadata)) {
  #   disaggregated = FALSE
  #   type = "A"
  # }
  
  if (disaggregated) {
    
    if (type == "A") {
      # fix order
      order <- borutadata |> 
        filter(group == "All") |>
        group_by(Feature) %>%
        summarise(mean = mean(MeanImportance)) %>%
        arrange(mean) %>%
        pull(Feature)
      borutadata$Feature <- factor(borutadata$Feature, levels = order)
      
      plot <- borutadata |> 
        mutate(MeanImportance = case_when(MeanImportance == -Inf ~ 0,
                                          MeanImportance != -Inf ~ MeanImportance),
               # group = factor(group, levels = c("All", "English Learners", "English-only Students"))
               group = factor(group, levels = c("All", "EL", "EO"))
        ) |> 
        ggplot(aes(x = Feature,
                   y = MeanImportance,
                   fill = group,
                   alpha = is.na(SelectionRate),
                   group = group)
               ) +
        geom_hline(aes(yintercept = 0), colour = "gray") +
        geom_col(position = "dodge") +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), #width = 0.2,
                      position = "dodge") +
        # facet_grid(cols = vars(group)) +
        labs(#x = "Feature",
             y = "Mean Importance",
             fill = "Group",
             # title = "Feature Importance from Boruta across Multiple Imputations",
             # subtitle = "Error bars represent 95% confidence intervals"
        ) + 
        scale_alpha_manual(values = c(1, .4)) +
        theme(legend.position = "none",
              axis.title.y = element_blank()
              ) + 
        guides(alpha = "none") +
        coord_flip() +
        scale_fill_manual(values = palette2)
      
    } else if (type == "B") {
      # fix order
      order <- borutadata |> 
        filter(group == "All") |> 
        group_by(Feature) %>%
        summarise(mean = mean(MeanImportance)) %>%
        arrange(desc(mean)) %>%
        pull(Feature)
      borutadata$Feature <- factor(borutadata$Feature, levels = order)
      
      plot <- borutadata |>
        ggplot(aes(x = group, y = MeanImportance, fill = group, alpha = is.na(SelectionRate))) +
        # ggplot(aes(x = MeanImportance, y = MeanImportance, fill = SelectionRate)) +
        geom_col() +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
        # coord_flip() +
        # scale_fill_viridis_c() +
        labs(#x = "Feature",
             y = "Mean Importance",
             fill = "Group",
             title = "Feature Importance from Boruta across Multiple Imputations",
             subtitle = "Error bars represent 95% confidence intervals") +
        scale_alpha_manual(values = c(1, .4)) + + 
        scale_fill_manual(values = palette2) +
        facet_grid(cols = vars(Feature)) +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
              ) +
        guides(alpha = "none")

    }
  } else if (!disaggregated) {
    # type A
    plot <- borutadata |> 
      filter(group == "All") |>
      ggplot(aes(x = reorder(Feature, MeanImportance),
                 y = MeanImportance,
                 alpha = is.na(SelectionRate))
      ) +
      geom_hline(aes(yintercept = 0), colour = "gray") +
      geom_col(fill = col.all) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
      coord_flip() +
      scale_alpha_manual(values = c(1, .4)) +
      # scale_fill_manual(values = palette1) +
      labs(x = "Feature", 
           y = "Mean Importance", fill = "Selection Rate",
           # title = "Feature Importance from Boruta across Multiple Imputations",
           # subtitle = "Error bars represent 95% confidence intervals"
           ) +
      theme(legend.position = "none")
  }
  return(plot)
}

# PREDICTION MODEL =====
build_prediction_model_cv <- function(data, grd, lg, mdl) {
  # Updated function to run a BOY to EOY prediction model (logistic regression)
  # and also identify optimal cut and method based on training data only
  # df: always the main df
  # grd: the grade for which to run the model
  # lg: the language for which to run the model
  # mdl: the model specification
  
  mdl <- paste0("risk ~ ", mdl)
  
  df.model0 <- data |> 
    filter(grade == grd,
           model_lang == lg,
           task_lang == lg,
           !is.na(risk)
    ) |> 
    select(c(student_id,
             el,
             group = district,
             risk,
             task,
             score
    )
    ) |> 
    unique() |> 
    pivot_wider(names_from = task, values_from = score, values_fn = mean) |> 
    ungroup()
  
  # run model with leave-one-group-out cross-validation (group = district)
  groups <- unique(df.model0$group)
  
  # df to store LOGO-CV results
  df.results <- data.frame(student_id = as.character(),
                           pred.logit = as.numeric(),
                           pred.prob = as.numeric(),
                           pred.risk = as.numeric()
  )
  opt.cuts <- c()
  # implement LOGO-CV ======
  for (g in groups) {
    
    # Skip problematic or small groups
    if (grd == "G2" & lg == "Spanish" & g == "district_15") next
    if (nrow(df.model0 %>% filter(group == g)) < 10) next
    
    # test-train split =====
    # Create 'ignore' flag for test group
    df.model <- df.model0 %>%
      mutate(ignore = if_else(group == g, TRUE, FALSE))
    
    # impute and fit model on training data oly ======
    # Impute missing data (train/test split maintained with `ignore`)
    imp <- mice(df.model, ignore = df.model$ignore, m = 5, print = FALSE, seed = 1)
    
    # Fit model using only training data
    fit.imp <- with(imp, glm(as.formula(mdl), family = "binomial"), subset = ignore == FALSE)
    
    # pool estimates for model output
    model = pool(fit.imp)
    
    # get pooled coefficients ====
    coef_vector <- summary(pool(fit.imp)) %>%
      select(term, estimate) %>%
      deframe()
    
    # determine optimal cutoff based on training data======
    # 1. Get training rows
    train_rows <- which(df.model$ignore == FALSE)
    # 2. Get and combine all completed datasets and models
    completed.data <- complete(imp, action = "all") 
    imputations <- tibble(
      imp = 1:5,
      data = completed.data,
      model = getfit(fit.imp)
    )
    # 3. Process and combine all training datasets
    df.train.all <- imputations %>%
      mutate(
        df.train = map(data, ~ .x[train_rows, ]),
        pred_prob = map2(model, df.train, ~ predict(.x, newdata = .y, type = "response")),
        df.train = map2(df.train, pred_prob, ~ mutate(.x, pred.prob = .y))
      ) %>%
      select(imp, df.train) %>%
      unnest(df.train)
    # 4. determine optimal cutoff based on training data
    opt.cut <- optimal.cutpoints(
      X = "pred.prob",
      status = "risk",
      tag.healthy = 0,
      method = "Youden",
      data = as.data.frame(df.train.all),
      ci.fit = FALSE,
      trace = FALSE
    )[["Youden"]]$Global$optimal.cutoff$cutoff
    # save optcut for later aggregation
    opt.cuts <- c(opt.cuts, opt.cut)
    
    # predictions for validation data =======
    # using pooled coefficients (from coef_vector) and fold-specific opt.cut
    # 1. Get validation rows
    val_rows <- which(df.model$ignore == TRUE)
    # 2. Get and combine all completed datasets and models
    completed.data <- complete(imp, action = "all") 
    # For each imputation, take test rows and predict with pooled coefficients
    df.test.preds <- map_dfr(completed.data, function(imp_data) {
  
      test_data <- imp_data[val_rows, ]
      # Add intercept column for matrix multiplication
      test_data <- test_data %>%
        mutate(`(Intercept)` = 1) %>%
        select(all_of(names(coef_vector)))  # select columns matching model terms including intercept
      # Calculate prob and logit
      pred.logit <- as.vector(as.matrix(test_data) %*% coef_vector)
      pred.prob <- plogis(pred.logit)
      tibble(
        student_id = imp_data[val_rows, "student_id"],
        pred.logit = pred.logit,
        pred.prob = pred.prob
      )
    }, .id = "imp")
    
    # Average predicted probabilities across imputations for each student
    df.results.temp <- df.test.preds %>%
      group_by(student_id) %>%
      summarize(
        pred.logit = mean(pred.logit),
        pred.prob = mean(pred.prob)
      ) %>%
      ungroup() %>% 
    # Add true risk from original df.model
      left_join(df.model %>% filter(ignore == TRUE) %>% select(student_id, risk), by = "student_id") %>%
      mutate(
        pred.risk = if_else(pred.prob > opt.cut, 1, 0)
      )
    
    # Append results
    df.results <- rbind(df.results, df.results.temp)
  }
  
  if (nrow(df.results) == 0) {
    warning("No test groups met criteria for modeling. Check data filtering.")
    return(NULL)
  }
  
  opt.cut <- mean(opt.cuts, na.rm = TRUE)
  return(list(results = df.results, model = model, opt.cut = opt.cut))
}

# PREDICTION MODEL =====
build_prediction_model_bootstrap_cv <- function(data, grd, lg, mdl, n_boot = 10, el_ratio = 0.2) {
  # Updated function to run a BOY to EOY prediction model (logistic regression)
  # and also identify optimal cut and method based on training data only
  # df: always the main df
  # grd: the grade for which to run the model
  # lg: the language for which to run the model
  # mdl: the model specification
  
  mdl <- paste0("risk ~ ", mdl)
  
  df.model0 <- data |> 
    filter(grade == grd,
           model_lang == lg,
           task_lang == lg,
           !is.na(risk)
    ) |> 
    select(c(student_id,
             el,
             group = district,
             risk,
             task,
             score
    )
    ) |> 
    unique() |> 
    pivot_wider(names_from = task, values_from = score, values_fn = mean) |> 
    ungroup()
  
  # run model with leave-one-group-out cross-validation (group = district)
  groups <- unique(df.model0$group)
  
  # df to store LOGO-CV results
  df.results <- data.frame(student_id = as.character(),
                           pred.logit = as.numeric(),
                           pred.prob = as.numeric(),
                           pred.risk = as.numeric()
  )
  opt.cuts <- c()
  # implement LOGO-CV ======
  for (g in groups) {
    
    # Skip problematic or small groups
    if (grd == "G2" & lg == "Spanish" & g == "district_15") next
    if (nrow(df.model0 %>% filter(group == g)) < 10) next
    
    # test-train split =====
    # Create 'ignore' flag for test group
    df.model <- df.model0 %>%
      mutate(ignore = if_else(group == g, TRUE, FALSE))
    
    # === stratified bootstrapping on training data only ===
    train_data_all <- list()
    df.train_base <- df.model %>% filter(!ignore)
    
    for (b in 1:n_boot) {
      n_total <- nrow(df.train_base)
      n_el <- round(n_total * el_ratio)
      n_eo <- n_total - n_el
      
      df_el <- df.train_base %>% filter(el == "EL") %>% sample_n(min(n_el, n()), replace = TRUE)
      df_eo <- df.train_base %>% filter(el == "EO") %>% sample_n(min(n_eo, n()), replace = TRUE)
      
      df_train_boot <- bind_rows(df_el, df_eo) %>% sample_frac(1)
      train_data_all[[b]] <- df_train_boot
    }
    
    df.train.all.boot <- bind_rows(train_data_all, .id = "boot")
    
    # Combine with validation
    df.model <- bind_rows(
      df.train.all.boot %>% select(-boot),
      df.model %>% filter(ignore)
    ) %>% mutate(ignore = student_id %in% df.model$student_id[df.model$ignore])
    
    # impute and fit model on training data only ======
    # Impute missing data (train/test split maintained with `ignore`)
    imp <- mice(df.model, ignore = df.model$ignore, m = 5, print = FALSE)
    
    # Fit model using only training data
    fit.imp <- with(imp, glm(as.formula(mdl), family = "binomial"), subset = ignore == FALSE)
    
    # pool estimates for model output
    model = pool(fit.imp)
    
    # get pooled coefficients ====
    coef_vector <- summary(pool(fit.imp)) %>%
      select(term, estimate) %>%
      deframe()
    
    # determine optimal cutoff based on training data======
    # 1. Get training rows
    train_rows <- which(df.model$ignore == FALSE)
    # 2. Get and combine all completed datasets and models
    completed.data <- complete(imp, action = "all") 
    imputations <- tibble(
      imp = 1:5,
      data = completed.data,
      model = getfit(fit.imp)
    )
    # 3. Process and combine all training datasets
    df.train.all.imputed <- imputations %>% 
      mutate(
        df.train = map(data, ~ .x[train_rows, ]),
        pred_prob = map2(model, df.train, ~ predict(.x, newdata = .y, type = "response")),
        df.train = map2(df.train, pred_prob, ~ mutate(.x, pred.prob = .y))
      ) %>%
      select(imp, df.train) %>%
      unnest(df.train)
    # 4. determine optimal cutoff based on training data
    opt.cut <- optimal.cutpoints(
      X = "pred.prob",
      status = "risk",
      tag.healthy = 0,
      method = "Youden",
      data = as.data.frame(df.train.all.imputed),
      ci.fit = FALSE,
      trace = FALSE
    )[["Youden"]]$Global$optimal.cutoff$cutoff
    # save optcut for later aggregation
    opt.cuts <- c(opt.cuts, opt.cut)
    
    # predictions for validation data =======
    # using pooled coefficients (from coef_vector) and fold-specific opt.cut
    # 1. Get validation rows
    val_rows <- which(df.model$ignore == TRUE)
    # 2. Get and combine all completed datasets and models
    completed.data <- complete(imp, action = "all") 
    # For each imputation, take test rows and predict with pooled coefficients
    df.test.preds <- map_dfr(completed.data, function(imp_data) {
      
      test_data <- imp_data[val_rows, ]
      # Add intercept column for matrix multiplication
      test_data <- test_data %>%
        mutate(`(Intercept)` = 1) %>%
        select(all_of(names(coef_vector)))  # select columns matching model terms including intercept
      # Calculate prob and logit
      pred.logit <- as.vector(as.matrix(test_data) %*% coef_vector)
      pred.prob <- plogis(pred.logit)
      tibble(
        student_id = imp_data[val_rows, "student_id"],
        pred.logit = pred.logit,
        pred.prob = pred.prob
      )
    }, .id = "imp")
    
    # Average predicted probabilities across imputations for each student
    df.results.temp <- df.test.preds %>%
      group_by(student_id) %>%
      summarize(
        pred.logit = mean(pred.logit),
        pred.prob = mean(pred.prob)
      ) %>%
      ungroup() %>% 
      # Add true risk from original df.model
      left_join(df.model %>% filter(ignore == TRUE) %>% select(student_id, risk), by = "student_id") %>%
      mutate(
        pred.risk = if_else(pred.prob > opt.cut, 1, 0)
      )
    
    # Append results
    df.results <- rbind(df.results, df.results.temp)
  }
  
  if (nrow(df.results) == 0) {
    warning("No test groups met criteria for modeling. Check data filtering.")
    return(NULL)
  }
  
  opt.cut <- mean(opt.cuts, na.rm = TRUE)
  return(list(results = df.results, model = model, opt.cut = opt.cut))
}

train_model_final <- function(data, grd, lg, mdl) {
  # df: always the main df
  # grd: the grade for which to run the model
  # lg: the language for which to run the model
  # mdl: the model specification
  
  mdl <- paste0("risk ~ ", mdl)
  
  df.model <- data |> 
    filter(grade == grd,
           model_lang == lg,
           task_lang == lg,
           !is.na(risk)
    ) |> 
    select(c(student_id,
             group = district,
             risk,
             task,
             score
    )
    ) |> 
    unique() |> 
    pivot_wider(names_from = task, values_from = score, values_fn = mean) |> 
    ungroup()

  # impute missing data ====
  imp <- mice(df.model, m = 5, print = FALSE, seed = 1)
  
  # fit model ====
  fit.imp <- with(imp, glm(as.formula(mdl), family = "binomial"))
  # pool estimates for model output
  model = pool(fit.imp)
  # get pooled coefficients ====
  coef_vector <- summary(pool(fit.imp)) %>%
    select(term, estimate) %>%
    deframe()
  
  # determine optimal cutoff ======
  # 1. Get and combine all completed datasets and models
  completed.data <- complete(imp, action = "all") 
  imputations <- tibble(
    imp = 1:5,
    data = completed.data,
    model = getfit(fit.imp)
  )
  # 2. Process and combine all training datasets
  df.all <- imputations %>%
    mutate(
      df.train = map(data, ~ .x),
      pred_prob = map2(model, df.train, ~ predict(.x, newdata = .y, type = "response")),
      df.train = map2(df.train, pred_prob, ~ mutate(.x, pred.prob = .y))
    ) %>%
    select(imp, df.train) %>%
    unnest(df.train)
  # 3. determine optimal cutoff
  opt.cut <- optimal.cutpoints(
    X = "pred.prob",
    status = "risk",
    tag.healthy = 0,
    method = "Youden",
    data = as.data.frame(df.all),
    ci.fit = FALSE,
    trace = FALSE
  )[["Youden"]]$Global$optimal.cutoff$cutoff

  return(list(model, opt.cut))
}


extract_metrics <- function(cm) {
  # Function to extract metrics and create a single-row dataframe
  
  # Extract overall statistics
  overall_stats <- cm$overall
  
  # Extract statistics by class (assuming binary classification)
  by_class_stats <- cm$byClass
  
  # Combine all metrics
  all_metrics <- c(overall_stats, by_class_stats)
  
  # Create a dataframe
  df_metrics <- data.frame(t(all_metrics))
  
  # Optionally, you can round the values to a specific number of decimal places
  df_metrics <- df_metrics %>% mutate(across(where(is.numeric), ~round(., 4)))
  
  return(df_metrics)
}

# (P)ROC PLOTTING =======

# Function to calculate ROC data
calculate_roc_information <- function(data, group) {
  
  actuals = data$risk
  preds = data$pred.prob
  
  # actuals = df.results.k.en$risk
  # preds = df.results.k.en$pred.prob
  # df = df.results.k.en
  # group = "All"
  
  roc <- roc(actuals, preds,
             direction = "<",
             levels = c(0, 1),
             ci = TRUE,
             n.thresholds = length(unique(preds)) + 1
  )
  
  proc <- pr.curve(scores.class0 = preds[actuals == 1],
                   scores.class1 = preds[actuals == 0],
                   curve = TRUE)
  
  
  optcut.youden <- coords(roc, "best", best.method = "youden")
  optcut.topleft <- coords(roc, "best", best.method = "closest.topleft")
  optcut.sens80 <- optimal.cutpoints(X = "pred.prob",
                                     status = "risk",
                                     tag.healthy = 0,
                                     methods = "MinValueSe",
                                     data = as.data.frame(data),
                                     pop.prev = NULL, 
                                     control = control.cutpoints(valueSe = 0.80),
                                     ci.fit = TRUE,
                                     conf.level = 0.95,
                                     trace = FALSE
  )$MinValueSe$Global$optimal.cutoff[[1]]
  optcut.spec80 <- optimal.cutpoints(X = "pred.prob",
                                     status = "risk",
                                     tag.healthy = 0,
                                     methods = "MinValueSe",
                                     data = as.data.frame(data),
                                     pop.prev = NULL, 
                                     control = control.cutpoints(valueSp = 0.80),
                                     ci.fit = TRUE,
                                     conf.level = 0.95,
                                     trace = FALSE
  )$MinValueSe$Global$optimal.cutoff[[1]]
  
  df.out <- data.frame(coord.FPR = 1 - roc$specificities,
                       coord.TPR = roc$sensitivities,
                       AUROC = auc(roc)[1],
                       CI.lower = roc$ci[1],
                       CI.upper = roc$ci[2],
                       group = group,
                       optcut.youden = optcut.youden[[1]],
                       optcut.topleft = optcut.topleft[[1]],
                       optcut.sens80 = optcut.sens80[[1]],
                       optcut.spec80 = optcut.spec80[[1]],
                       recall = proc$curve[, 1][1:length(roc$specificities)],
                       precision = proc$curve[, 2][1:length(roc$specificities)],
                       # recall = proc$curve[, 1][1:(length(proc$curve[, 1]) - 1)],
                       # precision = proc$curve[, 2][1:(length(proc$curve[, 2]) - 1)],
                       AUPROC = proc$auc.integral
  ) |> 
    # mutate(coord.FPR = 1 - coord.FPR) |> 
    unique()
  
  return(df.out)
}

plot_roc_curves_within_model <- function(data, disaggregated = TRUE, by = "EL") {
  
  # Calculate ROC data for entire sample
  df.roc.all <- calculate_roc_information(data, "All")
  
  if (disaggregated & by == "EL") {
    # check whether enough data exists
    skipEL = (nrow(df.plots %>% filter(el == "EL", risk == 1)) == 0 | nrow(df.plots %>% filter(el == "EL", risk == 0)) == 0)
    skipEO = (nrow(df.plots %>% filter(el == "EO", risk == 1)) == 0 | nrow(df.plots %>% filter(el == "EO", risk == 0)) == 0)
    
    # Calculate ROC data for subgroups and combine as applicable
    if (skipEL & skipEO) {
      df.roc <- df.roc.all
    } else if (skipEL) {
      df.roc.eo <- calculate_roc_information(data |> filter(el == "EO"), "EO")
      df.roc <- df.roc.all |> 
        rbind(df.roc.eo)
    } else if (skipEO) {
      df.roc.el <- calculate_roc_information(data |> filter(el == "EL"), "EL")
      df.roc <- df.roc.all |> 
        rbind(df.roc.el)
    } else {
      df.roc.el <- calculate_roc_information(data |> filter(el == "EL"), "EL")
      df.roc.eo <- calculate_roc_information(data |> filter(el == "EO"), "EO")
      df.roc <- df.roc.all |> 
        rbind(df.roc.el) %>% 
        rbind(df.roc.eo)
    }

    # ever_disability and black students are not shown due to low n
    # df.roc.iep <- calculate_roc_information(df |> filter(ever_disability == "IEP/504"), "Ever IEP/504")
    # df.roc.black <- calculate_roc_information(df |> filter(race == "Black/African American"), "African American")
    
    # Combine 
    df.roc <- df.roc |>
      group_by(group) |> 
      mutate(n = n(),
             legend.roc = paste0(group, " (AUC = ", round(mean(AUROC), 2), ")"),
             # legend.roc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUROC), 2), ")"),
             legend.proc = paste0(group, " (AUC = ", round(mean(AUPROC), 2), ")"),
             # legend.proc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUPROC), 2), ")"),
      ) |> 
      ungroup()
    
    youden_x <- 1 - mean(df.roc.all$optcut.youden.specificity)
    youden_y <- mean(df.roc.all$optcut.youden.sensitivity)
  
    # Create the ROC gplot
    plot.roc <- ggplot(df.roc, aes(x = coord.FPR, y = coord.TPR, color = legend.roc)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      geom_hline(aes(yintercept = youden_y), colour = "gray") +
      geom_vline(aes(xintercept = youden_x), colour = "gray") +
      geom_line() +
      labs(x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)",
           colour = "Group",
           #caption = "Notes. Where subgroup sample sizes do not add up to the n for All, this is due to \n students classified as English-proficient or with missing English proficiency information."
      ) +
      theme(legend.position = c(.75, .2)) +
      scale_colour_manual(values = palette1)
    
    # Create the PROC plot
    plot.proc <- ggplot(df.roc, aes(x = recall, y = precision, colour = legend.proc)) +
      geom_line() +
      # geom_area() +  # Add shaded area under curve
      xlim(0, 1) +
      ylim(0, 1) +
      labs(#title = "Precision-Recall Curve",
        x = "Recall",
        y = "Precision",
        colour = "Groups and AUCs") +
      # Add AUC annotations
      # annotate(x = 0.3, y = .3,  geom = "text", label = "Groups and AUCs") +
      # geom_text(data = df.roc %>% group_by(group) %>% slice(1),
      #           aes(x = 0.2, y = 0.25 - as.numeric(factor(group)) * 0.07, 
      #               label = sprintf("%s (n = %d), AUC = %.2f", group, n, AUPROC)),
      #           hjust = 0, vjust = 1, size = 3) +
      theme(legend.position = c(.25, .2)) +
      scale_colour_manual(values = palette1)
    
  } else if (disaggregated & by == "instruction") {
    
    # Calculate ROC data for subgroups
    df.roc.en <- calculate_roc_information(data |> filter(lg_inst == "English"), "English")
    df.roc.dual <- calculate_roc_information(data |> filter(lg_inst == "Dual Lg."), "Dual Lg.")
    # ever_disability and black students are not shown due to low n
    # df.roc.iep <- calculate_roc_information(df |> filter(ever_disability == "IEP/504"), "Ever IEP/504")
    # df.roc.black <- calculate_roc_information(df |> filter(race == "Black/African American"), "African American")
    
    # Combine 
    df.roc <- df.roc.all |> 
      rbind(df.roc.en) |> 
      rbind(df.roc.dual) |> 
      # rbind(df.roc.iep) |> 
      # rbind(df.roc.black) |>
      group_by(group) |> 
      mutate(n = n(),
             legend.roc = paste0(group, " (AUC = ", round(mean(AUROC), 2), ")"),
             # legend.roc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUROC), 2), ")"),
             legend.proc = paste0(group, " (AUC = ", round(mean(AUPROC), 2), ")"),
             # legend.proc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUPROC), 2), ")"),
      ) |> 
      ungroup()
    
    youden_x <- 1 - mean(df.roc.all$optcut.youden.specificity)
    youden_y <- mean(df.roc.all$optcut.youden.sensitivity)
    
    # Create the ROC gplot
    plot.roc <- ggplot(df.roc, aes(x = coord.FPR, y = coord.TPR, color = legend.roc)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      geom_hline(aes(yintercept = youden_y), colour = "gray") +
      geom_vline(aes(xintercept = youden_x), colour = "gray") +
      geom_line() +
      labs(x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)",
           colour = "Group",
           #caption = "Notes. Where subgroup sample sizes do not add up to the n for All, this is due to \n students classified as English-proficient or with missing English proficiency information."
      ) +
      theme(legend.position = c(.75, .2)) +
      scale_colour_manual(values = palette1)
    
    # Create the PROC plot
    plot.proc <- ggplot(df.roc, aes(x = recall, y = precision, colour = legend.proc)) +
      geom_line() +
      # geom_area() +  # Add shaded area under curve
      xlim(0, 1) +
      ylim(0, 1) +
      labs(#title = "Precision-Recall Curve",
        x = "Recall",
        y = "Precision",
        colour = "Groups and AUCs") +
      # Add AUC annotations
      # annotate(x = 0.3, y = .3,  geom = "text", label = "Groups and AUCs") +
      # geom_text(data = df.roc %>% group_by(group) %>% slice(1),
      #           aes(x = 0.2, y = 0.25 - as.numeric(factor(group)) * 0.07, 
      #               label = sprintf("%s (n = %d), AUC = %.2f", group, n, AUPROC)),
      #           hjust = 0, vjust = 1, size = 3) +
      theme(legend.position = c(.25, .2)) +
      scale_colour_manual(values = palette1)
    
  } else if (!disaggregated) {
    # Calculate ROC data for entire sample
    df.roc <- df.roc.all |> 
      group_by(group) |> 
      mutate(n = n(),
             # legend.roc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUROC), 2), ")")
             legend.roc = paste0(group, " (AUC = ", round(mean(AUROC), 2), ")")
      ) |> 
      ungroup()
    
    youden_x <- 1 - mean(df.roc$optcut.youden.specificity)
    youden_y <- mean(df.roc$optcut.youden.sensitivity)
    auc <-  mean(df.roc$AUROC)
      
    # Create the ROC gplot
    plot.roc <- ggplot(df.roc, aes(x = coord.FPR, y = coord.TPR, colour = legend.roc)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      geom_hline(aes(yintercept = youden_y), colour = "gray") +
      geom_vline(aes(xintercept = youden_x), colour = "gray") +
      geom_line(colour = col.all) +
      labs(x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)",
           colour = "Group and AUC",
           #caption = "Notes. Where subgroup sample sizes do not add up to the n for All, this is due to \n students classified as English-proficient or with missing English proficiency information."
      ) +
      theme(legend.position = c(.75, .2)) +
      scale_colour_manual(values = palette1)
    
    # Create the PROC plot
    plot.proc <- ggplot(df.roc, aes(x = recall, y = precision, colour = legend.proc)) +
      geom_line(colour = col.all) +
      # geom_area() +  # Add shaded area under curve
      xlim(0, 1) +
      ylim(0, 1) +
      labs(#title = "Precision-Recall Curve",
        x = "Recall",
        y = "Precision") +
      theme(legend.position = c(.25, .2)) +
      scale_colour_manual(values = palette1)
  }

  return(list(plot.roc, plot.proc))
}

plot_roc_curves_across_models <- function(dfs, group = "All", output = NA) {
  # takes in a list of dfs, from different models; only returns ROC curve for one group
  
  col <- palette2[group]
    
  if (group == "All") {
    dfs <- dfs |> 
      map(~ {
        .x  |> 
          mutate(el = "All")
      })
  } else if (group %in% c("EO", "EL")) {
    dfs <- dfs |> 
      map(~ {
        .x  |> 
          filter(el == group)
      })
  } else if (group %in% c("English", "Dual Lg.")) {
    dfs <- dfs |> 
      map(~ {
        .x  |> 
          filter(lg_inst == group)
      })
  }
  
  # Extract models
  extract_models <- function(df) {
    unique(df$model)
  }
  
  # Extract information from each dataframe
  models <- dfs |> 
    map(extract_models)
  
  # Calculate ROC data for each df and add model name
  dfs.roc <- map2(dfs, models, function(df, models) {
    roc_info <- calculate_roc_information(df, group)
    roc_info |> 
      mutate(group = models) # Store models as a list column
  }
  )
  
  # combine dfs
  df.roc <- dfs.roc |> 
    bind_rows() |> 
    group_by(group) |> 
    mutate(n = n(),
           legend.roc = paste0(group, " (AUC = ", round(mean(AUROC), 2), ")"),
           # legend.roc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUROC), 2), ")"),
           legend.proc = paste0(group, " (AUC = ", round(mean(AUPROC), 2), ")"),
           # legend.proc = paste0(group, " (n = ", n, ", AUC = ", round(mean(AUPROC), 2), ")"),
           ) |> 
    ungroup()
  
  # youden_x <- 1 - mean(df.roc$optcut.youden.specificity)
  # youden_y <- mean(df.roc$optcut.youden.sensitivity)
  
  # Create the ROC gplot
  plot.roc <- ggplot(df.roc, aes(x = coord.FPR,
                                 y = coord.TPR,
                                 # linetype = legend.roc
                                 colour = legend.roc
                                 )
                     ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_hline(aes(yintercept = youden_y), colour = "gray") +
    # geom_vline(aes(xintercept = youden_x), colour = "gray") +
    # geom_line(colour = col) +
    geom_line() +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         linetype = "Models and AUCs",
         colour = "Models and AUCs"
         #caption = "Notes. Where subgroup sample sizes do not add up to the n for All, this is due to \n students classified as English-proficient or with missing English proficiency information."
    ) +
    theme(legend.position = c(.65, .25)) +
    scale_linetype_manual(values = scale3)
  
  # Create the PROC plot
  plot.proc <- ggplot(df.roc, aes(x = recall,
                                  y = precision,
                                  # linetype = legend.proc
                                  colour = legend.proc
  )
                      ) +
    # geom_line(colour = col) +
    geom_line() +
    geom_abline(intercept = 1, slope = -1, linetype = "dashed", color = "gray") +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(#title = "Precision-Recall Curve",
      x = "Recall",
      y = "Precision",
      linetype = "Models and AUCs",
      colour = "Models and AUCs"
      ) +
    theme(legend.position = c(.45, .25))
  
  if (output == "roc") {
    return(plot.roc)
  } else if (output == "proc"){
    return(plot.proc)
  } else {
  return(list(plot.roc, plot.proc))
  }
}

# OTHERS ======
# Define the color formatting function for numeric columns
color_format <- function(x) {
  cell_spec(x, 'latex',
            color = case_when(x >= 0.8 ~ "#006600",
                              x >= 0.7 ~ "black",
                              x >= 0.6 ~ "orange",
                              !is.na(x) ~ "red",
                              is.na(x) ~ "gray"
            )
  )
}
# color_format <- function(x) {
#   cell_spec(x, 'latex', 
#             color = case_when(x >= 0.8 ~ "#006600",
#                                    x >= 0.7 ~ "#000000",  # Changed from "black"
#                                    x >= 0.6 ~ "#FFA500",  # Changed from "orange" 
#                                    !is.na(x) ~ "#FF0000"  # Changed from "red"
#             )
#   )
# }
background_format <- function(x) {
  cell_spec(x, 'latex', 
            background = case_when(x >= 0.8 ~ "#006600",
                              x >= 0.7 ~ "#000000",  # Changed from "black"
                              x >= 0.6 ~ "#FFA500",  # Changed from "orange" 
                              !is.na(x) ~ "#FF0000"  # Changed from "red"
            ),
            escape = FALSE,
            color = "white"  # Add text color for contrast
  )
}

create_dfs <- function() {
  # prepare df to store results ====
  df.results <- data.frame(grade = as.character(),
                           language = as.character(),
                           riskdef = as.character(),
                           outcometask = as.character(),
                           perc = as.integer(),
                           model = as.character(),
                           auc = as.integer(),
                           set.sens = as.integer(), # pre-set sensitivity
                           opt.cut = as.integer(), # optimal cut-off for entire sample
                           sens.all = as.integer(), # observed sensitivity (entire sample)
                           sens.eo = as.integer(), #                       (EO)
                           sens.ep = as.integer(), #                       (EP)
                           sens.el = as.integer(), #                       (EL)
                           sens.aa = as.integer(), #                       (AA)
                           sens.dis = as.integer(), #                      (DIS)
                           spec.all = as.integer(), # observed specificity (entire sample)
                           spec.eo = as.integer(), #                       (EO)
                           spec.ep = as.integer(), #                       (EP)
                           spec.el = as.integer(), #                       (EL)
                           spec.aa = as.integer(), #                       (AA)
                           spec.dis = as.integer(), #                      (DIS)
                           TP.all = as.integer(), # true positives (entire sample)
                           TP.eo = as.integer(), #                 (EO)
                           TP.ep = as.integer(), #                 (EP)
                           TP.el = as.integer(), #                 (EL)
                           TP.aa = as.integer(), #                 (AA)
                           TP.dis = as.integer(), #                (DIS)
                           FP.all = as.integer(), # false positives (entire sample)
                           FP.eo = as.integer(), #                  (EO)
                           FP.ep = as.integer(), #                  (EP)
                           FP.el = as.integer(), #                  (EL)
                           FP.aa = as.integer(), #                  (AA)
                           FP.dis = as.integer(), #                 (DIS)
                           FN.all = as.integer(), # false negatives  (entire sample)
                           FN.eo = as.integer(), #                   (EO)
                           FN.ep = as.integer(), #                   (EP)
                           FN.el = as.integer(), #                   (EL)
                           FN.aa = as.integer(), #                   (AA)
                           FN.dis = as.integer(), #                  (DIS)
                           TN.all = as.integer(), # true negatives (entire sample)
                           TN.eo = as.integer(), #                (EO)
                           TN.ep = as.integer(), #                (EP)
                           TN.el = as.integer(), #                (EL),
                           TN.aa = as.integer(), #                (AA),
                           TN.dis = as.integer(), #               (DIS),
                           PPV.all = as.integer(),
                           PPV.eo = as.integer(),
                           PPV.ep = as.integer(),
                           PPV.el = as.integer(),
                           PPV.aa = as.integer(),
                           PPV.dis = as.integer(),
                           NPV.all = as.integer(),
                           NPV.eo = as.integer(),
                           NPV.ep = as.integer(),
                           NPV.el = as.integer(),
                           NPV.aa = as.integer(),
                           NPV.dis = as.integer(),
                           BA.all = as.integer(),
                           BA.eo = as.integer(),
                           BA.ep = as.integer(),
                           BA.el = as.integer(),
                           BA.aa = as.integer(),
                           BA.dis = as.integer()
  )
  
  df.preds <- data.frame(student_id = as.character(),
                         pred.risk.optimised80 = as.integer(),
                         pred.risk.optimised85 = as.integer(),
                         pred.risk.optimised90 = as.integer(),
                         pred.risk.optimised95 = as.integer()
  )
  out = list(df.results, df.preds)
}