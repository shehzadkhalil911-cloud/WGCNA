# =============================================================================
# modeling_functions.R
# Functions for training base classifiers
# =============================================================================

#' Train a single model with given parameters
train_base_model <- function(
    formula,
    train_data,
    algorithm,
    tune_grid = NULL,
    train_control = NULL,
    seed = GLOBAL_SEED
) {
    set.seed(seed)
    if (is.null(train_control)) {
        train_control <- caret::trainControl(
            method = "repeatedcv",
            number = CV_FOLDS,
            repeats = CV_REPEATS,
            classProbs = TRUE,
            summaryFunction = twoClassSummary,
            savePredictions = "final",
            sampling = "up",
            allowParallel = USE_PARALLEL
        )
    }

    model <- caret::train(
        form = formula,
        data = train_data,
        method = algorithm,
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE,
        preProcess = NULL,
        verbose = FALSE
    )
    return(model)
}

#' Train Random Forest with ranger
train_rf_model <- function(train_data, train_control, tune_grid = RF_TUNE_GRID) {
    set.seed(GLOBAL_SEED)
    caret::train(
        Class ~ .,
        data = train_data,
        method = "ranger",
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE,
        importance = "permutation",
        verbose = FALSE
    )
}

#' Train Gradient Boosting with gbm
train_gbm_model <- function(train_data, train_control, tune_grid = GB_TUNE_GRID) {
    set.seed(GLOBAL_SEED)
    caret::train(
        Class ~ .,
        data = train_data,
        method = "gbm",
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE,
        verbose = FALSE
    )
}

#' Train SVM with linear kernel
train_svm_model <- function(train_data, train_control, tune_grid = SVM_TUNE_GRID) {
    set.seed(GLOBAL_SEED)
    caret::train(
        Class ~ .,
        data = train_data,
        method = "svmLinear",
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE,
        verbose = FALSE
    )
}

#' Train KNN
train_knn_model <- function(train_data, train_control, tune_grid = KNN_TUNE_GRID) {
    set.seed(GLOBAL_SEED)
    caret::train(
        Class ~ .,
        data = train_data,
        method = "knn",
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE
    )
}

#' Train Neural Network
train_nn_model <- function(train_data, train_control, tune_grid = NN_TUNE_GRID) {
    set.seed(GLOBAL_SEED)
    # PreProcess inside train using training data only
    caret::train(
        Class ~ .,
        data = train_data,
        method = "nnet",
        trControl = train_control,
        tuneGrid = tune_grid,
        metric = "ROC",
        maximize = TRUE,
        preProcess = c("center", "scale"),
        trace = FALSE,
        verbose = FALSE,
        MaxNWts = 10000
    )
}

#' Train all five base models
train_all_base_models <- function(train_data, train_control) {
    models <- list()
    models[["rf"]]   <- train_rf_model(train_data, train_control)
    models[["gbm"]]  <- train_gbm_model(train_data, train_control)
    models[["svm"]]  <- train_svm_model(train_data, train_control)
    models[["knn"]]  <- train_knn_model(train_data, train_control)
    models[["nnet"]] <- train_nn_model(train_data, train_control)
    return(models)
}
