#' dna_coding
#'
#' Encode DNA sequence strings into position feature matrices
#' @name dna_encoding
#' @param dna_strings DNA sequence vector
#' @return The encoded data frame, with each column representing a nucleotide position.
#' @export
#' @examples
#' # encoding mutiple dna sequences
#' sequences <- c("ATCG", "GCTA", "TTAA", "CCGG")
#' encoded <- dna_encoding(sequences)
#'
#' # to see the results
#' print(encoded)
library(randomForest)
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}




#' predict multiple sequences for m6A sites
#'
#' Predicting m6A sites across multiple sequences in a dataset
#'
#' @param ml_fit Fitted machine learning model
#' @param feature_df Data frame containing features
#' @param positive_threshold Positive class decision threshold, default value 0.5
#' @return Enhanced data frame containing prediction results
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' #load the model
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' # load the example data
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                      package = "m6APrediction"))
#'
#' #Perform batch forecasting

#'results <- prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#'print(results)

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  # Set factor levels to match those used during training.
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  # Encode DNA sequences
  dna_encoded <- dna_encoding(feature_df$DNA_5mer)

  # Only use the DNA position features required by the model (based on model output; do not use nt_pos3).
  dna_features <- dna_encoded[, c("nt_pos1", "nt_pos2", "nt_pos4", "nt_pos5")]

  # Merge all features
  prediction_features <- cbind(
    feature_df[c("gc_content", "RNA_type", "RNA_region", "exon_length",
                 "distance_to_junction", "evolutionary_conservation")],
    dna_features
  )

  # Ensure feature names match those used during model training.
  colnames(prediction_features) <- make.names(colnames(prediction_features))

  # Prediction
  predicted_prob <- predict(ml_fit, newdata = prediction_features, type = "prob")

  # Add predicted results to the original data frame
  feature_df$predicted_m6A_prob <- predicted_prob[, "Positive"]
  feature_df$predicted_m6A_status <- ifelse(feature_df$predicted_m6A_prob > positive_threshold,
                                            "Positive", "Negative")
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}







#' predict single sequence for m6A sites
#'
#' @param ml_fit Fitted machine learning model
#' @param gc_content GC content
#' @param RNA_type RNAt ype
#' @param RNA_region RNA region
#' @param exon_length Exon length
#' @param distance_to_junction distance to junction
#' @param evolutionary_conservation evolutionary conservation
#' @param DNA_5mer DNA 5mer sequence
#' @param positive_threshold Positive class decision threshold, default value 0.5
#' @return Named vectors containing prediction probabilities and states
#' @export
#' @examples
#' #load the model
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#'  #Perform prediction on a single sequence
#' result <- prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGAT",
#'   positive_threshold = 0.5
#' )
#'   print(result)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_row_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  # Set factor levels to match those used during training.
  single_row_df$RNA_type <- factor(single_row_df$RNA_type,
                                   levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  single_row_df$RNA_region <- factor(single_row_df$RNA_region,
                                     levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  # Use the prediction_multiple function for prediction
  result_df <- prediction_multiple(ml_fit, single_row_df, positive_threshold)

  # Extract prediction results and create named vectors
  returned_vector <- c(
    result_df$predicted_m6A_prob,
    as.character(result_df$predicted_m6A_status)
  )

  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")#Complete this function by writing code at the `___`
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
