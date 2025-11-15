README.md
markdown
# m6APrediction: m6A Site Prediction Using Machine Learning

An R package for predicting m6A modification sites in RNA sequences using machine learning models.

## Installation

You can install the development version of m6APrediction from GitHub using:

```r
# Install devtools if not already installed
if (!require("devtools")) {
  install.packages("devtools")
}

# Install m6APrediction from GitHub
devtools::install_github("imrealKen/m6APrediction")


Alternatively, using the remotes package:

r
# Install remotes if not already installed
if (!require("remotes")) {
  install.packages("remotes")
}

# Install m6APrediction from GitHub
remotes::install_github("imrealKen/m6APrediction")
Features
Single Sequence Prediction: Predict m6A sites for individual RNA sequences

Batch Prediction: Process multiple sequences efficiently

Random Forest Model: Utilizes trained machine learning model for accurate predictions

Comprehensive Features: Incorporates various biological features including:

GC content

RNA type and region

Exon length

Evolutionary conservation

Sequence motifs (5-mer)

Quick Start
Single Sequence Prediction
r
library(m6APrediction)

# Load the trained model
rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Predict m6A site for a single sequence
result <- prediction_single(
  ml_fit = rf_model,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "ATCGAT",
  positive_threshold = 0.5
)

print(result)
Batch Prediction
r
library(m6APrediction)

# Load the trained model
rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Load example data
example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", 
                                    package = "m6APrediction"))

# Perform batch prediction
results <- prediction_multiple(
  ml_model = rf_model,
  data = example_data,
  positive_threshold = 0.6
)

head(results)
Model Performance
Our machine learning model demonstrates strong performance in m6A site prediction:

ROC Curve
https://images/ROC_curve.png

The Receiver Operating Characteristic (ROC) curve shows excellent discriminatory power with high AUC values across different cross-validation folds.

Precision-Recall Curve
https://images/PRC_curve.png

The Precision-Recall Curve (PRC) demonstrates robust performance in class-imbalanced scenarios, which is common in biological sequence analysis.

Input Data Format
For batch prediction, your input data should be a CSV file with the following columns:

sequence: RNA sequence

position: Position in the sequence

Additional feature columns as required by the model

Output
The package returns prediction results including:

Prediction: Binary classification (0 = non-m6A site, 1 = m6A site)

Probability: Confidence score for the prediction

Sequence information: Original sequence and position data

Dependencies
randomForest

e1071

caret

Biostrings

seqinr

Citation
If you use m6APrediction in your research, please cite:

Liu, S. (2025). m6APrediction: m6A Site Prediction Using Machine Learning. R package version 1.0.1.

License
This package is licensed under the MIT License. See the LICENSE file for details.

Author
Siyuan Liu - Siyuan.Liu2302@student.xjtlu.edu.cn

Issues and Contributions
For bug reports, feature requests, or contributions, please open an issue on the GitHub repository.


