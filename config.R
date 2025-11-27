# Configuration file for Silica Spatial Controls Analysis
#
# SETUP INSTRUCTIONS:
# 1. Create a working directory for your analysis
# 2. Download source data from Zenodo (see README)
# 3. Place raw data in DATA_DIR/data/
# 4. Update DATA_DIR below to point to your working directory
# 5. Run scripts in order: Step1 → Step2 → Step3
#
# The workflow will create this structure:
#   DATA_DIR/
#   ├── data/                   # You put Zenodo downloads here
#   ├── harmonization_files/    # Step 1 creates intermediate data
#   ├── models/                 # Step 2 creates RF models & SHAP values
#   └── figures/                # Step 3 creates publication figures
#       ├── png/
#       └── pdf/

# ==============================================================================
# USER CONFIGURATION - EDIT THIS SECTION
# ==============================================================================

# Path to your analysis working directory
# Example: "/Users/yourname/my_sisyn_analysis"
# Windows example: "C:/Users/yourname/Documents/my_sisyn_analysis"
DATA_DIR <- "/path/to/your/analysis/directory"

# ==============================================================================
# DERIVED PATHS - DO NOT EDIT BELOW THIS LINE
# ==============================================================================

# Input data directories
RAW_DATA_DIR <- file.path(DATA_DIR, "data")
HARMONIZATION_DIR <- file.path(DATA_DIR, "harmonization_files", "inputs")
MODEL_DIR <- file.path(DATA_DIR, "models")
DRIVERS_SPLIT_FILE <- file.path(HARMONIZATION_DIR, "AllDrivers_recent30_split.csv")

# Output directories
FIGURES_DIR <- file.path(DATA_DIR, "figures")
OUTPUT_PNG_DIR <- file.path(FIGURES_DIR, "png")
OUTPUT_PDF_DIR <- file.path(FIGURES_DIR, "pdf")

# Create output directories if they don't exist
if (!dir.exists(OUTPUT_PNG_DIR)) dir.create(OUTPUT_PNG_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_PDF_DIR)) dir.create(OUTPUT_PDF_DIR, recursive = TRUE)

# Validate configuration
if (!dir.exists(DATA_DIR)) {
  stop(paste0(
    "\nERROR: DATA_DIR does not exist: ", DATA_DIR,
    "\n\nPlease:\n",
    "1. Create your working directory\n",
    "2. Update DATA_DIR in config.R to point to it\n"
  ))
}

# Check workflow progress (informational only)
step1_complete <- dir.exists(HARMONIZATION_DIR)
step2_complete <- dir.exists(MODEL_DIR)

message("Configuration loaded successfully!")
message("Data directory: ", DATA_DIR)
if (!step1_complete) {
  message("Note: Run Step 1 scripts to create harmonization_files/")
}
if (!step2_complete) {
  message("Note: Run Step 2 scripts to create models/")
}
