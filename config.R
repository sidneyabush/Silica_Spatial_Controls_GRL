# Configuration file for Silica Spatial Controls Analysis
#
# SETUP INSTRUCTIONS:
# 1. Create a working directory for your analysis
# 2. Download the source data (see README for Zenodo/USGS links)
# 3. Update DATA_DIR below to point to your working directory
# 4. Run scripts in order: Step1 → Step2 → Step3
#
# After running all steps, your directory will contain:
#   - harmonization_files/inputs/  (created by Step 1)
#   - Final_Models/                (created by Step 2)
#   - GRL_Materials/Final_Figures/ (created by Step 3)

# ==============================================================================
# USER CONFIGURATION - EDIT THIS SECTION
# ==============================================================================

# Path to the root directory where you extracted the Zenodo data
# Example: "/Users/yourname/Downloads/SiSyn_Data"
# Windows example: "C:/Users/yourname/Documents/SiSyn_Data"
DATA_DIR <- "/path/to/your/data/directory"

# ==============================================================================
# DERIVED PATHS - DO NOT EDIT BELOW THIS LINE
# ==============================================================================

# Input data directories
HARMONIZATION_DIR <- file.path(DATA_DIR, "harmonization_files", "inputs")
MODEL_OUTPUT_DIR <- file.path(DATA_DIR, "Final_Models")
MODEL_PERFORMANCE_DIR <- file.path(DATA_DIR, "model_performance")
DRIVERS_SPLIT_FILE <- file.path(HARMONIZATION_DIR, "AllDrivers_recent30_split.csv")

# Output directories
OUTPUT_PNG_DIR <- file.path(DATA_DIR, "GRL_Materials", "Final_Figures", "PNG")
OUTPUT_PDF_DIR <- file.path(DATA_DIR, "GRL_Materials", "Final_Figures", "PDF")

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
step2_complete <- dir.exists(MODEL_OUTPUT_DIR)

message("Configuration loaded successfully!")
message("Data directory: ", DATA_DIR)
if (!step1_complete) {
  message("Note: Run Step 1 scripts to create harmonization_files/")
}
if (!step2_complete) {
  message("Note: Run Step 2 scripts to create Final_Models/")
}
