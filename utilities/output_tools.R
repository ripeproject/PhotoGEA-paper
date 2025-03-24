# Define output directories
OUTPUT_DIR <- 'output'
TABLE_DIR <- 'tables'
FIGURE_DIR <- 'figures'
RDATA_DIR <- 'rdata'
SPREADSHEET_DIR <- 'for_spreadsheet_fits'

# Make sure the output directories exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

if (!dir.exists(file.path(OUTPUT_DIR, TABLE_DIR))) {
  dir.create(file.path(OUTPUT_DIR, TABLE_DIR))
}

if (!dir.exists(file.path(OUTPUT_DIR, FIGURE_DIR))) {
  dir.create(file.path(OUTPUT_DIR, FIGURE_DIR))
}

if (!dir.exists(file.path(OUTPUT_DIR, RDATA_DIR))) {
  dir.create(file.path(OUTPUT_DIR, RDATA_DIR))
}

if (!dir.exists(file.path(OUTPUT_DIR, SPREADSHEET_DIR))) {
  dir.create(file.path(OUTPUT_DIR, SPREADSHEET_DIR))
}
