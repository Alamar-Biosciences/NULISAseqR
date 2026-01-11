# Build the NULISAseqR bookdown guide
# Run this script from the package root directory

cat("Building NULISAseqR Guide Book...\n\n")

# Check if we're in the right directory
if (!dir.exists("vignettes/book")) {
  stop("Error: vignettes/book/ directory not found!\n",
       "Please check your working directory.")
}

# Load package in development mode
cat("Step 1: Loading package...\n")
devtools::load_all()

# Set working directory to vignettes/book
cat("Step 2: Navigating to book directory...\n")
original_wd <- getwd()
setwd("vignettes/book")

# Clean previous build
cat("Step 3: Cleaning previous build...\n")
if (file.exists("_main.Rmd")) file.remove("_main.Rmd")
bookdown::clean_book(TRUE)

# Build the book
cat("Step 4: Building book (this may take a few minutes)...\n")
tryCatch({
  bookdown::render_book(
    "index.Rmd",
    output_format = "bookdown::gitbook",
    clean = TRUE,
    quiet = FALSE
  )
  
  cat("\n✓ Book built successfully!\n\n")
  
  # Return to package root
  setwd(original_wd)
  
  # Display results
  cat("===============================================================================\n")
  cat("BUILD COMPLETE\n")
  cat("===============================================================================\n\n")
  
  cat("Output Location:\n")
  cat("  ", file.path(getwd(), "vignettes/book/_book/index.html"), "\n\n")
  
  # Open the book
  cat("Opening book in browser...\n")
  browseURL("vignettes/book/_book/index.html")
  
}, error = function(e) {
  setwd(original_wd)
  cat("\n✗ Build failed with error:\n")
  cat("  ", conditionMessage(e), "\n\n")
})
