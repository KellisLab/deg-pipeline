# Simple script to run all filter.design tests (10 focused tests)
library(tinytest)
library(caret)

# Source the functions and tests
source("newdeg.R")
source("test.R")

# Run all tests
cat("Running filter.design tests (10 focused tests)...\n")

# Basic functionality test
test_filter_design_basic_valid()

# Duplicate tests
test_filter_design_duplicates_error()
test_filter_design_duplicates_removal()
test_filter_design_multiple_duplicate_groups()

# Overspecified tests
test_filter_design_overspecified_error()
test_filter_design_overspecified_removal()

# Linear combinations tests
test_filter_design_linear_combos_error()
test_filter_design_linear_combos_removal()

# Combined issues tests
test_filter_design_duplicates_and_overspecified()
test_filter_design_all_removals()

cat("All 10 tests completed!\n") 