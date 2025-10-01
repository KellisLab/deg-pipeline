# Test file for filter.design function using tinytest
# Tests the design matrix filtering functionality - Focused on duplicates, linear combos, and overspecified

# Load required libraries
library(tinytest)
library(caret)

# Source the functions to test
source("newdeg.R")

# =============================================================================
# Test Suite for filter.design function - 10 Focused Tests
# =============================================================================

# Test 1: Basic functionality with valid design matrix
test_filter_design_basic_valid <- function() {
  # Create a simple valid design matrix
  design <- matrix(c(1,1,0,0, 0,0,1,1, 1,0,1,0), nrow=4, ncol=3)
  colnames(design) <- c("Intercept", "Group1", "Group2")
  
  # Test that function returns TRUE for valid design
  result <- filter.design(design)
  expect_true(result)
}

# Test 2: Duplicate columns detection - error case
test_filter_design_duplicates_error <- function() {
  # Create design matrix with duplicate columns
  design <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1), nrow=4, ncol=3)
  colnames(design) <- c("Group1", "Group1", "Group2")
  
  # Test error when remove_duplicates = FALSE
  expect_error(filter.design(design, remove_duplicates = FALSE))
}

# Test 3: Duplicate columns removal
test_filter_design_duplicates_removal <- function() {
  # Create design matrix with duplicate columns
  design <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1), nrow=4, ncol=3)
  colnames(design) <- c("Group1", "Group1", "Group2")
  
  # Test removal when remove_duplicates = TRUE
  result <- filter.design(design, remove_duplicates = TRUE)
  expect_equal(ncol(result), 2)  # Should have 2 columns after removing duplicate
  expect_equal(colnames(result), c("Group1", "Group2"))
}

# Test 4: Multiple duplicate groups
test_filter_design_multiple_duplicate_groups <- function() {
  # Create design matrix with multiple duplicate groups
  design <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1), nrow=4, ncol=4)
  colnames(design) <- c("Group1", "Group1", "Group2", "Group2")
  
  result <- filter.design(design, remove_duplicates = TRUE)
  expect_equal(ncol(result), 2)  # Should have 2 columns after removing duplicates
  expect_equal(colnames(result), c("Group1", "Group2"))
}

# Test 5: Overspecified design matrix - error case
test_filter_design_overspecified_error <- function() {
  # Create overspecified design (more columns than rows)
  design <- matrix(c(1,1, 0,0, 1,0, 0,1, 1,1), nrow=2, ncol=5)
  colnames(design) <- c("Intercept", "Group1", "Group2", "Group3", "Group4")
  
  # Test error when remove_overspecified = FALSE
  expect_error(filter.design(design, remove_overspecified = FALSE))
}

# Test 6: Overspecified design matrix - removal
test_filter_design_overspecified_removal <- function() {
  # Create overspecified design (more columns than rows)
  design <- matrix(c(1,1, 0,0, 1,0, 0,1, 1,1), nrow=2, ncol=5)
  colnames(design) <- c("Intercept", "Group1", "Group2", "Group3", "Group4")
  
  # Test removal when remove_overspecified = TRUE and remove_duplicates = TRUE
  result <- filter.design(design, remove_overspecified = TRUE, remove_duplicates = TRUE)
  expect_equal(ncol(result), 2)  # Should have 2 columns (same as rows)
  expect_equal(nrow(result), 2)
}

# Test 7: Linear combinations detection - error case
test_filter_design_linear_combos_error <- function() {
  # Create design matrix with linear combinations
  design <- matrix(c(1,1,1,1, 1,1,0,0, 0,0,1,1, 1,1,1,1), nrow=4, ncol=4)
  colnames(design) <- c("Intercept", "Group1", "Group2", "LinearCombo")
  
  # Test error when remove_linear_combos = FALSE
  expect_error(filter.design(design, remove_linear_combos = FALSE))
}

# Test 8: Linear combinations removal
test_filter_design_linear_combos_removal <- function() {
  # Create design matrix with linear combinations
  design <- matrix(c(1,1,1,1, 1,1,0,0, 0,0,1,1, 1,1,1,1), nrow=4, ncol=4)
  colnames(design) <- c("Intercept", "Group1", "Group2", "LinearCombo")
  
  # Test removal when remove_linear_combos = TRUE and remove_duplicates = TRUE
  result <- filter.design(design, remove_linear_combos = TRUE, remove_duplicates = TRUE)
  expect_true(ncol(result) < 4)  # Should have fewer columns
}

# Test 9: Combined issues - duplicates and overspecified
test_filter_design_duplicates_and_overspecified <- function() {
  # Create design with both duplicates and overspecified
  design <- matrix(c(1,1, 1,1, 0,0, 1,0, 0,1), nrow=2, ncol=5)
  colnames(design) <- c("Group1", "Group1", "Group2", "Group3", "Group4")
  
  result <- filter.design(design, remove_duplicates = TRUE, remove_overspecified = TRUE)
  expect_equal(ncol(result), 2)  # Should have 2 columns
}

# Test 10: All removal flags set to TRUE
test_filter_design_all_removals <- function() {
  # Create design with all types of issues
  design <- matrix(c(1,1, 1,1, 0,0, 1,0, 0,1), nrow=2, ncol=5)
  colnames(design) <- c("Group1", "Group1", "Group2", "Group3", "Group4")
  
  result <- filter.design(design, 
                         remove_duplicates = TRUE, 
                         remove_overspecified = TRUE, 
                         remove_linear_combos = TRUE)
  expect_equal(ncol(result), 2)  # Should have 2 columns
}

# =============================================================================
# Test Execution
# =============================================================================

# To run tests, use one of these approaches:
# 1. Run individual tests: test_filter_design_basic_valid()
# 2. Run all tests in environment: run_tests()
# 3. Run tests from file: run_test_file("test.R")

