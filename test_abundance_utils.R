source('abundanceFilteringUtils.R')

context('Most abundant genes')

sites <- data_frame(
    estAbund=sample(seq(1, 5), 10, replace=TRUE),
    # inconsistent with estAbund but enough for testing
    estAbundProp=sample(c(0.1, 0.2, 0.3,0.4, 0.5), 10, replace=TRUE), 
    nearest_refSeq_gene=sample(c('LMO2', 'TET2', 'BLA', 'XYZ'), 10, replace=TRUE)
) 

NUM_GENES <- 2
result <- getMostAbundantGenes(sites, NUM_GENES)

test_that("return list with 2 elements", {
    expect_is(result, 'list') 
    expect_named(result, c('abundanceCutoff', 'topGenes'))
})

test_that("return correct number of genes", {
    expect_is(result$topGenes, 'character') 
    expect_equal(length(result$topGenes), NUM_GENES)
    expect_true(all(result$topGenes %in% sites$nearest_refSeq_gene))
})

context('Mask low abundance genes')

white_list <- c('LMO2', 'XYZ')
sites_masked <- maskGenes(sites, white_list)

test_that("adds column to dataframe", {
   expect_is(sites_masked, 'data.frame')
   expect_true('maskedRefGeneName' %in% names(sites_masked))
})

test_that("only white list is kept", {
   expect_true(all(white_list %in% sites_masked$maskedRefGeneName))
   expect_false(all(setdiff(sites$nearest_refSeq_gene, white_list) 
        %in% sites_masked$maskedRefGeneName))
})
