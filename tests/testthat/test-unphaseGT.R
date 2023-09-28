gt_mat <- matrix(
    data = c("0|1", "1|1", "0/1"), nrow = 1, ncol = 3,
    dimnames = list("SNP1", c("sample1", "sample2", "sample3"))
)
gt_df <- data.frame(
    sample1 = "0|1",
    sample2 = "1|1",
    sample3 = "0/1",
    row.names = c("SNP1")
)
gt_unphased <- matrix(
    data = c("0/1", "1/1", "0/1"), nrow = 1, ncol = 3,
    dimnames = list("SNP1", c("sample1", "sample2", "sample3"))
)

test_that("unphaseGT returns correct output", {
    expect_equal(unphaseGT(gt_mat), gt_unphased)
    expect_equal(unphaseGT(gt_df), gt_unphased)
})
