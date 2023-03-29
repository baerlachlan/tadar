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
gt_out <- data.frame(
    sample1 = "0/1",
    sample2 = "1/1",
    sample3 = "0/1",
    row.names = c("SNP1")
)

test_that("correct output returned by unphaseGT", {
    expect_equal(unphaseGT(gt_mat), gt_out)
    expect_equal(unphaseGT(gt_df), gt_out)
})
