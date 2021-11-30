test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Pathway-gene table works", {
    pathway1 <- c("125", "3099", "126")
    pathway2 <- c("5232", "5230", "5162")
    pathList <- list("Path1" = pathway1, "Path2" = pathway2)

    # org.Hs.eg.db is required
    res <- pathwayGeneTab(pathwayList = pathList)
    expect_equal(nrow(res), 6)
})


test_that("Pathway summary generates a correct number of pathways", {

    pathTab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2 * (seq_len(12)), 4, 3)
    rownames(exprsMat) <- pathTab$ENSEMBL
    colnames(exprsMat) <- LETTERS[seq_len(3)]

    res <- pathwaySummary(exprsMat, pathTab, method = "x2")

    expect_equal(nrow(res), 2)
})


test_that("Pathway summary generates a correct number of samples", {

    pathTab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2 * seq_len(1:12), 4, 3)
    rownames(exprsMat) <- pathTab$ENSEMBL
    colnames(exprsMat) <- LETTERS[seq_len(3)]

    res <- pathwaySummary(exprsMat, pathTab, method = "x2")

    expect_equal(ncol(res), 3)
})


test_that("Pathway summary generates correct activity scores", {

    pathTab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2 * (1:12), 4, 3)
    rownames(exprsMat) <- pathTab$ENSEMBL
    colnames(exprsMat) <- LETTERS[1:3]

    res <- pathwaySummary(exprsMat, pathTab, method = "x2")

    expect_equal(res[1, 1], 2.5)
    expect_equal(res[2, 1], 12.5)
})
