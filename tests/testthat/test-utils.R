test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Pathway-gene table works", {
    Pathway1  <- c( "125","3099","126")
    Pathway2  <- c( "5232","5230","5162")
    path.list <- list("Path1" = Pathway1, "Path2" = Pathway2)

    # org.Hs.eg.db is required
    res <- pathwayGeneTab(pathwayList = path.list)
    expect_equal(nrow(res), 6)
})


test_that("Pathway summary generates a correct number of pathways", {

    path_tab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2*(seq_len(12)),4,3)
    rownames(exprsMat) <- path_tab$ENSEMBL
    colnames(exprsMat) <- LETTERS[seq_len(3)]

    res <- pathwaySummary(exprsMat, path_tab, method = "x2")

    expect_equal(nrow(res), 2)
})


test_that("Pathway summary generates a correct number of samples", {

    path_tab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2*seq_len(1:12),4,3)
    rownames(exprsMat) <- path_tab$ENSEMBL
    colnames(exprsMat) <- LETTERS[seq_len(3)]

    res <- pathwaySummary(exprsMat, path_tab, method = "x2")

    expect_equal(ncol(res), 3)
})


test_that("Pathway summary generates correct activity scores", {

    path_tab <- tibble::tribble(
        ~Pathway, ~ENTREZID,  ~ENSEMBL,
        "Path1",   "125",      "ENSG00000196616",
        "Path1",   "3099",     "ENSG00000159399",
        "Path2",   "5230",     "ENSG00000102144",
        "Path2",   "5162",     "ENSG00000168291"
    )

    exprsMat <- matrix(2*(1:12),4,3)
    rownames(exprsMat) <- path_tab$ENSEMBL
    colnames(exprsMat) <- LETTERS[1:3]

    res <- pathwaySummary(exprsMat, path_tab, method = "x2")

    expect_equal(res[1,1], 2.5)
    expect_equal(res[2,1], 12.5)
})
