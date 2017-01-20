
files <- list.files(path="out", pattern="*_fischertable.txt", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  IPR <- gsub('_fischertable.txt', '', x, ignore.case=FALSE, fixed=FALSE)
  IPR <- gsub('out/', '', IPR, ignore.case=FALSE, fixed=FALSE)

  fishertable <-data.frame()
  fishertable <- read.delim(x, header = FALSE, sep="\t")
  df2$V1 <- (fishertable$V2)
  df2$V2 <- (fishertable$V3)

  fishermatrix <-
  matrix(c(fishertable$V2[1], fishertable$V2[2], fishertable$V3[1], fishertable$V3[2]),
         nrow = 2,
         dimnames = list(Annotation = c("IPR", "Non-IPR"),
                          Contig_Set = c("LS", "Core")))

  results_2_way <- fisher.test(fishermatrix, y = NULL, workspace = 200000, hybrid = FALSE,
              control = list(), or = 1, alternative = "two.sided",
              conf.int = TRUE, conf.level = 0.95,
              simulate.p.value = FALSE, B = 2000)

  results_greater <- fisher.test(fishermatrix, y = NULL, workspace = 200000, hybrid = FALSE,
              control = list(), or = 1, alternative = "greater",
              conf.int = TRUE, conf.level = 0.95,
              simulate.p.value = FALSE, B = 2000)

  results_less <- fisher.test(fishermatrix, y = NULL, workspace = 200000, hybrid = FALSE,
              control = list(), or = 1, alternative = "less",
              conf.int = TRUE, conf.level = 0.95,
              simulate.p.value = FALSE, B = 2000)


  outline <- paste(IPR, results_2_way$p.value, results_greater$p.value, results_less$p.value,sep = "\t")
  print(IPR)
write(outline, file="results.txt", append=T)
})



#get config options
library(optparse)
library(VennDiagram, lib.loc="/home/armita/R-packages/")
opt_list = list(
    make_option("--inp", type="character", help="tab seperated contingency table"),
    make_option("--out", type="character", help="output file to write the results to")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out
