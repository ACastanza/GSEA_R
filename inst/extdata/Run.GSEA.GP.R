# GSEA 2.0 -- Gene Set Enrichment Analysis / Broad Institute Executable R script
# to run GSEA Analysis
suppressMessages(suppressWarnings(install.packages("getopt", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("getopt")))

suppressMessages(suppressWarnings(install.packages("optparse", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("optparse")))

suppressMessages(suppressWarnings(install.packages("devtools", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("devtools")))

suppressMessages(suppressWarnings(install.packages("dplyr", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("dplyr")))

suppressMessages(suppressWarnings(install.packages("BiocManager", repos = "https://cloud.r-project.org/", 
quiet = TRUE)))
suppressMessages(suppressWarnings(BiocManager::install("DESeq2", quiet = TRUE)))
suppressMessages(suppressWarnings(library("DESeq2")))

arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(suppressWarnings(library("utils")))
suppressMessages(suppressWarnings(library("tools")))

suppressMessages(suppressWarnings(install_github("ACastanza/GSEA_R@nwGSEA")))
suppressMessages(suppressWarnings(library("GSEA")))

option_list <- list(
make_option("--res", dest = "expression.dataset"),
make_option("--gmx_list", dest = "gene.sets.database"),
make_option("--nperm", dest = "number.of.permutations"),
make_option("--cls", dest = "phenotype.labels", default = NULL),
make_option("--permute", dest = "permutation.type"),
make_option("--collapse", dest = "collapse.dataset"),
make_option("--chip", dest = "chip.platform.file"),
make_option("--output_file_name", dest = "output.file.name", default = "gsea_result"),
make_option("--scoring_scheme", dest = "scoring.scheme"),
make_option("--metric", dest = "metric.for.ranking.genes"),
make_option("--set_max", dest = "max.gene.set.size"),
make_option("--set_min", dest = "min.gene.set.size"),
make_option("--mode", dest = "collapsing.mode.for.probe.sets.with.more.than.one.match"),
make_option("--rnd_type", dest = "randomization.mode"),
make_option("--plot_top_x", dest = "plot.graphs.for.the.top.sets.of.each.phenotype"),
make_option("--rnd_seed", dest = "random.seed"),
make_option("--save_rnd_lists", dest = "save.random.ranked.lists"),
make_option("--zip_report", dest = "create.zip", default = TRUE),
make_option("--network", dest = "usenetwork", default = TRUE),
make_option("--weight_mode", dest = "scoretype", default = "strength"),
make_option("--symbol_versions", dest = "msigdbversion", default = "7.2")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE, 
 args = arguments)$options

if (as.logical(opt$collapse) == TRUE) {
 inputchip <- opt$chip.platform.file
  collapsemode <- opt$collapsing.mode.for.probe.sets.with.more.than.one.match
} else {
 inputchip <- "NOCHIP"
 collapsemode <- "NOCOLLAPSE"
}

if (file_ext(opt$expression.dataset) == "rnk") {
 permutation <- "gene.labels"
 rankmethod <- "preranked"
} else {
 permutation <- opt$permutation.type
 rankmethod <- "GSEA"
 opt$phenotype.labels
}

nperms <- opt$number.of.permutations
maxsize <- opt$max.gene.set.size
minsize <- opt$min.gene.set.size

 usenetwork <- askYesNo("Use Network Weighted Ranking? ")

# if (usenetwork == TRUE) { nettype <- menu(c('Strength', 'Eigen Centrality'),
# graphics = FALSE, title = 'Weight Metric for Network Score') if (nettype == 1)
# { scoretype <- 'strength' } else if (nettype == 2) { scoretype <- 'centrality'
# } }

cat("\n")

outdir <- paste0(opt$output.file.name, "_reports")
dir.create(outdir)

if(!is.na(as.numeric(opt$random.seed))) {
rand.seed <- opt$random.seed
} else if (opt$random.seed == "timestamp"){
rand.seed <- as.POSIXct(Sys.time())
} else {
message("Invalid random seed value, defaulting to system timestamp")
rand.seed <- as.POSIXct(Sys.time())
}

db.source <- as.character(readLines(opt$gene.sets.database))
db.content <- unname(sapply(db.source, readLines))
gs.db <- unlist(db.content, use.names=FALSE)
names(gs.db)[1] <- paste0(basename(db.source), sep=" ", collapse=" ")

GSEA(
# Input/Output Files :-------------------------------------------------------------------------------
 input.ds = opt$expression.dataset,                    # Input gene expression dataset file in GCT format
 input.cls = opt$phenotype.labels,                  # Input class vector (phenotype) file in CLS format
 gs.db = gs.db,                          # Gene set database in GMT format
 input.chip = inputchip,               # CHIP File
 output.directory      = outdir,        # Directory where to store output and results (default: "")
#  Program parameters :-------------------------------------------------------------------------------
 doc.string            = opt$output.file.name,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
 reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
 nperm                 = as.integer(nperms),            # Number of random permutations (default: 1000)
 weighted.score.type   =  as.numeric(opt$scoring.scheme),              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
 topgs                 = as.integer(opt$plot.graphs.for.the.top.sets.of.each.phenotype),              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
 gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
 reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
 preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
 random.seed           = as.integer(rand.seed),            # Random number generator seed. (default: 123456)
 perm.type             = as.integer(opt$randomization.mode),               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
 fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
 replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
 collapse.dataset      = as.logical(opt$collapse.dataset), # Collapse dataset to gene symbols using a user provided chip file (default: F)
 collapse.mode         = collapsemode,
 save.intermediate.results = as.logical(opt$save.random.ranked.lists),           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
 use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
 gsea.type = rankmethod,                     # Select Standard GSEA (default) or preranked
 rank.metric = opt$metric.for.ranking.genes,
 network = opt$usenetwork,
 score.type = opt$scoretype,
 msigdbversion = opt$msigdbversion
 )
#----------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
   directory           = outdir,        # Directory where to store output and results (default: "")
   topgs = as.integer(opt$plot.graphs.for.the.top.sets.of.each.phenotype), # number of top scoring gene sets used for analysis
   height = 16,
   width = 16,
   gsea.type = rankmethod,
   doc.string = opt$output.file.name
)

if(as.logical(opt$create.zip) == TRUE) {
files2zip <- dir(outdir, full.names = TRUE)
zip(zipfile = outdir, files = files2zip)
}
