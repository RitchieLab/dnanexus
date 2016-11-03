#!/usr/bin/env Rscript

library(optparse, quietly=T)
library(SNPRelate, quietly=T)
library(reshape2, quietly=T)
library(parallel, quietly=T)
library(foreach, quietly=T)
library(doMC, quietly=TRUE)

registerDoMC()
options(cores=detectCores())


# given a matrix and a connected component (sorted by missingness), drop
# the minimal vertex cover for the connected component.
exact_drop <- function(mat, cc){
	for(n in 1:(length(cc) -1) ){
		all_comb <- combn(cc, n)
		for(j in 1:ncol(all_comb)){
			# find if this is a vertex cover
			rem_samp <- setdiff(cc, all_comb[,j])
			#print(rem_samp)
			if(sum(mat[rem_samp, rem_samp]) < 1){
				#print(all_comb[,j])
				return(all_comb[,j])
			}
		}
	}
}	

# given a matrix and a connected component (sorted by missingness), drop
# a vertex cover for the connected component found by a greedy algorithm
greedy_drop <- function(mat, cc){

	rem_samp <- cc
	while(sum(mat[rem_samp, rem_samp]) > 0){
		# find the index of the vertex of max. degree
		idx <- which.max(apply(mat[rem_samp, rem_samp],1,sum))
		rem_samp <- setdiff(rem_samp, idx)
	}
	
	return(setdiff(cc, rem_samp))
}

opt_list <- list(
	make_option(c("--bed"), action="store", default=NULL,
		help="BED file"),
	make_option(c("--bim"), action="store", default=NULL,
		help="BIM file"),
	make_option(c("--fam"), action="store", default=NULL,
		help="FAM file"),
	make_option(c("--missing"), action="store", default=NULL,
		help="Sample Missingness file"),
	make_option(c("--thresh"), action="store", default=0.125, type="numeric",
		help="Kinship coefficient threshold"),
	make_option(c("--rsq"), action="store", default=0.1, type="numeric",
		help="LD pruning R-squared threshold"),
	make_option(c("--out"), action="store", default="drop_list",
		help="Output file for dropped samples"),
	make_option(c("--exact"), action="store_true", default=F,
		help="Use an exact algorithm (very SLOW!)"))
		
op <- OptionParser(option_list=opt_list)

opts <- parse_args(op)

# TODO: Do some sanity checking here


# Get the appropriate vertex cover function
vc_fn <- greedy_drop
if(opts$exact){
	vc_fn <- exact_drop
}

# we need a file for the GDS file
gdsfn <- tempfile()

# convert the bed/bim/fam to gds format
snpgdsBED2GDS(opts$bed, opts$fam, opts$bim, gdsfn, compress.annotation="")

genofile <- snpgdsOpen(gdsfn, readonly=F)
# convert the IDs to "FID IID"
famd <- read.table("ACTG_phaseI_all.fam", header=F, colClasses=c('character','character'))
sample_data <- paste(famd[,1], famd[,2], sep=" ")
add.gdsn(genofile, "sample.id", sample_data, closezip=T, replace=T)

# Now, let's do a lttle LD pruning to make sure everything is independent (and speed up IBD)

snpset <- snpgdsLDpruning(genofile, ld.threshold=sqrt(opts$rsq),method="r")

# calculate the kinship coefficient for each pair
ibd <- snpgdsIBDMoM(genofile, snp.id=unlist(snpset), maf=0.05, missing.rate=0.05, num.thread=detectCores())
ibd_coeff <- snpgdsIBDSelection(ibd)

# generate an adjacency matrix for related pairs
dat <- acast(ibd_coeff, ID1 ~ ID2)
dat[is.na(dat)] = 0
dat <- dat + t(dat)
dat[dat < opts$thresh] = 0
dat[dat > 0] = 1

# sort everything by missingness here!
if(!is.null(opts$missing)){
	miss_data <- read.table(opts$missing, header=F, colClasses=c('character','character','numeric')	
	
	miss_data$ID <- paste(miss_data[,1], miss_data[,2], sep=" ")
	miss_data <- miss_data[order(miss_data[,3]),]
	
	# OK, we're now ordered by missingness, let's add anything not having a missing
	# rate to the end of the list
	id_list <- c(as.character(miss_data$ID), setdiff(rownames(dat), miss_dat$ID))
	
	# now, rearrange the "dat" rows and columns
	o <- match(rownames(dat), id_list)
	dat <- dat[o,o]	
}

#print(dat)

# OK, we now have an adjacency matrix, let's find connected components!
conn_comp <- list()
curr_comp_num = 0
comp_num <- rep(0, nrow(dat))
comp_queue <- c()
for (idx in 1:nrow(dat)){
	curr_comp <- c()
	if(comp_num[idx] == 0){
		curr_comp_num = curr_comp_num + 1
		comp_queue <- append(comp_queue, idx)

		while(length(comp_queue) > 0){
			curr_idx = comp_queue[1]
			comp_queue <- comp_queue[-1]
			if(comp_num[curr_idx] == 0){
				curr_comp <- append(curr_comp, curr_idx)
			}
			comp_num[curr_idx] = curr_comp_num

			for(adj_idx in 1:nrow(dat)){
				if(comp_num[adj_idx] == 0 & dat[adj_idx,curr_idx] != 0){
					comp_queue <- append(comp_queue, adj_idx)
				}
			}
		}
		
		conn_comp[[length(conn_comp) + 1]] = sort(curr_comp)
	}
}

#print(conn_comp)
#print(comp_num)

to_drop <- c()

# OK, now I have connected components, now time to get a minimal vertex cover
# We'll do this the easy (and naieve way)
todrop <- foreach(cc in conn_comp[lapply(conn_comp, length)>1], .combine="c") %dopar%{
	return(vc_fn(dat, cc))
}

write.table(rownames(dat)[to_drop], file=opts$out, row.names=F, col.names=F, quote=F)

