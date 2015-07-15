#!/usr/bin/env Rscript

library(optparse, quietly=T)
library(reshape2, quietly=T)
library(parallel, quietly=T)
library(foreach, quietly=T)
library(doMC, quietly=TRUE)

registerDoMC()
options(cores=2)


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
	make_option(c("--order"), action="store", default=NULL,
		help="Preference of Dropping"),
	make_option(c("--pairs"), action="store", default=NULL,
		help="File of pairs of relateds (one pair per line)"),
	make_option(c("--out"), action="store", default="drop_list",
		help="Output file for dropped samples"),
	make_option(c("--exact"), action="store_true", default=F,
		help="Use an exact algorithm (very SLOW!)"),
	make_option(c("--drop-before"), action="store_true", default=F,
		help="Drop samples not found in --order first"))
		
op <- OptionParser(option_list=opt_list)

opts <- parse_args(op)

# Get the appropriate vertex cover function
vc_fn <- greedy_drop
if(opts$exact){
	vc_fn <- exact_drop
}

# get a list of IDs for our adjancency matrix
dat <- read.table(opts$pairs, header=F)
miss <- c()
if(!is.null(opts$order)){
	miss <- read.table(opts$order, header=F)
}

names <- unique(c(miss[,1], dat[,1], dat[,2]))
o <- match(names, miss[,1])
allNames <- c(miss[,1], names[is.na(o)])
if(opts$"drop-before"){
	allNames <- c(names[is.na(o)], miss[,1])
} 

# generate the adjacency matrix
adjMat <- matrix(0, nrow=length(allNames), ncol=length(allNames))
rownames(adjMat) = allNames
rownames(adjMat) = allNames

for(i=1:nrow(dat)){
	adjMat[dat[i,1], dat[i,2]] = 1
	adjMat[dat[i,2], dat[i,1]] = 1
}


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

