library(pcalg)
library(bnlearn)
library(igraph)
library(GoodmanKruskal)

setwd("G:/Mandar/NCSU/Independent Study/Second chapter/Experiments/Real-world data")
network <- "alarm"
sample_size <- 1000
r_data_file <- paste(network, '_', sample_size, '.RData', sep='')
load(paste(network, '.rda', sep=''))
load(r_data_file)

disc_data <- data[[1]]
num_vars <- ncol(disc_data)
bn_amat <- amat(bn)
colnames(bn_amat) <- row.names(bn_amat) <- 1:num_vars
bn_graph <- igraph.to.graphNEL(graph_from_adjacency_matrix(bn_amat, mode="directed"))
unique_vals <- apply(disc_data, 2, function(x) length(unique(x)))

# Get top-K associated variables and refine them
assoc_obj <- getAssocMeasures(disc_data)
cv_mat <- assoc_obj[[1]]
gk_mat <- assoc_obj[[2]]

top_k_list <- getTopK(cv_mat)
top_k_cv <- top_k_list[[1]]
G_cv <- top_k_list[[2]]

top_k_list <- getTopK(gk_mat)
top_k_gk <- top_k_list[[1]]
G_gk <- top_k_list[[2]]

# skel_cv <- buildCausalGraph(alarm_data, top_k_cv)
# skel_gk <- buildCausalGraph(alarm_data, top_k_gk)

pc_obj <- getPC(disc_data, top_k_cv, G_cv)
G_cv <- pc_obj[[1]]
new_top_k_cv <- pc_obj[[2]]

G_skel_cv <- mergeTopk(new_top_k_cv)
idx <- which(G_skel_cv==TRUE)
cv_amat <- matrix(0, num_vars, num_vars)
cv_amat[idx] <- 1
cv_graph <- igraph.to.graphNEL(graph_from_adjacency_matrix(cv_amat, mode="undirected"))

compareGraphs(cv_graph, dag[[1]])

pc_obj <- getPC(alarm_data, top_k_gk, G_gk)
G_gk <- pc_obj[[1]]
new_top_k_gk <- pc_obj[[2]]
G_skel_gk <- mergeTopk(new_top_k_gk)
idx <- which(G_skel_gk==TRUE)
gk_amat <- matrix(0, num_vars, num_vars)
gk_amat[idx] <- 1
gk_graph <- igraph.to.graphNEL(graph_from_adjacency_matrix(gk_amat, mode="undirected"))

pc_obj[[j]] <- pc(suffStat=list(dm=data[[j]], nlev=unique_vals, adaptDF=FALSE), indepTest=disCItest, alpha=alpha, p=ncol(data[[j]]), skel.method="original", verbose=F)
pc_stable_obj[[j]] <- pc(suffStat=list(dm=data[[j]], nlev=unique_vals, adaptDF=FALSE), indepTest=disCItest, alpha=alpha, p=ncol(data[[j]]), skel.method="stable", verbose=F)

# Get the true parent-child set
true_pc <- list()
vars <- 1:num_vars

amat <- wgtMatrix(dag[[1]])
for(i in 1:num_vars){
  true_pc[[i]] <- c(which(amat[i, ]!=0), which(amat[, i]!=0))
}

# Perturb the true parent-child set with different threshold of noise
noise_th <- c(0.1, 0.2, 0.3, 0.4, 0.5)
num_iter <- 10
pc_size <- unlist(lapply(true_pc, length))
num_pc_vars <- sum(pc_size)
pc_dist <- table(pc_size[pc_size!=0])
pc_size_frac <- array(0, length(pc_dist))
for(i in 1:length(pc_dist)){
    pc_size_frac[i] <- round(i*pc_dist[i]/num_pc_vars, 2)
}

for(n in noise_th){
    noise_var_dist <- round(pc_dist*1:length(pc_dist)*n)
    for(itr in 1:num_iter){
        for(s in 1:length(pc_dist)){
            vars_pc_s <- which(pc_size==s)
            if(length(vars_pc_s >= noise_var_dist[s])){
                vars_1 <- sample(vars_pc_s, noise_var_dist[s])
                true_pc <- perturbPC(vars_1, true_pc, noise_var_dist[s], order=s)
            }else{

            }
        }
    }
}

# Pick top-K variables
getTopK <- function(mat){
    top_k <- list()
    G <- matrix(FALSE, nrow=nrow(mat), ncol=nrow(mat))
    for(i in 1:ncol(mat)){
        sort_vals <- sort(mat[i, ], decreasing=TRUE, index.return=T)
        idx <- which(sort_vals$x!=0)
        sort_vals$x <- sort_vals$x[idx]
        sort_vals$ix <- sort_vals$ix[idx]
        k <- find.maximum.distance.point(sort_vals$x)
        k <- length(sort_vals$x)
        top_k[[i]] <- sort_vals$ix[1:k]
        G[i, sort_vals$ix[1:k]] <- TRUE
    }
    return(list(top_k, G))
}

find.maximum.distance.point <- function(y, x=1:length(y)){   
    allCoord = rbind(y, x)

    firstPoint = allCoord[,1]
    lineVec = allCoord[,length(y)] - firstPoint
    lineVecN = lineVec / sqrt(sum(lineVec^2))

    vecFromFirst = allCoord - firstPoint
    scalarProduct = lineVecN %*% vecFromFirst

    vecFromFirstParallel = t(scalarProduct) %*% lineVecN
    vecToLine = t(vecFromFirst) - vecFromFirstParallel
    distToLine = sqrt(rowSums(vecToLine^2,2))
    which.max(distToLine)
}

getAssocMeasures <- function(org_data){
    data <- data.frame(org_data)
    num_vars <- ncol(data)
    bool_mat <- matrix(1, nrow=num_vars, ncol=num_vars)
    diag(bool_mat) <- 0
    rem_edges <- which(lower.tri(bool_mat)==1, arr.ind=TRUE)
    for(i in 1:num_vars){
        data[, i] <- as.factor(data[, i])
    }

    # Get Cramer's V and Goodman-Krusal association measures
    cv_mat <- matrix(0, num_vars, num_vars)
    gk_mat <- matrix(0, num_vars, num_vars)

    for(e in 1:nrow(rem_edges)){
      i <- rem_edges[e, 1]
      j <- rem_edges[e, 2]
      chi_obj <- chisq.test(data[, i], data[, j])
      # if(chi_obj$p.value < 0.05){
        cv_mat[i, j] <- cv_mat[j, i] <- sqrt(chi_obj$statistic/(length(data[, i])*(min(length(unique(data[, i])), length(unique(data[, j])))-1)))
      # }
      gk_val <- GKtau(data[, i], data[, j])
      gk_mat[i, j] <- gk_mat[j, i] <- max(gk_val$tauxy, gk_val$tauyx)

    }
    return(list(cv_mat, gk_mat))
}

cv.test <- function(x,y) {
  CV <- sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(as.numeric(CV))
}

# Refine top-K variables
getPC <- function(data, top_k, G, cor_mat){
    if(is.null(cor_mat)){
        unique_vals <- apply(data, 2, function(x) length(unique(x)))
    }
    rem_edges <- which(G==TRUE, arr.ind=TRUE)
    num_vars <- ncol(data)
    seq_p <- seq(num_vars)
    pVals <- matrix(-Inf, nrow=num_vars, ncol=num_vars)
    sepset <- lapply(seq_p, function(.) vector("list", num_vars))
    l <- 1
    for(i in 18:18){
        adj_i <- pc_i <- top_k[[i]]
        j_vars <- rev(pc_i)
        l <- 0
        non_adj <- c()
        repeat{
            pc_i <- top_k[[i]]
            sep_set <- matrix(0, nrow=l, ncol=num_vars)
            for(j in j_vars){
                pc_j <- top_k[[j]]
                # S_1 <- intersect(pc_i, pc_j)
                # S_2 <- setdiff(pc_i, c(S_1, j))
                # print(cat('pc_i', pc_i, ' l: ', l, ' S:', union(S_1, S_2)))
                if(length(pc_i)-1>=l){
                    if(l==1){
                        # S <- matrix(union(S_1, S_2), nrow=1)
                        S <- matrix(pc_i, nrow=1)
                    }else{
                        # S <- combn(union(S_1, S_2), l)
                        S <- combn(pc_i, l)
                    }
                    for(s in 1:ncol(S)){
                        if(j %in% S[, s]) next
                        if(is.null(cor_mat)){
                            suffStat <- list(dm=data, nlev=unique_vals, adaptDF=FALSE)
                            p_val <- disCItest(i, j, S=S[, s], suffStat)                            
                        }else{
                            suffStat <- list(C=cor(data), n=nrow(data))
                            p_val <- gaussCItest(i, j, S[, s], suffStat)                            
                        }
                        if(p_val > 0.05){
                            print(cat('First test, i: ',i, ' j:', j, ' S:', S[, s], ' p_val: ', p_val))
                            pdsep_res <- checkDSep(data, sepset, i, j, S[, s], non_adj)
                            pdsep <- pdsep_res[[1]]
                            sepset <-pdsep_res[[2]]
                            if(pdsep){
                                non_adj <- c(non_adj, j)
                                pVals[i, j] <- pVals[j, i] <- p_val
                                sepset[[i]][[j]] <- S[, s]
                                pc_i <- pc_i[!pc_i %in% j]
                                G[i, j] <- G[j, i] <- FALSE
                                break
                            }
                        }
                    }
                }
            }
            # Check for robust d-separators
            # if(length(non_adj)>0){
            #     for(j in non_adj){
            #         # Check if sep_set[, j] d-separates i from its top-k var set
            #         pdsep <- checkDSep(data, sepset, i, j, non_adj)
            #         # best case when the true parent-child set is in sep_set
            #         print(cat('i: ', i, ' j:', j, ' pdsep:', pdsep))
            #         if(all(sepset[[i]][[j]] %in% adj_i) && pdsep){
            #             # print(cat('Best case, i: ',i, ' j:', j, ' S:', sep_set[, j]))
            #             pc_i <- pc_i[!pc_i %in% j]
            #             G[i, j] <- G[j, i] <- FALSE
            #         }else{
            #             pc_j <- top_k[[j]]
            #             # S_1 <- intersect(adj_i, pc_j)
            #             # S_2 <- setdiff(adj_i, c(S_1, j, sep_set[, j]))
            #             S_1 <- setdiff(adj_i, c(j, sepset[[i]][[j]]))
            #             # print(cat('l: ', l, ' S:', union(S_1, S_2)))
            #             if(length(S_1)>=l){
            #                 if(l==1){
            #                     # S <- matrix(union(S_1, S_2), nrow=1)
            #                     S <- matrix(S_1, nrow=1)
            #                 }else{
            #                     # S <- combn(union(S_1, S_2), l)
            #                     S <- combn(S_1, l)
            #                 }
            #                 for(s in 1:ncol(S)){
            #                     suffStat <- list(dm=data, nlev=unique_vals, adaptDF=FALSE)
            #                     p_val <- disCItest(i, j, S=S[, s], suffStat)
            #                     if(p_val > 0.05){
            #                       sepset[[i]][[j]] <- S[, s]
            #                       print(cat('Check edge, i: ',i, ' j:', j, ' S:', S[, s]))
            #                       pc_i <- pc_i[!pc_i %in% j]
            #                       # pc_vars <- pc_vars[!pc_vars %in% j]
            #                       G[i, j] <- G[j, i] <- FALSE
            #                       break
            #                     }
            #                 }                            
            #             }
            #         }
            #     }
            # }
            # print(cat('i: ', i, ' pc_i:', pc_i))
            top_k[[i]] <- pc_i
            j_vars <- rev(pc_i)
            l <- l + 1
            if(l > max(unlist(lapply(top_k, function(x) length(x))))-1) break
        }
    }
    return(list(G, top_k))
}

checkDSep <- function(data, sepset, i, j, S, non_adj){
    pdsep_var <- S
    unique_vals <- apply(data, 2, function(x) length(unique(x)))
    idx <- NULL
    non_adj_sub <- c()
    for(s in non_adj){
        if(j %in% sepset[[i]][[s]]){
            d_sep <- sepset[[i]][[s]]
            d_sep <- union(d_sep[!d_sep %in% j], setdiff(pdsep_var, s))
            print(cat('s: ', s, 'd_sep: ', d_sep, ' pdsep_var:', pdsep_var))
            p_val <- disCItest(i, s, S=d_sep, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
            non_adj_sub <- c(non_adj_sub, s)
            if(p_val<0.05)
                return(list(FALSE, sepset))
        }
    }
    count <- length(non_adj_sub)
    if(count>0){
        for(c in non_adj_sub){
            d_sep <- sepset[[i]][[c]]
            sepset[[i]][[c]] <- union(d_sep[!d_sep %in% j], setdiff(pdsep_var, c))
        }
    }
    return(list(TRUE, sepset))
}

# Merge all top-k and build skeleton
mergeTopk <- function(top_k){
    num_vars <- length(top_k)
    skel <- matrix(FALSE, num_vars, num_vars)
    for(i in 1:num_vars){
        pc_i <- top_k[[i]]
        for(p in pc_i){
            if(i %in% top_k[[p]])
               skel[i, p] <- skel[p, i] <- TRUE
        }
    }
    return(skel)
}

# Build causal graph using PC-stable by prioritizing true_pc in conditioning set
buildCausalGraph <- function(data, top_k_pc){
    num_vars <- ncol(data)
    vars <- 1:num_vars
    G <- matrix(TRUE, nrow=num_vars, ncol=num_vars)
    diag(G) <- FALSE
    row.names(G) <- colnames(G) <- 1:num_vars
    l <- 1
    unique_vals <- apply(data, 2, function(x) length(unique(x)))
    suffStat <- list(dm=data, nlev=unique_vals, adaptDF=FALSE)
    for(i in 1:num_vars){
        j_vars <- vars[!vars %in% c(i, top_k_pc[[i]])]
        for(j in j_vars){
            p_val <- disCItest(x=i, y=j, S=NULL, suffStat=suffStat)
            if(p_val>0.05){
                print(p_val)
                G[i, j] <- G[j, i] <- FALSE
            }
        }
    }

    repeat{
        for(i in 1:num_vars){
            j_vars <- which(G[i, ]!=0, arr.ind=T)
            j_vars <- j_vars[!j_vars %in% c(top_k_pc[[i]])]
            if(l<=length(top_k_pc[[i]])){
                S_i <- combn(top_k_pc[[i]], l)
            }else{
                next
            }
            for(j in j_vars){
                for(s in 1:ncol(S_i)){
                    p_val <- disCItest(x=i, y=j, S=S_i[, s], suffStat=suffStat)
                    if(p_val>0.05){
                        G[i, j] <- G[j, i] <- FALSE
                    }            
                }
            }
        }
        l <- l + 1
        if(l>max(unlist(lapply(top_k_pc, length))))  break
    }
    return(G)
}

new_data <- data.frame(alarm_data)
for(i in 1:ncol(new_data)){
  new_data[, i] <- factor(new_data[, i])
}

gk_mat <- GKtauDataframe(new_data)
for(i in 1:num_vars){
  sort_gk_vals <- sort(gk_mat[i, ], decreasing=TRUE, index.return=TRUE)
  nbrs[i, ] <- sort_gk_vals$ix[!sort_gk_vals$ix %in% i][1:K]
  unique_vals <- apply(alarm_data[, c(i, nbrs[i, ])], 2, function(x) length(unique(x)))
  V <- as.character(c(i, nbrs[i, ]))
  suffStat <- list(dm=alarm_data[, c(i, nbrs[i, ])], nlev=unique_vals, adaptDF=FALSE)
  pc_fit[[i]] <- pc(suffStat, indepTest = disCItest, alpha=0.05, labels=V)
  amat[[i]] <- wgtMatrix(pc_fit[[i]])
  amat[[i]] <- amat[[i]] | t(amat[[i]])
}

comb_amat <- matrix(0, num_vars, num_vars)
for(i in 1:num_vars){
  amat_i <- amat[[i]]
  nbrs_i <- as.numeric(names(which(amat_i[1, ]!=0)))
  for(n in nbrs_i){
    amat_n <- amat[[n]]
    nbrs_n <- as.numeric(names(which(amat_n[1, ]!=0)))
    if(i %in% nbrs_n){
      comb_amat[i, n] <- comb_amat[n, i] <- 1
    }
  }
}

comb_graph <- igraph.to.graphNEL(graph_from_adjacency_matrix(comb_amat, mode="undirected"))
for(i in 1:num_iter){
  rand_idx[i, ] <- sample(1:num_vars, num_sub_vars)
  unique_vals <- apply(alarm_data[, rand_idx[i, ]], 2, function(x) length(unique(x)))
  V <- as.character(rand_idx[i, ])
  suffStat <- list(dm=alarm_data[, rand_idx[i, ]], nlev=unique_vals, adaptDF=FALSE)
  pc_fit[[i]] <- pc(suffStat, indepTest = disCItest, alpha=0.05, labels=V)
  amat[[i]] <- wgtMatrix(pc_fit[[i]])
}

undir_edges <- dir_edges <- data.frame(matrix(0, nrow=1, ncol=3))

for(i in 1:num_iter){
  var_names <- as.numeric(colnames(amat[[i]]))
  all_edges <- showEdgeList(pc_fit[[i]])
  
  edges <- updateEdgeStat(var_names, all_edges, undir_edges, dir_edges)
  undir_edges <- edges[[1]]
  dir_edges <- edges[[2]]
}
undir_edges <- undir_edges[-1, ]
dir_edges <- dir_edges[-1, ]


updateEdgeStat <- function(var_names, all_edges, undir_edges, dir_edges){
  for(e in 1:nrow(all_edges$undir)){
    x <- var_names[all_edges$undir[e, 1]]
    y <- var_names[all_edges$undir[e, 2]]
    edge <- data.frame(X1=x, X2=y)
    rev_edge <- data.frame(X1=y, X2=x)
    if(nrow(merge(edge, undir_edges[, -3]))>0 || nrow(merge(rev_edge, undir_edges[, -3]))>0){
      bool_vec <- apply(undir_edges[, -3], 1, function(x) all(x==edge))
      row_idx <- which(bool_vec==TRUE)
      undir_edges[row_idx, 3] <- undir_edges[row_idx, 3] + 1
    }else{
      edge <- data.frame(X1=x, X2=y, X3=1)
      undir_edges <- rbind(undir_edges, edge)
    } 
  }
  
  for(e in 1:nrow(all_edges$direct)){
    x <- var_names[all_edges$direct[e, 1]]
    y <- var_names[all_edges$direct[e, 2]]
    edge <- data.frame(X1=x, X2=y)
    if(nrow(merge(edge, dir_edges[, -3]))>0){
      bool_vec <- apply(dir_edges[, -3], 1, function(x) all(x==edge))
      row_idx <- which(bool_vec==TRUE)
      dir_edges[row_idx, 3] <- dir_edges[row_idx, 3] + 1
    }else{
      edge <- data.frame(X1=x, X2=y, X3=1)
      dir_edges <- rbind(dir_edges, edge)
    } 
  }
  return(list(undir_edges, dir_edges))
}

perturbPC <- function(vars, true_pc, num_noise_vars, order){
    quotient <- num_noise_vars/vars
    rem <- num_noise_vars%%vars
    num_vars <- 1:length(true_pc)
    if(quotient==0){
        vars_1 <- sample(vars, num_noise_vars)
        for(i in 1:length(vars_1)){
            org_pc <- true_pc[[vars_1[i]]]
            noise_var <- sample(num_vars[-c(org_pc)], 1)
            true_pc[[vars_1[i]]][1] <- noise_var
        }
    }else{

    }
}