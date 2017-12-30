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

getPC <- function(disc_data, top_k, G, unique_vals){
    if(is.null(unique_vals)){
        cor_mat <- abs(cor(disc_data))
    }
    # Remove the following two lines
    top_k <- top_k_cv
    G <- G_cv

    num_vars <- ncol(disc_data)
    seq_p <- seq(num_vars)
    pVals <- matrix(-Inf, nrow=num_vars, ncol=num_vars)
    sepset <- lapply(seq_p, function(.) vector("list", num_vars))
    l <- 0
    repeat{
    	sepset_discard <- lapply(seq_p, function(.) vector("list", num_vars))
        non_adj <- vector("list", num_vars)
        adj <- top_k
        stableAdj <- FALSE
        while(!stableAdj){
	        count <- 0
	    	for(i in 1:num_vars){
	            pc_i <- top_k[[i]]
				j_vars <- rev(adj[[i]])
	            for(j in j_vars){
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
	                        if(l>0 && (j %in% S[, s] || invalidSepset(S[, s], l, sepset_discard[[i]][[j]]))) next
	                        if(is.null(unique_vals)){
 	                            suffStat <- list(C=cor(disc_data), n=nrow(disc_data))
                            	p_val <- gaussCItest(i, j, S[, s], suffStat)                      
	                        }else{     
	                            suffStat <- list(dm=disc_data, nlev=unique_vals, adaptDF=FALSE)
	                            p_val <- disCItest(i, j, S=S[, s], suffStat)                      
	                        }
	                        if(p_val > 0.05){
	                            print(cat('First test, i: ',i, ' j:', j, ' S:', S[, s], ' p_val: ', p_val))
	                            pdsep_res <- checkDSep(disc_data, sepset, i, j, S[, s], non_adj[[i]])
	                            pdsep <- pdsep_res[[1]]
	                            sepset <-pdsep_res[[2]]
	                            sepset_old <- pdsep_res[[3]]
	                            if(pdsep){
	                                non_adj[[i]] <- c(non_adj[[i]], j)
	                                pVals[i, j] <- pVals[j, i] <- p_val
	                                sepset[[i]][[j]] <- S[, s]
	                                pc_i <- pc_i[!pc_i %in% j]
	                                G[i, j] <- G[j, i] <- FALSE
	                                break
	                            }else{
	                            	sepset_ij <- sepset_discard[[i]][[j]]
	                            	sepset_ij[[length(sepset_ij)+1]] <- S[, s]
	                            	sepset_discard[[i]][[j]] <- sepset_ij
	                            }
	                        }
	                    }
	                }
	            }
	            top_k[[i]] <- pc_i
	        }
	        # Check for the necessary path condition
	        if(l>0){
		        for(i in 1:num_vars){
		        	non_adj_cpy <- non_adj[[i]]
		        	res <- checkNPC(G, top_k, sepset, sepset_old, sepset_discard, i, non_adj[[i]])
		        	G <- res[[1]]
		        	top_k <- res[[2]]
		        	non_adj_sub <- res[[3]]
		        	sepset <- res[[4]]
		        	sepset_old <- res[[5]]
		        	sepset_discard <- res[[6]]
		        	if(is.null(non_adj_sub)){
		        		non_adj[i] <- list(NULL)
		        	}else{
		        		non_adj[[i]] <- res[[3]]
		        	} 

		        	adj_set <- setdiff(non_adj_cpy, non_adj[[i]])
		        	if(is.null(adj_set)){
	        			adj[i] <- list(NULL)
		        	}else{
		        		adj[[i]] <- adj_set
		        	}
		        	count <- count + length(adj[[i]])
		        	# if(length(adj[[i]])>0)
		        		# print(cat('i: ', i,' Edges added: ', adj[[i]]))
		        }
	        	stableAdj <- ifelse(count>0, FALSE, TRUE)	        	
	        }else{
	        	stableAdj <- TRUE
	        }
        }
        l <- l + 1
        if(l > max(unlist(lapply(top_k, function(x) length(x))))-1) break
    }
    return(list(G, top_k))
}

performCITest <- function(data, G, cor_mat, sepset, pc_i){

}

invalidSepset <- function(s, l, sepset_discard){
	if(length(sepset_discard)>0){
		n_sepset <- length(sepset_discard)
		for(i in 1:n_sepset){
			if(all(s %in% sepset_discard[[i]])){
				return(TRUE)
			}
		}		
	}
	return(FALSE)
}

checkDSep <- function(data, sepset, i, j, S, non_adj){
    pdsep_var <- S
    unique_vals <- apply(data, 2, function(x) length(unique(x)))
    idx <- NULL
    non_adj_sub <- c()
    sepset_old <- sepset
    for(s in non_adj){
        if(j %in% sepset[[i]][[s]] && length(sepset[[i]][[s]]>0)){
            d_sep <- sepset[[i]][[s]]
            d_sep <- union(d_sep[!d_sep %in% j], setdiff(pdsep_var, s))
            print(cat('s: ', s, 'd_sep: ', d_sep, ' pdsep_var:', pdsep_var))
            p_val <- disCItest(i, s, S=d_sep, suffStat=list(dm=data, nlev=unique_vals, adaptDF=FALSE))
            non_adj_sub <- c(non_adj_sub, s)
            if(p_val<0.05)
                return(list(FALSE, sepset, sepset))
        }
    }
    count <- length(non_adj_sub)
    if(count>0){
        for(c in non_adj_sub){
            d_sep <- sepset[[i]][[c]]
            sepset_old[[i]][[c]] <- d_sep
            sepset[[i]][[c]] <- union(d_sep[!d_sep %in% j], setdiff(pdsep_var, c))
        }
    }
    return(list(TRUE, sepset, sepset_old))
}

checkNPC <- function(G, top_k, sepset, sepset_old, sepset_discard, i, non_adj){
	pc_i <- top_k[[i]]
	non_adj_cpy <- non_adj
	for(j in non_adj){
		s_new <- sepset[[i]][[j]]
		s_old <- sepset_old[[i]][[j]]
		flag_new <- bfs(s_new, i, j, top_k)
		if(length(s_new)!=length(s_old) || !all(s_new %in% s_old)){
			flag_old <- bfs(s_old, i, j, top_k)
		}else{
			flag_old <- flag_new
			# if(i==21){
			# 	print(cat('Same Path i:', i, ' j: ', j, ' s_new:', s_new, ' s_old:', s_old))
			# 	print(cat('top_k: ', top_k[[i]]))
			# }
		}
		if(flag_new && !flag_old){
			# print(cat('New Path i:', i, ' j: ', j, ' s_new:', s_new, ' s_old:', s_old))
			sepset_old[[i]][[j]] <- s_new
		}
		if(!flag_new && flag_old){
			# print(cat('Old Path i:', i, ' j: ', j, ' s_new:', s_new, ' s_old:', s_old))
			sepset[[i]][[j]] <- s_old
		}
		if(!flag_new && !flag_old){
			# print(cat('Path not found i:', i, ' j: ', j, ' s_new:', s_new, ' s_old:', s_old))
			if(is.null(sepset_discard[[i]][[j]])){
				sepset_discard[[i]][[j]] <- list(s_new, s_old)
			}else{
				sepset_ij <- sepset_discard[[i]][[j]]
				sepset_ij[[length(sepset_ij)+1]] <- s_new
				sepset_ij[[length(sepset_ij)+1]] <- s_old
				sepset_discard[[i]][[j]] <- sepset_ij
			}
			# sepset_discard[[i]][[j]] <- c(sepset_discard[[i]][[j]], s_new, s_old)
			sepset[[i]][[j]] <- sepset_old[[i]][[j]] <- integer(0)
			top_k[[i]] <- union(top_k[[i]], j)
			G[i, j] <- G[j, i] <- TRUE
			non_adj_cpy <- non_adj_cpy[!non_adj_cpy %in% j]
		}
	}
	non_adj <- non_adj_cpy
	return(list(G, top_k, non_adj, sepset, sepset_old, sepset_discard))
}

bfs <- function(sep_set, i, j, top_k){
	if(length(sep_set)==0)
		return(FALSE)
	count <- array(0, length(sep_set))
	k <- 1
	for(s in sep_set){
		open_queue <- s
		visited_vars <- i
		while(length(open_queue)>0){
			v <- open_queue[1]
			open_queue <- open_queue[-1]
			if(v == j){
				count[k] <- 1
				k <- k + 1
			}
			v_child <- top_k[[v]]
			v_child <- v_child[!v_child %in% intersect(v_child, visited_vars)]
			if(length(open_queue)>0){
				v_child <- v_child[!v_child %in% open_queue]
			}
			open_queue <- union(open_queue, v_child)
			visited_vars <- union(visited_vars, v)
			# print(cat('open_queue: ', open_queue, ' v_child:', v_child))
		}		
	}
	if(sum(count)==length(sep_set))	return(TRUE)
	return(FALSE)
}