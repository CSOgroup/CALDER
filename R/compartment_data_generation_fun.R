	
	compress_mat_fast_tmp = function(input_mat, compress_size)
	  {
	         mat_block <- function(n, r) suppressWarnings( matrix(c(rep(1, r), rep(0, n)), n, n/r) )
	                n = ncol(input_mat)
	                mat2prod = mat_block(n, compress_size)
	                # return(t(input_mat%*%mat2prod / compress_size))
	                return(t(matrix_multiplication_cpp(input_mat, mat2prod)))
	        }

	mat_10to40kb = function(mat2compress, bin_size, bin_size_initial) ## 10kb to 40 kb
	{
		compress_size = bin_size / bin_size_initial
		len = nrow(mat2compress) - nrow(mat2compress)%%compress_size
		mat2compress = mat2compress[1:len, 1:len]
		c_mat_compressed = compress_mat_fast_tmp( as.matrix(mat2compress), compress_size=compress_size )
		mat_compressed = compress_mat_fast_tmp( c_mat_compressed, compress_size=compress_size )
		rownames(mat_compressed) = colnames(mat_compressed) = as.character( 1:nrow(mat_compressed) )
		return(mat_compressed)
	}

	# if(!file.exists(contact_mat_file)) contact_mat_file = paste0('/mnt/ndata/Yuanlong/2.Results/1.Juicer/', CELL_LINE, '/contact_mat/mat_chr', chr, '_', bin_size_initial_kb, 'kb_ob.txt.gz')

	## bin_size_initial is the binsize of your input matrix, can be different from the bin_size of your planned analysis
	contact_mat_processing = function(contact_mat_file, bin_size, bin_size_initial=bin_size)
	{	
		
		compress_size = ifelse(bin_size < 40E3, 1, 1)
		zero_ratio = 0.01

		combined_xk_oe_raw = data.table::fread(contact_mat_file)

		## this code generates the compartment domains

		combined_xk_oe_raw = subset(combined_xk_oe_raw, !is.na(V3))
		combined_xk_oe_raw[,1] = combined_xk_oe_raw[,1]/bin_size_initial
		combined_xk_oe_raw[,2] = combined_xk_oe_raw[,2]/bin_size_initial
		combined_xk_oe = combined_xk_oe_raw

		colnames(combined_xk_oe) = c('pos_1', 'pos_2', 'val')	
		if(!all(combined_xk_oe[[2]] >= combined_xk_oe[[1]])) stop('\nYou provided matrix does not represent an upper triangular matrix!\n\n')

		oe_size = max(max(combined_xk_oe[[1]]), max(combined_xk_oe[[2]])) + 1 ## should +1 because combined_xk_oe index starts from 0 (bin 0 represents: 0-10E3, checked by looking at the juicebox map, 2018-11-19)
		mat_oe_sparse = Matrix::Matrix(0, nrow=oe_size, ncol=oe_size)
		mat_oe_sparse[cbind(combined_xk_oe[[1]]+1, combined_xk_oe[[2]]+1)] <- combined_xk_oe[[3]]

		rownames(mat_oe_sparse) = colnames(mat_oe_sparse) = as.character( 1:nrow(mat_oe_sparse) )
		
		mat_oe_sparse = Matrix::forceSymmetric(mat_oe_sparse, uplo='U')
		if(bin_size!=bin_size_initial) mat_oe_sparse = mat_10to40kb( mat_oe_sparse, bin_size, bin_size_initial )
		A_oe = remove_blank_cols(mat_oe_sparse, sparse=TRUE, ratio=zero_ratio) ## has the same rows/cols as A
		if(nrow(A_oe) < 100) A_oe = remove_blank_cols(mat_oe_sparse, sparse=TRUE, ratio=0) ## when all are dense
		while(min(apply(A_oe, 1, sd))==0) ## sometimes after removing the cols / rows, the remained part will all be 0
		{
			A_oe = remove_blank_cols(A_oe, sparse=TRUE, ratio=1E-7) ## has the same rows/cols as A
			if(nrow(A_oe) < 1) stop('ERROR IN GENERATING MEANINGFUL A_oe at the data generating step')
		}

		##########################################################

		len = nrow(A_oe) - nrow(A_oe)%%compress_size
		A_oe_2_compress = A_oe[, 1:len]

		bin_names = rownames(A_oe)

		A_oe_compressed = compress_mat_fast( as.matrix(A_oe_2_compress), compress_size=compress_size )
		colnames(A_oe_compressed) = bin_names
		rm(A_oe_2_compress); gc()
		
		range(A_oe_compressed)
		# # sum(A_oe_compressed > 1000)
		# # A_oe_compressed[A_oe_compressed > 1000] = 1000
		A_oe_compressed_sparse = A_oe_compressed
		A_oe_compressed = as.matrix(A_oe_compressed)
		A_oe_compressed_log = log2(A_oe_compressed + 1)
		
		# #########################################################
		# cat('compute correlation matrix ... ')

		cA_oe_compressed_log = fast_cor(A_oe_compressed_log)
		ccA_oe_compressed_log = fast_cor(cA_oe_compressed_log)

		# cat('compute correlation matrix done ... ')

		# #########################################################
		# # ccA_oe_compressed_atanh = atanh(ccA_oe_compressed - 1E-7)
		ccA_oe_compressed_log_atanh = atanh(ccA_oe_compressed_log / (1+1E-7))

		# rm(A_oe_compressed, A_oe_compressed_sparse, cA_oe_compressed_log)
		gc()
		# #########################################################
		# cat('ready to compute compartment domains\n')

		out = list(A_oe=A_oe, atanh_score=ccA_oe_compressed_log_atanh)

		return(out)
	}

