############################################################
CALDER_main = function(contact_mat_file, 
					   chr, 
					   bin_size, 
					   out_dir, 
					   sub_domains = TRUE, 
					   save_intermediate_data = TRUE,
					   genome='hg19') {
	required_packages = c('doParallel', 'GenomicRanges', 'R.utils', 'factoextra', 'maptools')
	sapply(required_packages, require, character.only = TRUE, quietly = TRUE)

	dir.create(out_dir, showWarnings = FALSE)
	intermediate_data = CALDER_CD_hierarchy(contact_mat_file, 
											chr, 
											bin_size, 
											out_dir, 
											save_intermediate_data,
											genome)

	if(sub_domains==TRUE) {
		CALDER_sub_domains(intermediate_data=intermediate_data, 
						   chr=chr, 
						   out_dir=out_dir, 
						   bin_size=bin_size)
	}
}

