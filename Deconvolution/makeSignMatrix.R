function(st_data, expression_values=c("normalized"), logbase=2, cluster_column="leiden_clus",
                          name=NULL, return_gobject= TRUE)
  
{
  expr_values = merfish_st@norm_expr
  nolog_expr = logbase^(expr_values) - 1
  cell_metadata = pDataDT(merfish_st)
  if (!cluster_column %in% colnames(cell_metadata)) {
    stop("\n cluster column not found \n")
  }
  cluster = cell_metadata[[cluster_column]]
  Sig_exp = as.matrix(Sig_exp)
  intersect_gene = intersect(rownames(Sig_exp), rownames(nolog_expr))
  filter_Sig = Sig_exp[intersect_gene, ]
  #filter_expr = nolog_expr[intersect_gene, ]
  #filter_log_expr = expr_values[intersect_gene, ]
  
  #signature matrix for PAGE and Hypergeomtric
  filter_Sig = filter_Sig[rowSums(filter_Sig) > 0, ]
  enrich_matrix = matrix(0, nrow = dim(filter_Sig)[1], ncol = dim(filter_Sig)[2])
  rowmax_col = Rfast::rowMaxs(filter_Sig)
  for (i in 1:length(rowmax_col)) {
    enrich_matrix[i, rowmax_col[i]] = 1
  }
  colsum_ct_binary <- colSums(enrich_matrix)
  for (i in 1:length(colsum_ct_binary)) {
    if (colsum_ct_binary[i] <= 2) {
      rank = rank(-ct_exp[, i])
      enrich_matrix[rank <= 2, i] = 1
    }
  }
  rownames(enrich_matrix) = rownames(filter_Sig)
  colnames(enrich_matrix) = colnames(filter_Sig)
  return(enrich_matrix)
}

