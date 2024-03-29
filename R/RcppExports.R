# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

determine_lower_tri_indices <- function(input, noc_dist_mat) {
    .Call(`_DoD_determine_lower_tri_indices`, input, noc_dist_mat)
}

trimmed_quantile_diff <- function(sorted_distx, sorted_disty, beta, p) {
    .Call(`_DoD_trimmed_quantile_diff`, sorted_distx, sorted_disty, beta, p)
}

Bootstrap <- function(number, beta, dist_sample, m, sorteddist, p) {
    .Call(`_DoD_Bootstrap`, number, beta, dist_sample, m, sorteddist, p)
}

