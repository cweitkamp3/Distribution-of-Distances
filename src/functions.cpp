#include <RcppArmadillo.h>
#include <stdexcept>
#include <RcppArmadilloExtensions/sample.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec determine_lower_tri_indices(const arma::uvec& input, int noc_dist_mat){
  size_t boot_m = input.size();
  arma::uvec result(boot_m*(boot_m-1)/2);
  int counter =0;
  for(int i = 0; i < boot_m; i++){
    for(int j = i+1; j<boot_m; j++){
      result[counter]= input[i]*noc_dist_mat+input[j];
      counter += 1;
    }
  }
  return result;
}



// L_2 distance between two empirical quantile functions
// [[Rcpp::export]]
double trimmed_quantile_diff(const arma::vec sorted_distx,const arma::vec sorted_disty, double beta, double p){
  if(beta >=0.5){
    throw std::invalid_argument( "Error: Beta has to be chosen in [0,0.5)");
  }
  double cost = 0;
  const size_t n = sorted_distx.size();
  const size_t m = sorted_disty.size();
  //std::cout <<n<<std::endl;
  //std::cout <<m<<std::endl;

  int start_index_u = floor(beta*n);
  int end_index_u =   ceil((1.-beta)*n);
  //std::cout << (beta)*n<<std::endl;
  //std::cout << (1.-beta)*n<<std::endl;
  //std::cout << start_index_u<<std::endl;
  //std::cout << end_index_u<<std::endl;

  if((start_index_u+1)== end_index_u){
      throw std::invalid_argument("Beta has been chosen too large in comparison to the sample size.");
    }
    int start_index_v = floor(beta*m);
    int end_index_v = ceil((1.-beta)*m);
    //std::cout << (beta)*m<<std::endl;
    //std::cout << (1.-beta)*m<<std::endl;
    //std::cout << start_index_v<<std::endl;
    //std::cout << end_index_v<<std::endl;

    if((start_index_v+1)== end_index_v){
      throw std::invalid_argument("Beta has been chosen too large in comparison to the sample size.");
    }
    int i = start_index_u;
    double w_i =  (start_index_u +1)/(double)n-beta;
    int j = start_index_v;
    double  w_j = (start_index_v +1)/(double)m-beta;
    //std::cout << w_i<<std::endl;
    //std::cout << w_j<<std::endl;
    double m_ij = 0;
    while(TRUE){
      m_ij = pow(fabs(sorted_distx[i] - sorted_disty[j]),p);
      if (w_i < w_j || j == m-1){
        cost = cost+m_ij * w_i;
        /*std::cout << "First"<<std::endl;
        std::cout << m_ij<<std::endl;
        std::cout << w_i<<std::endl;
        std::cout << m_ij * w_i<<std::endl;
        std::cout << cost<<std::endl;*/

        i = i+1;
        if (i == end_index_u){
          break;
        }
        w_j = w_j-w_i;
        if(i<end_index_u-1){
          w_i = 1/(double)n;
        }
        else{
          w_i=  (1-beta)-(end_index_u-1)/(double)n;
        }
      }
      else{
        cost = cost+ m_ij * w_j;
        /*std::cout << "Second"<<std::endl;
        std::cout << m_ij<<std::endl;
        std::cout << w_i<<std::endl;
        std::cout << m_ij * w_j<<std::endl;
        std::cout << cost<<std::endl;*/

        j = j+1;
        if (j == end_index_v){
          break;
        }
        w_i =w_i- w_j;
        if(j <end_index_v-1){
          w_j = 1/(double)m;
        } else {
          w_j=  (1-beta)-(end_index_v-1)/(double)m;
          }
        }
      }
    return(cost);
  }

// [[Rcpp::export]]
arma::vec Bootstrap(int number, const double beta, const arma::mat& dist_sample, const int m, const arma::vec& sorteddist, const int p){
  //std::cout << "0"<<std::endl;
  arma::vec result(number);
  int n = dist_sample.n_rows;
  //std::cout << "0.5"<<std::endl;
  for (int i =0; i<number; i++){
    R_CheckUserInterrupt();
    arma::uvec sequence = arma::linspace<arma::uvec>(0, n-1, n);
    arma::uvec sampled_numbers=Rcpp::RcppArmadillo::sample(sequence, m, true);

    arma::uvec indices= determine_lower_tri_indices(arma::sort(sampled_numbers),n);

    arma::vec distbootstrapped = dist_sample.elem(indices);
    //std::cout << distbootstrapped <<std::endl;

    result[i] = ((double)m)*trimmed_quantile_diff(arma::sort(distbootstrapped),sorteddist,beta,p);
  }
  return(result);
}

