// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <algorithm>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_hello_world() {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Identity(3, 3);
    // Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);
    // Do not use Random() here to not promote use of a non-R RNG
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Zero(3, 3);
    for (auto i=0; i<m2.rows(); i++)
        for (auto j=0; j<m2.cols(); j++)
            m2(i,j) = R::rnorm(0, 1);

    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_outerproduct(const Eigen::VectorXd & x) {
    Eigen::MatrixXd m = x * x.transpose();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcppeigen_innerproduct(const Eigen::VectorXd & x) {
    double v = x.transpose() * x;
    return v;
}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcppeigen_bothproducts(const Eigen::VectorXd & x) {
    Eigen::MatrixXd op = x * x.transpose();
    double          ip = x.transpose() * x;
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}


// [[Rcpp::export]]
Eigen::MatrixXd get_dist_zhu_eigen(Eigen::Map<Eigen::MatrixXd> A) {
  int N = A.rows();
  Eigen::MatrixXd dist_matrix(N, N);
  Eigen::MatrixXd A_sq = (A * A) / N;
  
  for (int i = 0; i < N - 1; i++) {
    for (int ip = i + 1; ip < N; ip++) {
      Eigen::VectorXd tgtvec = (A_sq.row(i) - A_sq.row(ip)).cwiseAbs();
      tgtvec[i] = 0;
      tgtvec[ip] = 0;
      double max_d = tgtvec.maxCoeff();
      dist_matrix(i, ip) = max_d;
      dist_matrix(ip, i) = max_d; // Mirror the matrix
    }
  }
  return dist_matrix;
}





// [[Rcpp::export]]
Eigen::MatrixXd calculate_p_hat_matrix(Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXi> neibors_matrix) {
  int N = A.rows();
  Eigen::MatrixXd p_hat_matrix = Eigen::MatrixXd::Zero(N, N);
  
  for (int a = 0; a < N; ++a) {
    for (int b = a; b < N; ++b) { // Loop starts from 'a' to include diagonal
      
      const Eigen::VectorXi index1 = neibors_matrix.row(a);
      const Eigen::VectorXi index2 = neibors_matrix.row(b);
      
      double sum = 0.0;
      int num1 = index1.size();
      int num2 = index2.size();
      
      for (int i = 0; i < num1; ++i) {
        for (int j = 0; j < num2; ++j) {
          int row = index1[i];
          int col = index2[j];
          sum += A(row, col);
        }
      }
      
      // Calculate mean
      double mean_val = sum/(num1*num2);
      p_hat_matrix(a, b) = mean_val;
      if (a != b) {
        p_hat_matrix(b, a) = mean_val; // Mirror to lower triangular part
      }
    }
  }
  
  return p_hat_matrix;
}


// [[Rcpp::export]]
Eigen::MatrixXd sample_from_p_cpp(Eigen::Map<Eigen::MatrixXd> p_hat_matrix,
                              Eigen::Map<Eigen::VectorXi> blist, bool no_loop) {
  int N = p_hat_matrix.rows();
  
  Eigen::MatrixXd p_hat_matrix_b(N, N);
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      p_hat_matrix_b(i, j) = p_hat_matrix(blist(i), blist(j));
    }
  }
  
  Eigen::MatrixXd random_matrix = Eigen::MatrixXd::Random(N, N);
  // Adjust to range [0, 1]
  random_matrix = (random_matrix + Eigen::MatrixXd::Constant(N, N, 1.0)) / 2.0;
  
  Eigen::MatrixXd g_adj_nb(N, N);
  
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      g_adj_nb(i, j) = 1 * (random_matrix(i, j) < p_hat_matrix_b(i, j));
      g_adj_nb(j, i) = g_adj_nb(i, j); // Mirror to lower triangular part
    }
  }
  
  if (!no_loop) {
    for (int i = 0; i < N; i++) {
      g_adj_nb(i, i) = 1 * (random_matrix(i, i) < p_hat_matrix_b(i, i));
    }
  }
  
  return g_adj_nb;
  
}