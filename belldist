#include <Rcpp.h>
#include <cmath>
#include <vector>
// Lambert W function for the principal branch (W_0)
// Using Halley's method for numerical approximation
double lambertW0(double z) {
  if (z == 0.0) {
    return 0.0;
  }
  if (z == -1.0 / std::exp(1.0)) {
    return -1.0;
  }
  if (z < -1.0 / std::exp(1.0)) {
    return NAN; // No real solution exists
  }
  
  double w = z;
  double tolerance = 1e-10;
  for (int i = 0; i < 100; ++i) {
    double ew = std::exp(w);
    double wew = w * ew;
    double wewz = wew - z;
    double w1 = w + 1.0;
    double delta_w = wewz / (ew * w1 - (w + 2.0) * wewz / (2.0 * w1));
    w -= delta_w;
    if (std::abs(delta_w) < tolerance) {
      return w;
    }
  }
  return w; // Return the last approximation if convergence not achieved
}

// Expose the function to R
// [[Rcpp::export]]
double lambertW(double z) {
  return lambertW0(z);
}
// [[Rcpp::export]]
double bell(int n) {
  std::vector<double> bell_numbers(n + 1, 0.0);
  bell_numbers[0] = 1;
  for (int i = 1; i <= n; ++i) {
    for (int j = 0; j < i; ++j) {
      bell_numbers[i] += bell_numbers[j] * std::pow(i - 1, j) / std::tgamma(j + 1);
    }
  }
  return bell_numbers[n];
}

// [[Rcpp::export]]
Rcpp::NumericVector dbell(Rcpp::IntegerVector x, double theta, bool log_prob = false) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    double Bx = bell(x[i]);
    double lf = x[i] * std::log(theta) - std::exp(theta) + 1 + std::log(Bx) - std::lgamma(x[i] + 1);
    result[i] = log_prob ? lf : std::exp(lf);
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector pbell(Rcpp::IntegerVector q, double theta, bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    double prob = q[i] < 0 ? 0 : Rcpp::sum(dbell(Rcpp::seq(0, q[i]), theta, false));
    if (!lower_tail) {
      prob = 1 - prob;
    }
    result[i] = log_p ? std::log(prob) : prob;
  }
  return result;
}
// Function to compute the quantile function of the Bell distribution
Rcpp::IntegerVector qbell(Rcpp::NumericVector p, double theta, bool log_p = false) {
  int n = p.size();
  Rcpp::IntegerVector result(n);
  for (int i = 0; i < n; ++i) {
    double target_p = log_p ? std::exp(p[i]) : p[i];
    double cumulative_p = 0;
    int k = 0;
    while (cumulative_p <= target_p) {
      cumulative_p = pbell(Rcpp::IntegerVector::create(k), theta, true, false)[0];
      ++k;
    }
    result[i] = k - 1;
  }
  return result;
}
// [[Rcpp::export]]
Rcpp::IntegerVector rbell(int n, double theta) {
  Rcpp::NumericVector uniform_randoms = Rcpp::runif(n);
  return qbell(uniform_randoms, theta);
}



// Function to compute quantile residuals for a Bell regression model
Rcpp::NumericVector qresiduals(Rcpp::List object) {
  Rcpp::NumericVector y = object["y"];
  Rcpp::NumericVector mu = object["mu"];
  int n = y.size();
  Rcpp::NumericVector residuals(n);
  for (int i = 0; i < n; ++i) {
    double theta = lambertW0(mu[i]);
    double a = pbell(Rcpp::IntegerVector::create(y[i] - 1), theta)(0);
    double b = pbell(Rcpp::IntegerVector::create(y[i]), theta)(0);
    double u = R::runif(a, b);
    residuals[i] = R::qnorm(u, 0.0, 1.0, 1, 0);
  }
  return residuals;
}

// Rcpp export for the qresiduals function
// [[Rcpp::export]]
Rcpp::NumericVector cpp_qresiduals(Rcpp::List object) {
  return qresiduals(object);
}
// Function to compute the log-likelihood of the Bell distribution
double log_likelihood_bell(Rcpp::IntegerVector x, double theta) {
  int n = x.size();
  double log_likelihood = 0.0;
  for (int i = 0; i < n; ++i) {
    double Bx = bell(x[i]);
    double lf = x[i] * std::log(theta) - std::exp(theta) + 1 + std::log(Bx) - std::lgamma(x[i] + 1);
    log_likelihood += lf;
  }
  return log_likelihood;
}

// Rcpp export for the log_likelihood_bell function
// [[Rcpp::export]]
double cpp_log_likelihood_bell(Rcpp::IntegerVector x, double theta) {
  return log_likelihood_bell(x, theta);
}

// Function to compute vmu variance function of mu
double vmu(double mu) {
  return mu * (1 + lambertW0(mu));
}
// Rcpp export for vmu function
// [[Rcpp::export]]
double cpp_vmu(double mu) {
  return vmu(mu);
}
