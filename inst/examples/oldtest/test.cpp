#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector test1(NumericVector p) {
  NumericVector q(p.size());
  q=p*2;
  return(q);
}

// [[Rcpp::export]]
NumericVector test2() {
  NumericVector m (10);
  NumericVector two(2);
  for (int i=0;i<m.size();i++)  m[i]=i;
  two = m[Range(2,3)];
  return  two[ifelse(runif(1)>0.5,1,0)];
}


#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector useTimer() {
    int n = 1000000;

    // start the timer
    Timer timer;
    timer.step("start");        // record the starting point

    for (int i=0; i<n; i++) {
        GetRNGstate();
        PutRNGstate();
    }
    timer.step("get/put");      // record the first step

    for (int i=0; i<n; i++) {
        GetRNGstate();
        rnorm(10, 0.0, 1.0);
        PutRNGstate();
    }
    timer.step("g/p+rnorm()");  // record the second step

    for (int i=0; i<n; i++) {
        // empty loop
    }
    timer.step("empty loop");   // record the final step

    NumericVector res(timer);   // 
    for (int i=0; i<res.size(); i++) {
        res[i] = res[i] / n;
    }
    return res;
}
