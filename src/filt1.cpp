#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

//' filt1_dip: transforming the phred score to relative genotype likelihood
//'
//' @param n an integer, the number of individuals contain in the vcf
//' @param input a string, the dir to cleaned vcf file
//' @param output a string, the dir to the output file
//' @export
// [[Rcpp::export]]
void filt1_dip(int n, std::string input, std::string output) {
  int a, pos, g1, g2, g3, cov;
  double l1, l2, l3;
  char sname[20];

  Rcout<<n<<" individuals in the input"<<std::endl;

  FILE* fin = fopen(input.c_str(), "r");
  FILE* fout = fopen(output.c_str(), "w");
  Rcout<<"Reading in input..."<<std::endl;
  while (std::fscanf(fin,"%s %d ",&sname,&pos) != EOF){
    std::fprintf(fout,"%s %d ", sname, pos);
    for (a=0; a<n; ++a) {
      fscanf (fin,"%d\n", &cov);
      if (cov == 999) {
        fscanf (fin,"%*s %*s %*s %*s %*s ");
        fprintf (fout,"-1. -1. -1. ");
      }
      else {
        fscanf (fin,"%*s %*s %*s %*s %d,%d,%d ", &g1, &g2, &g3);
        l1 = pow (10., -g1/10.);
        l2 = pow (10., -g2/10.);
        l3 = pow (10., -g3/10.);
        fprintf (fout,"%lf %lf %lf ", l1/(l1+l2+l3), l2/(l1+l2+l3), l3/(l1+l2+l3));
      }
    }
    fprintf(fout,"\n");
  }
  Rcout<<"Calculation finished."<<std::endl;
  fclose(fin);
  fclose(fout);
}

//' filt1_hap: transforming the phred score to relative genotype likelihood
//'
//' @param n an integer, the number of individuals contain in the vcf
//' @param input a string, the dir to cleaned vcf file
//' @param output a string, the dir to the output file
//' @export
// [[Rcpp::export]]
void filt1_hap(int n, std::string input, std::string output) {
  int a, pos, g1, g2, g3, cov;
  double l1, l2, l3;
  char sname[20];

  Rcout<<n<<"individuals in the input"<<std::endl;

  FILE* fin = fopen(input.c_str(), "r");
  FILE* fout = fopen(output.c_str(), "w");
  Rcout<<"Reading in input..."<<std::endl;
  while (std::fscanf(fin,"%s %d ",&sname,&pos) != EOF){
    std::fprintf(fout,"%s %d ", sname, pos);
    for (a=0; a<n; ++a) {
      fscanf (fin,"%d\n", &cov);
      if (cov == 999) {
        fscanf (fin,"%*s %*s %*s %*s ");
        fprintf (fout,"-1. -1. ");
      }
      else {
        fscanf (fin,"%*s %*s %*s %d,%d ", &g1, &g3);
        l1 = pow (10., -g1/10.);
        l3 = pow (10., -g3/10.);
        fprintf (fout,"%lf %lf", l1/(l1+l3),l3/(l1+l3));
      }
    }
    fprintf(fout,"\n");
  }
  Rcout<<"Calculation finished."<<std::endl;
  fclose(fin);
  fclose(fout);
}
