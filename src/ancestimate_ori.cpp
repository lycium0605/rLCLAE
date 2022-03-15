#include <Rcpp.h>
using namespace Rcpp;

//' @title anccall_c_ori: Generating ancestry call along the chromosome
//'
//' @param deltaf The cut-off for minimum delta f between two reference populations
//' @param window Size of the sliding window
//' @param SMAX The number of all snp sites
//' @param anclikdir the dir to anclik file
//' @param output the dir of the output anccall file
// [[Rcpp::export]]
void anccall_c_ori(double deltaf, int window, int SMAX,
              std::string anclikdir,
              std::string output,
              std::string chrom, std::string indiv,
              double mode, int n) {

  int a, c1, c2, c0, b, c, s, nsum, *anc, *anc2, *majmode;
  double *pos, *lik2, *lik1, *lik0, p0, p1, p2, d, delta, t, posit;
  FILE *anclik, *out;
  t=window;
  d=deltaf;
  Rcout<<"Finding "<<SMAX<<" snps in the file. Scanning ancestral likelihood..."<<std::endl;
  SMAX=SMAX+200;
  lik2 = (double *) malloc (SMAX * sizeof (double));
  lik1 = (double *) malloc (SMAX * sizeof (double));
  lik0 = (double *) malloc (SMAX * sizeof (double));
  pos = (double *) malloc (SMAX * sizeof (double));
  anc = (int *) malloc (SMAX * sizeof (int));
  anc2 = (int *) malloc (SMAX * sizeof (int));
  majmode = (int *) malloc (5 * sizeof (int));
  anclik=fopen(anclikdir.c_str(),"r");
  out=fopen(output.c_str(),"w");
  s=0;
  c=0;
  while (fscanf(anclik,"%lf\t%lf\t%lf\t%lf\t%lf\n", &posit, &delta, &p2, &p1, &p0) != EOF) {
    if ((delta > (d - 1.E-10)) && (p0>-1.E-10 && p2>-1.E-10)) { //for hap p1==0
      if(p0>p1 && p0>p2){
        c=0;
        fprintf(out,"%f\t%d\n",posit, c);
      }
      else if(p2>p0 && p2>p1){
        c=2;
        fprintf(out,"%f\t%d\n",posit, c);
      }
      else{
        c=-1;
        fprintf(out,"%f\t%d\n",posit, c);
      }
    }
    //else{
      //++s;
    //}
  }
  Rcout<<"Scanning finished, "<<s<<" informative sites found."<<std::endl;

  fclose(anclik);
  fclose(out);
  }
