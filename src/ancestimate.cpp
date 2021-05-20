#include <Rcpp.h>
using namespace Rcpp;

//' ancest_c: Generating ancestry call along the chromosome
//'
//' @param deltaf The cut-off for minimum delta f between two reference populations
//' @param window Size of the sliding window
//' @param SMAX The number of all snp sites
//' @param anclikdir the dir to anclik file
//' @param output the dir of the output ancfreq file
//' @export
// [[Rcpp::export]]
void anccall_c(double deltaf, int window, int SMAX,
              std::string anclikdir, std::string output) {

  int a, c1, c2, c0, b, c, s, *anc;
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
  anclik=fopen(anclikdir.c_str(),"r");
  out=fopen(output.c_str(),"w");
  s=0;
  while (fscanf(anclik,"%lf\t%lf\t%lf\t%lf\t%lf\n", &posit, &delta, &p2, &p1, &p0) != EOF) {
    if ((delta > (d - 1.E-10)) && p0>1.E-10 && p1>1.E-10 && p2>1.E-10) {
      lik2[s] = log (p2);
      lik1[s] = log (p1);
      lik0[s] = log (p0);
      pos[s] = posit;
      ++s;
    }
  }
  Rcout<<"Scanning finished"<<std::endl;

  for (a=0; a<s; ++a) {
    p0 = p1 = p2 = 0.;

    for (b=0; b<s && pos[b]<(pos[a]+t*.5); ++b){
      if (pos[b]>(pos[a]-t*.5)) {
        p0 += lik0[b];
        p1 += lik1[b];
        p2 += lik2[b];
      }
    }
      if (p0>p1 && p0>p2)
        anc[a] = 0;
      else if (p1>p0 && p1>p2)
        anc[a] = 1;
      else if (p2>p0 && p2>p1)
        anc[a] = 2;
  }
  Rcout<<"First round of ancestry call finished. Starting smoothing with the sliding window..."<<std::endl;
  for (a=0; a<s; ++a) {
    c1 = c2 = c0 = 0;
    for (b=0; b<s; ++b){
      if ((pos[a]-pos[b]) < t*.5 && (pos[a]-pos[b])>(t*(-.5))) {
        if (anc[b] == 0){++c0;}

        else if (anc[b] == 1){++c1;}

        else if (anc[b] == 2){++c2;}

      }
    }
      if (c0>c1 && c0>c2)
        {c = 0;}
      else if (c1>c0 && c1>c2)
        {c = 1;}
      else if (c2>c0 && c2>c1)
        {c = 2;}
      else
        {c = -1;}

      fprintf(out,"%f\t%d\n", pos[a], c);
  }
  Rcout<<"Smoothing finshed.\nOutput stored in "<<output<<std::endl;
  fclose(anclik);
  fclose(out);
  }
