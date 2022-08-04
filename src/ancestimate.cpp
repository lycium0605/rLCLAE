#include <Rcpp.h>
using namespace Rcpp;

//' @title anccall_c: Generating ancestry call along the chromosome
//'
//' @param deltaf The cut-off for minimum delta f between two reference populations
//' @param window Size of the sliding window
//' @param SMAX The number of all snp sites
//' @param anclikdir the dir to anclik file
//' @param output the dir of the output ancfreq file
// [[Rcpp::export]]
void anccall_c(double deltaf, int window, int SMAX,
              std::string anclikdir, std::string output,
              std::string chrom, std::string indiv,
              double mode, int n,
              double round1, double round2, double round3,
              double zero_value) {

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
  Rcout<<"Adding small values "<<zero_value<<" to avoid zeroes."<<std::endl;
  while (fscanf(anclik,"%lf\t%lf\t%lf\t%lf\t%lf\n", &posit, &delta, &p2, &p1, &p0) != EOF) {
    if ((delta > (d - 1.E-10)) && (p0>1.E-10 || p2>1.E-10)) { //for hap p1==0
      lik2[s] = log (p2+zero_value); //avoid 0
      lik1[s] = log (p1+zero_value);
      lik0[s] = log (p0+zero_value);
      pos[s] = posit;
      ++s;
    }
    //else{
      //++s;
    //}
  }
  Rcout<<"Scanning finished, "<<s<<" informative sites found."<<std::endl;

  //round1=round1/s; adjust window size based on snp density, recover if needed

  Rcout<<"Using "<<round1*t<<" as window size for round1."<<std::endl;
  for (a=0; a<s; ++a) {
    p0 = p1 = p2 = 0.;

    for (b=0; b<s && pos[b]<(pos[a]+t*.5*round1); ++b){
      if (pos[b]>(pos[a]-t*.5*round1)) {
        p0 += lik0[b];
        p1 += lik1[b];
        p2 += lik2[b];
      }
    }
      if (p0>p1 && p0>p2){
        anc[a] = 0;
      }
      else if (p1>p0 && p1>p2){
        anc[a] = 1;
        //Rcout<<"Finding het, "<<p0<<" "<<p1<<" "<<p2<<" "<<std::endl;
      }
      else if (p2>p0 && p2>p1){
        anc[a] = 2;
      }
      else{
        anc[a] = -1;
        Rcout<<"Finding outlier, "<<p0<<" "<<p1<<" "<<p2<<" "<<std::endl;
        Rcout<<"Warning, warning, outlier warning, what is happening"<<std::endl;
      }
  }
  Rcout<<"First round of ancestry call finished. Starting smoothing with the sliding window..."<<std::endl;
    for (a=0; a<s; ++a) {
    c1 = c2 = c0 = 0;
    for (b=0; b<s; ++b){
      if ((pos[a]-pos[b]) < t*.5*round2 && (pos[a]-pos[b])>(t*(-.5)*round2)) {
        if (anc[b] == 0){++c0;}

        else if (anc[b] == 1){++c1;
          //Rcout<<"Finding het in 1, "<<lik0[b]<<" "<<lik1[b]<<" "<<lik2[b]<<" "<<std::endl;
          }

        else if (anc[b] == 2){++c2;}

      }
    }
      if (c0>c1 && c0>c2)
        {//c = 0;
        anc2[a] = 0;
        }
      else if (c1>c0 && c1>c2)
        {//c = 1;
        anc2[a] = 1;
        //Rcout<<"Findint het, "<<c0<<" "<<c1<<" "<<c2<<" "<<std::endl;
        }
      else if (c2>c0 && c2>c1)
        {//c = 2;
        anc2[a] = 2;
        }
      else
        {//c = -1;
        anc2[a] = -1;
        }
      fprintf(out,"%f\t%d\n", pos[a], anc2[a]);
  }
  Rcout<<"Second round of smoothing finshed."<<std::endl;
  /*
   for (a=0; a<s; ++a) {
    c1 = c2 = c0 = 0;
    for (b=0; b<s; ++b){
      //using double length window
      if ((pos[a]-pos[b]) < t*.5*round3 && (pos[a]-pos[b])>(t*(-.5)*round3)) {
        if (anc2[b] == 0){++c0;}

        else if (anc2[b] == 1){++c1;}

        else if (anc2[b] == 2){++c2;}

      }
    }
    if (c0>c1 && c0>c2)
    {c = 0;
      //anc[a] = 0;
    }
    else if (c1>c0 && c1>c2)
    {c = 1;
      //anc[a] = 1;
    }
    else if (c2>c0 && c2>c1)
    {c = 2;
      //anc[a] = 2;
    }
    else
    {c = -1;
      //anc[a] = -1;
    }
    majmode[0] = c0;
    majmode[1] = c1;
    majmode[2] = c2;
    nsum = c0+c1+c2;
    if (c!=-1 && nsum>=n && majmode[c]/nsum >= mode){
      //snp,call,chrom,indiv,n,mode,n_mode,perc
      //fprintf(out,"%f\t%d\t%s\t%s\t%d\t%d\t%d\t%f\n",
              //pos[a], anc2[a],chrom.c_str(),indiv.c_str(),nsum,c,majmode[c],double(majmode[c])/double(nsum));
      fprintf(out,"%f\t%d\n",pos[a], c);
    }
  }
  */
  Rcout<<"No round 3.\n Filter used:\n"<<
    "Minimum required snp in a window:"<< n <<
    "\n Minimum percentage of the majority call:"<< mode <<
    "\n Output stored in "<<output<<std::endl;
  fclose(anclik);
  fclose(out);
  }
