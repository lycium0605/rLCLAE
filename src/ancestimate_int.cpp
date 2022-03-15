#include <Rcpp.h>
using namespace Rcpp;

//' @title anccall_c_int: Generating ancestry call along the chromosome
//'
//' @param deltaf The cut-off for minimum delta f between two reference populations
//' @param window Size of the sliding window
//' @param SMAX The number of all snp sites
//' @param anclikdir the dir to anclik file
//' @param int1 output intermediate file 1 call
//' @param int2 output intermediate file 2 after 1st smoothing
//' @param output the dir of the output anccall file
// [[Rcpp::export]]
void anccall_c_int(double deltaf, int window, int SMAX,
              std::string anclikdir,
              std::string int1, std::string int2,
              std::string output,
              std::string chrom, std::string indiv,
              double mode, int n) {

  int a, c1, c2, c0, b, c, s, nsum, *anc, *anc2, *majmode;
  double *pos, *lik2, *lik1, *lik0, p0, p1, p2, d, delta, t, posit;
  FILE *anclik, *out, *inter1, *inter2;
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
  inter1=fopen(int1.c_str(),"w");
  inter2=fopen(int2.c_str(),"w");
  s=0;
  while (fscanf(anclik,"%lf\t%lf\t%lf\t%lf\t%lf\n", &posit, &delta, &p2, &p1, &p0) != EOF) {
    if ((delta > (d - 1.E-10)) && (p0>-1.E-10 && p2>-1.E-10)) { //for hap p1==0
      lik2[s] = log (p2+1.E-6); //avoid 0
      lik1[s] = log (p1+1.E-6);
      lik0[s] = log (p0+1.E-6);
      pos[s] = posit;
      ++s;
    }
    //else{
      //++s;
    //}
  }
  Rcout<<"Scanning finished, "<<s<<" informative sites found."<<std::endl;

  for (a=0; a<s; ++a) {
    p0 = p1 = p2 = 0.;

    for (b=0; b<s && pos[b]<(pos[a]+t*.5); ++b){
      if (pos[b]>(pos[a]-t*.5)) {
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
      fprintf(inter1,"%f\t%d\n",pos[a], anc[a]);
  }
  Rcout<<"First round of ancestry call finished. Starting smoothing with the sliding window..."<<std::endl;
    for (a=0; a<s; ++a) {
    c1 = c2 = c0 = 0;
    for (b=0; b<s; ++b){
      if ((pos[a]-pos[b]) < t*.5 && (pos[a]-pos[b])>(t*(-.5))) {
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
      fprintf(inter2,"%f\t%d\n",pos[a], anc2[a]);
  }
  Rcout<<"Second round of smoothing finshed."<<std::endl;
  for (a=0; a<s; ++a) {
    c1 = c2 = c0 = 0;
    for (b=0; b<s; ++b){
      //using double length window
      if ((pos[a]-pos[b]) < t && (pos[a]-pos[b])>(t*(-1))) {
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
    //if (c!=-1 && nsum>=n && majmode[c]/nsum >= mode){
    if (nsum>=n && majmode[c]/nsum >= mode){
      //snp,call,chrom,indiv,n,mode,n_mode,perc
      //fprintf(out,"%f\t%d\t%s\t%s\t%d\t%d\t%d\t%f\n",
              //pos[a], anc2[a],chrom.c_str(),indiv.c_str(),nsum,c,majmode[c],double(majmode[c])/double(nsum));
      fprintf(out,"%f\t%d\n",pos[a], c);
    }
  }
  Rcout<<"Third round of smoothing / filtering finshed.\n Filter used:\n"<<
    "Minimum required snp in a window:"<< n <<
    "\n Minimum percentage of the majority call:"<< mode <<
    "\n Output stored in "<<output<<std::endl;
  fclose(anclik);
  fclose(out);
  fclose(inter1);
  fclose(inter2);
  }
