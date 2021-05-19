#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

//' ancfreq_c: Calculating the allele frequency in two different ancestral population
//'
//' @param n an integer, the number of individuals contain in the input
//' @param type 1 or 2, specifying diploid/haploid data.
//' @param pop1 the dir of a file containing individuals belonging to reference population 1
//' @param pop2 the dir of a file containing individuals belonging to reference population 2
//' @param input the dir of the input genolik file
//' @param output the dir of the output ancfreq file
//' @export
// [[Rcpp::export]]
void ancfreq_c(int n, int type, std::string pop1, std::string pop2, std::string input, std::string output) {
  int a, b, test1, pos, flag, *Pop1, *Pop2, c, anc1, anc2;
  double L1, L2, L3, f1, f2, **lik, deltaf, anc1_num, anc2_num;

  Rcout<<"Finding "<<n<<" individuals in the input"<<std::endl;
  FILE *pop1_, *pop2_, *geno, *out;
  Pop1 = (int *) malloc ((n+1) * sizeof (int));
  Pop2 = (int *) malloc ((n+1) * sizeof (int));
  for (a=1; a<=n; ++a) {
    Pop1[a] = 0;
    Pop2[a] = 0;
  }
  lik = (double **) malloc ((n+1) * sizeof (double *));
  for (a=1; a<=n; ++a)
    lik[a] = (double *) malloc (3 * sizeof (double *));

  pop1_=fopen(pop1.c_str(),"r");
  std::fscanf(pop1_, "%d ", &b);

  Rcout<<"Finding "<<b<<" individuals for reference population 1"<<std::endl;
  for (a=0; a<b; ++a) {
    fscanf (pop1_, "%d ", &c);
    Pop1[c] = 1;
  }
  fclose (pop1_);

  pop2_=fopen(pop2.c_str(),"r");
  std::fscanf(pop2_, "%d ", &b);
  Rcout<<"Finding "<<b<<" individuals for reference population 2"<<std::endl;
  for (a=0; a<b; ++a) {
    fscanf (pop2_, "%d ", &c);
    Pop2[c] = 1;
  }
  fclose (pop2_);

  geno=fopen(input.c_str(),"r");
  out=fopen(output.c_str(),"w");

  Rcout<<"Starting ancestral allele frequency calculation"<<std::endl;

  if(type==2){
    while (std::fscanf (geno,"%*s %d ", &pos) != EOF) {
      flag = 0;
      for (a=1; a<=n; ++a)
        fscanf (geno,"%lf %lf %lf", &lik[a][0], &lik[a][1], &lik[a][2]);

      f1 = f2 = 0.;
      L1 = L2 = L3 = 0.;
      for (a=1; a<=n; ++a)
        if (Pop1[a]>0) {
          if (lik[a][0] > -.1)
            L1 += lik[a][0];
          if (lik[a][1] > -.1)
            L2 += lik[a][1];
          if (lik[a][2] > -.1)
            L3 += lik[a][2];
        }
        if ((L1+L2+L3)>0.){
          f1 = (L1+.5*L2) / (L1+L2+L3);
          anc1_num = L1+L2+L3;}
        else
          flag = 1;

        L1 = L2 = L3 = 0.;
        for (a=1; a<=n; ++a)
          if (Pop2[a]>0) {
            if (lik[a][0] > -.1)
              L1 += lik[a][0];
            if (lik[a][1] > -.1)
              L2 += lik[a][1];
            if (lik[a][2] > -.1)
              L3 += lik[a][2];
          }
          if ((L1+L2+L3)>0.){
            f2 = (L1+.5*L2) / (L1+L2+L3);
            anc2_num = L1+L2+L3;}
          else
            flag = 1;

          if (flag==0) {
            deltaf = (f1 > f2) ? (f1-f2) : (f2-f1) ;
            anc1 = int(anc1_num);
            anc2 = int(anc2_num);
            std::fprintf(out,"%d\t%.8lf\t%.8lf\t%.6lf\t%d\t%d\n", pos, f1, f2,deltaf, anc1, anc2);
          }
    }
  }
  else if (type==1){
    while (std::fscanf (geno,"%*s %d ", &pos) != EOF) {
      flag = 0;
      for (a=1; a<=n; ++a)
        fscanf (geno,"%lf %lf ", &lik[a][0], &lik[a][1]);

      f1 = f2 = 0.;
      L1 = L2 = L3 = 0.;
      for (a=1; a<=n; ++a)
        if (Pop1[a]>0) {
          if (lik[a][0] > -.1)
            L1 += lik[a][0];
          if (lik[a][1] > -.1)
            L3 += lik[a][1];
        }
        if ((L1+L2+L3)>0.){
          f1 = (L1+.5*L2) / (L1+L2+L3);
          anc1_num = L1+L3;}
        else
          flag = 1;

        L1 = L2 = L3 = 0.;
        for (a=1; a<=n; ++a)
          if (Pop2[a]>0) {
            if (lik[a][0] > -.1)
              L1 += lik[a][0];
            if (lik[a][1] > -.1)
              L3 += lik[a][1];
          }
          if ((L1+L2+L3)>0.){
            f2 = (L1+.5*L2) / (L1+L2+L3);
            anc2_num = L1+L3;}
          else
            flag = 1;

          if (flag==0) {
            deltaf = (f1 > f2) ? (f1-f2) : (f2-f1) ;
            anc1 = int(anc1_num);
            anc2 = int(anc2_num);
            std::fprintf(out,"%d\t%.8lf\t%.8lf\t%.6lf\t%d\t%d\n", pos, f1, f2,deltaf, anc1, anc2);
          }
    }
  }

  fclose(geno);
  fclose(out);
  Rcout<<"Ancestral allele frequency calculation finished."<<std::endl;
}

//' ancfreq_merge: Merging haploid & dipoid dataset
//'
//' @param hapfreq The dir to file containing ancfreq for haploid data
//' @param dipfreq The dir to file containing ancfreq for diploid data
//' @param type An integer, 0-intersection 1-keep all haploid snp 2-keep all diploid snp 3-union
//' @param outputdir The output dir
//' @export
// [[Rcpp::export]]
void ancfreq_merge(std::string hapfreq, std::string dipfreq, std::string outputdir, int type) {
  int pos_m, pos_f, pos, anc1_m, anc2_m,anc1_f, anc2_f,anc1,anc2;
  double f1_m, f2_m, f1_f, f2_f,f1, f2, deltaf;
  FILE *hap, *dip , *out;
  //pop1_=fopen(pop1.c_str(),"r");
  hap=fopen(hapfreq.c_str(),"r");
  dip=fopen(dipfreq.c_str(),"r");
  out=fopen(outputdir.c_str(),"w");
  while (std::fscanf (hap,"%d\t%lf\t%lf\t%d\t%d",
                      &pos_m, &f1_m, &f2_m, &anc1_m, &anc2_m)!= EOF){
    fscanf(dip,"%d\t%lf\t%lf\t%d\t%d", &pos_f, &f1_f, &f2_f, &anc1_f, &anc2_f);
    //For data-specifec snps
    while (pos_m != pos_f) {
      //For haploid specific snps
      if (pos_m < pos_f){
        if(type==1||type==3){
          fprintf(out,"%d\t%lf\t%lf\t%d\t%d",pos_m, f1_m, f2_m, anc1_m, anc2_m);
        }
        fscanf (hap,"%d\t%lf\t%lf\t%d\t%d",&pos_m, &f1_m, &f2_m, &anc1_m, &anc2_m);
      }
      //For diploid specific snps
      else if (pos_m > pos_f){
        if(type==2||type==3){
          fprintf(out,"%d\t%lf\t%lf\t%d\t%d",pos_f, f1_f, f2_f, anc1_f, anc2_f);
        }
        fscanf (dip,"%d\t%lf\t%lf\t%d\t%d",&pos_f, &f1_f, &f2_f, &anc1_f, &anc2_f);
      }
    }
    //For shared snps
    pos = pos_f;
    f1 = (f1_m * anc1_m + 2 * (f1_f * anc1_f))/(2 * anc1_f + anc1_m);
    f2 = (f2_m * anc2_m + 2 * (f2_f * anc2_f))/(2 * anc2_f + anc2_m);
    deltaf = (f1 > f2) ? (f1-f2) : (f2-f1) ;
    anc1 = anc1_f+anc1_m;
    anc2 = anc2_f+anc2_m;
    fprintf(out,"%d\t%.8lf\t%.8lf\t%.6lf\t%d\t%d\n", pos, f1, f2, deltaf,anc1,anc2);
  }
  fclose(out);
  fclose(hap);
  fclose(dip);
}

