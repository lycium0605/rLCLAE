#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

//' anclik_c: Calculating the ancestral likelihood of a test individual
//'
//' @param sum an integer, the number of individuals contain in the input
//' @param type 1 or 2, specifying diploid/haploid data.
//' @param testid the individual to be estimate
//' @param genolik the dir to genolik file
//' @param ancfreq the dir to ancfreq file
//' @param output the dir of the output ancfreq file
//' @export
// [[Rcpp::export]]
void anclik_c(int sum, int type, int testid,
               std::string genolik, std::string ancfreq, std::string output) {

  Rcout<<"Finding "<<sum<<" individuals in the input"<<std::endl;

  int a, pos, pos_anc, flag;
  double gl1, gl2, gl3, L1, L2, L3, f1, f2, deltaf;
  FILE *geno, *freq, *out;

  geno=fopen(genolik.c_str(),"r");
  freq=fopen(ancfreq.c_str(),"r");
  out=fopen(output.c_str(),"w");
  flag=0;

  Rcout<<"Starting ancestral likelihood estimation for individual "<<testid<<std::endl;

  while (fscanf (geno,"%*s %d ", &pos) != EOF){
    //Rcout<<pos<<std::endl;
    if(fscanf(freq,"%d\t%lf\t%lf\t%lf\t%*d\t%*d", &pos_anc, &f1, &f2, &deltaf)==EOF){
      Rcout<<"Fscanf meets EOF in ancestral frequncy file."<<std::endl;
      flag=1;
    }
    while(pos != pos_anc){
      if(pos<pos_anc){
        if(fscanf(geno,"%*[^\n]%*c")==EOF){
          Rcout<<"Reaching the end to genolikelihood file."<<std::endl;
          flag=1;
          break;
        }
        if(fscanf(geno,"%*s %d ", &pos)==EOF){
          Rcout<<"Reaching the end to genolikelihood file."<<std::endl;
          flag=1;
          break;
        }
        /*
        for (a=1;a<=sum;++a){
          if(type==2){
            if(fscanf(geno,"%*lf %*lf %*lf ")==EOF){
              Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
            }
                 }
          if(type==1){
            if(fscanf(geno,"%*lf %*lf ")==EOF){
              Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
            }
                 }
          }
        */
      }
      else if(pos>pos_anc){
        if(fscanf(freq,"%d\t%lf\t%lf\t%lf\t%*d\t%*d", &pos_anc, &f1, &f2, &deltaf)==EOF){
          Rcout<<"Reaching the end of the ancestral frequency file."<<std::endl;
          flag=1;
          break;
        }
        //Rcout<<"The snp in ancestral frequency file is missing in the genotype likelihood file, please check."<<std::endl;
        //flag=1;
      }
      }
    if(flag==1){
      break;
    }
    else if(flag==0){
      for (a=1;a<testid;++a){
        if(type==2){
          if(fscanf(geno,"%*lf %*lf %*lf ")==EOF){
            Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
          }
               }
        if(type==1){
          if(fscanf(geno,"%*lf %*lf ")==EOF){
            Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
          }
            }
        }
      gl2=0;
      if(type==2){
        if(fscanf(geno,"%lf %lf %lf ", &gl1, &gl2, &gl3)==EOF){
          Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
        }
             }
      if(type==1){
        if(fscanf(geno,"%lf %lf ", &gl1, &gl3)==EOF){
          Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
        }
             }
      for (a=testid+1;a<=sum;++a){
        if(type==2){
          if(fscanf(geno,"%*lf %*lf %*lf ")==EOF){
            Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
          }
               }
        if(type==1){
          if(fscanf(geno,"%*lf %*lf ")==EOF){
            Rcout<<"Fscanf meets EOF, please recheck input!"<<std::endl;
          }
          }
        }
      if(gl1 != -1 && gl2 != -1 && gl3 != -1){
        if(type==2){
          L1 = gl1*f1*f1 + gl2*2.*f1*(1.-f1) + gl3*(1.-f1)*(1.-f1);
          L2 = gl1*f1*f2 + gl2*(f1*(1.-f2) + f2*(1.-f1)) +gl3*(1.-f1)*(1.-f2);
          L3 = gl1*f2*f2 + gl2*2.*f2*(1.-f2) + gl3*(1.-f2)*(1.-f2);
        }
        else if(type==1){
          L1=gl1*f1+gl3*(1.-f1);
          L2=0;
          L3=gl1*f2+gl3*(1.-f2);
        }
        fprintf(out,"%d\t%lf\t%lf\t%lf\t%lf\n",pos, deltaf, L1, L2, L3);
      }
    }
  }
  fclose(geno);
  fclose(freq);
  fclose(out);
  Rcout<<"Ancestral likelihood estimation finished."<<std::endl;
}
