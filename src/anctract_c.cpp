#include <Rcpp.h>
using namespace Rcpp;

//' @title anctract_c: Generating ancestry tract from call.
//'
//' @param chrlen Length of the chromosome
//' @param excludelen length to be excluded
//' @param input input file dir
//' @param output output file dir
// [[Rcpp::export]]
void anctract_c(int chrlen, int excludelen,
                std::string input, std::string output){
  int pos0, pos, state0,state,flag;
  double brk0,brk1, start, end, tractlength;
  char *chrom, *chrom0, *name, *name0;
  FILE *call, *out;
  
  chrom = new char[50];
  chrom0 = new char[50];                                                                                                                             
  name = new char[50];
  name0 = new char[50];
  
  flag=0;
  
  //open input and output file
  call=fopen(input.c_str(),"r");
  out=fopen(output.c_str(),"w");
  
  //Set starting point at the starting point
  start=excludelen;
  pos0=0;
  
  while(pos0<=start){
    if(fscanf(call,"%d\t%*d\t%s\t%s\t%*d\t%d\t%*f\t%*f",
              &pos0,chrom0,name0,&state0)==EOF){
      flag++;
      break;
    }
  }
  brk0=pos0;
  
  //Find the first tract endpoint
  while(fscanf(call,"%d\t%*d\t%s\t%s\t%*d\t%d\t%*f\t%*f",
               &pos,chrom,name,&state)!=EOF){
    if(pos<=chrlen-excludelen && state0!=state){
      brk1=pos0;
      end=(brk1-brk0)/2;
      tractlength=end-start;
      fprintf(out,"%s\t%f\t%f\t%d\t%f\t%s\n",
              chrom0,start,end,state0,tractlength,name0);
      start = end;
      chrom0 =  chrom;
      name0 = name;
      state0 = state;
      pos0 = pos;
    }
    else if(pos>chrlen-excludelen){
      end = chrlen-excludelen;
      tractlength=end-start;
      fprintf(out,"%s\t%f\t%f\t%d\t%f\t%s\n",
              chrom0,start,end,state0,tractlength,name0);
    }
  }
  
  //close file
  fclose(call);
  fclose(out);
  
}
