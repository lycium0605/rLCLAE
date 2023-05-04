#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>

using namespace Rcpp;

//' get_genolik_c: Get genolik directly from vcf
//'
//' @param input_file a string, the dir to cleaned vcf file
//' @param output_file a string, the dir to the output file
// [[Rcpp::export]]
void get_genolik_c(std::string input_file, std::string output_file) {
  std::ifstream file(input_file);
  std::ofstream outfile(output_file);
  std::string line;
  double l1, l2, l3;

  if (file.is_open() && outfile.is_open()) {
    while (std::getline(file, line)) {
      // Skip lines that start with a '#'
      if (line[0] == '#') {
        continue;
      }

      std::stringstream ss(line);
      std::string field1, field2, field9;
      std::vector<std::string> field9_parts;

      std::getline(ss, field1, '\t');
      std::getline(ss, field2, '\t');

      // Skip fields 3-8
      for (int i = 0; i < 6; i++) {
        std::string temp;
        std::getline(ss, temp, '\t');
      }

      std::getline(ss, field9, '\t');

      // Split field 9 by ':'
      std::stringstream field9_ss(field9);
      std::string part;
      while (std::getline(field9_ss, part, ':')) {
        field9_parts.push_back(part);
      }

      int n = field9_parts.size();

      // Check if "PL" is in field9_parts
      bool has_pl = false;
      int pl_idx;
      for (int i = 0; i < n; i++) {
        if (field9_parts[i] == "PL") {
          has_pl = true;
          pl_idx = i;
          break;
        }
      }

      if(has_pl){
        outfile << field1 << " " << field2;
        // Get the remaining fields until the end of the line
        std::vector<std::string> remaining_fields;
        while (std::getline(ss, field9, '\t')) {
          remaining_fields.push_back(field9);
        }

        // Process each field9 in remaining_fields
        for (int i = 0; i < remaining_fields.size(); i++) {
          std::stringstream field9_ss(remaining_fields[i]);
          std::vector<std::string> field9_parts;
          while (std::getline(field9_ss, part, ':')) {
            field9_parts.push_back(part);
          }
          if (field9_parts.size() > pl_idx ) {
            std::stringstream pl_ss(field9_parts[pl_idx]);
            std::vector<std::string> pl_parts;
            while (std::getline(pl_ss, part, ',')) {
              pl_parts.push_back(part);
            }
            try{
              // std::cout<<pl_parts[0]<<std::endl;
              int gl1 = std::stoi(pl_parts[0]);
              // std::cout<<pl_parts[1]<<std::endl;
              int gl2 = std::stod(pl_parts[1]);
              // std::cout<<pl_parts[2]<<std::endl;
              int gl3 = std::stod(pl_parts[2]);
              l1 = pow (10., -gl1/10.);
              l2 = pow (10., -gl2/10.);
              l3 = pow (10., -gl3/10.);

              if(l1+l2+l3==3){
                outfile <<" -1. -1. -1.";
              }else{
                outfile << std::fixed << std::setprecision(6) << " " << l1/(l1+l2+l3) << " "<< l2/(l1+l2+l3) << " "<< l3/(l1+l2+l3);
                              }
            }catch(const std::invalid_argument& e){
              std::cerr << "Error: " << e.what() << std::endl;
            }




          } else {
            perror("Didn't find field expected");
          }
        }
        outfile << std::endl;

      }else{
        perror("PL field not found");
      }

    }

    file.close();
    outfile.close();
  }
}
