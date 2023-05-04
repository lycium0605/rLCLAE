#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>

using namespace Rcpp;

//' get_genolik_extended_c: Get genolik directly from vcf
//'
//' @param input_file a string, the dir to cleaned vcf file
//' @param output_file a string, the dir to the output file
// [[Rcpp::export]]
void get_genolik_extended_c(std::string input_file, std::string output_file, std::string format_type) {
  std::ifstream file(input_file);
  std::ofstream outfile(output_file);
  std::string line;

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
        if (field9_parts[i] == format_type) {
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
              double gl1 = std::stod(pl_parts[0]);
              // std::cout<<pl_parts[1]<<std::endl;
              double gl2 = std::stod(pl_parts[1]);
              // std::cout<<pl_parts[2]<<std::endl;
              double gl3 = std::stod(pl_parts[2]);

              if(std::fabs(gl1+gl2+gl3)<0.1){
                outfile <<" -1. -1. -1.";
              }else{
                outfile << std::fixed << std::setprecision(6) << " " << gl1 << " "<< gl2 << " "<< gl3;
              }
            }catch(const std::invalid_argument& e){
              std::cerr << "Error: " << e.what() << std::endl;
            }




          } else {
            perror("Didn't find number of field expected");
          }
        }
        outfile << std::endl;

      }else{
        perror("Input field not found");
      }

    }

    file.close();
    outfile.close();
  }
}
