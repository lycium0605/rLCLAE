

#'datacheck
#'A function to check if the file contains certain wrong character types
#'
#' @param input the input file dir
#' @param character the character types regular expression i.e. '[^0-9]'
#' @param field the field to be kept, default for all
#'
#' @return 0 pass 1 error
# @export
#'
#' @examples
#' datacheck(genolik,'[^0-9[:space:]]','2-')
#'
datacheck<-function(input,character="[^0-9[:space:]]",field='1-'){

  #Check characters
  cut = ' | cut -f '
  grep = ' | grep \"'
  out =  "\" > /dev/null && echo 'Unexpected character, please double check your data' || echo 'Pass character test'"
  command_ = paste(sep='', 'cat ',input,cut,field,grep,character,out)
  #print(command_)
  system(command = command_,intern = TRUE) -> charcheck
  print(charcheck)
  if(charcheck=='Pass character test'){
    return (0)
  }
  else{
    return (1)
  }

}
