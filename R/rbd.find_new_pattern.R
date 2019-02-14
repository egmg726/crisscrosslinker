

rbd.find_new_pattern <- function(pattern_num,list_of_cleaves,cleave_offset,direction=c('less','more')){

  new_pattern_num <- pattern_num

  #make a boolean statement up here before the rest of the
  #will depend on the direction statement
  #return NULL?
  #can turn character (min or max)
  if((direction == "less") && (new_pattern_num > max(list_of_cleaves))){
    return('max') #if the end pattern is beyond the last letter --> return 'max'
  } else if((direction == "more") && (new_pattern_num < min(list_of_cleaves))){
    return('min') #if the start pattern is beyond the last letter --> return 'min'
  }

  if(!(pattern_num %in% list_of_cleaves)){

    #list_of_cleaves_sub <- list_of_cleaves[pattern_num > list_of_cleaves]
    pattern_num2 <- pattern_num
    cleave_site_found <- FALSE
    first_round <- TRUE


    #need to first check if the difference between pattern_num2 and

    while(cleave_site_found == FALSE){
      if(direction == "more"){
        list_of_cleaves_sub <- list_of_cleaves[pattern_num2 > list_of_cleaves]
        last_cleave <- tail(list_of_cleaves_sub,n=1) #last integer
        last_cleave2 <- tail(list_of_cleaves_sub,n=2)[1] #second to last integer

        # if(((pattern_num2-last_cleave) < cleave_offset) && first_round == TRUE){
        #   #if the difference is more, use the last one
        #   cleave_site_found <- TRUE
        #   new_pattern_num <- last_cleave
        #   break
        # }

      } else if(direction == "less"){
        list_of_cleaves_sub <- list_of_cleaves[pattern_num2 < list_of_cleaves]
        last_cleave2 <- head(list_of_cleaves_sub,n=1) #last integer
        last_cleave <- head(list_of_cleaves_sub,n=2)[2] #second to last integer
      } else {
        warning("No direction chosen --> Returning NULL")
        return(NULL)
      }

      if(is.na(last_cleave) || is.na(last_cleave2)){
        last_cleave <- pattern_num2
        cleave_site_found <- TRUE
        if(direction == "more"){
          new_pattern_num <- min(list_of_cleaves)
        } else if(direction == "less"){
          new_pattern_num <- max(list_of_cleaves)
        }
        break
      }

      if((last_cleave-last_cleave2) <= cleave_offset){
        #if TRUE, the difference is less, so the start pattern should be changed to the last_cleave
        pattern_num2 <- last_cleave
        first_round <- FALSE
        next
        #need to account for the fact that need to take last_cleave2 if on more than the second iteration
        #so if this
      } else {#end if((last_cleave-last_cleave2) > cleave_offset){

        cleave_site_found <- TRUE
        if(direction == "more"){
          new_pattern_num <- list_of_cleaves[match(last_cleave,list_of_cleaves)]+1
        } else if(direction == "less"){
          new_pattern_num <- list_of_cleaves[match(last_cleave,list_of_cleaves)-1]
        }


      }#end else to if((last_cleave-last_cleave2) > cleave_offset){

    } #end while(cleave_site_found == FALSE){
  } #end if(!(pattern_num %in% list_of_cleaves)){

  if(direction == "more"){
    if(pattern_num < new_pattern_num){
      new_pattern_num <- pattern_num
    }
  } else if(direction == "less"){
    if(pattern_num > new_pattern_num){
      new_pattern_num <- pattern_num
    }
  } else {
    warning("No direction chosen --> Returning NULL")
    return(NULL)
  }


  return(new_pattern_num)

} #end function rbd.find_new_pattern

