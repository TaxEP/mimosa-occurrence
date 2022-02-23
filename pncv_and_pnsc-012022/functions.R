#===========#
# FUNCTIONS #
#===========#

#=======================================================================================#

#===============#
# replace.names #
#===============#

replace.names <- function(x, top, bottom, check.by, replace.by,
                          not.replace = replace.by, replace = replace.by,
                          print.bottom = TRUE,
                          print.replaced = TRUE,
                          return = TRUE){
  require(RecordLinkage)
  y <- x
  z <- c()
  for(i in 1:length(x)){
    for(j in 1:length(check.by))
      if(jarowinkler(x[i], check.by[j]) >= top & (x[i] %in% not.replace) == 0){
        y[i] <- replace.by
      } else if(x[i] %in% replace){
        y[i] <- replace.by
      } else if(jarowinkler(x[i], check.by[j]) >= bottom){
        z <- append(z, x[i])
      }
  } 
  if(print.bottom){
    print(list("Above bottom names" = unique(z[z %in% y & z != replace.by])))
  }
  if(print.replaced){
    print(list("Replaced names" = unique(x[c(x %in% y) == 0])))
  }
  if(return){
    return(y)
  }
}


#=======================================================================================#

#===============#
# generate.strg #
#===============#

#generate names to input in the argument 'check.by' from 'replace.names' function
generate.names <- function(strg1, strg2){
  s1 <- paste(strg1, strg2, sep = " ")
  s2 <- paste(strg2, strg1, sep = ", ")
  s3 <- paste(strg1, toupper(strg2), sep = " ")
  s4 <- paste(toupper(strg2), strg1, sep = ", ")
  return(c(s1, s2, s3, s4))
}

#=======================================================================================#

#=========#
# titling #
#=========#

titling <- function(string, exceptions = c("de", "da", "das", "do", "dos")){
  ex_title <- str_to_title(exceptions)
  ex_low <- str_to_lower(exceptions)
  ex_up <- str_to_upper(exceptions)
  exceptions <- c(ex_title, ex_low, ex_up)
  str_vec <- unlist(str_split(string, " "))
  for(i in 1:length(str_vec)){
    if(str_vec[i] %in% exceptions){
      str_vec[i] <- str_to_lower(str_vec[i])
    } else {
      str_vec[i] <- str_to_title(str_vec[i])
    }
  }
  string <- paste(str_vec, collapse = " ")
  return(string)
}

#=======================================================================================#

#=============#
# correct.mun #
#=============#

correct.mun <- function(x){
  require(stringi)
  chars <- c("Ã¢" = "â", "Ã?" = "á", "Ã¡" = "á", "Ã§" = "ç",
             "Ã#?#" = "ô", "Ã³" = "ó", "Ã©" = "é", "Ã£" = "ã", 
             "Ãº" = "ú", "Ãª" = "ê", "Ãµ" = "õ", "-" = " ",
             "." = "")
  for(i in 1:length(chars)){ #Replacing junk characters
    x <- gsub(pattern = names(chars[i]),
              x = x,
              replacement = unname(chars[i]), fixed = TRUE)
  }
  x <- trimws(x) #Removing white spaces
  for(i in 1:length(x)){
    if(!is.na(x[i])){
      x[i] <- titling(x[i]) #Converting into title format
    } 
  }
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT") #Removing accents
  for(i in 1:length(x)){ #Merging canonically equivalent strings
    if(!is.na(x[i])){
      for(j in 1:length(x)){
        if(!is.na(x[j])){
          if(stri_compare(x[i], x[j]) == 0){
            x[j] <- x[i]
          }
        }
      }
    }
  }
  x <- gsub("[^[:alnum:] ]", "", x) #Removing non-alphanumeric characters
  return(x)
}

