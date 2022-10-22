.onLoad <- function(libname, pkgname){

  package.folder <- system.file(package = 'ngvb')
  folder  <- paste0(package.folder,'/cfiles/')

  file.exists(paste0(folder,'cgeneric-ngvb.so'))

  if(!  file.exists(paste0(folder,'cgeneric-ngvb.so'))) {

    ## generate the so
    system(paste0('gcc -c -fpic -I', shQuote(folder), ' ',
                  shQuote(paste0(folder, 'cgeneric-ngvb.c')),' -o ',
                  shQuote(paste0(folder, 'cgeneric-ngvb.o'))))
    system(paste0('gcc -shared ',
                  shQuote(paste0(folder,'cgeneric-ngvb.o')),
                  ' -o ',
                  shQuote(paste0(folder, 'cgeneric-ngvb.so'))))
  }


}
