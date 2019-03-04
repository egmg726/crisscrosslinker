#'PPI: Load Data
#'
#'This function loads the data from your XL-MS experiments to a format that is readable by ppi.combineData
#'
#'@param xldata Name of files or loaded data.frames. If using pLink2 data, it is recommended to use the name of file.
#'@param fasta_file FASTA file as loaded by seqinr::read.fasta()
#'@param file_directory Directory where XL-MS data is located. Defaults to getwd().
#'@param datatype Type of data being loaded. Will automatically detect otherwise. Defaults to NULL.
#'@export

ppi.loadData <- function(xldata, fasta_file, file_directory = getwd(), datatype = NULL){

  #should also load fasta_file here if not already loaded
  list_of_protein_names <- names(fasta_file)
  #typeof


  xldata_loaded <- list()
  #need to initially keep a list of names and dataframes and then go through them?
  #replace the dataframes with the real data structure afterwards?

  if(typeof(xldata) == 'character'){
    #if character --> need to load the files
    #should allow the user to add a directory

    for(file_name in xldata){

      #need to choose the command based on the final characters of the file_name
      #endsWith --> if statement, turn into a

      #will need to make this a separate function --> need to deal with
      #this --> should implement a tryerror loop since there is an issue with
      #reading CSV files from the filtered
      #should immediately make into the xlink format then if this is the case

      #can check here first to see if the file name ends with the plink suffixes


      if(!(file_name %in% list.files(file_directory))){
        warning(paste('File not found in specified directory, skipping',file_name,'\n'))
        next
      }

      if(endsWith(file_name,'csv')){

        #can turn this part into a function


        #should check if file actually exists up here before attempting to load it??


        #can have boolean up here to see what kind of data it is if it has not
        #been already pre-defined

        loaded_data <- tryCatch(read.csv(paste0(file_directory,'/',file_name)),
                                error = function(err){
                                  cat('Loaded as readline\n')
                                  return(readLines(file_name))})

        #instead of boolean, can just use typeof to see if it is
        #character => readLines
        #list => dataframe
        if(typeof(loaded_data) == 'character'){
          #if TRUE, this has been loaded in via readLines

          #if the file has the phrase 'filtered_cross-linked_peptides' in it OR
          if((grepl('filtered_cross-linked_peptides',file_name)) || (!(FALSE %in% (loaded_data[1:2] == c("Peptide_Order,Peptide,Peptide_Mass,Modifications,Proteins,Protein_Type",
                                                                                                         ",Spectrum_Order,Title,Charge,Precursor_Mass,Evalue,Score,Precursor_Mass_Error(Da),Precursor_Mass_Error(ppm),Alpha_Evalue,Beta_Evalue"))))){
            #If TRUE, then it is a filtered_cross-linked_peptides file and should be loaded via the
            #function


            #need the filtered crosslinked peptides function here

            xlink_list <- make_plink2_master_list(plink2_peptides = loaded_data,
                                                  list_of_protein_names = list_of_protein_names,
                                                  plink2_type_of_file = 'cross-linked',
                                                  master_list_index = 0)


            #add to bigger list of files


          } else if(grepl('filtered_loop-linked_peptides',file_name)){

            xlink_list <- make_plink2_master_list(plink2_peptides = loaded_data,
                                                  list_of_protein_names = list_of_protein_names,
                                                  plink2_type_of_file = 'loop-linked',
                                                  master_list_index = 0)

          } #then check if is a looplinked file (also readLines??)


          #grepl('filtered_loop-linked_peptides',file_name)

        } else if(typeof(loaded_data) == 'list'){ #end if(typeof(loaded_data) == 'character')
          #if TRUE, this has been loaded in via read.csv

          #will need to make this into a function so that other datatypes can be added to it
          #as time goes on

          if(paste(as.character(loaded_data[1,]),collapse = ',') == "NA,Order,Spectrum,Score,Calc_M,Delta,ppm,Modification,SampleID,Engine,Rank,Proteins"  &&
             paste(colnames(loaded_data[1:2]),collapse=',') == 'Order,Sequence'){
            #if TRUE, the file structure matches that of the regular plink1 file
            #can also have or statement if it has been designated as a plink1 file

            xlink_list <- load_plink_file(loaded_data)

            #if it has been loaded as a xlsx file that means it is not a plink2 file

          } else {

            #will need the other file types here as well that it will check and load as needed


            #if else --> currently unsupported file type
            #should give warning
            warning(paste('Unsupported xlsx file for',file_name))
            #all warnings should be in log file??

          }


        } else { #end else if(typeof(loaded_data) == 'list')
          warning('Unknown datatype loaded')
          next
        } #end else to else if(typeof(loaded_data) == 'list')


        #if the file ends with this
        #cross-linked_peptides
        #OR the column names match (first 2 lines of loaded_data)

        #check the first couple of lines??
        #have the lines already pre-loaded to check directly against them
        #can check first if it has the suffix in it first

        #so if it is
        #how to check the



        #loaded_data <- read.csv(paste0(file_directory,'/',file_name))
        # crosslinked_peptides <- readLines(crosslinked_peptides_file_name)
        # #looplinked_peptides <- readLines(looplinked_peptides_file_name)
        #
        # #can probably make the crosslinked master list not have to import the list of protein names
        # #since it doesn't use it
        # crosslinked_master_list <- make_plink2_master_list(plink2_peptides = crosslinked_peptides,
        #                                                    list_of_protein_names = list_of_protein_names,
        #                                                    plink2_type_of_file = 'cross-linked',
        #                                                    master_list_index = 0)




      } else if(endsWith(file_name,'xlsx')){
        loaded_data <- data.frame(read.xlsx(paste0(file_directory,'/',file_name)))


        #once the data has been loaded --> should just have 1 function?
        #can take the file_name as input as well?

        #for regular plink file
        if(paste(as.character(loaded_data[1,]),collapse = ',') == "NA,Order,Spectrum,Score,Calc_M,Delta,ppm,Modification,SampleID,Engine,Rank,Proteins"  &&
           paste(colnames(loaded_data[1:2]),collapse=',') == 'Order,Sequence'){
          #if TRUE, the file structure matches that of the regular plink1 file
          #can also have or statement if it has been designated as a plink1 file

          xlink_list <- load_plink_file(loaded_data)

          #if it has been loaded as a xlsx file that means it is not a plink2 file

        } else {

          #will need the other file types here as well that it will check and load as needed


          #if else --> currently unsupported file type
          #should give warning
          warning(paste('Unsupported xlsx file for',file_name))
          #all warnings should be in log file??

        }
        #have else statement for unsupported files

        #should have entire workflow (that needs to be repeated) in a function


        #check what kind of data it is (right now just plink)
        #should have way of validating if it is plink1
        #before the final loading occurs
        #can check the column names

      } else {
        warning('Unknown file type, skipping this upload')
        next
      }

      #once the xlink_list has been loaded in, it will be added to the master list
      xldata_loaded[[file_name]] <- xlink_list


    } #end for(file_name in xldata){


  } else if(typeof(xldata) == 'list'){#end if(typeof(xldata) == 'character')

    #first check if it's 1 dataframe
    #if yes, just load that
    #if no, it's a list()
    #the go through each of those (names()?) and see if it's a character or list


    if(!is.null(nrow(xldata))){
      #if TRUE, it is just the 1 dataframe

      xldata_loaded[['dataset']] <- xldata

    } else { #end if(is.null(nrow(xldata)))
      #if TRUE, it is a list()

      for(xld_num in 1:length(xldata)){

        xld <- xldata[[xld_num]]
        #then check what type xld is
        if(typeof(xld) == 'character'){


          if(endsWith(xld,'csv')){

            loaded_data <- read.csv(paste0(file_directory,'/',xld))

          } else if(endsWith(xld,'xlsx')){
            loaded_data <- data.frame(read.xlsx(paste0(file_directory,'/',xld)))

          } else {
            warning('Unknown file type, skipping this upload')
            next
          }

          #if character, need to load as done before

        } else if(typeof(xld) == 'list'){ #end if(typeof(xld) == 'character'){

          #double make sure that it is really a dataframe
          if(!is.null(nrow(xld))){

            #if TRUE --> there are rows in the list
            loaded_data <- xld

          } else { #end if(!is.null(nrow(xld))){
            warning('List within a list not permitted')
          } #end else to if(!is.null(nrow(xld))){


        } else { #end else if(typeof(xld) = 'list'){
          warning('Unknown type dectected in list')


        } #end else to else if(typeof(xld) = 'list'){

        #establish data_name based on whether or not
        #data_name

        if(is.null(names(xldata))){
          #if NULL, just use either the file names
          #or name the dataset based on the number
          data_name <- paste0('dataset',xld_num)

        } else { #end if(is.null(names(xldata))){
          #if TRUE, they do have names and should be used
          data_name <- names(xldata)[[xld_num]]

        } #end else to if(is.null(names(xldata))){



        xldata_loaded[[data_name]] <- loaded_data

      } #end for(xld_num in length(xldata))

      #might not have names
      #should go by length
      #and then check if names

    } #end else if(typeof(xldata) == 'list')
  }
  return(xldata_loaded)

  #if type of is list, should first check if it is a dataframe
  #if not, then will need to cycle through
  #and then check if it

  #users should be able to use alias for their datasets for better readability
  #if just a list (or list with no names), will use the names that have been given


  #should load them all as a list with the list names as the file names

  #should also export the names of the datasets

  #will need to keep track of the names of the datasets
  #if it's just one dataset, should just name it?
  #can also remove it

  #need to specify which type of data?
  #can check if there is a type specified, otherwise will try to automatically detect
  #can check the column names


  #step 1: check what kind of data is it
  #if it's a list, will need to iterate through it
  #if it is a dataframe

  #step 2: based on the type of data, it will be loaded so that will be

  #split up the xlink_list so that each dataset will still be assigned to the name
  #of the file that it comes from

  #return the xlink_list
} #end function ppi.loadData
