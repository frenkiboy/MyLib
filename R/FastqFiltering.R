### INFO: Filters the fastq file on the given parameters	
### DATE: 25.11.2011.
### AUTHOR: Vedran Franke

# {0} TEST DATA
infile=
#/{0} TEST DATA


# {0} PARSING ARGUMENTS
library(optparse)
option_list = list(
				make_option(c("-i", "--infile"), type='character', action="store"),
				make_option(c("-o", "--outpath"), type='character', action="store"),
				make_option(c("-w", "--width"), type="integer", action='store', default=200),
				make_option(c("-u", "--uniq"), type="logical", action='store', default=TRUE)
				)
opt = parse_args(OptionParser(option_list = option_list))
#/{0} PARSING ARGUMENTS


# {1} LIBRARIES
#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	#/{{3}} MAIN
#/{2} CODE



