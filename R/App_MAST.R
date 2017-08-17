run_MAST = function(motif_file=NULL, sequence_file=NULL,  mast_path = 'mast', strand=TRUE){

  # -w show wak matches
  # -remcoor remove correlated motifs
  # -dblist specified

  if(is.null(motif_file))
    stop('Please specify the motif file')

  if(is.null(sequence_file))
    stop('Please specify the sequence file')

  command = paste(
    mast,
    motif_file,
    sequence_file,
    '-oc',
    '-remcorr',
    '-w'
    )

  if(!strand)
    command = paste(command, '-norc')

  syste(command)
}
