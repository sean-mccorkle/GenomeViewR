# Program:     viewer1.R
# Programmer:  Sean R. McCorkle
#              Biology Dept, Brookhaven National Laboratory
# Language:    R
#
# Description: prototype TopHat read alignment plot system
#         
# 

source( "genomevieweR.R" )


                            ################
                            # Main Program #
                            ################


viewer1 <- function( gene, xrange=NULL )
   {
    genefile <- paste( "gene", gene, "tab", sep="." )

    genedata <- read.table( genefile )
    names( genedata ) <- c( "chrom","caller","type","start","stop","score", 
                            "d1", "d2", "descrip" )

    setfile <- paste( "cut", gene, "sam", sep="." )

    setdata <- sort_by_position( read_sam( setfile ) )

    if ( is.null( xrange ) )  
       {  r <- getrange( c( genedata$start, genedata$stop ) ) }
    else
       {  r <- xrange }
    chr <- genedata$chrom[1]

    # Initialize with an empty plot

    plot( 0, 0, xlim=r, ylim=c(1,100), type="n",
           main = paste( gene, "     chromosome", chr ),
           xlab="position (nt)",
           ylab=" "
         )

    plot_gene( genedata )

    rug( setdata$pos )

    plot_reads( setdata, 15 )
   }


multi_view <- function( genes )
   {
    for ( gene in genes )
       viewer1( gene )
   }

multi_pdf <- function( genes )
   {
    for ( gene in genes )
      {
       pdffile <- paste( gene, ".pdf", sep="") 
       cat( "doing ", pdffile, "\n" )
       pdf( pdffile, height=5,width=10 )
       viewer1( gene )
       dev.off()
      }
   }
