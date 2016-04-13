# Module:      genomeviewer.R
# Programmer:  Sean R. McCorkle
#              Biology Dept, Brookhaven National Laboratory
# Language:    R
#
# Description: prototype TopHat read alignment plot system
#         
# 

                              ###########
                              # globals #
                              ###########

range_extend <- 1000    # 1 kb
y_gene_trans <- 1
exon_thick <- 4
gene_color <- c("blue")   # add more for mutiple genes

read_y_delta = 1
read_color <- "black"
read_break_color <- "violet"

                             #############
                             # functions #
                             #############

                ########################################
                # Lower level functions (non-plotting) #
                ########################################

getrange <- function( x )    # finds the min & max, extends both
   {
    c( min( x ) - range_extend, max( x ) + range_extend )
   }

# get_gene_id( desc ) extracts a gene identifier from the mRNA description
# string.   I'm setting this up to handle multiple transcripts, so desc is 
# presumed to be a vector of strings 

get_gene_id <- function( desc )
   {
    sapply( strsplit( desc, ";" ), function( v ) substring( v[1], 4 ) )
   }

# read the useful subset of columns from the given SAM file
# Since SAM files have variable number of columns, we can't use
# read.table()

read_sam <- function( filename )
   {
    str_v <- strsplit( readLines( filename ), "\t" )

    df <- data.frame(  flag = sapply( str_v, function(r) as.numeric(r[2])),
                       chrom = sapply( str_v, function(r) r[3] ),
                       pos = sapply( str_v, function(r) as.numeric(r[4])),
                       mapq = sapply( str_v, function(r) as.numeric(r[5])),
                       cigar = sapply( str_v, function(r) r[6] ),
                       rnext = sapply( str_v, function(r) r[7] ),
                       pnext = sapply( str_v, function(r) as.numeric(r[8])),
                       tlen = sapply( str_v, function(r) as.numeric(r[9]))
                     )
    df
   }

# this returns the above style data frame sorted by position

sort_by_position <- function( df )
   {
    df[order( df$pos ),]
   }

#
# this reads one CIGAR string, ala "27M2D10M75N16M", and returns
# a list of specifiers, each a sublist of "Type" and number:
#    ( ("M", 27), ("D", 2), ("M", 10), ("N", 75), ("M", 16) )
#

parse_cigar <- function( s )
   {
    res <- list()
    while ( nchar(s) > 0 )
       {
        x1 <- regexpr( "[0-9]+[A-Z]", s )
        l <- attr( x1, "match.length")
        x2 <- x1 + l - 1
        if ( x1 != 1 )
           { stop( paste( "could not find spec in ", s ) ) }
        spec <- substr( s, x1, x2 )
        k <- substr( spec, 1, l-1 )
        type <- substring( spec, l )
        s <- substring( s, x2+1 )
        res <- c( res, list( list( type, as.numeric(k) ) ) )
       }
    return( res )
   }


                          #####################
                          # Ploting functions #
                          #####################

# this plots the thin transcript line across the whole transcript,
# and prints the ID above it, centered

plot_gene_transcript <- function( genedata )
   {
    trans_ind <- genedata$type == "mRNA"
    segments( genedata$start[trans_ind], y_gene_trans, 
              genedata$stop[trans_ind], y_gene_trans, 
              lwd=1, col=gene_color )
    text( (genedata$start[trans_ind] + genedata$stop[trans_ind]) / 2, 
           y_gene_trans, 
           #labels = as.character( genedata$symb ),
           labels = get_gene_id( as.character(genedata$desc[trans_ind]) ),
           col=gene_color,
           cex=0.6,
           adj=c(0.5,-1) )

   }

# this puts thick boxes for every exon in the gene

plot_gene_exons <- function( genedata )
   {
    exon_ind <- genedata$type == "CDS"
    segments( genedata$start[exon_ind], y_gene_trans, 
              genedata$stop[exon_ind], y_gene_trans,
              lwd=exon_thick, col=gene_color )
   }


#  plot the entire gene graphic

plot_gene <- function( genedata )
   {
    plot_gene_transcript( genedata )
    plot_gene_exons( genedata )
   }

# parse_segments( setdata ) scans through positions and cigar strings and
# returns a list of stops and starts with read codes indicating exon/intron
# for used by plot_reads() and coverage calculation

parse_segments <- function( setdata )
   {
    segs <- list()
    for ( i in 1:length( setdata$pos ) )    # for each sam read
       {
        cspecs <- parse_cigar( as.character(setdata$cigar[i]) )
        x <- setdata$pos[i]
        for ( spec in cspecs )    # for each specifier in the cigar string
           {
            if ( spec[[1]] == "M" ) {    # match N positions => exon?
                segs <- c( segs, 
                      list( list( "E", x, x + spec[[2]]  ) ) )
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "N" ) {  # skip N positions => intron?
                segs <- c( segs, 
                           list( list( "I", x, x + spec[[2]] ) ) )
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "D" ) {  # delete N positions?
                x <- x - spec[[2]]            # try backing up
            } else if ( spec[[1]] == "I" ) {  # insert N positions?
                x <- x + spec[[2]]            # try moving forward
            } else {
               stop( paste( "CIGAR spec", spec[[1]], "unknown.  S is ", as.character( setdata$cigar[i]) ) )
            }
            # print( segs[[ length(segs) ]] )
           }
       }
    segs       # okay we're done, return list
   }

# plot reads

plot_reads <- function( setdata, ystart )
   {
    last_bottom <- 0
    xstarts <- c()
    xstops <- c()
    colors <- c()
    ys <- c()
    y <- ystart

    for ( i in 1:length( setdata$pos ) )    # for each sam read
       {
        cspecs <- parse_cigar( as.character(setdata$cigar[i]) )
        x <- setdata$pos[i]
        for ( spec in cspecs )    # for each specifier in the cigar string
           {
            if ( spec[[1]] == "M" ) {    # match N positions
                ys <- c( ys,  y )
                xstarts <- c( xstarts, x )
                xstops <- c( xstops, x + spec[[2]] )
                colors <- c( colors, read_color )
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "N" ) {  # skip N positions
                ys <- c( ys,  y )
                xstarts <- c( xstarts, x )
                xstops <- c( xstops, x + spec[[2]] )
                colors <- c( colors, read_break_color )  # different color
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "D" ) {  # delete N positions?
                x <- x - spec[[2]]            # try backing up
            } else if ( spec[[1]] == "I" ) {  # insert N positions?
                x <- x + spec[[2]]            # try moving forward
            } else {
               stop( paste( "CIGAR spec", spec[[1]], "unknown.  S is ", as.character( setdata$cigar[i]) ) )
            }
            #cat( i, "  ", length(xstarts), ":", spec[[1]], "  ", spec[[2]], "  ", 
            #     xstarts[length(xstarts)], ", ", xstops[length(xstops)], 
            #     ", ", ys[length(ys)], "\n" )

           }

        y <- y + read_y_delta

        if ( xstarts[length(xstarts)] > last_bottom )
           { 
            y <- ystart 
            last_bottom <- xstops[length(xstops)] 
           }
       }
    segments( xstarts, ys, xstops, ys, col=colors )
   }
