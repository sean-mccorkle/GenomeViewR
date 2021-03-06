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
y_gene_trans <- 0
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


                           ###################
                           # File operations #
                           ###################

#
# read index file (generated by build_index.pl for the samfile) and 
# return it as a data table.
#
load_sam_index <- function( filename )
   {
    cat( "load sam index ", filename, "\n" )
    read.table( filename,
                colClasses=c("character","integer","integer","numeric"), 
                col.names = c("chrom","binstart","line_num","filepos" ) 
              )
   }



# given SAMindex table, find index of lower bound of chromosome block 
# containing position pos.  Returns index i for samdindex[i,]

find_range_index <- function( samindex, chrom, pos )
   {
    # TODO - handle last
    #n <- length( samindex$chrom )
    #achrom <- samindex$chrom[2:n]        # vector of chroms offset by one
    #abinstart <-samindex$binstart[2:n]   # vector of binstarts offset by o

    #cat( "find_range_index n =", n, "\n" )
    #(1:(n-1))[    samindex$chrom[1:(n-1)] == chrom 
    #            & samindex$binstart[1:(n-1)] <= pos
    #            & abinstart >=pos  
    #         ]

    cat( "find range index for ", chrom, " ", pos, "\n" )
    n <- length( samindex$binstart )

    cind <- (1:n)[samindex$chrom == chrom]
    cat( "chrom indices are ", cind, "\n" )
    i <- min(cind)
    k <- max(cind)
    while ( i <= k && samindex$binstart[i+1] <= pos )
        i <- i + 1
    cat( "got: ", i, "\n" )
    print( samindex[i,] )
    i
   }

#
# open SAMfile for reading and seeking.  Returns an R connection
#
open_samfile_read <- function( filename )
   {
    cat( "open samfile read ", filename, "\n" )
    con <- file( filename )
    open( con, "r" )
    con
   }

#
# given a chromosome and position, return all reads from the SAMfile within
# the 1000000 block containing that read
#
#  Note to self:  extend with two blocks if too close to edge
#
extract_sam_record_block <- function( samindex, samfilecon, chrom, pos )
   {
    k <- find_range_index( samindex, chrom, pos )
    cat( "isSeekable? ", isSeekable(samfilecon), "\n" )

    if ( ! isSeekable(samfilecon) )
       {  stop( "samfile is not seekable\n" ); }

    nlines <- samindex$line_num[k+1] - samindex$line_num[k]
    cat( "k is ", k, " nlines is ", nlines, "\n" );

    seek( samfilecon, samindex$filepos[k] )

    readLines( samfilecon, n = nlines ) 
   }


# read the useful subset of columns from the given SAM file
# Since SAM files have variable number of columns, we can't use
# read.table()

parse_sam_lines <- function( samlines, range = c(0, 1e10) )
   {
    str_v <- strsplit( samlines, "\t" )

    d <- data.frame(  flag = sapply( str_v, function(r) as.numeric(r[2])),
                      chrom = sapply( str_v, function(r) r[3] ),
                      pos = sapply( str_v, function(r) as.numeric(r[4])),
                      mapq = sapply( str_v, function(r) as.numeric(r[5])),
                      cigar = sapply( str_v, function(r) r[6] ),
                      rnext = sapply( str_v, function(r) r[7] ),
                      pnext = sapply( str_v, function(r) as.numeric(r[8])),
                      tlen = sapply( str_v, function(r) as.numeric(r[9]))
                   )
    cat( "parse_sam_lines, main length ", length( d$pos ), " in-range ",
         length( d$pos[ d$pos >= range[1] & d$pos <= range[2] ] ),
         "\n" )
    d[ d$pos >= range[1] & d$pos <= range[2], ]
   }

read_sam <- function( filename )
   {
    parse_sam_lines( readLines( filename ) )
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
    if ( is.na( s ) ) 
        return( NA )
    res <- list()
    while ( nchar(s) > 0 )
       {
        x1 <- regexpr( "[0-9]+[A-Z]", s )
        l <- attr( x1, "match.length")
        x2 <- x1 + l - 1
        if ( (! is.numeric(x1)) || x1 != 1 )
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
              lwd=2, col=gene_color )
    strands <- as.character( genedata$d1[trans_ind] )
    cat( "strands are ", strands, "\n" )

    # label top strand ("+") transcripts above the transcript line
    if ( sum( strands=="+" ) > 0 )
       text( (genedata$start[trans_ind] + genedata$stop[trans_ind])[strands=="+"] / 2, 
           y_gene_trans + 5, 
           #labels = as.character( genedata$symb ),
           labels = get_gene_id( as.character(genedata$desc[trans_ind][strands=="+"]) ),
           col=gene_color,
           cex=0.6,
           adj=c(0.5,0) )

    # label bottom strand ("-") transcripts below the transcript line
    if ( sum( strands=="-" ) > 0 )
      text( (genedata$start[trans_ind] + genedata$stop[trans_ind])[strands=="-"] / 2, 
           y_gene_trans - 5, 
           #labels = as.character( genedata$symb ),
           labels = get_gene_id( as.character(genedata$desc[trans_ind][strands=="-"]) ),
           col=gene_color,
           cex=0.6,
           adj=c(0.5,1) )

   }


# this puts thick boxes for every exon in the gene

plot_gene_exons <- function( genedata )
   {
    if ( is.data.frame( genedata ) && nrow( genedata ) == 0 ) 
        return()
    exon_ind <- genedata$type == "CDS"
    if ( sum( exon_ind ) <= 0 )
        return()
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


plot_blastx_hsps <- function( blastxdata )
   {
    if ( is.data.frame( blastxdata ) && nrow( blastxdata ) == 0 ) 
        return()
    #cat( "in plot_blastx_hsps(), blastxdata is:\n" )
    #print( str( blastxdata ) )
    #print( blastxdata )
    #cat( "nrow( blastxdata ) == ", nrow( blastxdata ), "\n" )
    segments( blastxdata$start + blastxdata$q_start, y_gene_trans,
              blastxdata$start + blastxdata$q_end, y_gene_trans,
              lwd = 2, 
              col = gene_color )
   }

#
# make a legend box with some blastx hit results
#
list_blastx_hsps <- function( blastxdata )
   {
    if ( is.data.frame( blastxdata ) && nrow( blastxdata ) == 0 ) 
        return()
    save_family <- par()$family
    par(family="mono")
    legend( "topright", 
            legend = paste( blastxdata$gene, 
                            format( blastxdata$percent_id, width=7, justify="right"), 
                            format( blastxdata$alignment_length, width=4, justify="right" ),
                            sep = ""
                         ),
            text.col = gene_color,
            cex = 0.6,
            bty = "n"
          )
    par(family=save_family)
   }


# parse_segments( setdata ) scans through positions and cigar strings and
# returns a data frame, a "segs" table, of stops and starts with read codes 
# indicating exon/intron for used by plot_reads() and coverage calculation

# returns NA if no data in setdata

parse_segments <- function( setdata )
   {
    n_reads <- length( setdata$pos ) 
    if ( n_reads < 1 )
        return( NA )

    n_init <- 4 * n_reads                  # make initial guess (overestimate) for num rows

    cat( "creating matrix ", n_init, "\n" )
    segs <- matrix( data = NA, nrow = n_init, ncol = 4 )
    cat( "created matrix\n" )

    j <- 1
    for ( i in 1:n_reads )                         # for each sam read
       {
        if ( i %% 100 == 0 ) cat( "read ", i, "\n" );
        cspecs <- parse_cigar( as.character( setdata$cigar[i] ) )
        x <- setdata$pos[i]
        for ( spec in cspecs )    # for each specifier in the cigar string
           {
            if ( spec[[1]] == "M" ) {    # match N positions => exon?
                #segs[j,] <- c( "E", i, x, x + spec[[2]] )
                segs[j,] <- c( 1, i, x, x + spec[[2]] )   # 1 means "E"
                j <- j + 1
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "N" ) {  # skip N positions => intron?
                #segs[j,] <- c( "I", i, x, x + spec[[2]] )
                segs[j,] <- c( 2, i, x, x + spec[[2]] )    # 2 means "I"
                j <- j + 1
                x <- x + spec[[2]] + 1
            } else if ( spec[[1]] == "D" ) {  # delete N positions?
                x <- x - spec[[2]]            # try backing up
            } else if ( spec[[1]] == "I" ) {  # insert N positions?
                x <- x + spec[[2]]            # try moving forward
            } else {
               stop( paste( "CIGAR spec", spec[[1]], "unknown.  S is ", as.character( setdata$cigar[i]) ) )
            }
            # print( segs[[ length(segs) ]] )
            if ( j > n_init )
               stop( paste( "overran", n_init, "initial records in segs\n" ) )
           }
       }
    clip <- 1:(j-1)
    cds <- rep( "E", j-1 )         # okay we're done, return clipped data frame
    cds[ segs[clip,1] == 2 ] <- "I"
    data.frame( code = cds,
                readnum = segs[clip,2],
                start = segs[clip,3],
                stop = segs[clip,4] )
   }


# coverage plots

# get_setdata_xrange() returns a vector of length two holding the maximum and minimum
# of start & stop positions in the given segs table

get_setdata_xrange <- function( segs )
   {
    return( c( min( segs$start, segs$stop ), max( segs$start, segs$stop ) ) )
   }

#
# given a "segs" table, compute and return a single vector of the coverage.
# If the xrange not supplied, it is calculated from segs
#
# returns NA if segs is NA
#
compute_coverage <- function( segs, xrange=NA )
   {
    if ( is.na( segs ) )
       return( NA )
       
    if ( is.na( xrange ) )
        xrange <- get_setdata_xrange( segs )

    xoff <- xrange[1] - 1                    # what if zero?
    len <- xrange[2] - xrange[1] + 1
 
    cov <- rep( 0, len )
    ex_ind <-  segs$code == "E"
    starts <- segs$start[ex_ind] - xoff;
    stops <- segs$stop[ex_ind] - xoff;

    n_ex <- length( starts )
    cat( n_ex, "exons\n" )
    for ( i in 1:n_ex )
       {
        #if ( i %% 100 == 0 )
        #   {
        #    cat( i, ":", starts[i], " ", stops[i], "\n" )
        #    #print( cov[ starts[i]:stops[i]
        #   }
        cov[ (starts[i]):(stops[i]) ] <- cov[ (starts[i]):(stops[i]) ] + 1
       }
    cov
   }

plot_coverage <- function( segs, yoffset, col="black", cov=NA )
   {
    if ( is.na( segs ) )
        return
    print( head(segs) )
    xrange <- get_setdata_xrange( segs )

    if ( is.na( cov ) )
        cov <- compute_coverage( segs, xrange )
    # add in extra endpoint zeros to spoof the plot program into 
    # extending the plot to the windo boundaries
    points( c( 0, xrange[1], (xrange[1]):(xrange[2]), xrange[2], 1e9 ), 
            c( 0, 0,         cov,                     0, 0 ) + yoffset,
            type="s", col=col )
   }

# plot reads

plot_reads <- function( setdata, segs, ystart )
   {
    print( head(segs) )
    #nreads <- length( setdata$pos )
    nreads <- 100
    y <- ystart

    xstarts <- c()
    xstops <- c()
    ys <- c()
    colors <- c()

    for ( i in 1:nreads )    # for each sam read
       {
        if ( i %% 100 == 0 ) cat( "read ", i, "\n" );
        rsegs <- segs[ segs$readnum == i, ]    # make subtable of rows for this one read

        m <- length( rsegs$code )

        rcols <- rep( read_color, m )
        rcols[ rsegs$code == "u" ] <- read_break_color;
         
        xstarts <- c( xstarts, rsegs$start )
        xstops  <- c( xstops,  rsegs$stop )
        ys      <- c( ys, rep( y, m ) )
        colors  <- c( colors, rcols )
         
        y <- y + read_y_delta

        # if ( xstarts[length(xstarts)] > last_bottom )
        #   { 
        #    y <- ystart 
        #    last_bottom <- xstops[length(xstops)] 
        #   }
       }
    print( xstarts )
    print( xstops ) 
    print( ys )
    print( colors )
    segments( xstarts, ys, xstops, ys, col=colors )
   }


old_plot_reads <- function( setdata, ystart )
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
