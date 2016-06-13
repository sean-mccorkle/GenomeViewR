# Module:      multiview.R
# Programmer:  Sean R. McCorkle
# Language:    R
#
# Description: routines for creating multi-track coverage plots
#
#              There are three routines intended for the user here:
#
#           1) multiview_by_transcript( transcript_name, extent )
#
#                      optional extent may be string percentage "20%" or
#                      numeric length in nucleotides (both indicate the 
#                      extention beyond transcript start and stop, default 
#                      is "10%")
#
#              examples:
#
#              multiview_by_transcript( "BnaC09g12820D" )
#              multiview_by_transcript( "BnaC09g12820D", extent )
#
#           2) multiview_region( chromsome, start, stop )
#
#                      plot the region of given chromosome between start and stop
#
#              example:
#
#              multiview_region( "chrA01", 15000000, 15075000 )
#
#
#           3) multiview_by_position( chromosome, target_position, extent )
#
#                      plot the region +/- extent on either side of 
#                      target_position on chromosome.
#
#              example:
#
#              multiview_by_position( "chrC09", 7345679, 20000 )
#
#

source( "/Projects/Plants/Brassica_napus/Graphics/GenomeViewR/genomeviewer.R" )
    

# multiview_by_pos( chromosome, target_position, extent )

#
# this only reads genefile into table genedata if it doesn't exist,
# or force=TRUE
#
load_genefile <- function( genefile, force = FALSE )
   {
    if ( force || ! exists( "genedata" ) )
       {
        cat( "loading genefile", genefile, "\n" )
        genedata <<- read.table( genefile,
                                 col.names = c( "chrom","caller","type","start","stop",
                                                "score", "d1", "d2", "descrip" )
                               )
       }
   }

#
#  (In this routine, extent must be supplied and must be in nt)
# 
multiview_by_pos <- function( chromosome, target_position, extent )
   {
    pos_range <- target_position + c( -extent, extent )

    #pos_range <- c( floor(target_position/1000000), ceiling(target_position/1000000) ) * 1000000
    #
    # gene annotation first if not loaded already
    #
    load_genefile( "/Projects/Plants/Brassica_napus/Brassica_napus.annotation_v5.gff3" )

    #
    # six unguided directories
    #
    unguided_dirs <- 
       paste( "/Projects/Plants/Brassica_napus/Tophat_hpc1/tophat_out_allchrom_", 
             # c( "580-1", "580-2", "580-3", "581-1", "581-2", "581-3" ), 
              c( "580-1", "581-2", "581-3", "581-1", "580-2", "580-3" ),   # regrouped
             sep="" )

    #
    # six guided directories
    #
    guided_dirs <- 
       paste( "/Projects/Plants/Brassica_napus/Tophat_hpc1/tophat_out_guided_allchrom_",
              # as.character( 1:6 ),
              as.character( c( 1, 5, 6, 4, 2, 3 ) ),  # regrouped
              sep="" )

    #
    # read unguided data 
    #
    unguided_setdata <- as.list( rep( "", 6 ) )
  
    for ( i in 1:6 )
       {
        samindex <- load_sam_index( paste( unguided_dirs[i], "sorted_accepted_hits.index", sep="/" ) )
        con <- open_samfile_read( paste( unguided_dirs[i], "sorted_accepted_hits.sam", sep="/" ) )
        lines <- extract_sam_record_block( samindex, con, chromosome, target_position )
        cat( lines[1], "\n" )
        cat( lines[length(lines)], "\n" )
        unguided_setdata[[i]] <- parse_sam_lines( lines, range = pos_range )
        close( con )
        print( head( unguided_setdata[[i]] ) )
       }
    
    unguided_segs <- as.list( rep( NA, 6 ) )
    
    for ( i in 1:6 )
       {
        unguided_segs[[i]] <- parse_segments( unguided_setdata[[i]] )
        cat( "unguided_segs[[", i,"]]\n" )
        print( head( unguided_segs[[i]] ) )
       }
    
    #
    # read guided data
    #
    guided_setdata <- as.list( rep( NA, 6 ) )
    
    for ( i in 1:6 )
       {
        samindex <- load_sam_index( paste( guided_dirs[i], "sorted_accepted_hits.index", sep="/" ) )
        con <- open_samfile_read( paste( guided_dirs[i], "sorted_accepted_hits.sam", sep="/" ) )
        lines <- extract_sam_record_block( samindex, con, chromosome, target_position )
        cat( lines[1], "\n" )
        cat( lines[length(lines)], "\n" )
        guided_setdata[[i]] <- parse_sam_lines( lines, range = pos_range )
        close( con )
        print( head( guided_setdata[[i]] ) )
       }
    
    guided_segs <- as.list( rep( NA, 6 ) )
    
    for ( i in 1:6 )
       {
        guided_segs[[i]] <- parse_segments( guided_setdata[[i]] )
        cat( "guided_segs[[", i,"]]\n" )
        print( head( guided_segs[[i]] ) )
       }
    
    # First, pre-compute coverage arrays for all tracks, so we can ascertain
    # maxima and minima
    
    unguided_cov <- as.list( rep( NA, 6 ) )
    guided_cov <- as.list( rep( NA, 6 ) )
    
    maxcov <- rep( 0, 6 )
    
    for ( i in 1:6 ) 
       {
        unguided_cov[[i]] <- compute_coverage( unguided_segs[[i]] )
        guided_cov[[i]]   <- compute_coverage( guided_segs[[i]]   )
        maxcov[i] <- max( unguided_cov[[i]], guided_cov[[i]] )
        cat( "maxcov[",i,"] is ", maxcov[i], "\n" )
       }
    
    gene <- ""
    ysep <- round( max( maxcov ) * 1.05 )
    cat( "ysep is ", ysep, "\n" )
    ylevs = ((1:6) - 0.75) * ysep
    ymax <- ylevs[6] + ysep
    
    plot( 0, 0, xlim=pos_range, ylim=c(1,ymax), type="n",                # inital empty plot
          main = paste( gene, "     chromosome", chromosome ),
          xlab = "position (nt)",
          ylab = "coverage"
         )
    
    #
    # Plot the coverage histograms, use pre-computed coverage arrays
    #
    for ( i in 1:6 )
       {
        plot_coverage( unguided_segs[[i]], 
                       ylevs[i], 
                       col="red",
                       cov = unguided_cov[[i]]
                     )
        plot_coverage( guided_segs[[i]], 
                       ylevs[i], 
                       col="black",
                       cov = guided_cov[[i]]
                     )
       }
    
    text( pos_range[1], ylevs, 
          # labels = c("580-1","580-2","580-3","581-1","581-2","581-3"),
          labels = c( "580-1", "581-2", "581-3", "581-1", "580-2", "580-3" ),  # regrouped
          adj=c(0,0),
     )

    plot_gene( genedata[   genedata$chrom == chromosome 
                         & genedata$start >= pos_range[1] 
                         & genedata$stop  <= pos_range[2], 
                       ] 
             )
    print( maxcov )

   }   # end multiview_by_position

#
#  this makes the plot given chromesome and start,stop
# 
multiview_region <- function( chromosome, start, stop )
   {
    center <- (start + stop) / 2
    ext <- (stop - start) / 2
    multiview_by_pos( chromosome, center, ext )
   }



#
# here, extent is the plot extension beyond both ends of transcript
# start and stop
#

lookup_transcript <- function( trans )
   {
    #maybe we should precompute some of these big tables?
    # or service multiple transcript ids?

    transdata <- genedata[genedata$type=="mRNA",c("chrom","start","stop","descrip")]

    # this is taking ~ 4 secs:
    transnames <- sapply( strsplit( as.character( transdata$descrip ), ";" ),
                          function( r ) substring( r[1], 4 ) 
                        )
    transdata[ transnames == trans, ]
   }

#
# In this routine, extent (as in extension) means extend beyond each end of the transcript.  
# It can either be a numeric length in nt, or be a string percentage, ie "10%",
# indicating fraction of transcript length on either side.
# Default is "%10%".
#

multiview_by_transcript <- function( transcript, extent="10%" )
   {
    #
    # gene annotation first if not loaded already
    #
    load_genefile( "/Projects/Plants/Brassica_napus/Brassica_napus.annotation_v5.gff3" )

    locus <- lookup_transcript( transcript )

    print( locus )

    nt_extent <- 0  # extent in nucleotides

    if ( is.numeric( extent ) ) {
        nt_extent <- extent
    } else {
        if ( ! is.character( extent ) )
            stop( "multiview_by_transcript: extent must be numeric or string percentage" )
        n <- nchar( extent )
        if ( substr( extent, n, n ) == "%" )
            n <- n - 1
        nt_extent <- ( locus$stop - locus$start ) * (as.numeric( substr( extent, 1, n ) ) / 100 )
    }
    
    cat( "multiview_by_pos", as.character(locus$chrom), " ",
                      (locus$stop + locus$start) / 2, " ", 
                      (locus$stop - locus$start) / 2 + nt_extent , "\n" )

    multiview_by_pos( as.character(locus$chrom),
                      (locus$stop + locus$start) / 2, 
                      (locus$stop - locus$start) / 2 + nt_extent 
                    )
   }
    

