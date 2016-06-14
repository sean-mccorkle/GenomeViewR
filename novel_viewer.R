# Module:      novel_viewer.R
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
#           2) novel_view_region( chromsome, start, stop )
#
#                      plot the region of given chromosome between start and stop
#
#              example:
#
#              multiview_region( "chrA01", 15000000, 15075000 )
#
#
#           3) novel_view_by_position( chromosome, target_position, extent )
#
#                      plot the region +/- extent on either side of 
#                      target_position on chromosome.
#
#              example:
#
#              multiview_by_position( "chrC09", 7345679, 20000 )
#
#

source( "genomeviewer.R" )
    
# multiview_by_pos( chromosome, target_position, extent )

#
# this only reads genefile into table genedata if it doesn't exist,
#  or force=TRUE
#
load_blastxfile <- function( blastxfile, force = FALSE )
   {
    if ( force || ! exists( "blastxdata" ) )
       {
        cat( "loading blastxfile", blastxfile, "\n" )
        blastxdata <<- read.table( blastxfile,
                                 colClasses = c( "factor",
                                                 "integer",
                                                 "integer",
                                                 "character",
                                                 "numeric",
                                                 "integer",
                                                 "integer",
                                                 "integer",
                                                 "integer",
                                                 "integer",
                                                 "integer",
                                                 "integer",
                                                 "numeric",
                                                 "numeric"
                                                ),
                                 col.names = c( "chrom",
                                                "start",
                                                "stop",
                                                "gene",
                                                "percent_id",
                                                "alignment_length", 
                                                "mismatches", 
                                                "gap_opens", 
                                                "q_start", 
                                                "q_end", 
                                                "s_start", 
                                                "s_end", 
                                                "evalue", 
                                                "bit_score"
                                              )
                            )
       }
   }




#
#
#  (In this routine, extent must be supplied and must be in nt)
# 
novel_view_by_pos <- function( chromosome, target_position, extent, title="" )
   {
    pos_range <- target_position + c( -extent, extent )

    #pos_range <- c( floor(target_position/1000000), ceiling(target_position/1000000) ) * 1000000
    #
    # gene annotation first if not loaded already
    #
    #load_genefile( "/Projects/Plants/Brassica_napus/Brassica_napus.annotation_v5.gff3" )

    #
    # six novel files
    #
    novel_files <- 
       paste( "../../no_trans_overlaps",
             # c( "580-1", "580-2", "580-3", "581-1", "581-2", "581-3" ), 
             # c( "580-1", "581-2", "581-3", "581-1", "580-2", "580-3" ),   # regrouped
             c( "1", "5", "6", "4", "2", "3" ), #regroupd 
             sep="." )
    
    #
    # read novel (unguided, no transcript overlap) data 
    #
    novel_setdata <- as.list( rep( "", 6 ) )
  
    for ( i in 1:6 )
       {
        samindex <- load_sam_index( paste( novel_files[i], "index", sep="." ) )
        con <- open_samfile_read( paste( novel_files[i], "sam", sep="." ) )
        lines <- extract_sam_record_block( samindex, con, chromosome, target_position )
        cat( lines[1], "\n" )
        cat( lines[length(lines)], "\n" )
        novel_setdata[[i]] <- parse_sam_lines( lines, range = pos_range )
        close( con )
        print( head( novel_setdata[[i]] ) )
       }
    
    novel_segs <- as.list( rep( NA, 6 ) )
    
    for ( i in 1:6 )
       {
        novel_segs[[i]] <- parse_segments( novel_setdata[[i]] )
        cat( "novel_segs[[", i,"]]\n" )
        print( head( novel_segs[[i]] ) )
       }
    
    # First, pre-compute coverage arrays for all tracks, so we can ascertain
    # maxima and minima
    
    novel_cov <- as.list( rep( NA, 6 ) )
    
    maxcov <- rep( 0, 6 )
    
    for ( i in 1:6 ) 
       {
        novel_cov[[i]] <- compute_coverage( novel_segs[[i]] )
        maxcov[i] <- max( novel_cov[[i]] )
        cat( "maxcov[",i,"] is ", maxcov[i], "\n" )
       }
    
    gene <- ""
    ysep <- round( max( maxcov ) * 1.05 )
    cat( "ysep is ", ysep, "\n" )
    ylevs = ((1:6) - 0.75) * ysep
    ymax <- ylevs[6] + ysep
    
    plot( 0, 0, xlim=pos_range, ylim=c(1,ymax), type="n",                # inital empty plot
          main = title,
          xlab = "position (nt)",
          ylab = "coverage"
         )
    #
    # Plot the coverage histograms, use pre-computed coverage arrays
    #
    for ( i in 1:6 )
       {
        plot_coverage( novel_segs[[i]], 
                       ylevs[i], 
                       col="red",
                       cov = novel_cov[[i]]
                     )
       }
    
    text( pos_range[1], ylevs, 
          # labels = c("580-1","580-2","580-3","581-1","581-2","581-3"),
          labels = c( "580-1", "581-2", "581-3", "581-1", "580-2", "580-3" ),  # regrouped
          adj=c(0,0),
     )

    blxdata <- blastxdata[   blastxdata$chrom == chromosome 
                           & blastxdata$start >= pos_range[1] 
                           & blastxdata$stop  <= pos_range[2], ]

    plot_blastx_hsps( blxdata )
    list_blastx_hsps( blxdata )

    print( maxcov )
   }   # end novel_view_by_pos

#
#  this makes the same plot given chromesome and start,stop
# 
novel_view_region <- function( chromosome, start, stop, title="" )
   {
    center <- (start + stop) / 2
    ext <- (stop - start) / 2
    novel_view_by_pos( chromosome, center, ext, title )
   }


plotone <- function( cl )
  {
   locus <- paste( cl$chrom[1], ":", cl$start[1], "-", cl$stop[1], sep="" )
   cat( paste( rep("=", 80), sep="", collapse=""), "\n" )
   cat( "new locus", locus, "  ", toString( cl$num[1] ), "reads", "\n" )
   cat( paste( rep("=", 80), sep="", collapse=""), "\n" )
   center <- (cl$start[1] + cl$stop[1] ) / 2
   novel_view_by_pos( cl$chrom, center, 7500, 
                      paste( "novel cluster ", locus,
                             "\n", toString( cl$num[1] ), " reads",
                             sep=""
                            ) 
                    )
  }

clusts <- read.table( "../../above_mean.1000.clusts", col.names=c("chrom","start","stop","num" ) )

for ( i in 1:100 )
   plotone( clusts[i,] )

