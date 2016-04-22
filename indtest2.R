source( "genomeviewer.R" )

loaded <- TRUE
#loaded <- FALSE

if ( ! loaded )
{

chromosome <- "chrC03"
target_position <- 7234532

#
# gene annotation first
#
genefile <- "../../Brassica_napus.annotation_v5.gff3"

genedata <- read.table( genefile )
names( genedata ) <- c( "chrom","caller","type","start","stop","score", 
                            "d1", "d2", "descrip" )

#
# six unguided directories
#
unguided_dirs <- 
   paste( "/Projects/Plants/Brassica_napus/Tophat_hpc1/tophat_out_allchrom_", 
          c( "580-1", "580-2", "580-3", "581-1", "581-2", "581-3" ), 
         sep="" )

#
# six guided directories
#
guided_dirs <- 
   paste( "/Projects/Plants/Brassica_napus/Tophat_hpc1/tophat_out_guided_allchrom_",
          as.character( 1:6 ),
          sep="" )

# get unguided data

unguided_setdata <- as.list( rep( "", 6 ) )

for ( i in 1:6 )
   {
    samindex <- load_sam_index( paste( unguided_dirs[i], "sorted_accepted_hits.index", sep="/" ) )
    con <- open_samfile_read( paste( unguided_dirs[i], "sorted_accepted_hits.sam", sep="/" ) )
    lines <- extract_sam_record_block( samindex, con, chromosome, target_position )
    cat( lines[1], "\n" )
    cat( lines[length(lines)], "\n" )
    unguided_setdata[[i]] <- parse_sam_lines( lines )
    close( con )
    print( head( unguided_setdata[[i]] ) )
   }

# get guided data

guided_setdata <- as.list( rep( "", 6 ) )

for ( i in 1:6 )
   {
    samindex <- load_sam_index( paste( guided_dirs[i], "sorted_accepted_hits.index", sep="/" ) )
    con <- open_samfile_read( paste( guided_dirs[i], "sorted_accepted_hits.sam", sep="/" ) )
    lines <- extract_sam_record_block( samindex, con, chromosome, target_position )
    cat( lines[1], "\n" )
    cat( lines[length(lines)], "\n" )
    guided_setdata[[i]] <- parse_sam_lines( lines )
    close( con )
    print( head( guided_setdata[[i]] ) )
   }

guided_segs <- as.list( rep( "", 6 ) )

for ( i in 1:6 )
   {
    guided_segs[[i]] <- parse_segments( guided_setdata[[i]] )
    cat( "guided_segs[[", i,"]]\n" )
    print(head( guided_segs[[i]] ) )
   }

}   # end load section


#####################################################################################



# Initialize with an empty plot

r <- getrange( setdata$pos )
#r <- c(7200000,7500000)
#r <- c(7360000,7380000)

# First, pre-compute coverage arrays for all tracks, so we can ascertain
# maxima and minima

unguided_cov <- as.list( rep( "", 6 ) )
guided_cov <- as.list( rep( "", 6 ) )

maxcov <- rep( 0, 6 )

for ( i in 1:6 ) 
   {
    unguided_cov[[i]] <- compute_coverage( unguided_segs[[i]] )
    guided_cov[[i]]   <- compute_coverage( guided_segs[[i]]   )
    maxcov[i] <- max( unguided_cov[[i]], guided_cov[[i]] )
   }

gene <- ""
ysep <- round( max( maxcov ) * 1.05 )
cat( "ysep is ", ysep, "\n" )
ylevs = ((1:6) - 0.75) * ysep
ymax <- ylevs[6] + ysep

plot( 0, 0, xlim=r, ylim=c(1,ymax), type="n",                # inital empty plot
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

text( r[1], ylevs, labels=c("580-1","580-2","580-3","581-1","581-2","581-3"),adj=c(1,0) )

plot_gene( genedata[ genedata$chrom == chromosome & genedata$start >= r[1] & genedata$stop <= r[2], ] )

print( maxcov )




