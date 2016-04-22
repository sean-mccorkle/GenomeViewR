
* first pass at multiplot function

     . parameters chromosome and position

     . ability to add or substract tracks

* check if desired range is close to 1,000,000 segment boundary in either 
  direction.  If so, include the previous or subsequent range as appropriate

     . then trim!!!
x abstract out segment compilation so that both plots and coverage histograms use it

* filter on matched pairs

x tracks

* program to read.table() genes gff3 file once, and then save as a binary
  so the binary can be hopfully quickly loaded.

* Coverage histograms

    x coverage
    
       x try R based first
        
    x hist plot

    * abstract out coverage array  calculation so maximum can be obtained for
      scale setting.

    * associated vertical scale

x Include all genes in window

x Faster retrieval from SAM files

    x indexing SAM files?

    * when reading SAM blocks, if pos is too close to edge, read two blocks.

    * BAM file support?

