#Great Ape Recombination (GAR) Maps 

The files in this repository are part of the GARMaps project (GAR = Great Ape Recombination Maps). To cite this repository, use: http://dx.doi.org/10.5281/zenodo.13975

###1. Rate estimation Pipeline   
  
####1A. Initial Filtering   

  i. Thin sites: Sites were filtered to thin SNPs within 15bp of each other. Although there is a function for this in VCFtools, it removes both SNPs if they are within 15bp of each other. We wanted to only remove one SNP, so we wrote a separate script to complete this filtering step. This script is located in a separate repository along with other vcf conversion scripts called [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/). It is called _thinVCF.pl_. Usage information can be found in the corresponding README file in that repository.   

 ii. All other initial filtering was done using VCFtools. Options used include: 

>--mac 1 --max-missing-count --hwe 0.001  
  
####1B. Reciprocal liftover.   

  i. Thinned VCF was converted to extended BED file for input into liftOver keeping original coordinates for later comparison. Each site was reciprocally lifted-over using the UCSC tool liftOver (options: minMatch=0.1; -bedPlus=3).    

  ii. Between liftover steps, _liftover_cleanup_new.pl_ was run to move add to extended BED file the converted coordinates from first liftover step for later comparison and add a distance column for original interval and converted coordinates.   

  iii. Filtering of the mapped sites following the reciprocal liftover process first cleaned output using _reciprocal_liftover_cleanup_new.pl_ script to add a distance column to the new coordinates. Then _filter_rl_intervals_new.pl_ script was run to retain only sites that mapped back to the within 100bp of original position and intervals within 20bp of the length of the original interval.   

####1C. Intersection between assemblies.
  
  i. First filtered sites from hg18 mapping were recorded using VCFtools (option: --filtered-sites). Then, the corresponding hg18 coordinates were extracted from non-human primate VCF from the info field using VCFtools (option: --get-INFO 'hg18'). This output was then converted to a format to similar to the '.kept.sites' output from VCFtools using the script _INFOtokept.sites.pl_.
    
  ii. Intersection was done using awk:   
    
```Shell
	$ awk 'NR==FNR{a[$0]=$0;next}a[$0]' converted-info-field-output filtered-sites-output-fromhg18 >OUTPUT
```  

  iii. Converted coordinates back: Using script _new_Intersection.pl_, the hg18 coordinates pulled from the non-human primate VCF were used to convert the intersected hg18 coordinates back to the non-human primate coordinates to reduce each corresponding VCF file using VCFtools.  
   
>Options: --positions --filtered-sites
  
####1D: Synteny
  
  i. Get coordinates in both reference genomes: Similar to intersection step, synteny requires coordinates in both assemblies. These are joined using _convert-coord-gorgor2hg18.pl_ and the same INFO field output from above and the filtered-sites output from step 1Ciii.    
    
  ii. Step 1: The output from 1Di is input into _synteny.pl_ to generate syntenic blocks between the two assemblies based on similar orientations and consistent interval lengths without gaps greater than 50kb.
    
    
  iii. Step 2-3: The first set of syntenic blocks are shunted through the program _synteny_step2.pl_ twice to collapse consecutive blocks with gaps less than 300 SNPs and matching orientations.    
  
####1E. Computational Phasing
    
  i. FastPHASE imputation: First an initial round of rough phasing was done to perform imputation of missing sites. This was done within synteny blocks to reduce the computational load. First, vcf files were converted to fastPHASE input using _vcf2fastPHASE.pl_ located in a separate repository along with other vcf conversion scripts called [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/). Usage information can be found in the corresponding README file in that repository. The program [fastPHASE](https://els.comotion.uw.edu/express_license_technologies/fastphase) was then run with option -K10. The output was converted back to vcf using _fastPHASE2VCF.pl_ also found in the [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/) repository.    
    
  ii. Minor allele frequency filtering: Due to a low level of contamination in the initial Great Ape dataset and the fact that minor alleles do not contribute much information to recombination rate estimation, we imposed a minor allele cutoff at 0.05 using [VCFtools](https://vcftools.github.io/index.html) (option: --maf 0.05). Blocks with less than 300 SNPs remaining were dropped as they would not be suitable for rate estimation.   
    
  iii. More accurate phasing using PHASE: To allow for more accurate rate estimation especially at fine-scales, we added an additional computational phasing step using the program [PHASE](http://stephenslab.uchicago.edu/phase/download.html). First, the data were split into 400 SNP blocks with 50 SNP overlap using 'CreateChrBedFromVCF2.pl' found in the 'vcf-conversion-tools' repository. Then, for each 400 SNP block the VCF file was converted to PHASE input using _vcf2PHASE.pl_ found in the [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/) repository. PHASE was run with the following options: 200 1 300 ‐MR ‐F.05 ‐l10 ‐x5 –X5. Then, phase blocks were joined by chromosome using _join_phase_blocks.pl_.   

####1F. Rate Estimation
    
  i. Convert PHASE output to VCF: Using initial minor allele filtered VCF file and PHASE output file, a VCF was generated using _PHASE2VCF.pl_ which can be found in the [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/) repository.   
    
  ii. Break chromosomes into 4000 SNP segments: Each block was further reduced to 4k SNP segments for rate estimation using _CreateChrBedFromVCF.p'_ found in the [vcf-conversion-tools](https://github.com/lstevison/vcf-conversion-tools/) repository.   
    
  iii. Rate estimation: Recombination rates for each SNP interval were estimated using LDhat 2.1 with options: -its 60000000 -bpen 5 -samp 40000. The [LDhat program](https://sourceforge.net/projects/ldhat/files/?source=navbar) **lkgen** was used to generate a corresponding lookup table for each dataset. Then the LDhat program **stat** was used to remove the burn in and create a "res" file with option -burn 500 corresponding to a 20 million iteration burn in. Then the program _first_column.pl_ was used to generate a convergence summary file from the rates and bounds outputs.    

####1G. Post-filtering
    
  i. Combine LDhat outputs: Using the program _make_rho_map_filtered4.p'_, the set of LDhat outputs for each chromosome were combined into a single chromosome map file removing overlapping sites and converting rate estimates for sites with high rho to maximum of 100 and nearby sites to 150. Also, R scripts for each chromosome to plot rate estimates and convergence summaries are generated. Then, corresponding hg18 coordinates were added to the map file using the program _new_convert-map_gorgor2hg18.pl_. Next, synteny information including block and orientation were added to the map file using _'filter_supercontigs3.pl'_ as are previous site filtering information.

###2. Analysis 
  
####2A. Get published maps.
    
  i. PanMap. Once we downloaded the file: _lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10CEU_10YRI_with_chimp_super_clean_genotype_map_and_HapMap_Pops.txt_ from the [PanMap website](http://panmap.uchicago.edu/), we used the program _split_maps.pl_ to split the file based on the corresponding columns for into each sub map and converted recombination rate to rho/kb . This program also worked to get published maps in a similar format to our files which include both start and stop coordinate for each recombination interval to account for the breaks in the map due to synteny blocks.    
    
  ii. [Human population specific recombination rates](https://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/) were downloaded. As will the PanMap data, we wanted to convert both CEU and YRI genetic map files to have start and stop coordinates for each recombination interval and the rates in cM/Mb. We used awk to do this, assuming an Ne=10000 for CEU and Ne=XX for YRI.    
  Example: 

```Shell
awk 'NR>1 && OFS="\t" {if(v)print "chr1",int(v),int($1),int((v+$1)/2),rate,map,(rate*4)/10;v=$1;rate=$2;map=$3}' genetic_map_chr1_CEU_b36.txt
```   

  iii. We combined all chromosomes into a full genome map file to use below in the multi-syntenic comparison. We also proceeded to generate genome files of 1kb regions with average and peak recombination rate for each map file as described later in 2Di for the other map files. To do this, we made a file with the first and last coordinates on each chromosome similar to a program described below. We then used _make_bins.pl_ to get 1kb bins across whole genome. Next, we used the unix command _split_ to the expanded BED file into separate files each containing only 5000 lines. We then used _rate_at_hotspots.pl_ to run a parallel jobs based on each split file to get average and peak recombination rate information from the genetic maps file for each 1kb region. We also used the program _extractseq.pl_ to simultaneously extract the FASTA sequence for each region using the unmasked version of each genome file and calculate the GC content for each region. We then wrote these individual rate and GC outputs to a main output file for each genome.

####2B. Multi-syntenic comparison.  
  
  i. Combined BED file. First, we used the program _map_block_boundaries.pl_ to get the boundaries of each syntenic block in the final map file using the final map file and a BED summary of the syntenic blocks across each genome. We then used liftover to convert the species-specific coordinates to hg18 for all genomes and sorted the output by chromosome. We then used the BEDtools program _multiIntersectBed_ to get the intersection between all non-human syntenic block coordinates (options: -header -names). We then removed all blocks below 1Mb. We then used the intersected set of coordinates to extract from each sorted hg18 bed file the orientation information for each syntenic block with the BEDtools program _intersectBed_ (options: -wba). Finally, we used the unix tool _paste_ to join the columns of each file into a combined output.   
  
  ii. Reduce and bin maps. Next, we generated a binning for the combined BED file based on various size bins we wanted to do using the program _thin_reduced_maps.pl_ and specifying the size (example 1000000 for 1Mb). Then, we used the program _reduce_maps.pl_ to reduce each genetic map file based on the combined BED file, inputting the column containing the chromosome information, rate information and format on the command line. We then sorted each reduced map file and binned using the program _bin_maps.pl_ and the output of _thin_reduced_maps.pl_ to bin each map into specific size bins. Finally, we used the unix command paste to combine the rate columns from all reduced/binned maps into a mutli-synteny file for all the data for each bin size.    

####2C. Hotspot Analysis. 
    
  i. To get a composite of recombination rates at +/- 20kb from all of the sorted hg18 hotspots, we used the program _hotspot_calc_new.pl_ using as input the map files for each population. The output is two columns, the first is the relative distance to the center of the hotspot ranging from -20 to +20. The second column contains the recombination rate from the map that corresponds to that distance. For each set of hotspots, the outputs for all six maps were plotted in R using a loess smoothing with enp.target set to 20 for each dataset to ensure the same scaling for all files. 

  ii. To create cumulative distribution of sequence and rate information, we used the program _plot_regression.pl_ and the corresponding map file for each species. This file also has two columns which start at zero for each and end at roughly 100 for each. The first column represents the cumulative physical distance fraction and the second column represents the commutative recombination rate fraction. It is the second column which is used in R to plot the Gini coefficient and Lorenz curve using the _ineq_ package. 

####2D. PRDM9 Analysis.

  i. First, we generated PWM for each sub motif using zf.princeton.edu (option: polynomial SVM) and the species-specific coordinates for hotspots. Also, we generated a genome-wide file with recombination rate and GC content for every 1kb region in the genome. To do this, we first got the coordinates for the syntenic blocks for the whole genomes of each map to ensure the regions we generated would have recombination information. We then used _make_bins.pl_ to split this BED file so that it would have an individual entry for each 1kb region within each set of coordinates. Next, we used the unix command 'split' to the expanded BED file into separate files each containing only 5000 lines. We then used _rate_at_hotspots.pl_ to run a parallel jobs based on each split file to get average and peak recombination rate information from the genetic maps file for each 1kb region. We also used the program _extractseq.pl_ to simultaneously extract the FASTA sequence for each region using the unmasked version of each genome file and calculate the GC content for each region. We then wrote these individual rate and GC outputs to a main output file for each genome. 

  ii. Next, we did a count of potential binding sites for each PRDM9 sub-motif in both hotspot regions and cold spot regions. To get matched cold spots, we sorted the set of hotspots using _sort_by_chr.sh_ from the repository [random-scripts](https://github.com/lstevison/random-scripts) using species-specific chromosome lists for each map. Then, we used the program _rate_at_hotspots.pl_ to calculate average and peak rate for each hotspot region based on the genetic map files. Next, we used _extractseq.pl_ to extract FASTA sequence from the unmasked genome file for all hotspots. We then used the program _match_cold_spots_final.pl_ to read in the FASTA of hotspots, the hotspot file with rates, the genome file and output matched regions based on an average and peak rate cutoff and a distance metric. Once we had a set of matched hotspot and coldspot regions, we generated corresponding bed files for the coordinates for each and used _extractseq.pl_ to extract the FASTA sequence for each region from the masked genome files. We then ran the program [fimo, from the meme suite](http://meme-suite.org/doc/fimo.html?man_type=web) on both the hotspot and coldspot sequences based on the relevant PRDM9 sub-motifs (option: --parse-genomic-coord). We then combined the fimo outputs for hotspot and coldspot regions back into the original matched set file to include a count of significant fimo hits using the program _panmap_motif_cleanup.pl_.

  iii. Third, we performed a genome-wide search for PRDM9 binding, regardless of hotspot locations. This also relied on the genome file of 1kb regions. For the purposes of running fimo on each region, we split the genome file into smaller files using the unix command _split_ that were each 10k lines. We then ran fimo on each subset of 10k lines in parallel. First, we extracted the FASTA sequence of each 1kb region using the masked version of the genome files with the program _extractseq.pl_. We then ran the program **fimo** on each 1kb region based on the relevant PRDM9 sub-motifs (option: --parse-genomic-coord). Finally, we combined the fimo outputs for each individual fimo run using the unix command _cat_ . We then summarized the individual fimo outputs for each sub-motif using the original genome file and the program _chr_motif_cleanup.pl_. 

###3. Folders with final map data

####3A. Recombination map files for bonobo, Nigerian chimpanzee and gorilla. 

Each genetic map is in one of the main folders and within the folder is a genetic map for each chromosome. See main readme.txt file for information on headers within each subfolder. _Coming soon, we will post the species-specific set of hotspots and the hg18 liftOver of these hotspots in each sub-directory._

####3B. Raw recombination hotspots, and hotspots with matched coldspot files. 

In main directory, hotspots identified were lifted over to hg18 using liftover and sorted. In sub-directory, recombination hotspots and matched coldspot for bonobo, western chimpanzee, Nigerian chimpanzee, and gorilla together with number of potential binding sites of each PRDM9 submotif searched. Both average and peak recombination rate as well as %GC and %N are provided for both hotspots and coldspots.
