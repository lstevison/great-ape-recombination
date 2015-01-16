The files in this repository are part of the GARMaps project (GAR = Great Ape Recombination Maps). 

1. Rate estimation Pipeline

 1A. Initial Filtering

  1Ai. Thin sites: Sites were filtered to thin SNPs within 15bp of each other. Although there is a function for this in VCFtools, it removes both SNPs if they are within 15bp of each other. We wanted to only remove one SNP, so we wrote a separate script to complete this filtering step. This script is located in a separate repository along with other vcf conversion scripts called 'vcf-conversion-tools'. It is called 'thinVCF.pl'. Usage information can be found in the corresponding README file in that repository.

  1Aii. All other initial filtering was done using VCFtools. Options used include: --mac 1 --max-missing-count --hwe 0.001  


 1B. Reciprocal liftover.

  1Bi. Thinned VCF was converted to extended BED file for input into liftOver keeping original coordinates for later comparison. Each site was reciprocally lifted-over using the UCSC tool liftOver (options: minMatch=0.1; -bedPlus=3). 

  1Bii. Between liftover steps, 'liftover_cleanup_new.pl' was run to move add to extended BED file the converted coordinates from first liftover step for later comparison and add a distance column for original interval and converted coordinates.

  1Biii. Filtering of the mapped sites following the reciprocal liftover process first cleaned output using 'reciprocal_liftover_cleanup_new.pl' script to add a distance column to the new coordinates. Then 'filter_rl_intervals_new.pl' script was run to retain only sites that mapped back to the within 100bp of original position and intervals within 20bp of the length of the original interval. 

 1C. Intersection between assemblies.

  1Ci. First filtered sites from hg18 mapping were recorded using VCFtools (option: --filtered-sites). Then, the corresponding hg18 coordinates were extracted from non-human primate VCF from the info field using VCFtools (option: --get-INFO 'hg18'). This output was then converted to a format to similar to the '.kept.sites' output from VCFtools using the script 'INFOtokept.sites.pl'. 

  1Cii. Intersection was done using awk:

	$ awk 'NR==FNR{a[$0]=$0;next}a[$0]' converted-info-field-output filtered-sites-output-fromhg18 >OUTPUT

  1Ciii. Converted coordinates back: Using script 'new_Intersection.pl', the hg18 coordinates pulled from the non-human primate VCF were used to convert the intersected hg18 coordinates back to the non-human primate coordinates to reduce each corresponding VCF file using VCFtools (options: --positions --filtered-sites).

 1D: Synteny

  1Di. Get coordinates in both reference genomes: Similar to intersection step, synteny requires coordinates in both assemblies. These are joined using 'convert-coord-gorgor2hg18.pl' and the same INFO field output from above and the filtered-sites output from step 1Ciii. 

  1Dii. Step 1: The output from 1Di is input into 'synteny.pl' to generate syntenic blocks between the two assemblies based on similar orientations and consistent interval lengths without gaps greater than 50kb. 

  1Diii. Step 2-3: The first set of syntenic blocks are shunted through the program 'synteny_step2.pl' twice to collapse consecutive blocks with gaps less than 300 SNPs and matching orientations. 

 1E. Computational Phasing

  1Ei. FastPHASE imputation: first an initial round of rough phasing was done to perform imputation of missing sites. This was done within synteny blocks to reduce the computational load. First, vcf files were converted to fastPHASE input using 'vcf2fastPHASE.pl' located in a separate repository along with other vcf conversion scripts called 'vcf-conversion-tools'. Usage information can be found in the corresponding README file in that repository. fastPHASE was then run with option -K10. The output was converted back to vcf using 'fastPHASE2VCF.pl' also found in the 'vcf-conversion-tools' repository. 

  1Eii. Minor allele frequency filtering: Due to a low level of contamination in the initial Great Ape dataset and the fact that minor alleles do not contribute much information to recombination rate estimation, we imposed a minor allele cutoff at 0.05 using VCFtools (option: --maf 0.05). Blocks with less than 300 SNPs remaining were dropped as they would not be suitable for rate estimation.

  1Eiii. More accurate phasing using PHASE: To allow for more accurate rate estimation especially at fine-scales, we added an additional computational phasing step using the program PHASE. First, the data were split into 400 SNP blocks with 50 SNP overlap using 'CreateChrBedFromVCF2.pl' found in the 'vcf-conversion-tools' repository. Then, for each 400 SNP block the VCF file was converted to PHASE input using 'vcf2PHASE.pl' found in the 'vcf-conversion-tools' repository. PHASE was run with the following options: 200 1 300 ‐MR ‐F.05 ‐l10 ‐x5 –X5. Then, phase blocks were joined by chromosome using 'join_phase_blocks.pl'. 

 1F. Rate Estimation

  1Fi. Convert PHASE output to VCF: Using initial minor allele filtered VCF file and PHASE output file, a VCF was generated using 'PHASE2VCF.pl' which can be found in the 'vcf-conversion-tools' repository.

  1Fii. Break chromosomes into 4000 SNP segments: Each block was further reduced to 4k SNP segments for rate estimation using 'CreateChrBedFromVCF.pl' found in the 'vcf-conversion-tools' repository.

  1Fiii. Rate estimation: Recombination rates for each SNP interval were estimated using LDhat 2.1 with options: -its 60000000 -bpen 5 -samp 40000. The LDhat program 'lkgen' was used to generate a corresponding lookup table for each dataset. Then the LDhat program 'stat' was used to remove the burn in and create a "res" file with option -burn 500 corresponding to a 20 million iteration burn in. Then the program 'first_column.pl' was used to generate a convergence summary file from the rates and bounds outputs. 

 1G. Post-filtering

  1Gi. Combine LDhat outputs: Using the program 'make_rho_map_filtered4.pl', the set of LDhat outputs for each chromosome were combined into a single chromosome map file removing overlapping sites and converting rate estimates for sites with high rho to maximum of 100 and nearby sites to 150. Also, R scripts for each chromosome to plot rate estimates and convergence summaries are generated. Then, corresponding hg18 coordinates were added to the map file using the program 'new_convert-map_gorgor2hg18.pl'. Next, synteny information including block and orientation were added to the map file using 'filter_supercontigs3.pl' as are previous site filtering information.

2. Analysis

 2A. Get published maps.

  2Ai. PanMap. Once we downloaded the file: 'lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10CEU_10YRI_with_chimp_super_clean_genotype_map_and_HapMap_Pops.txt' from the PanMap website, we used the program 'split_maps.pl' to split the file based on the corresponding columns for into each sub map and converted recombination rate to rho/kb . This program also worked to get published maps in a similar format to our files which include both start and stop coordinate for each recombination interval to account for the breaks in the map due to synteny blocks. 

  2Aii. Human population specific recombination rates were downloaded from: https://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/. As will the PanMap data, we wanted to convert both CEU and YRI genetic map files to have start and stop coordinates for each recombination interval and the rates in cM/Mb. We used awk to do this, assuming an Ne=10000 for CEU and Ne=XX for YRI (ex. awk 'NR>1 && OFS="\t" {if(v)print "chr1",int(v),int($1),int((v+$1)/2),rate,map,(rate*4)/10;v=$1;rate=$2;map=$3}' genetic_map_chr1_CEU_b36.txt). We combined all chromosomes into a full genome map file to use below in the multi-syntenic comparison. We also proceeded to generate genome files of 1kb regions with average and peak recombination rate for each map file as described later in 2Di for the other map files. To do this, we made a file with the first and last coordinates on each chromosome similar to a program described below. We then used 'make_bins.pl' to get 1kb bins across whole genome. Next, we used the unix command 'split' to the expanded BED file into separate files each containing only 5000 lines. We then used 'rate_at_hotspots.pl' to run a parallel jobs based on each split file to get average and peak recombination rate information from the genetic maps file for each 1kb region. We also used the program 'extractseq.pl' to simultaneously extract the FASTA sequence for each region using the unmasked version of each genome file and calculate the GC content for each region. We then wrote these individual rate and GC outputs to a main output file for each genome.

 2B. Multi-syntenic comparison.

  2Bi. Combined BED file. First, we used the program 'map_block_boundaries.pl' to get the boundaries of each syntenic block in the final map file using the final map file and a BED summary of the syntenic blocks across each genome. We then used liftover to convert the species-specific coordinates to hg18 for all genomes and sorted the output by chromosome. We then used the BEDtools program 'multiIntersectBed' to get the intersection between all non-human syntenic block coordinates (options: -header -names). We then removed all blocks below 1Mb. We then used the intersected set of coordinates to extract from each sorted hg18 bed file the orientation information for each syntenic block with the BEDtools program 'intersectBed' (options: -wba). Finally, we used the unix tool 'paste' to join the columns of each file into a combined output.

  2Bii. Reduce and bin maps. Next, we generated a binning for the combined BED file based on various size bins we wanted to do using the program 'thin_reduced_maps.pl' and specifying the size (example 1000000 for 1Mb). Then, we used the program 'reduce_maps.pl' to reduce each genetic map file based on the combined BED file, inputting the column containing the chromosome information, rate information and format on the command line. We then sorted each reduced map file and binned using the program 'bin_maps.pl' and the output of 'thin_reduced_maps.pl' to bin each map into specific size bins. Finally, we used the unix command paste to combine the rate columns from all reduced/binned maps into a mutli-synteny file for all the data for each bin size. 

 2C. Hotspot Analysis. 

  2Ci. To get a composite of recombination rates at +/- 20kb from all of the sorted hg18 hotspots, we used the program 'hotspot_calc_new.pl' using as input the map files for each population. The output is two columns, the first is the relative distance to the center of the hotspot ranging from -20 to +20. The second column contains the recombination rate from the map that corresponds to that distance. For each set of hotspots, the outputs for all six maps were plotted in R using a loess smoothing with enp.target set to 20 for each dataset to ensure the same scaling for all files. 

  2Cii. To create cumulative distribution of sequence and rate information, we used the program 'plot_regression.pl' and the corresponding map file for each species. This file also has two columns which start at zero for each and end at roughly 100 for each. The first column represents the cumulative physical distance fraction and the second column represents the commutative recombination rate fraction. It is the second column which is used in R to plot the Gini coefficient and Lorenz curve using the 'ineq' package. 

 2D. PRDM9 Analysis.

  2Di. First, we generated PWM for each sub motif using zf.princeton.edu (option: polynomial SVM) and the species-specific coordinates for hotspots. Also, we generated a genome-wide file with recombination rate and GC content for every 1kb region in the genome. To do this, we first got the coordinates for the syntenic blocks for the whole genomes of each map to ensure the regions we generated would have recombination information. We then used 'make_bins.pl' to split this BED file so that it would have an individual entry for each 1kb region within each set of coordinates. Next, we used the unix command 'split' to the expanded BED file into separate files each containing only 5000 lines. We then used 'rate_at_hotspots.pl' to run a parallel jobs based on each split file to get average and peak recombination rate information from the genetic maps file for each 1kb region. We also used the program 'extractseq.pl' to simultaneously extract the FASTA sequence for each region using the unmasked version of each genome file and calculate the GC content for each region. We then wrote these individual rate and GC outputs to a main output file for each genome. 

  2Dii. Next, we did a count of potential binding sites for each PRDM9 sub-motif in both hotspot regions and cold spot regions. To get matched cold spots, we sorted the set of hotspots using 'sort_by_chr.sh' from the repository 'random-scripts' using species-specific chromosome lists for each map. Then, we used the program 'rate_at_hotspots.pl' to calculate average and peak rate for each hotspot region based on the genetic map files. Next, we used 'extractseq.pl' to extract FASTA sequence from the unmasked genome file for all hotspots. We then used the program 'match_cold_spots_final.pl' to read in the FASTA of hotspots, the hotspot file with rates, the genome file and output matched regions based on an average and peak rate cutoff and a distance metric. Once we had a set of matched hotspot and coldspot regions, we generated corresponding bed files for the coordinates for each and used 'extractseq.pl' to extract the FASTA sequence for each region from the masked genome files. We then ran the program 'fimo' on both the hotspot and coldspot sequences based on the relevant PRDM9 sub-motifs (option: --parse-genomic-coord). We then combined the fimo outputs for hotspot and coldspot regions back into the original matched set file to include a count of significant fimo hits using the program 'panmap_motif_cleanup.pl'.

  2Diii. Third, we performed a genome-wide search for PRDM9 binding, regardless of hotspot locations. This also relied on the genome file of 1kb regions. For the purposes of running fimo on each region, we split the genome file into smaller files using the unix command 'split' that were each 10k lines. We then ran fimo on each subset of 10k lines in parallel. First, we extracted the FASTA sequence of each 1kb region using the masked version of the genome files with the program 'extractseq.pl'. We then ran the program 'fimo' on each 1kb region based on the relevant PRDM9 sub-motifs (option: --parse-genomic-coord). Finally, we combined the fimo outputs for each individual fimo run using the unix command 'cat' . We then summarized the individual fimo outputs for each sub-motif using the original genome file and the program 'chr_motif_cleanup.pl'. 

3. Links to data in Dryad

 3A. Recombination map files for bonobo, Nigerian chimpanzee and gorilla. Each genetic map is in one of the main folders and within the folder is a genetic map for each chromosome. See main readme.txt file for information on headers within each subfolder. Also, the species-specific set of hotspots and the hg18 liftOver of these hotspots is in each sub-directory.

 3B. Recombination hotspots and matched coldspots for bonobo, western chimpanzee, Nigerian chimpanzee, and gorilla together with number of potential binding sites of each PRDM9 submotif searched. Both average and peak recombination rate as well as %GC and %N are provided for both hotspots and coldspots. Hotspot regions that had multiple matches include the number of matches in a 100kb window.



