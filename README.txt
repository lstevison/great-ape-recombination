The files in this repository are part of the GARMaps project (GAR = Great Ape Recombination Maps). 

1. Rate estimation Pipeline

 1A. Initial Filtering

  1Ai. Thin sites: Sites were filtered to thin SNPs within 15bp of each other. Although there is a function for this in VCFtools, it removes both SNPs if they are within 15bp of each other. We wanted to only remove one SNP, so we wrote a separate script to complete this filtering step. This script is located in a separate repository along with other vcf conversion scripts called 'vcf-conversion-tools'. It is called 'thinVCF.pl'. Usage information can be found in the corresponding README file in that repository.

  1Aii. All other initial filtering was done using VCFtools. Options used include: --mac 1 --max-missing-count --hwe 0.001  


 1B. Reciprocal leftover.

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

 2A. Get published maps:

 2B. Multi-syntenic comparison:

 2C. Hotspot analysis: Running MLEHOT stuff can go here.

 2X. PRDM9 Analysis: Stuff running PRDM9 and such can go here.


3. Links to data in Dryad
