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

  1Ciii. Converted coordinates back: Using script 'new_Intersection.pl', the hg18 coordinates pulled from the non-human primate VCF were used to convert the intersected hg18 coordinates back to the non-human primate coordinates to reduce each corresponding VCF file using VCFtools (option: --positions).

 1D: Synteny

 1E. Computational Phasing

 1F. Rate Estimation

 1G. Post-filtering

2. Analysis

 2A. PRDM9 Analysis


3. Links to data in Dryad
