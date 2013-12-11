--------------------------------------------------------------------------------
Rover - ROVER-PCR Variant Caller: read-pair overlap considerate variant-calling
software for PCR-based massively parallel sequencing datasets
--------------------------------------------------------------------------------

Version: 1.1.0

Authors: Bernard J Pope (2,3), Tú Nguyen-Dumont (1), Fleur Hammet (1) and
         Daniel J Park (1)

         (1) Genetic Epidemiology Laboratory, Department of Pathology,
             The University of Melbourne.
         (2) Victorian Life Sciences Computation Initiative (VLSCI).
         (3) Department of Computing and Information Systems,
             The University of Melbourne.
         

Web:     https://github.com/bjpop/rover

License: BSD

Requirements: Python 2.7, and the PySam library
(http://code.google.com/p/pysam/).

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------

ROVER-PCR Variant Caller enables users to quickly and accurately identify
genetic variants from PCR-targeted, overlapping paired-end MPS datasets. The
open-source availability of the software and threshold tailorability enables
broad access for a range of PCR-MPS users.

ROVER-PCR Variant Caller is implemented in Python and runs on any Posix 
compliant system (Linux, OS X).

The software accepts a tab-delimited file (TSV) listing the coordinates of the
target-specific primers used for targeted enrichment based on a specified
genome-build. It also accepts aligned sequence files resulting from mapping
 to the same genome-build. ROVER-PCR Variant Caller identifies the amplicon
a given read-pair represents and removes the primer sequences by using the
 mapping co-ordinates and primer co-ordinates. It considers completely
 overlapping read-pairs with respect to primer-intervening sequence. Only 
when a variant is observed in both reads of a read-pair does the signal 
contribute to a tally of read-pairs containing or not containing the variant.
A user-defined threshold informs the minimum number of and proportion of 
read-pairs a variant must be observed in for a ‘call’ to be made. ROVER-PCR 
Variant Caller also reports the depth of coverage across amplicons to 
facilitate the identification of any regions that may require further
screening.

--------------------------------------------------------------------------------
Command line usage:
--------------------------------------------------------------------------------

usage: rover [-h] [--version] --primers PRIMERS [--overlap OVERLAP]
             [--log FILE] --out FILE [--proportionthresh N] [--absthresh N]
             [--qualthresh N] [--coverdir COVERDIR]
             bams [bams ...]

Consider mapped reads to amplicon sites

positional arguments:
  bams                  bam files containing mapped reads

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --primers PRIMERS     File name of primer coordinates in TSV format.
  --overlap OVERLAP     Minimum fraction overlap of read to block region.
                        Defaults to 0.9.
  --log FILE            Log progress in FILENAME, defaults to stdout.
  --out FILE            Name of output file containing called variants.
  --proportionthresh N  Keep variants which appear in this proportion of the
                        read pairs for a given target region, and bin
                        otherwise. Defaults to 0.05.
  --absthresh N         Only keep variants which appear in at least this many
                        read pairs. Defaults to 2.
  --qualthresh N        Minimum base quality score (phred).
  --coverdir COVERDIR   Directory to write coverage files, defaults to current
                        working directory.

Explanation of the arguments:

   -h

      Print a help message and exit.

   --version

      Print the version number of rover and exit.

   --primers PRIMERS

      Required.

      A tab-delimited list of the coordinates of the target-specific primers
      used for targeted enrichment based on a specified genome-build.
      The format is compatible with the BED file format:

      chromosome start end 

   --overlap OVERLAP

      Optional. Defaults to 0.9.

      Minimum fraction overlap of read to block region.

      0.5 means at least half of a read should overlap a given target region
      for the read to be considered for that region. 1.0 would mean the entire
      read must overlap the region.

   --log FILE

      Optional. Defaults to the standard output.

      Write a log of the program's progress to this file. 

   --out FILE

      Required.

      Name of the output file created by Rover. This file contains the
      variants called by the program.

   --proportionthresh N

      Optional. Defaults to 0.05.

      Only keep variants which appear in this proportion of the read pairs for
      a given target region, and bin otherwise. A variant must appear in
      both reads of a pair to be counted. The proportion is calculated as
      follows:

          N = number of pairs containing this variant in both reads
          T = number of read pairs overlapping the target region

          proportion = N/T


      Note: variants must pass BOTH the proportionthresh and absthresh thresholds
      to be kept. If they fail either test then they are binned.

      That is to say:

          if N/T  >= proportionthresh and N >= absthresh:
              keep the variant
          else:
              bin the variant

   --absthresh N

      Optional. Defaults to 2 read pairs.

      Only keep variants which appear in at least this many read pairs
      for a given target region.

      See comments above about proportionthresh.

   --qualthresh N

      Optional. Minimum phred quality score for bases appearing in SNVs and
      insertions. If this argument is set Rover will only consider SNVs and
      insertions where all DNA bases in those variants have a quality score
      greater than or equal to the argument. For example, setting this
      argument to 35 will cause Rover to discard any SNVs or insertions
      containing any bases with a score below 35. If the argument is not
      set then Rover will not consider quality scores in its decision
      to keep or discard a variant.
 
   --coverdir COVERDIR

      Optional. Defaults to current working directory.

      Directory to write the coverage files. If the directory does not
      exist Rover will not try to create it. You must create the directory
      yourself.

      The format of the coverage files is a TSV with:

      chr     block_start     block_end       num_pairs

   bams [bams ...] 

      One or more BAM files containing mapped reads.

--------------------------------------------------------------------------------
Example usage (should be all on one line)
--------------------------------------------------------------------------------

   rover --primers primer_coords.tsv --log rover_log --out variants
         --proportionthresh 0.15 --absthresh 2 --coverdir coverage_files
         sample1.bam sample2.bam sample3.bam

This assumes that the coordinates for your primer regions are in the file
primer_coords.tsv. The detected variants will appear in the output file
variants. The names of the samples will be taken from the prefix of
the bam file name, in this case "sample1" "sample2" and "sample3".
Coverage files containing the number of read pairs which mapped to
each region will be output in coverage_files/sample1.coverage
coverage_files/sample2.coverage and coverage_files/sample3.coverage.
A log file describing the actions taken by the program
will be stored in rover_log.

--------------------------------------------------------------------------------
