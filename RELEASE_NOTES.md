# kb_concoct release notes
=========================================

ToDo
-----
* Add kb_parallels support - see kbase-metabat example
* Add different ways to extract coverage (e.g. jgi_summarize_bam_contig_depths vs coverM, vs ?)
* Expose seed parameters (e.g. HISAT2)

1.3.1
-----
* Added HISAT2 read aligner as mapping option
* Large code de-linting effort
* Added tests for different read mapping tools
* Added bowtie2 very sensitive mode

1.3.0
-----
* Refactored most code for kb_concoctTest, streamlined variable and function usage
* Added new method (stats.sh) to populate summary files, exposed many other statistics for genome bin summary files
* Removed bash script dependency
* Added support for multiple read file inputs
* Added auto-detection of read file type for read mapping - single vs paired-end interleaved file

1.2.0
-----
* Added multiple read mapping software (bowtie2, minimap2, and bwa), in addition to existing bbmap as selectable options

1.1.0
-----
* Added parameter options for: no_cov_normalization, no_total_coverage, and total_percentage_pca

1.0.0
-----
* Initial release

0.0.0
-----
* Module created by kb-sdk init
