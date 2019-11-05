#!/bin/bash

echo "Running summary_utils.sh script"
# set defaults
BINS='final_bins'
BINSDIR='concoct_output_dir'

### ---------------------- ###
###        Functions       ###
testRun()
{
  command=$1
  errorCode=$2
  if [[ $errorCode != 0 ]]; then
    echo "Command failed: $command"
    echo "$errorCode"
    exit 1
  else
    echo "completed $command successfully"
  fi
}
### ---------------------- ###
set -x


echo "Test to see if working"
# The "binned_contigs_obj_ref" requires a file called *.summary so
# I'm creating a fake one.  The binned contigs object needs to be updated
# so this is not necessary.
cd $BINSDIR

echo -e "Bin name\tCompleteness\tGenome size\tGC content" > $BINS/notreal.summary
#testRun "create $BINS/notreal.summary" $?

# Unfortunately, the MetagenomeUtils function that loads the bins into "binned_contig_obj_ref"
# uses a regular expression that requires the bin names to be in the form <something>.123.fasta
# So I'm converting bin names here.


# if we have bins, change their names so that MetagenomeUtils will
# recognize them.
# Test that metabat produced some bins
BINSARRAY=($(find $BINS -name "*.fa"))
testRun "find $BINS -name '*.fa'"
COUNT=0
for i in ${BINSARRAY[@]}; do
    if [[ $i =~ "unbinned" ]]; then
      newID="unbinned"
    else
      COUNT=$((COUNT+1))
      binID=$(basename $i)
      binID=${binID%.fa}
      binID=${binID#bin.}
      newID=$(printf "%03d\n" $binID)
    fi

    mv $i $BINS/bin.$newID.fasta

	# lets populate the .summary file (this is a work around until MetagenomeUtils gets fixed)
	# The collumns are binID, read coverage, Total Contig Length, GC%
    echo -e "bin.$newID.fasta\t0\t0\t0" >> $BINS/notreal.summary
done

# cleanup large files
#rm -r $BINSDIR
