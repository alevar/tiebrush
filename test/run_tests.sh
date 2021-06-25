#!/bin/env bash
/bin/rm -f tst_* t?/tst_*

diff_check () {
 if [[ -z "$2" ]]; then 
  echo 'diff_check requires 2 files as parameters!'
  exit 1
 fi
 fa=$1
 fb=$2
 ext="${fb##*.}"
 if [[ $ext == "bam" ]]; then
   if diff -q <(samtools view $fa) <(samtools view $fb) &>/dev/null; then
      echo "OK ($fb)"
   else
      echo "Error: test failed ($fa != $fb)"
      exit 1
   fi
 else
   if diff -q $fa $fb &>/dev/null; then
      echo "OK ($fb)"
   else
      echo "Error: test failed ($fa != $fb)"
      exit 1
   fi
 fi
}

../tiebrush -o t1/tst_t1.bam t1/t1s[0-9].bam
diff_check t1/tst_t1.bam t1/t1.bam

../tiebrush -o t2/tst_t2.bam t2/t2s[0-9].bam
diff_check t2/tst_t2.bam t2/t2.bam

../tiebrush -o tst_t12.bam t1/tst_t1.bam t2/tst_t2.bam
diff_check tst_t12.bam t12.bam

../tiecov -s t1/tst_t1.sample -c t1/tst_t1.coverage -j t1/tst_t1.junctions t1/t1.bam
diff_check t1/tst_t1.sample.bedgraph t1/t1.sample.bedgraph
diff_check t1/tst_t1.coverage.bedgraph t1/t1.coverage.bedgraph
diff_check t1/tst_t1.junctions.bed t1/t1.junctions.bed

../tiecov -s t2/tst_t2.sample -c t2/tst_t2.coverage -j t2/tst_t2.junctions t2/t2.bam
diff_check t2/tst_t2.sample.bedgraph t2/t2.sample.bedgraph
diff_check t2/tst_t2.coverage.bedgraph t2/t2.coverage.bedgraph
diff_check t2/tst_t2.junctions.bed t2/t2.junctions.bed
