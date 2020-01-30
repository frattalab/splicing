
#Submit to the cluster, give it a unique name
#$ -S /bin/bash
# 2. prompt for name of file or directory
echo -n "file or directory name: "
# ...  and read it
read HANDLE
# mv {,new.}original.filename
# # 2. b - check if it exists and is readable
# if [ ! -r "$HANDLE" ]
# then
#     echo "$HANDLE is not readable";
#     # if not, exit with an exit code != 0
#     exit 2;
# fi

for f in $HANDLE*.Aligned.sorted.out.bam; do
  echo "$f"
  echo "$f.insert_size_metrics.txt"
  java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar \
  CollectInsertSizeMetrics I=$f O=$f.insert_size_metrics.txt \
  H=$f.insert_size_histogram.pdf
done
