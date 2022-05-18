for SAMPLE in barcode{01..08} 
do 
	sbatch map_ONT.sh $SAMPLE\_pass
	sbatch map_ONT.sh $SAMPLE\_fail
done