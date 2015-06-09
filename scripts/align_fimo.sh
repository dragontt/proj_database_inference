#! /bin/bash

# bash align_fimo.sh resources/scertf_names.txt resources/scertf_pwms output/scertf Saccharomyces_cerevisiae

# Input variables
IN_TF_LIST=$1 	# list of tf names
IN_TF_PWM=$2	# directory of tf pwm
OUT_FIMO=$3		# directory of fimo alignment output 
SPECIES=$4	# species, used for background freq file
# OUT_RANK=$4		# direcotry of the ranked lists of fimo output

# Algin scertf pfms individually 
PROJ_DIR=$HOME/proj_database_inference

mkdir -p $PROJ_DIR/$OUT_FIMO
mkdir -p $PROJ_DIR/$OUT_RANK

counter=0

while read -a line
do
	counter=$[$counter +1]
	motif=${line[0]}
	echo  "*** Processing $motif ... $counter"
	cd $HOME/usr/meme/bin
	species_bmf=$PROJ_DIR/resources/cisbp_all_species_bg_freq/$SPECIES.bmf
	species_bmf_uniform=$PROJ_DIR/resources/cisbp_all_species_bg_freq/_uniform_species.bmf
	if [ -e $species_bmf ]; then
		echo "Using background freq file $species_bmf"
		./fimo -o $PROJ_DIR/$OUT_FIMO/$motif --thresh 5e-3 --bgfile $species_bmf $PROJ_DIR/$IN_TF_PWM/$motif $PROJ_DIR/resources/s_cerevisiae.promoters.fasta
	else
		echo "Using uniform background freq file"
		./fimo -o $PROJ_DIR/$OUT_FIMO/$motif --thresh 5e-3 --bgfile $species_bmf_uniform $PROJ_DIR/$IN_TF_PWM/$motif $PROJ_DIR/resources/s_cerevisiae.promoters.fasta 
	fi
	sed ' 1d ' $PROJ_DIR/$OUT_FIMO/$motif/fimo.txt | cut -f 1,2,7 > $PROJ_DIR/$OUT_FIMO/$motif/temp.txt
	ruby $PROJ_DIR/scripts/estimate_affinity.rb -i $PROJ_DIR/$OUT_FIMO/$motif/temp.txt > $PROJ_DIR/$OUT_FIMO/$motif.summary
	mv $PROJ_DIR/$OUT_FIMO/$motif/fimo.txt $PROJ_DIR/$OUT_FIMO/$motif.fimo
	rm -r $PROJ_DIR/$OUT_FIMO/$motif
	echo "*** Done"
done < $PROJ_DIR/$IN_TF_LIST

# # Compute the rankded lists
# echo "*** Computing rankings ... "
# python $PROJ_DIR/scripts/rank_fimo.py $PROJ_DIR/$OUT_FIMO -t $PROJ_DIR/resources/target_names.txt -o $PROJ_DIR/$OUT_RANK
# echo "*** Done"

echo "*** ALL DONE! ***"
