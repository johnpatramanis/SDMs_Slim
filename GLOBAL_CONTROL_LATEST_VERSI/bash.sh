
for KAY in 1 2 3 5 8
do

for SIGMA in 2
do

slim -d seed=1000 -d sigma=$SIGMA -d K=$KAY SDM_Version_2_exper.slim

done
done

python3 SDM_analysis.py
python3 compare_matrices.py 
