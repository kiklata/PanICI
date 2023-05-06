path=/home/shpc_100839/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/pre
outpath=/home/shpc_100839/PanCancerICI/TumorMetabolism/Data/scData/AntiPD1_BC/compass
file=BIOKEY_7

mkdir $outpath/$file

compass --data $path/$file.tsv --num-processes 24 --species homo_sapiens --output-dir $outpath/$file --microcluster-size 10 

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
