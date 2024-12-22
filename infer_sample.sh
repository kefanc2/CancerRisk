

# echo "Enter the barcode: "  
# read barcode
# echo "Enter maf file name"
# read mafname

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 arg1 arg2"
    exit 1
fi

# Access command-line arguments
project_name="$1"
barcode="$2"

# Prepare ssm file
python3 ssm_parser.py data/$project_name/$barcode.maf ./result/$barcode.txt

echo "running ssm complete"

# cleaning up
mv result/$barcode.txt ./phylowgs/input/
cd ./phylowgs
rm -rf chains
mkdir chains

source activate py2
echo "environment activated, running multievolve..."

# Run the multievolve script
python multievolve.py --num-chains 4 --ssms ./input/$barcode.txt --cnvs empty_cnv.txt

echo "running complete, writing results..."

# collect results
mkdir $barcode
python write_results.py $barcode ./chains/trees.zip $barcode/$barcode.summ.json.gz $barcode/$barcode.muts.json.gz $barcode/$barcode.mutass.zip

echo "unzipping results..."

# clean up results
cd $barcode
unzip $barcode.mutass.zip
mkdir $barcode.mutass
mv *.json ./$barcode.mutass
gzip -d $barcode.summ.json.gz
rm $barcode.mutass.zip
rm $barcode.muts.json.gz
cd ..
mv input/$barcode.txt $barcode
mv $barcode ../result/$project_name
cd ..

conda deactivate

# post processingq
echo "parsing tree to mutation csv..."

python3 tree_parser.py $project_name $barcode result/$project_name/$barcode/ project_data/$project_name/gene_list.txt result/$project_name/$barcode/$barcode.txt

# mapping to pathways
cd pathway_map
python3 pathway_map_sample.py ../result/$project_name/$barcode/$barcode.csv ../result/$project_name/$barcode/pathway_counts.csv $project_name/gene_to_id.csv $project_name/pantherGeneList.txt

# cleaning up temp files
cd ..
cd result/$project_name/$barcode
rm -rf $barcode.mutass
rm $barcode.summ.json
