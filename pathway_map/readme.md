## Structure of this folder:
```
/your-repo
├── TCGA-x/
│   ├── gene_ids.txt
│   ├── gene_to_id.csv
│   ├── pantherGeneList.txt (this should be downloaded from patherdb using gene_ids.txt as input)
│   └── pathway_counts.csv
├── pathway_map_project.py
├── pathway_map_sample.py
├── save_gene_id.py
└── readme.md
```
## Usage
pathway_map_sample.py is called by the script for processing sample of a single patient. 

pathway_map_project.py maps the mutations in a cohort to the corresponding pathways and saves it under the corresponding folder.

save_gene_id.py writes the gene_ids.txt and gene_to_id.csv under the corresponding folder. Then use the gene ids as an input to Panther database.