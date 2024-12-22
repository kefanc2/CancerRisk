## Usage
Download the data and put it under a folder with the corresponding name. Then run
```
python get_gene_list.py <project_name>
python join_clinical.py <project_name>
```
## Structure:
```
count_processor.py <project name>
|
|-<project name>
  |-clinical.csv
  |-mutations.csv
  |-(counts.csv)
  |-(joined_counts.csv)
```