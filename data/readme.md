Structure of this folder:
```
/your-repo
├── cohort_data/
├── TCGA-x/
├── data_splitter.py
└── readme.md
```

When working on a new cohort, put the aggregated data under cohort_data. This should be a single .maf file. Then run
```
python data_splitter.py <project_name>
```
If .maf files is already available for each patient, directly put those under the corresponding cohort.