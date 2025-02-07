# CancerRisk
Computationally reconstructing the evolution of cancer progression risk

## Abstract
Understanding the evolution of cancer in its early stages is critical to identifying key drivers of cancer progression and developing better early diagnostics or prophylactic treatments. Early cancer is difficult to observe, though, since it is generally asymptomatic until extensive genetic damage has accumulated.  In this study, we develop a computational approach to infer how once-healthy cells enter into and become committed to a pathway of aggressive cancer.  We accomplish this through a strategy of using tumor phylogenetics to look backwards in time to earlier stages of tumor development combined with machine learning to infer how progression risk changes over those stages.  We apply this paradigm to point mutation data from a set of cohorts from the Cancer Genome Atlas (TCGA) to formulate models of how progression risk evolves from the earliest stages of tumor growth, as well as how this evolution varies within and between cohorts.  The results suggest general mechanisms by which risk develops as a cell population commits to aggressive cancer, but with significant variability between cohorts and individuals.  These results imply limits to the potential for earlier diagnosis and intervention while also providing grounds for hope in extending these beyond current practice.


## Workflow
![plot](pipeline.png)

## Requirements
### 1. Create a virtual environment
Create a virtual environment for this repo and PhyloWGS, respectively. Note that PhyloWGS only runs in Python2.

### 2. Install dependencies
```
joblib==1.4.2
matplotlib==3.8.4
mygene==3.2.2
numpy==2.2.1
pandas==2.2.3
scikit_learn==1.4.2
scikit_survival==0.23.0
scipy==1.14.1
shap==0.46.0
tqdm==4.66.5
tqdm==4.66.1
```

### 3. Clone the PhyloWGS repository into the root directory of your repository and install it
See <https://github.com/morrislab/phylowgs>
The repo should now look like the follows:
```
/your-repo
├── phylowgs/
├── other-files-in-repo
└── README.md
```

More detailed instructions are included in the readme.md's in each subdirectory.
Weights and other data is available at <https://drive.google.com/drive/folders/1Nvst5XrAbCaG4Su1Zj0cqaay_Hci2OrI?usp=sharing>