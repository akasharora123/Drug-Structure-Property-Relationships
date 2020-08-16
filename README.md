# Quantitative Drug Structure-Property Relationships (QSPR)
> Featurization/Fingerprinting of chemical molecules, and linear regression to predict property as a function of chemical structure. 

## Overview
Quantitative structure-property relationship (QSPR) modeling is a regression or classification modeling method used in drug-discovery and Cheminformatics fields to identify potential drug molecules for targeted diseases. QSPR modeling is based on a fundamental assumption that similar molecules have similar properties/activities. With this assumption, one of the goals of QSPR modeling is to identify new molecules that have property/substructure similar to the property/substructure of a molecule that is proven or strongly-hypothesized to cure or prevent the disease.

QSPR modeling involves two steps:
- **Featurization/Fingerprinting:** This step transforms the chemical structures into molecular descriptors represented as numerical vectors which are required as inputs to machine-learning models (second step). This is a crucial step in QSPR modeling as the generated descriptors should contain all the necessary information about substructures/chemical groups present in the drug molecule which are ultimately responsible for imparting the specific property to the drug. There are several schemes to do fingerprinting, ranging from as simple as a dictionary of molecular properties to as complex as the use of Graph Neural Networks. Here we use four schemes (all generated using RDKit): 
	- Morgan (or Circular or ECFP) fingerprints    
	- RDKit topological fingerprints
	- MACCS Keys (dictionary based)
	- RDKit molecular descriptors (molecular properties)

- **Machine-learning modeling:** The second step is to fit an appropriate regression (or classification) model to correlate the molecular descriptors/features with the property of the drug molecules. Here, we demonstrate this step for an important physicochemical property called aqueous solubility (S) (typically measured as log(S)). We compare various regression models and examine their performance using common metrics. 

## Files Description
- qsar: main notebook for analysis
- featurizer.py: class (wrapper over RDKit) for generating different feature types
- utilities.py: some utilitiy function for plot aesthetics and model fitting
- best_regularized_models.py: class to find hyperparameters-optimized models using Gridsearch and Cross-validation.

## Prerequisites
- Python
- RDKit
- Pandas
- Seaborn
- Matplotlib
- Yellowbrick
- Statsmodels
- Scikit-learn

## Data and References:
- Data: MoleculeNet Datasets (http://moleculenet.ai/datasets-1)
- ESOL:  Estimating Aqueous Solubility Directly from Molecular Structure (https://pubs.acs.org/doi/abs/10.1021/ci034243x) 
- Molecular Descriptors and Fingerprints: (http://datascience.unm.edu/biomed505/Course/Cheminformatics/basic/descs_fingers/molec_descs_fingerprints.htm)

## License
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)
- **[MIT license](http://opensource.org/licenses/mit-license.php)**

