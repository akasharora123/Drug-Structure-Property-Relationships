# Quantitative Drug Structure-Property Relationships (QSPR)
> Featurization/Fingerprinting of chemical molecules, and linear regression to predict solubility as a function of chemical structure of drugs  

## Overview
Quantitative structure-property relationship (QSPR) modeling is a regression or classification modeling method used in drug-discovery or Cheminformatics field to identify potential drug molecules for targeted diseases. QSPR modeling is based on a fundamental assumption that similar molecules have similar activities. With this assumption, one of the goals of QSPR modeling is to identify new molecules that have property (or activity) similar to the property of a molecule that is proven or strongly-hypothesized to cure or prevent the disease.

QSPR modeling involves two steps:
- Featurization/Fingerprinting: This step transforms the chemical structures into molecular descriptors or numerical vectors that are further as inputs to machine-learning models. This is a crucial step in QSPR modeling as the generated descriptors are expected to contain all the necessary information about substructures/chemical groups present in the drug molecule, and these substructures/chemical groups are highly correlated the property or activity of the drugs. There are several schemes to do fingerprinting ranging from as simples as dictionary of molecular properties to as complex as the use of Graph Neural Network. Here we use the fours schemes: 
	-- Morgan (or Circular or ECFP) and topological fingerprints.    
	-- RDKit topological fingerprints
	-- MACCS Keys
	-- RDKit Molecular Descriptors

- Machine-learning modeling: The second step is to fit an appropriate regression (or classification) model to correlate the molecular descriptors with the property/activity of the drug molecules. Here, we demonstrate this step for an important physico-chemical property of drug molecules called aqueous solubility (S) (typically measured as log(S)). We compare various regression models and examine their performance using common metrics. 

## Data and References:
- Data: MoleculeNet/Datasets (http://moleculenet.ai/datasets-1)
- ESOL:  Estimating Aqueous Solubility Directly from Molecular Structure (https://pubs.acs.org/doi/abs/10.1021/ci034243x) 
- Fingerprinting: Molecular Descriptors and Fingerprints (http://datascience.unm.edu/biomed505/Course/Cheminformatics/basic/descs_fingers/molec_descs_fingerprints.htm)
 

## Prerequisites
- Python
- RDKit
- Pandas
- Seaborn
- MatplotLib
- Scikit-learn

## License
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)
- **[MIT license](http://opensource.org/licenses/mit-license.php)**

