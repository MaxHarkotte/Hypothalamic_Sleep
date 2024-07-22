# Hypothalamic Sleep

Scripts for the analysis of ephys recordings from the rats' lateral hyothalamus combined with surface EEG and EMG. 

## Data 
All data is stored and documented on the born_animal file server in 'born_animal\Hypothalamic_Sleep\data'. Details about the recording parameters are stored under \docs.

## Repository structure 
- **src/** contains all source code for analysis in Matlab and Python 
- **notebooks/** contains jupyter notebooks for exploratory data analysis and interactive results 
- **results/** contains outputs from the analysis 
- **docs/** contains all documentation for the project. References (papers, etc.) are shared in the Zotero Group.
- **environment/** contains files specifying the project environment for both matlab and python based anaylsis 

### Preprocessing 
Raw data is stored as .ncs (Neuralynx file format) and preprocessed in matlab using fieldtrip. 

## Cluster access
To run the analysis on the CIN cluster you need your CIN Account details and connect via ssh with the following commands: 

ssh cin_account@172.25.250.112 -p 60222 -X

ssh born1 -X

Matlab is installed here: 

/usr/local/MATLAB/R2020b/bin/matlab -&

Fieldtrip is installed here: 

'born_animal\Hypothalamic_Sleep\fieldtrip'

