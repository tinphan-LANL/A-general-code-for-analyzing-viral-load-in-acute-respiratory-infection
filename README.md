# A-general-code-for-analyzing-viral-load-in-acute-respiratory-infection
LANL O# O4798

Brief Summary: Code developed to fit and analyze viral load data for acute respiratory infections such as SARS-CoV-2. The code can also implement different treatments and generate in-silico cohorts to examine the effects of treatments.

Detailed Description: The code includes a mathematical model and statical method for viral dynamics. The mathematical model and statistical method are built on previously published modeling studies of within-host viral dynamics in acute respiratory infection from LANL. The code is implemented in Monolix and MATLAB. The code fits the mathematical model to viral load data of a cohort of individuals using a non-linear mixed effect approach (via Monolix). THe distributions of the model parameters are used to generate cohorts of in-silico patients. Hypothetical treatments are then applied to these in-silico patients under different modeling assumptions. Standard statistical tests are then applied to detemrine statistical significance between different outcomes.

Order of code execution: 
  (1) Download the de-identified viral load data, model text file, Monolix and MATLAB files.
  (2) Run Monolix file for fitting.
  (3) Run MATLAB files in folder #1 to replot the model fit.
  (4) Run MATLAB files in folder #2 to generate in-silico cohort.
  (5) Run MATLAB files in folder #3 to plot the results of different treatments (subfolders of folder #3).

  This set of code will generate main text figures 2 and 3 in https://doi.org/10.1101/2024.09.13.613000.  

  Straightforward modifications of these codes will also give the figures in the SM. If you have any questions, email ttphan (at) lanl (dot) gov.
