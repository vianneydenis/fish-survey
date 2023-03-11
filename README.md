<p align="center">
  <img alt="logo" src="./.github/logo.png" width="25%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img alt="fish" src="./.github/fish.png" width="25%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img alt="fish" src="./.github/dna.png" width="25%">
</p>

# Data from: Navigating the scales of diversity in subtropical and coastal fish assemblages ascertained by eDNA and visual surveys

---

This Hsuetal_2022_dataset_README.txt file was generated on 2022-08-05 by Tsai-Hsuan Tony Hsu (tonytsaihsuanhsu@gmail.com).



1. **Author Information**

	First author
		Name: Tsai-Hsuan Tony Hsu
		Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan

	Second author
		Name: Dr. Wei-Jen Chen
		Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
 		
	Third & corresponding author 
		Name: Dr Vianney Denis
		Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
		Email: vianney.denis@gmail.com

2. **Date of data collection**: 2021

3. **Geographic location of data collection**: northern Taiwan, West Pacific

4. **Funding sources that supported the collection of the data**: Ocean Conservation Administration of Taiwan, Ministry of Science and Technology of Taiwan 

5. **Recommended citation for this dataset**
Hsu, Tsai-Hsuan Tony; Chen, Wei-Jen; Denis, Vianney (2023), Data from: Navigating the scales of diversity in subtropical and coastal fish assemblages ascertained by eDNA and visual surveys, Dryad, Dataset, https://doi.org/10.5061/dryad.3xsj3txk4

[Data and R script are also available through the GitHub repository https://github.com/vianneydenis/fish-survey.git in order to replicate our analyses]


## DESCRIPTION OF THE DATA AND FILE STRUCTURE

### DATA & FILE OVERVIEW

1. **Description of dataset**

These data were generated to navigate different scales of fish diversity and compare results of eDNA with those from visual surveys in northern Taiwan.

2. **File List**

   File 1: Hsuetal_dataset_Site.csv
   File 1 description: Information of the 21 sampling sites

   File 2: Hsuetal_dataset_DOV.csv
   File 2 description: Fish abundance and total length data from diver-operated video surveys among the 21 sampling sites
  
   File 3: Hsuetal_dataset_UVC.csv
   File 3 description: Fish abundance and total length data from underwater visual census surveys among the 21 sampling sites

   File 4: Hsuetal_dataset_eDNA.xlsx
   File 4 description: Species annotation and read coverage for fish eDNA data among the 21 sampling sites

   File 5: Hsuetal_dataset_BV.csv
   File 5 description: Count of benthic variables identified by the CoralNet among the 21 sampling sites

   File 6: Hsuetal_dataset_level.csv
   File 6 description: Species' vertical position in the water column

   File 7: Hsuetal_dataset_Lw_ab.csv
   File 7 description: Length-weight equation coefficients a and b retreived from FishBase to calculate biomass

   File 8: Hsuetal_dataset_flow.csv
   File 8 description: Flow data recorded by a current meter, Marotte HS, at each site

   File 9: Hsuetal_dataset_gd.csv
   File 9 description: Pairwise gepgraphic distance (km) between sites

   File 10: Hsuetal_dataset_NGS.zip
   File 10 description: A zip file containing paired-end fastq files obtaining from NGS for the 21 sites and the filtration (FB) & extraction blanks (EB)

### METHODOLOGICAL INFORMATION

A detailed description of data acquisition and processing can be found in the published manuscript in the Ecological Indicators (https://doi.org/10.1016/j.ecolind.2023.110044).


#### DATA-SPECIFIC INFORMATION


##### **Hsuetal_dataset_Site.csv**

1. Number of variables/columns: 6

2. Number of cases/rows: 22

3. Missing data codes

    None

4. Variable List

    + Column A - Date: sampling date
    + Column B - Site: survey sites
    + Column C - Code: abbreviation for each site
    + Column D – Lattitude
    + Column E – Longitude
    + Column F – Area: Northern coast and outlying islands
    

5. Abbreviations used
  
    + for Area: 
      + C: Coastal sites along the northern coast
      + I: Outlying islands

##### **Hsuetal_dataset_DOV.csv**

1. Number of variables/columns: 6

2. Number of cases/rows: 429

3. Missing data codes 

    None

4. Variable List

    + Column A – Method
    + Column B – Site: abbreviation for each site
    + Column C – Family
    + Column D – Species
    + Column E – Length: total length of fish individuals
    + Column F – Number: the number of fish individuals

##### **Hsuetal_dataset_UVC.csv**

1. Number of variables/columns: 6

2. Number of cases/rows: 591

3. Missing data codes: 

    None

4. Variable List 

    + Column A – Method
    + Column B – Site: abbreviation for each site
    + Column C – Family
    + Column D – Species
    + Column E – Length: total length of fish individuals (cm)
    + Column F – Number: the number of fish individuals


##### **Hsuetal_dataset_eDNA.csv**

1. Number of variables/columns: 6

2. Number of cases/rows: 1704

3. Missing data codes

    None

4. Variable List

    + Column A – Site: abbreviation for each site
    + Column B – Family
    + Column C – Species: species name automatically assigned by the MiFish pipeline
    + Column D – Ratio: the coverage of the sequence reads among total reads in each site
    + Column E – Valid_as: confirmed species name
    + Column F – Chinese: Chinese name of the species

##### **Hsuetal_dataset_BV.csv**

1. Number of variables/columns: 11

2. Number of cases/rows: 22

3. Missing data codes 

    None

4. Variable List

    + Column A – Site: abbreviation for each site
    + Column B – CCA: crustose oralline algae
    + Column C – hard_coral: Scleractinain corals
    + Column D – macroalgae: macroalgae
    + Column E – other_life: other living organisms
    + Column F – soft_coral: soft corals
    + Column G – sponge: sponge
    + Column H – sub_stable: stable substrates
    + Column I – sub_unstable: unstable substrates
    + Column J – turf: turf algae
    + Column K – zoanthid: zoanthids
    + Column L – other: other abiotic items

##### **Hsuetal_dataset_level.csv**

1. Number of variables/columns: 2

2. Number of cases/rows: 483

3. Missing data codes 

    None

4. Variable List 

    + Column A – Species
    + Column B – level: species' vertical positions in the water column

5. Abbreviations used: 

    + for level
      + B: benthic
      + P: pelagic
      +BP: bentho-pelgic
      
##### **Hsuetal_dataset_Lw_ab.csv**    
      
1. Number of variables/columns: 3

2. Number of cases/rows: 129

3. Missing data codes 

    None

4. Variable List 

    + Column A – Species
    + Column B – mean_a: averaged coefficient a
    + Column C – mean_b: averaged coefficient b

##### **Hsuetal_dataset_flow.csv**  

1. Number of variables/columns: 4

2. Number of cases/rows: 22

3. Missing data codes 

    None

4. Variable List 

    + Column A – Date: sampling date
    + Column B – Site: survey sites
    + Column C - Code: abbreviation for each site
    + Column D – Speed: flow speed in m/s
    
##### **Hsuetal_dataset_gd.csv**  

1. Number of variables/columns: 21

2. Number of cases/rows: 22

3. Missing data codes 

    None

4. Variable List

    + Column A – Site
    + Column B to V – Code of the 21 survey sites




