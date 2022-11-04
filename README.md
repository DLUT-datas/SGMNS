SGMNS
==
SGMNS User Manual
--
* *1. Introduction and Installation<br>

* *2. Input files<br>

	* *2.1 Background database<br>

	* *2.2 Files for structure annotation<br>

* *3. Software running instructions<br>

	* *3.1 Parameter settings for MGMN<br>

	* *3.2 Metabolite annotation<br>
		* *3.2.1 Using the default background database
		* *3.2.2 Building a new background network


 
# 1. Introduction and Installation

* Structure-Guided Molecular Network Strategy (SGMNS) is a free network-based and spectral library-independent tool for deep annotation of untargeted ultra-performance liquid chromatography-high resolution mass spectrometry (UPLC-HRMS) metabolomics data. The source program is in the Python language, and the development environment is Spyder 3.7.

# 2. Input files

## 2.1 Background database
* If you want to use the default background network, you can skip this step.<br>
* If you want to build a new background network instead of using the default one for the subsequent metabolite annotation, you need to provide a background database form. This form should be stored in csv format file and be named ‘Background database_POS.csv’ for positive ion mode and ‘Background database_NEG.csv’ for negative ion mode (Figure 1). The background database contains the ID, SMILES, exact mass, retention time (RT/min), formula, and Name. 

 ![image](https://user-images.githubusercontent.com/49905529/165009662-69977fed-3277-405f-bca5-c048ed669343.png)<br>
Figure 1. Background database

## 2.2 Files for structure annotation
* Peak table: The original data obtained from the biological sample after UPLC-HRMS analysis was processed by the peak picking software (such as Compound Discoverer, MS-DIAL, MZmine 2, XCMS, etc.) to obtain the peak table, including MS1, RT, exact mass, and calibration RT (Figure 2). The experimental RT should be corrected to the chromatographic condition used in the background database, which is called calibration RT. Then the ion fusion tools were used to remove the redundancy of metabolic features, such as CAMERA, mz.unity, and so on. The peak table should be named ‘Peak table_POS.csv’ for positive ion mode and ‘Peak table_NEG.csv’ for negative ion mode.

 ![image](https://user-images.githubusercontent.com/49905529/165009685-967df104-2e96-46ea-a639-2e43b5e96d6e.png)<br>
Figure 2. Peak table

* MGF files: The raw data was translated to the MGF files with the open-source software MSConvert, which was known as MS/MS spectra data (Figure 3). The MGF files should be named ‘MSMS1.mgf’, ‘MSMS2.mgf’, …, ‘MSMSn.mgf’.

 ![image](https://user-images.githubusercontent.com/49905529/165009723-b8cd79d4-ecd3-49e9-adca-6fa4b4f9b909.png)<br>
Figure 3. MGF files

* Seed metabolites: A list of seed metabolites is needed, in which the seed metabolites were identified by searching the MS/MS spectral library (Figure 4). The list of seed metabolites should be named ‘Seed metabolites_POS.csv’ for positive ion mode and ‘Seed metabolites_NEG.csv’ for negative ion mode.

 ![image](https://user-images.githubusercontent.com/49905529/165009744-b3a8f903-2ee0-4473-a968-92bb4834249e.png)<br>
Figure 4. Seed metabolites

# 3. Software running instructions
* The executable file is published in a folder, as shown in Figure 5. Double-click SGMNS.exe in the folder ‘SGMNS’ to run the executable file.

 ![image](https://user-images.githubusercontent.com/49905529/165009757-6060d4af-14cf-4597-9e53-6e8b8fc58500.png)<br>
Figure 5. The folder where the software is published

## 3.1 Parameter settings for SGMNS
* Enter parameters in turn according to the software prompts, and then press ‘Enter’ to continue. The default parameters are shown in Figure 6.

 ![image](https://user-images.githubusercontent.com/49905529/165009772-fde6511b-4ec7-4b4a-8248-6ef762ff716c.png)<br>
Figure 6. Software interface for parameter settings

## 3.2 Metabolite annotation
### 3.2.1 Using the default background database
* When you decide to use the default background network for metabolite annotation, the software interface is shown in Figure 7.
* If you use SGMNS for the first time, you should set parameters as shown in Figure 8.

 ![image](https://user-images.githubusercontent.com/49905529/165009783-7992d029-3086-4532-9096-6b6a8e2e81e8.png)<br>
Figure 7. Software interface for Metabolite annotation using the default background network

### 3.2.2 Building a new background network
* When you decide to build a new background network for the subsequent metabolite annotation, the software interface is shown in Figure 8.

  ![image](https://user-images.githubusercontent.com/49905529/165009789-dff81596-c6fe-4cd8-a162-8a5d960330d1.png)<br>
Figure 8. Software interface for Metabolite annotation using a new background network
