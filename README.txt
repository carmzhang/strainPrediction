README - Carolyn Zhang 20190506

Code-20190503 - sequence analysis (folder):

(1) Scripts used for sequence analysis

Data Files (folder):

(1) blanked data - folder contains all of the growth dynamics data that have been blanked and stored in .mat formatted files

(2) data dependency files - all data files required to run the program files

(3) environmentalIsolates_data&code - the data files 

(4) raw data files - all raw data files from plate readers, file name includes date of experiment, isolate numbers, and growth condition

(5) raw data keio - all raw data files from Kulp et al. which contain the growth dynamics of the keio collection

(6) 16s&SNP_data_20190504 - all .fasta files for generating phylogenic trees for the environmental isolates (by taxonomic order)

(7) all_strains.xlsx - list of strain IDs for the clinical isolates

Program files (folder):

(1) Functions sourced from online (MATHWORKS): 
	bramila_mantel.m
	ciplot.m
	densityplot
 	histf
	drawCurve.m
	svmdecision.m
	subaxis.m

(2) Files to generate figures:
	Run each file to recreate the figures. More details can be found in the comments at the top of each file.
	
	Main Text Figures:
		Figure_visualizations.m
		Figure2c.m
		Figure3_fourStrainExample.m
		Figure3c_confidenceInterval.m
		Figure4_correlation_envIsolates.m
		Figure4_correlations_clinicalIsolates.m
		Figure5a_wgs_phenotype_arcs.m
		Figure5b_resistance_prediction_GC_roc.m
	
	Supplementary Figures:
		Supp_correlationsClinicalIsolates.m
		SuppFigure_visualizeData.m
		predictStrain_clinicalIsolates_changetime.m - Supplementary Figure 1.1
		plot_phylotree_envIsolates.m (for Supplementary section 3.1 - not shown)

(3) Files to recreate predictions:
	Run each file to generate predictions of genetic identity or antibiotic resistance. More details can be found in the comments at the top of each file.

	predictResistance_clinicalIsolates.m
	predictResistance_WGS.m
	predictStrain_clinicalIsolates.m
	predictStrain_EnvironmentalIsolates.m
	predictStrain_keio_allPlates.m
	predictStrain_keio_singlePlate.m

(4) Files for figures and predictions using growth rate (alternative features):
	Run each file to generate either predictions or figures using growth rate as the features. More details can be found in the comments at the top of each file.

	Figure4_correlation_envIsolates_traditionalGR.m
	Figure4_correlations_clinicalIsolates_traditionalGR.m
	Figure5b_resistance_prediction_GC_roc_traditionalGR.m
	importClinicalIsolates_traditionalGR.m
	importEnvironmentalIsolates_traditionalGR.m
	predictResistance_clinicalIsolates_traditionalGR.m
	predictStrain_clinicalIsolates_traditionalGR.m
	predictStrain_EnvironmentalIsolates_traditionalGR.m
	predictStrain_keio_allPlates_traditionalGR.m
	predictStrain_keio_singlePlate_traditionalGR.m
	Supp_correlationsClinicalIsolates_traditionalGR.m

(5) Other files are functions used to run the underlying computations for other scripts.
	
*2-3 use time derivative of the growth curves as the features*



