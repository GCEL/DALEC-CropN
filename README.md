# DALEC-CropN
DALEC-Crop model developments for simulating the impacts of N dilution on the ACM equation used to determine the leaf N and temperature limited rate of photosynthesis.

Code used for inferring crop leaf nitrogen per area (LNA) based on the calibration of a semi-empirical critical N dilution function derived from the ATEC project experimental field trials. Code adapted from DALEC_GSI_DFOL_FR_CROP.f90 - key changes are the use of N dilution slope and intercept parameters applied in lines 376 to 389. From a calibration of parameters 1-34 from a prior calibration step the CARDAMOM MDF framework optimises the slope and intercept of the N dilusion (parameters 35 and 36) in order to find an N dilution fit consistent with an LAI observational time-series.

![TRIALS_Github](https://user-images.githubusercontent.com/43847496/97121237-9dd99800-1714-11eb-8e95-5672db870db2.png)

![N_fit_Github](https://user-images.githubusercontent.com/43847496/97121388-bdbd8b80-1715-11eb-9fd9-b6ebfeb1e580.png)
