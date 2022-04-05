# GCA_Code
This is code for paper: Heteroscedastic Growth Curve Modeling with Shape-Restricted Splines

## chickweight_analysis.R
End-to-end code to fit chickweight data with the proposed iteratively reweighted methods for clustered data.

## pancrea_analysis.R
End-to-end code to fit pancrea data with the proposed iteratively reweighted methods for independent data. I-splines and restrictions were applied to fit monotone increasing pattern.

## pancrea_cspline.R
End-to-end code to fit pancrea data with the proposed iteratively reweighted methods for independent data. C-splines and restrictions were applied to fit monotone increasing and concave pattern.

## lmerIter.R and lmerMLE.R
Self-defined functions for the proposed models and algorithm. They will be used in the chickweight and pancrea data analysis code.

Other .R scripts are supplement functions to summarize the fitting results and get the fitted value, such as the error variance. 
