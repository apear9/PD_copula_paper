# READ ME

This repository contains R code that can be used to reproduce the figures and analysis in Pearse, A.R. &amp; Bondell, H. (2025). Power-divergence copulas: A new class of Archimedean copulas. arXiv:2510.06177

The subfolders in this repository correspond to version 1 (V1) and version 2 (V2) of the paper, as shown on arXiv.

Every script has been named according to the following template:

[Sec#]\_[Description]\_[Fig#].R

where

* [Sec#] corresponds to the Section of the manuscript or its supplement where the results are presented.
* [Description] is a short description of what is produced by the script.
* [Fig#] is optional, and it indicates which figure(s) from the manuscript or its supplement is (are) produced by the script.

All scripts in the subfolder V1 do not require any external data, and they only rely on data provided by R packages.

Some scripts in the subfolder V2 require a download of freely and publicly available data from FRED:

* S&P Dow Jones Indices LLC, Dow Jones Industrial Average [DJIA], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/DJIA, November 21, 2026. 
* S&P Dow Jones Indices LLC, S&P 500 [SP500], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/SP500, November 21, 2026. 