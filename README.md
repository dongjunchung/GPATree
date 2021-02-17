# GPATree
GPA-Tree is a statistical approach to integrate GWAS summary statistics and functional annotation information within a unified framework. Specifically, by combining a decision tree algorithm with a hierarchical modeling framework, GPA-Tree simultaneously implements association mapping and identifies key combinations of functional annotations related to disease risk-associated SNPs. 

# Installtion
To install the development version of GPATree, use the following commands.

```{r}
#install.packages("devtools")
library(devtools)
install_github("asthakhatiwada/GPATree")
```

# Usage
The 'GPATree' vignette provides a framework for the step-by-step data analysis using the 'GPATree' package. The following two help pages provide a good starting point for the genetic analysis using the 'GPATree' package, including the overview of 'GPATree' package and the example command lines:

```{r}
library(GPATree)
package?GPATree
```

# ShinyGPATree
The ShinyGPATree app provides visualization of the GPA-Tree model, identifies risk-associated SNPs, and characterizes the combinations of functional annotations that can describe the risk-associated SNPs. The app can also be utilized to improve the visualization of the GPA-Tree model fit to collate or separate layers of the model (add or remove leaves) by changing the $cp$ parameter. The number of non-risk-associated and risk-associated SNPs that can be characterized by combinations of functional annotations are also automatically updated based on different user selected input options. The following commands can be used to initialize the ShinyGPATree app:

```{r}
fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
ShinyGPATree(fit)
```
