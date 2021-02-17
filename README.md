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
The 'GPATree-vignette.pdf' provides a framework for the step-by-step genetic analysis using the 'GPATree' package. Help file generated using the following code also provides a good starting point to use the 'GPATree' package, including some example command lines:

```{r}
library(GPATree)
package?GPATree
```

# ShinyGPATree
The ShinyGPATree app provides visualization of the GPA-Tree model, identifies risk-associated SNPs, and characterizes the combinations of functional annotations that can describe the risk-associated SNPs. The app can also be utilized to improve the visualization of the GPA-Tree model fit to collate or separate layers of the model (add or remove leaves). The number of non-risk-associated and risk-associated SNPs that can be characterized by combinations of functional annotations are also automatically updated based on different user selected input options. The following commands can be used to initialize the ShinyGPATree app:

```{r}
fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
ShinyGPATree(fit)
```
initial desc
