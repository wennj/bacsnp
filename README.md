 <img src="https://github.com/wennj/bacsnp/blob/master/vignettes/bacsnp_logo.png" height="100"> 

# bacsnp

Bacsnp is a tool that can help decipher populations of baculoviruses. It is written in R programming language and uses files of the Variant Call Format, which were created with Mpileup and is able to assign specificities to detected SNP positions.

How the bacsnp can be integrated into a workflow can be found in detail in the following publication:

*Wennmann, Jörg T.; Fan, Jiangbin; Jehle, Johannes A. **2020**. ["Bacsnp: Using Single Nucleotide Polymorphism (SNP) Specificities and Frequencies to Identify Genotype Composition in Baculoviruses."](https://www.mdpi.com/1999-4915/12/6/625) Viruses 12, no. 6: 625.*

## Getting Started

In the following section I will try to explain how to install bacsnp for R on your computer.

### Prerequisites

The tool was written in the R programming language and therefore requires the latest version of R to be installed on the computer. In addition, I recommend using RStudio.
The following packages for R should also be installed:

* devtools
* ggplot2
* VariantAnnotation (Bioconductor package)

To install the **devtool** and **ggplot2** package, the following steps can be carried out in the R console:

```
install.packages("devtools")
install.packages("ggplot2")
```

Make sure you have the latest version of **Bioconductor** installed. You can install it for R by following the instructions on the [Bioconductor homepage](https://www.bioconductor.org/install/).
Then you can install the **VariantAnnotation** package.

```
BiocManager::install("VariantAnnotation")
```

Check if everything is installed correctly by starting the packages.

```
library(devtools)
library(ggplot2)
library(VariantAnnotation)
```

The packages should start without any error messages.

### How to install bacsnp from Github

Next we look at how we install bacsnp directly from Github in R. This step requires the devtools package, which is available via CRAN. It must be installed and initiated.

```
library(devtools)
```

The bacsnp package is locaetd in my repository wennj/bacsnp. For its installation the you type

```
install_github("wennj/bacsnp")
```

After installation, you should be able to start it as usual to use it straight away.

```
library(bacsnp)
```

## Running the tests

To be continued.

### Break down into end to end tests

To be continued.

### And coding style tests

To be continued.

## Contributing

To be continued.

## Versioning

To be continued.

## Authors

To be continued.

## License

To be continued.

## Acknowledgments

To be continued.
