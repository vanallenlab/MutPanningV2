# MutPanningV2
This repository contains the source code of the revised version of MutPanning. MutPanning is publicly available under the BSD3-Clause open source license.

MutPanning is designed to detect rare cancer driver genes from aggregated whole-exome sequencing data. Most approaches detect cancer genes based on their mutational excess, i.e. they search for genes with an increased number of nonsynonymous mutations above the background mutation rate. MutPanning further accounts for the nucleotide context around mutations and searches for genes with an excess of mutations in unusual sequence contexts that deviate from the characteristic sequence context around passenger mutations.

INTRODUCTION

MutPanning analyzes aggregated DNA sequencing data of tumor patients to identify genes that are likely to be functionally relevant, based on their abundance of nonsynonymous mutations or their increased number of mutations in unusual nucleotide contexts that deviate from the background mutational process. 
 
The name MutPanning is inspired by the words "mutation" and "panning". The goal of the MutPanning algorithm is to discover new tumor genes in aggregated sequencing data, i.e. to "pan" the few tumor-relevant driver mutations from the abundance of functionally neutral passenger mutations in the background. Previous approaches for cancer gene discovery were mostly based on mutational recurrence, i.e. they detected cancer genes based on their excess of nonsynonymous mutation above the local background mutation rate.  Further, they search for mutations that occur in functionally important genomic positions, as predicted by bioinformatical scores). These approaches are highly effective in tumor types, for which the average background mutation rate (i.e., the total mutational burden) is low or moderate.
 
The ability to detect driver genes can be increased by considering the nucleotide context around mutations in the statistical model. MutPanning utilizes the observation that most passenger mutations are surrounded by characteristic nucleotide sequence contexts, reflecting the background mutational process active in a given tumor. In contrast, driver mutations are localized towards functionally important positions, which are not necessarily surrounded by the same nucleotide contexts as passenger mutations. Hence, in addition to mutational excess, MutPanning searches for genes with an excess of mutations in unusual sequence contexts that deviate from the characteristic sequence context around passenger mutations. That way, MutPanning actively suppresses mutations in its test statistics that are likely to be passenger mutations based on their surrounding nucleotide contexts. Considering the nucleotide context is particularly useful in tumor types with high background mutation rates and high nucleotide context specificity (e.g., melanoma, bladder, endometrial, or colorectal cancer).

ALGORITHM
 
Most passenger mutations occur in characteristic nucleotide contexts that reflect the mutational process active in a given tumor. MutPanning searches for mutations in “unusual” nucleotide contexts that deviate from this background mutational process. In these positions, passenger mutations are rare and mutations are thus a strong indicator of the shift of driver mutations towards functionally important positions.
 
The main steps of MutPanning are as follows (adopted from Dietlein et al.): 
(i) Model the mutation probability of each genomic position in the human exome depending on its surrounding nucleotide context and the regional background mutation rate. 
(ii) Given a gene with n nonsynonymous mutations, use a Monte Carlo simulation approach to simulate a large number of random “scenarios” in which n or more nonsynonymous mutations are randomly distributed along the same gene . 
(iii) Compare the number and positions of mutations in each random scenario with the observed mutations in gene . Based on these comparisons, derive a p-value for the gene. 
(iv) Combine this p-value with additional statistical components that account for insertions and deletions, the abundance of deleterious mutations, and mutational clustering.

USAGE
 
You can run the algorithm for multiple cancer types at the same time. All you need to run MutPanning is a mutation file (*.maf) and a sample file (*.txt):
- Mutation File (*.maf): mutations that you would like to include in the analysis 
- Sample File (*.txt): contains sample IDs and assocates them with cancer types 
Unless you are familiar with the MutPanning algorithm, we recommend running MutPanning with standard parameters. MutPanning was only tested to run with standard parameters. A precompiled desktop version of MutPanning is available on cancer-genes.org. In addition, MutPanning can be run as a module on genepattern.org

COMPILATION AND EXECUTION OF THE SOURCE CODE IN THIS REPOSITORY

Precompiled files of the interactive desktop version of MutPanning can be downloaded from www.cancer-genes.org (section "Downloads") as well as the supplement of our paper (Dietlein et al., cf. below). Furthermore, MutPanning can be run as an online module on the GenePattern platform. If you prefer to compile and execute MutPanning from scratch (e.g., to make changes to the source code, to change the memory usage of MutPanning, or to run MutPanning via the command line), please follow the instructions in this section.
In order to compile the source code on your computer, the Java Development Kit (JDK) needs to be installed. If you are unsure whether JDK is installed on your computer, you can check the java version with the command
java -version

1.) Download the source code on your local computer. Click on "Clone or download" (green button at the top right of this page) and click on "Download ZIP" in the dropdown menu. Unzip the zip file, which creates a folder named "MutPanningV2-master" on your computer.

2.) On your local computer, open the command line (e.g., Terminal on macOS). Store the folder name where you downloaded the source code in the variable $ff via the command
ff="<folder name>"
For instance, on my computer:
ff="/Users/fdietlein/Downloads/MutPanningV2-master/"

3.) Compile your source code via the command
javac -classpath $ff/commons-math3-3.6.1.jar:$ff/jdistlib-0.4.5-bin.jar $ff/*.java
Alternatively, you can replace all occurrences of $ff with the folder name. For instance, on my computer:
javac -classpath /Users/fdietlein/Downloads/MutPanningV2-master/commons-math3-3.6.1.jar:/Users/fdietlein/Downloads/MutPanningV2-master/jdistlib-0.4.5-bin.jar /Users/fdietlein/Downloads/MutPanningV2-master/*.java
This command will compile all Java source code files in the folder $ff (e.g., on my computer /Users/fdietlein/Downloads/MutPanningV2-master/) using the libraries commons-math3-3.6.1.jar and jdistlib-0.4.5-bin.jar, which are passed to the compiler through the classpath argument.

4.) After the compilation of the source code, you can choose between two versions of MutPanning to execute. A) the interactive desktop version, which is equivalent to the version provided on cancer-genes.org and the paper supplement; B) a command-line version, which allows you to pass all the arguments through the command line directly. These arguments need to be manually entered through the dialog window in version A.

A) To execute the desktop version:
java -Xmx8G -classpath $ff/commons-math3-3.6.1.jar:$ff/jdistlib-0.4.5-bin.jar:$ff MutPanningDesktop
Variable $ff contains the folder where you downloaded the source code and was defined in step 2.
Alternatively, you can replace all occurrences of $ff by the folder name. For instance, on my computer:
java -Xmx8G -classpath /Users/fdietlein/Downloads/MutPanningV2-master/commons-math3-3.6.1.jar:/Users/fdietlein/Downloads/MutPanningV2-master/jdistlib-0.4.5-bin.jar:/Users/fdietlein/Downloads/MutPanningV2-master/ MutPanningDesktop

This command executes the class MutPanningDesktop and uses the libraries commons-math3-3.6.1.jar and jdistlib-0.4.5-bin.jar that are passed through the classpath argument. The argument -Xmx8G regulates the memory usage of MutPanning (maximum 8GB). In our tests, this argument worked for all maf files, including large maf files with >10,000 samples. If you process small maf files, you can reduce the memory usage of MutPanning (e.g., -Xmx4G to allocate 4GB). If you use larger maf files than used in our study, you may consider increasing the memory allocated to MutPanning (e.g., -Xmx16G to allocate 16GB). Upon execution of this command, a dialog window will open that guides you through the execution of MutPanning. Upon completion, the dialog window will automatically redirect you to the final results files.

B) To execute the command line version:
The command-line version expects the arguments below in the following order. 
1) root file, where to write intermediate and result files
2) maf file (standard value: root file/MutationsComplete.maf)
3) sample annotation file (standard value: root file/SamplesComplete.txt)
4) path to Hg19 folder (standard value root file/Hg19/)
5) minimal no. samples per cluster (standard value 3)
6) minimal no. mutations per cluster (standard value 1000)
7) min no. samples for CBASE (standard value 100)
8) min no. mutations for CBASE (standard value 5000)

Since all arguments are passed as a String to MutPanning, they have to be surrounded by quotation marks, including the integer values in arguments 5-8. Arguments should be separated by a space (no comma, tab, etc.). Only the first argument is mandatory. If maf and sample file, as well as the Hg19 folder, are in the root file (cf. standard values above), you do not have to specify them explicitly. 
If you would like to execute MutPanning with the standard parameters (recommended), you only need to specify arguments 1-4.
For instance, to execute the command line version on my computer with root folder in /Users/fdietlein/Desktop/Test/, maf file in /Users/fdietlein/Documents/DockerRevised/mutpanning_testfiles/MutationsSkin.maf, sample file in /Users/fdietlein/Documents/DockerRevised/mutpanning_testfiles/SamplesSkin.txt and the Hg19 folder in /Users/fdietlein/Downloads/FilesMutPanningResourse/Hg19/


java -Xmx8G -classpath $ff/commons-math3-3.6.1.jar:$ff/jdistlib-0.4.5-bin.jar:$ff MutPanning "/Users/fdietlein/Desktop/Test/" "/Users/fdietlein/Documents/DockerRevised/mutpanning_testfiles/MutationsSkin.maf" "/Users/fdietlein/Documents/DockerRevised/mutpanning_testfiles/SamplesSkin.txt" "/Users/fdietlein/Downloads/FilesMutPanningResourse/Hg19/"

Variable $ff contains the folder where you downloaded the source code and was defined in step 2. Alternatively, you can replace all occurrences of $ff by the folder name. For instance, on my computer:

java -Xmx8G -classpath /Users/fdietlein/Downloads/MutPanningV2-master/commons-math3-3.6.1.jar:/Users/fdietlein/Downloads/MutPanningV2-master/jdistlib-0.4.5-bin.jar:/Users/fdietlein/Downloads/MutPanningV2-master/ MutPanning "/Users/fdietlein/Documents/MutPanningTest/" "/Users/fdietlein/Documents/mutpanning_testfiles/MutationsSkin.maf" "/Users/fdietlein/Documents/mutpanning_testfiles/SamplesSkin.txt" "/Users/fdietlein/Documents/Hg19/"

This command executes the class MutPanning with the specified argument and uses the libraries commons-math3-3.6.1.jar and jdistlib-0.4.5-bin.jar that are passed through the classpath argument. The argument -Xmx8G regulates the memory usage of MutPanning (maximum 8GB). In our tests, this argument worked for all maf files, including large maf files with >10,000 samples. If you process only small maf files, you can reduce the memory usage of MutPanning (e.g., -Xmx2G to allocate 2GB). If you use larger maf files than processed in our study, you may consider increasing the memory allocated to MutPanning (e.g., -Xmx16G to allocate 16GB). Upon execution of this command, MutPanning launches with the specified arguments. You can follow the progress in the command line window. Upon completion, the results are in the root file specified in argument 1 (subfolder "Significance Filtered").


REFERENCES
 
Dietlein F, Weghorn D, Taylor-Weiner A, Richters A, et al. Identification of cancer driver genes based on nucleotide context. Under review. (preprint available on biorxiv)
