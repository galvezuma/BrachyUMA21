# BrachyUMA21
Short Java programs to deal with Brachypodium distachyon ecotypes.

This project contains three programs to extract information from data files extracted from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). In particular:
* **es.uma.keyword.LookForKeyword**: Allows iterating through the fasta and gff3 files associated to several species and extract information related to a family of genes. Before working with this package, you have to download some files from Phytozome, configure some lines of source Java files and recompile.
* **es.uma.html.GenerateSVG**: Takes as input the previously extracted data and generates an HTML file with the distribution of the genes. This file is specifically prepared to work with the 10 dehydrins found in *Brachypodium distachyon*. This program can be executed as is.
* **es.uma.motif.HTMLDecorator**: Iterates over a fasta file and converts it into an HTML file with highlighted motives. This program can be executed as is.

These programs are not intented to simply execute them but, instead, the code must be modified slightly if any developer when to adapt to their needs. In the next sections are explained how to work with them.

## es.uma.keyword.LookForKeyword
This program uses as input an external structure of files that the user has to download from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). This has been tested with an use case in which 54 ecotypes of *Brachypodium distachyon* have been downloaded. Before compilation, you have to change the line:
```Java
    private static final String DIR_NAME = "genomes/Phytozome/PhytozomeV13";
```
to adapt to your own directory that contains the files taken from Phytozome. In particular, you have to download the files with extension `.gene.gff3.gz`, `.protein_primaryTranscriptOnly.fa.gz` and `genes.fasta`. In addition, the file `datafiles/Order.txt` contains the order in which you want to generate the data on the extracted genes. Finally, the genes extracted are those that verify the sentences:
```Java
if (line.toUpperCase().contains("PF00257") || 
    line.toUpperCase().contains("DEHYDRIN")) { ...
```
By default this refers to the LEA2 gene family of Dehydrins whose PFAM code is PF00257. You have to modify these lines to adapt them to your needs.

Once these changes have been carried out, compilation is carried by executing the next commands:
```
cd scr
javac es/uma/keyword/LookForKeyword.java
```
and executed using 
```
java es.uma.keyword.LookForKeyword
```
The file `datafiles/dehydrins_ecotypes.data` contains an example of the output of this program once executed over 54 ecotypes of *Brachypodium distachyon* downloaded from Phytozome. In addition, it is also given a file `datafiles/Order.txt` that refers to these 54 ecotypes.

## es.uma.html.GenerateSVG
This program uses the data extracted from the previous one (es.uma.keyword.LookForKeyword). As an example, a file `datafiles/dehydrins_ecotypes.data` is provided so you could run it from scratch. Compilation is carried by executing the next commands:
```
cd scr
javac es/uma/html/GenerateSVG.java
```
and executed using 
```
cd ..
java -cp ./src es.uma.html.GenerateSVG
```
This will display on console a trace of the genes processed and, in daddition, a file `result.html` will be generated in the root directory. Such a file has the next appearance:


