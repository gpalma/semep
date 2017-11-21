# semEP

## 1.  About

semEP is an edge partitioning approach that combines
a data mining framework for link prediction, semantic knowledge
(similarities) from ontologies, and an algorithmic approach
to partition the edges of a heterogeneous graph.

For more information about semEP see:

Guillermo Palma, Maria-Esther Vidal and Louiqa Raschid.
*Drug-Target Interaction Prediction Using Semantic Similarity and Edge
Partitioning*. Proceedings of the 12th International Semantic Web
Conference (ISWC 2013). Italy. 2014. [(PDF)](http://ldc.usb.ve/~gpalma/papers/semEP-ISWC14.pdf)

## 2. Content

* AUTHORS: list of contributors of the semEP project.
* doc: Documentation about semEP. 
* LICENSE: GPL version 2.
* Makefile: builds semEP.
* README.md: this file.
* src: source code.
* test: datasets to test semEP
  * test/yamanishi: link prediction dataset proposed by Yamanishi et al. [1].
  * test/pawels:   drug side-effect dataset proposed by Pauwels et al. [2].
* VERSION: software version.
  
[1] K. Bleakley and Y. Yamanishi.
*Supervised prediction of drug-target interactions
using bipartite local models*. Bioinformatics, 25(18):2397-2403, 2009.

[2] E. Pauwels, V. Stoven, and Y. Yamanishi.
*Predicting drug side-effect profiles: a chemical fragment-based approach*.
BMC Bioinformatics, 12:169, 2011.

## 3. License

GNU GENERAL PUBLIC LICENSE Version 2.

## 4. Requirements

* GNU Colecction Compiler (GCC) or Clang.
* GNU make (make).

semEP has been tested on FreeBSD, GNU/Linux, and OS X.

## 5. Installation

Clean and generate necessary files:

`$>make clean`

`$>make`

## 6. Usage

semEP has 7 command line arguments and all them
are mandatory. Here, mandatory means that without specifying this
argument, the program won't work.

semEP command synopsis:

`semEP  <-l left threshold> <-r right threshold> <left matrix>  <left vertices> <right matrix> <right vertices> <bipartite graph>`

where:

* __bipartite grap__: file with the bipartite graph (BPG) to study.
* __left threshold__: threshold of similarity between the left side vertices of the BPG.
* __right threshold__: threshold of similarity between the right side vertices of the BPG.
* __left matrix__: similarity matrix file with the similarities between the left side vertices of the BPG.
* __left vertices__: file with the list of vertices of the left side of the BPG. The order of the vertices corresponds to the __left matrix__. 
* __right matrix__: similarity matrix file with the similarities between the right side vertices of the BPG.
* __right vertices__: file with the list of vertices of the right side of the BPG. The order of the vertices corresponds to the __right matrix__. 

## 7. Running some samples:

* Compute the partitions in a drug-target bipartite graph and get predicted links on nuclear receptor dataset:

`./semEP -l 0.3061 -r 0.1614 test/yamanishi/nr/nr_matrix_drugs.txt test/yamanishi/nr/nr_drugs.txt test/yamanishi/nr/nr_matrix_targets.txt test/yamanishi/nr/nr_targets.txt test/yamanishi/nr/nr_drug-target_graph.txt`

* Compute the partitions in a drug-target bipartite graph and get predicted links on enzyme dataset:

`./semEP -l 0.2174 -r 0.0198 test/yamanishi/e/e_matrix_drugs.txt test/yamanishi/e/e_drugs.txt test/yamanishi/e/e_matrix_targets.txt test/yamanishi/e/e_targets.txt test/yamanishi/e/e_drug-target_graph.txt`

## 8 semEP input

### 8.1. File format of the bipartite graph

A bipartite graph is a graph whose vertices can be divided into
two disjoint sets *L* and *R* such that every edge connects a
vertex in *L* to one in *R*. We can denote a bipartite graph
as *G=(L,R,E)*, with *E* denoting the edges of the graph. We call
the set *L* as *the left side* of the bipartite graph, and the set *R*
as *the right side* of the bipartite graph. 

The first line of the bipartite graph file, contains the number of edges.
The following lines in the file shows the edges data, where
each line corresponds to a edge and contains the node of the left
side, the node of the right side, and the cost of the edge.

	[number of edges n]
	[right-node edge-1][TAB][lef-tnode edge-1][TAB][edge-cost-1]
	..
	..
	..
	[right-node edge-n][TAB][left-node edge-n][TAB][edge-cost-n]

### 8.2. File format of the similarity matrix

Files with the matrices must contain the similarities between the vertices of the left side,
or the right side. The file format is as follows:

	[number of rows and columns]
	[sim vertex-1 vertex-1][SPC]...[SPC][sim vertex-1 vertex-n]
	..
	..
	..
	[sim vertex-n vertex-1][SPC]...[SPC][sim vertex-n vertex-n]

### 8.3. File format of the list vertices

The file format is as follows:

	[number of vertices n]
	[vertex 1]
	..
	..
	..
	[vertex n]

Where *[vertex x]* is the identifier of vertex in the left side or right 
side of the bipartite graph. The order of each vertex correspond to the
position of the vertex in the similarity matrix.

## 9. semEP ouput

semEP produces two outputs:

1. a directory with suffix *-Clusters* that contains the clusters with the
edge partitioning of the bipartite graph. Each cluster correspond to a file on the directory;
2. a file with suffix *-Predictions.txt* that contains the semEP predictions for each cluster
and the probability of each prediction. 

## 10. Contact

I hope you find semEP an useful tool. Please, let me know
any comment, problem, bug, or suggestion.


[Guillermo Palma](https://gpalma.github.io/)

[palma at l3s dot de ](mailto:palma@l3s.de)




