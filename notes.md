# Module 04 - Transmission Genomics notes

## 13 Jul 2020
### Video 01 - Part 1
#### New Data - A Big Opportunity

Infections evolve over time (e.g. TB and antibiotic resistance).

Sequencing revolution ~2007.
NextGen sequencing tools developed and cost of sequencing plummeted.

Acquired resistance vs. primary resistance (side bar).
"[A] Index case for two-person outbreak" (A infects B)

#### It's hard to control outbreaks.
We don't see exactly when infection events happen.
We don't know who's infecting who and when.
We can use data to understand & control transmission -- core task of transmission epidemiology.

microreact.org/project/west-african-ebola-epidemic?tt=rc

#### Sequencing less expensive now than ever.
Can now sequence thousands or tens of thousands of virus or bacteria samples in a single study.
Unprecedented, high resolution view of how organisms changing over time & space, from person to person in outbreaks.
Understand transmission in a better way by using standard techniques.

#### Bugs acquire small genetic variation as they spread.
nextflu.org
Downsample to be as globally representative as possible.
Shows phylogenetic structure.
Clade designation derived from phylogenetic structure.
As time goes on, clades grow & diversify because viruses pick up small amounts of genetic variation.

#### Bugs even vary from person to person.
e.g. the game of "telephone"; acttga -> actcga -> actcga -> actcgc -> actcac

#### Can we use sequences to understand transmission?
Can we learn about who infected whom and when from sequence data?

Why do we want to know who infected whom?
* we don't really need to know at the level of specific individual infections
* we do want to know when, where and how transmission takes place
* drilling into the details of transmission can help us understand this
* and it can help to know when infection _did not_ occur

#### Data on genetic variation can help
* help with who infected whom - improve outbreak control
* bigger scale: help w/ choosing best vaccines, best antibiotics

### Video 01 - Part 2
#### Challenges in genomic epidemiology & beyond.
Goal: try to understanda transmission -- "who infected whom" & how confident are we in that?

Outbreak questions:
* how fast do we have to find cases?
* how do we find missing cases? (case finding)
* what are early signs of a big outbreak? could be local. look for hidden/cryptic transmission

Large-scale questions
* how to choose what goes in a vaccine? (diversity)
* how to use antibiotics most wisely? (drug resistance -- acquired or transmitted?)

These two layers connect -- think about who infected whom, estimate it, and determine who has primary vs. required resistance in course of treatment.
Gap in understanding these questions & using "big fancy" sequence data of how pathogens evolve.

#### Relatedness is key
The answers to our questions aren't just in the data, but in the connections among data points.

The sequence `AACCATAGGT` doesn't mean much for transmission on its own.
But with two, `AACCATAGGT` and `GACCATAGGT`, we know we have two very similar things.
Patterns of relatedness used to understand transmission.

Why we can't just throw neural networks at the problem.
We have complex relationships.

#### Three roles for sequences in public health
1. Infer global routes of movement
* nextstrain.org, microreact.org
* challenges from different sampling and sequencing in different places (e.g. perhaps Europe sequenced more than other parts of the world)

2. Infer population dynamics back through time
* field of phylodynamcis - see Julia Palacios' module!
* limited in terms of _direct public health action_; more ecological

3. Analysis of transmission in localized outbreaks
* potentially directly useful for public health in the short term

#### Transmission in localized outbreaks
The first step is to broadly find out what you have
Clustering: put similar sequences into groups

#### How to put sequences into groups: clustering
A cluster: more similar to each other than to objects outside of the group
* SNP - single nucelotide polymorphism. e.g. change from A -> C at a single site
* Simplest clustering method: Sequences from A & B are placed in the same cluster if they differ by _k_ SNPs or fewer
* More sophisticated clustering methods:
    - Account for time, rate variation, selection. (e.g. our own 'transcluster' (Stimson et al MBE 2019 Beyond the SNP threshold))
    - Phylogenetic methods: cluster picker, cluster picker II, cluster matcher (developed for HIV primarily)
    - Integrated methods -- combine diagnosis, genotyping, and genetic similarity (Poon et al Lancet HIV)

#### Sequence cluster in public health 1
Exclude suspected transmission events:
1. for ex, 2 cases that are epidemiologically linked
2. if their viral or bacterial sequences are very different, there was liekly no direct transmission (they were infected somewhere else)
3. the apparent link was spurious ('false positive' epi link)

#### Sequence clusters in public health 2
Identify sources of infection that did not have epidemiological links
1. sequences firmly place B in a cluster w/ A & E (for ex, sequences are 1 SNP away)
2. B has no epidemiological links to A, E, or their contacts
3. this can help identify previously unknown exposures: true links were missing ('false negative' epi link)

#### Ex. COVID-19 genomic epidemiology in Australia
Shows epi curve in Victoria

#### SARS-CoV-2 data from Australia
Brief data description
* 1388 lab-confirmed cases in Victoria
* 62% travellers
* 27% known contacts
* 10% unknown source of exposure
* 1242 sequenced
* 1085 passed QC
* Max 15 SNPs compared to Wuhan 1

#### Clusters in Australian (Victoria) SARS-CoV-2 data
* authors used ClusterPicker to divide genomes into clusters
* 737 of 1085 were in any cluster
* 76 clusters: median 5, median duration 13 days. suggests good control (and could provide a serial interval estimate, too)
* 34 clusters were entirely overseas travellers
* 34 mixed, typically the first case was a traveller
* 81 seq with unknown exposure (from epi point of view) land in 24 clusters

Genomics provide much higher resolution picture over a histogram.

#### Special clusters
* epi clusters: groups of cases thought to be linked together on the basis of epidemiological (not genome) data, e.g. where people live, timing, health care, suspected exposure
* genomic clusters: similar genomes grouped together on the basis of (sort of) genetic distance
* four distinct epidemiological clustesr -> one gemoic cluster -> find links they didn't know about
* one big epi cluster separated into 4 distinct genomic clusters: exclude links they thought they knew about

Data required for this nice work: SARS-CoV-2 sequences, suspected exposure times and sources from epidemiology


### Video 01 - Part 3
#### Early work: transmission with sequences
* Cottam et al: Integrating genetic and epidemiological data to determine transmission pathways of foot-and-mouth disease virus. Proc. Biol. Sci. 2008
    - Phylogenetic trees constrain transmission events
* Jombart T, Eggo RM, Dodd PJ, Balloux F (2011) Reconstructing disease outbreaks from genetic data: a graph approach. Heredity 106: 383390
    - Parsimonious: transmission tree w/ fewest mutations

#### Early work cont'd
* 2012 Ypma RJ, Bataille AM, Stegeman A, Koch G, Wallinga J, et al. Unravelling transmission trees of infectious diseases by combining genetic and epi. data. Proc Biol Sci 279
* 2012 Morelli MJ, Thebaud G, Chadoeuf J, King DP, Haydon DT, et al. A Bayesian inference framework to reconstruct transmission trees using epi. and genetic data. PLoS Comput Biol 8
    - Both use a unified likelihood of genetic and epidemiological data - needs full sampling
* 2014 Outbreaker: Jombart T ...

#### Outbreaker and SARS
The data
* sequences
* times of sample collection
* number of cases
Additional inputs (Derived from data):
* generation time - time from infection to symptom onset
* time from infection to sample collection

#### The math behind Outbreaker
Outbreaker is a Bayesian MCMC method, and it uses _augmentation_.
"Augmented" quantities for each case _i_:
* the infector for each case, or most recently sampled infector
* the number of unknokwn intermediates b/w _i_ and _i_'s ancestor
* The date when _i_ was infected

The likelihood relies on these augmented data.

#### Outbreaker's likelihood in brief
**Likelihood for case** _i_ _alpha_

Model for how sequences change over time: P(sequence | ancestor's sequence, intermediates, evolution mode) x
Model for public health system: P(sampled time | infection time) x
Model for epidemiology of infection: P(infection time | parent's infection time, intermediates)

**A 60-second MCMC briefing**:
* propose augmented data: transmission tree, intermediates, times of infection
* compute the above likelihood
* accept or reject (in a principled way, depending on proposal method & likelihood)
* run for a long time. collect "posterior samples" of these quantities

#### Outbreaker results: simulated data

#### Application: SARS (The first SARS!)
* 13 SARS "genomes" but note, sequencing was very different then
* Ruan et al, Lancet, Volume 361, Issue 9371, 24 May 2003

#### Discussion
You'll hear more about the advent of sequencing - much more data now
Outbreaker is now in the broader package 'outbreaker2'
there was more in these papers not described
both these methods have limitations
* they do not capture the shared ancestry of the pathogens
* for that, we need phylo trees
* they do not accommodate variation w/in hosts

#### What's next?
* reconstructing transmission trees
* intro to genomics for genomic epi
* non-phylo outbreak reconstructions in outbreaker
* phylo trees: theory & practice
* transphylo - genomic epi w/ trees
* research forefronts: bringing in more data
* research forefronts: SARS-CoV-2 and COVID-19


### Video 02 - Part 1
#### Introduction to genomics for transmission inference
From genomes to phylogenies

Molecular info provides new ways to combat disesaes!
Why are they so important? What are they made of? Where are they?

What is a genome?
An organism's complete set of chromosomes (genetic info).
Humans have 23 pairs of chromosomes.

Section of DNA we call a gene.
Each genome contains all of the info needed to build and maintain an organism.

DNA passed from parent => descendants.

How can we look at mutations in DNA to infer tree structure?

#### DNA replication
Two strands are separated like two sides of a zipper
The enzyme DNA polymerase moves along the exposed DNA strand, joining newly arrived nucleotides into a new DNA strand.

DNA replication is very accurate but sometimes errors happen.
During DNA replication, an incorrect base may be added to the growing chain.
These changes are called mutations.

Other types of damage to DNA cause mutations.
e.g. exposure to radiation, carcinogens

#### Mutations
Changes in the DNA can affect one chromosome (inversion) or comprise multichromosomal events (translocation)

Overview of some chromosomal transolcations involved in different cancers and implicated in other conditions, e.g. szichophrenia

Mutations can be duplication, insertion, deletion (can be small sections of DNA)

Sometimes mutations can be beneficial for the organism that carries them, but not for the host

We need to update flu vaccines every year b/c the flu mutates so fast!

Mutations are the raw materials of evolution, e.g. adaptations to hot and cold climates.

Accumulated mutations in the DNA can help us combat diseases. How can we use it?

We need to read DNA.
How to sequence DNA? Chain termination method

1. Isolate DNA.
  A single stranded DNA fragment is isolated for which the base sequence is to be determined (the template).
2. Create "special" basis.
  Force DNA to replicate w/ DNA Polymerase.
  Provide source of free bases.
  Smaller set of some kind of special bases (colored w/ fluorescent dye).
  When they attach onto DNA fragment, stop replication process (cut DNA).
3. Put everything together: unknown DNA, normal bases, special bases.
4. Let the DNA be replicated.
  DNA polymerase stops replication when a special base is added.
5. Repeat replication a lot of times!
  Sequences always start at the same point, but they end in random points when a special base is added.
  If we do this enough times, a special base could be attached to each position in the DNA sequence at least once.
6. Separate and read syntehsized fragments (all of different lengths depending on where in sequence special base added).
  Lay them out in "conveyor belt" that sends them towards a laser.
  The color at the end of each fragment is detected by a laser beam.
  Shorter fragments move faster than longer ones.

The chain termination method takes too long when DNA sequences are big.
The human genome contains approximately 3 billion of these base pairs.

Shotgun sequencing
* DNA is broken up randomly into numerous small segments, which are sequenced using the chain termination method to obtain reads (fragments).
* Con: We have to draw enough random segments that we are confident entirety of genome has been sampled.

The last two decades have seen a revolution in genome sequencing.

Recommendation: Life, the science of biology, covers genome sequencing.


### Video 02 - Part 2
#### How many random fragments are enough to obtain the original DNA sequence?

G = length of genome
L = read length
N = number of reads

What is the minimum number of reads required to cover all the sequence?
```
LN >= G
```
As the fragments are randomly distributed, it is unlikely that the minimum will be enough.

Let's look at the coverage per nucleotide instead.
The average coverage per nucleotide is:
```
c = (LN)/G
```
What is the probability of a read with size L covers 1 nucleotide in the sequence of size G?

```
L = 13
G = 26
p = L/(G-L + 1) = 13/14 = 0.928
```

In more general terms...
What is the probability of N(4) reads covers 1 nucleotide exactly x (2) times?
Binomial distribution

```
P(x reads covering one base) = (N x)p^x(1 - p)^(N - x)
```

p is small: size of the read (L=100 bases) is very small compared to the size of the sequence (G=23 billion bases).
Small read ~100 bases long.

The binomial distribution converges towards the Poisson distribution when p tends to zero:

```
P(x reads covering one base) = (N x)p^x(1 - p)^(N - x)
P(x reads covering one base) = (e^-c * c^x)/(x!)
p = L/(G-L+1)
c = (LN)/G
```

Poisson distribution w/ rate = quantity c.
Avg. number of reads that cover particular nucleotide.

How mnay random fragments are enough to obtain the original DNA sequence (i.e. obtain coverage at every nucleotide)?
```
c = (LN)/G
```

Generally want to increase N until % of genome not covered is low (< 25%).

#### How to obtain the original DNA sequence from its fragments?
Graphs!
Brides of Konigsberg problem.
Ideas for solving our problem, as well as the entire branch of mathematics, known today as graph theory, can be traced back 300 years ago to this city.
(Konigsberg, presnet-day Kaliningrad, Russia)

At the time, Konigsberg's residents enjoyed strolling through their city, and they wondered if every part of the city could be visited by walking across each of the seven bridges exactly once and returning to one's starting location.
City has 4 parts separated by the Pregel River.

Leonhard Euler made a conceptual breakthrough.
Euler's first insight was to reprsent each landmass as a point (called a node) and each brdige as a line segment (called an edge) connecting the two appropriate points.
This creates a graph, a network of nodes connected by edges.

It was not possible; when it is possible, we call it Eulerian cycle (or path).

200 years later, in 1946, Dutch mathematician Nicolaas de Bruijn became interested in the superstring problem.
Shortest possible word to make that contains all words within. e.g. angelfish (angel, gel, elf, fish)

ex.
Alphabet: 0 and 1
Input: all possible 3-mers 000, 001, 010, 011, 100, 101, 110, 111
Output: 0001110100

De Bruijn answered this question by borrowing Euler's solution to the Bridge problem.
Every possible (k-1)-mer is assigned to a node
prefix, suffix
110
11 -> 10

Two nodes are connected by a directed edge if there is a word whose prefix is the former and whose suffix is the latter.

Find word that forms Eulerian circuit.

#### 50 years later...
Obtain original DNA sequence from its fragments with De Bruijn graphs.

#### Can you assemble your DNA with De Bruijn graphs?
Our alphabet: ATCG
Input: reads (k-mers)
ATCAG, TCAGT

Output: DNA sequence (Eulerian circuit)

#### Genome assembly software
SOAPdenovo, Velvet, ALLPATHS, ABySS


### Video 03 - Part 1
#### Non-phylogenetic transmission reconstruction

Epidemiological outbreak reconstruction
outbreak datat alone can be used for outbreak reconstruction, but genetic data offer a high-resolution source of extra info.

What can genomic data offer?
* extra detail
* resolve transmission where epi data are hard to get/have gaps
* genomic data now much easier, cheaper, and faster to get than ever before (real-time sequencing even becoming possible in the field)

# Challenge: create a single framework/likelihood incorporating genomic + epi data
Imagine we have 3 people infected in an outbreak... A, B, and C.

We want to combine our genomic info and our epi info to best narrow down which possible path the infection took.

Infeasible to test all transmission trees.

Early approaches included:
Maximum likelihood approach
* first, restrict to all transmission trees which are consistent with known infections
* then, use genomic data to constrain the set of possible transmission trees (parsimony)
* then, calculate the likelihood of each remaining tree based on the epi information -- e.g. the chance each individual (a farm) was infected on a given day/able to infect others on a given day (consensus tree)

Graph theory approach to find 'genetically parsimonious' transmission trees.
Algorithm SeqTrack finds the optimum branching in a directed graph.
* Create a distance matrix w/ genetic data (i.e. number of mutations between each strain)
* Create a connected, directed graph with weights w_ij equal to the genetic distance
* Remove edge _ij_ if t_j < t_i
* Find the spanning directed tree optimizing (i.e. minimizing) sum of weights along nodes, w_ij

Consider sample collection dates (t).
Assume anyone who was sampled before someone else can't have been infected by them.
Some limitations:
* all cases come from single index case e.g. a single sampled ancestor
* all cases are known and sampled

SeqTrack also:
* Assumes that individuals became infectious in the order they are sampled
* Has no uncertainty in the output transmission tree or probabilistic parameters
