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


### Video 03 - Part 2
#### A quick primer 1: generation time and sampling time

**generation time** = the time interval b/w the infection of an individual and their seeding of new secondary cases

**sampling time** = the time interval b/w infection and collection of an isolate

#### A quick primer 2: Markov Chain Monte Carlo (MCMC)
A popular method for exploring complex and/or high-dimensional spaces - e.g. transmission trees.

The main idea:
posterior - the probability of our model parameters given the data
likelihood
prior

When this quantity is hard to maximize directly, we instead form a Markov chain w/ equilibrium distribution equal to the posterior distribution, and then take many samples from this chain.

Essentially, we approx. the posterior distribution by random sampling from a probabilistic space (of all possible transmission trees).

Data-augmented MCMC is a method for dealing w/ missing data w/in an MCMC algorithm.
As well as sampling from the parameter space at each step of the Markov chain, we also sample values for the missing data.

In transmission inference, missing data might be the time of infection of the cases (since typically we only know sampling times) or the number of unsampled cases, for example.


### Video 03 - Part 3
#### outbreaker and outbreaker2

Creates a unified likelihood for genetic & epi data within a Bayesian framework which allows more estimation and greater flexibility.

outbreaker2 a more customizable version of outbreaker.

We're going to focus on the core outbreaker model today.

Data:
* N sampled cases, each w/ genetic sequenc s_i and time of sampling t_i

Quantities:
* d(s_i, s_j) = number of mutations (distance) b/w sequences i and j
* l(s_i, s_j) = number of nucleotide positions which can be compared i and j
* w = distribution of the generation time
* f = distribution of the sampling time

Augmented data:
* α_i = index of the most recent sampled ancestor of i
* κ_i = number of (unsampled) generations b/w i and α_i
* (T^inf)_i = date of infection of i

Parameters:
* mu = mutation rate, per site per generation of infection
* pi = proportion of unsampled cases
are estimated as well as the transmission tree

Posterior distribution:
...

All cases are assumed to be conditionally independent, given the identity of their most recent sampled ancestor, so the likelihood decomposes to:
...

The pseudo-likelihood is further decomposed into genetic and epidemiological components.
For each case i = 1, ..., N:
...
genetic part, epi part, constant (α_i)


#### Genetic part
The outbreaker genetic model assumes no w/in host genetic diversity (i.e. an individual can only be infected once), and so mutations are direct features of transmission events.
All transmission events are assumed independent, and the genetic pseudo-likelihood is very fast to compute.

Genetic pseudo-likelihood of case i = probability of observing genetic distance d(s_i, s_α_i) b/w sequence s_i and the ancestral sequence s_α_i w/ i and α_i separated by κ_i generations.

As a method designed for shorter timescale outbreaks, reverse mutations are considered negligible.

#### Epidemiological part
Time of sampling given time of infection * time of infection given knowledge of inferior * number of missing cases given rate of missing cases

probability of obtaining one 'success' (sampling a case) after κ_i-1 'failures' (unobserved cases), w/ probability of success pi


In total, genetic part * epi part = overall likelihood
That forms the core of the outbreaker model.
The likelihood expressions introduced in the previous slides are combined w/ priors for the mutation rate mu and proportion of unsampled cases pi.

mu is given a uniform prior on [0,1] - corresponding to an assumption of scarce prior info on this.
pi is given a beta distributed prior w/ parameters controlled by the user of outbreaker.
This is flexible prior which can ...

Authors also introduce a method for detecting imported cases -- i.e. cases that are not descended from another case in the outbreak.
In an initial step of the model, genetic outliers are detected relative to the other samples in the dataset.
A 'global influence' GI_i is calculated for each sampled case, defined as
...

where GPL is the genetic pseudo-likelihood.
This is calculated over the first few samples of the MCMC, say 50.


## 14 Jul 2020
### Video 01 - Part 1
#### Intro to Phylogenetics

**Phylogenetics** is the study of the evolutionary history and relationships among individuals or groups of organisms.
Phylogenetics goes back to the origin of evolution.

Darwin's classic tree example "I think"

There are many applications of pylogenetics:
* identifying pathogens
* classifications
* answering biological questions
* bioinformatics
* forensics

There are a lot of learning materials online. EBI: ebi.ac.uk/training/online/course/introduction-phylogenetics

#### What are phylogenetic trees?
**Phylogenetic tree** is a tree structure (with nodes and edges but no cycles) in which tips correspond to organisms and internal nodes correspond to inferred common ancestors.
**Root** usually an inferred oldest node, or node separating your study organisms from an "outgroup".
There may not bea root.

```r
load("demotree.Rdata")
library(ape)
plot(tr, edge.width=2, edge.col="grey", tip.col=c("grey","grey","grey","blue","blue","blue",rep("grey",6)))
```

If we are studying the blue organisms, we might want to root the tree this way (but we don't have to -- it depends on the question).

Same tree w/ a different root. -- structurally the same tree

The horizontal axis, branch lengths, are in units of _genetic distance_, typically something like the number of substitutions per site in molecular (sequence) data.
Sometimes we create a _timed_ phylogenetic tree.
Its branch lengths are in units of time.

The vertical axis -- NOTHING!

In the above tree, Organism 11 is as closely related to Organism 8 as Organism 10 is to Organism 9.

#### Terminology
**sister groups** - two groups descending from the same ancestor
**clade** - an ancestor and all of its descendents. This is also sometimse called a _monophyletic clade_.
**common ancestor** (of a set of tips) - a node that is an ancestor of all tips in the set
**most recent common ancestor** (MRCA) of a set of tips - the common ancestor of the st of tips that is farthest from the root

#### Lots of trees
Tree shapes -- for a 45-leaf tree, there are 10^33 possible shapes, even if only considering trees w/ clades of 10, 10, 12, and 13 nodes each.

So is it practical to find the truly best tree among all the possible trees?
Not really, but it works quite well.
There is usually some uncertainty.


### Video 01 - Part 2
#### Pylogenetic analysis
Stages
1. start with a question
2. identify a model and parameters that could answer the question
3. collect sequence data that would help to answer the question
4. identify the orthologous sequences
5. align sequences (multiple sequence alignment)
6. estimate tree and other parameters given the data and model
7. estimate the error associated w/ the tree and/or parameter estimates
8. does it answer your question?
Unlock new biological insight!

We'll start from step 6 in this lecture.
1-3 are specific to your research question.
4-5 are other areas of bioinformatics.

#### What data do we need?
Inherited material (DNA, RNA).
Molecular (sequence) data are usually the key for phylogenetic tree reconstruction.

Quick, inexpensive, reliable.
Sequences are highly specific and rich in info.

We could use morphological or phenotypic data.

The _characters_ are the entries of the aligned data.
In molecular data, each charatcer is a genetic site or locus.
Ex. suppose the multiple sequence alignment looks like this:
```fasta
>A130
atgaaaccct cgcgccctta
```

#### Distance- and character- based tree building
1. Distance based methods
Compute distances b/w pairs of sequences and find a tree (e.g. using clustering) that best describe these.
Ex. neighbor joining

In fact the distances are still based on characters, but hte characters aren't considered one at a time --
they are pooled together to create pairwise distances b/w all the sequences in the dataset.

2. Character-based methods
Consider sequences of DNA, RNA, morphological characters _one at a time_, and seek the tree with an optimal "score": likelihood or parsimony;
or use a Bayesian method to compute a posterior collection of phylo trees.

#### Can we get the best tree?

Either way, we need to understand how to compare sequences to each other.
We need to either compute a distance b/w all the pairs of sequences or to figure out a likelihood or score for the tree.

To do these comparisons, we use models that describe how sequences evolve.

#### Models of sequence evolution
What is this and how does it describe distance?

How can we compute the distance b/w two sequences?
A common approach: the # of differences per site.
E.g. for the two sequences `AACCTATCGG` and `AATCTAACGG` differ in two places, so they would be 2/10 = 0.2 substitutions per site apart.

In other words, ifd the # of differences is _m_ and the total number of sites is _n_, then this basic distance is
```
d_s = m/n
```

However: what if there were actually 2 changes at the same site?
This simple computation would actually underestimate the true distance.

There are many models for sequence evolution.
They define the probability that a character at a site wil evolve into a diff. character at that site over a specified time period.
They can take various complexities into account.
Models of sequence evolution can be used to define distances b/w two sequences, and they can be used to simulate evolving sequences.

A,C,T,G: adenine, cytosine, guanine, thymine

A and G: 2-ring structures, purine
C and T: 1-ring structures, pyrimidine

Transition: a purine-purine change, or a pyrimidine-pyrimidine change (C-T or A-G)
Transversion: a pyrimidine-purine change (or vice versa)

The Jukes-Cantor model treats these as having the same probability.
But it's reasonalbe to think that transitions might be more probable than transversions.
Also, A, C, T, and G might not all be equally frequent.

#### Assumptions and parameters in several common methods
The Jukes-Cantor model and other models take this into account.
In the JC model, any possible change is equally likely.
It is the simplest model of sequence evolution.

Conceptually the link is this: The JC model (and higher-complexity variations) define the probability that one sequence evolves into another.
The "distance" is related to the amt. of time that this would take in the model.

The Jukes-Cantor model is the simplest, bu they are all based on homogeneous continuous-time Markov chains.
These are beyond teh scope of this course.
The wiki page on these models is quite good as a quick reference: en.wikipedia.org/wiki/Models_of_DNA_evolution.
"Markov Processes for Stochastic Modeling" by Masaaki Kijima is a math reference.

If the substitution rate is _mu_ substitutions per time, then the mean number of substitutions in time _t_ is _mu*t_.
The CTMC has a rate matrix _Q_ whose off-diagonal entries are the instantaneous rates of transition A -> C,
A -> T, C -> G, etc. and these are all equal to _mu_/4.
The transition matrix is P(t) = e^(Qt).

Consider the probability that a change from state i to state j is observed in time t, if the substitution rate is mu.
We have a probability e^(-mu*t) that no event happens in time t.
The rate of leaving each state is 3*mu/4 since 3 of the 4 "changes" (A to A, C, T, G) arrive at a different base.
So the expected number of changes in time t is
```
v = 3/(4*mu*t)
```

We can write `mu * t = 4/3v`.

The probability of any change at the site is 3 times this, becausethere are 3 choices of base that are different from the current one.
```
P(change) = 3/4 - (3/4)*e^(mu * t)
```
People use this relationship to relate a distance b/w two sequences (v := d) to the empirical fraction of the sequence's total length where there is a change b/w the two.

The Jukes-Cantor model's distance is
...

where d_s is the same as above (# of differences / sequence length).

JC: ACTG all equally frequent and all changes equally likely.
K80: (Kimura 1980) ACTG all equally frequent, but transitions more likely than transversions
F81: as with JC but ACTG are not all equally frequent
HKY85: Transitions and transversions are not equally likely and ACTG are not all equally frequent

In pratice, these models are implemented w/in software packages used to build trees.

#### Building phylo trees: the central task of phylogenetics
Neighbor joining
- distances are defined using sequence evolution models
- idea: join closest tips, make new node, repeat
- details: differing clustering algorithms
- single linkage, complete linkage, UPGMA
- software: bionj, nj in R
- ref: Saitou and Nei 1987

The neighbor joining algorithm:
Intuition: start w/ all the tips separate.
Create a new node that joins the two that are closest to each other.
Then define the distances b/w these two tips and the new node.
Then define the distances b/w the new node and all the other tips.
Forget about those two tips, but keep the new node.
Repeat until everything is joined up.

NJ strengths: can use realistic substitution model;
computationall efficient: no need to compare lots of trees for an optimality condition;
great for large datasets, short distances

Tree-like distances: if the distances are perfect (and reflect the true tree), the tree is guaranteed to be correct.

Weaknesses: if the distances aren't perfect then the tree may well be wrong.
In some circumstances the branch lengths can be negative!
The method is sensitive to gaps in the sequences.
The model does not allow for diff sites to have diff "matches" to the tree -- they all just contribute to the distance.
Distance measures may be inaccurate for large distances.
