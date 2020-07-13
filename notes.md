# Module 04 - Transmission Genomics notes

## 13 Jul 2020
### Video 01
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

### Video 02
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


### Video 03
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
