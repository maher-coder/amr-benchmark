---
title: "Reading bacterial genomes to predict antibiotic resistance"
subtitle: "A clinical companion to a cross-cohort study"
author:
  - name: Maher el Ouahabi
email: "maher.elouk2@gmail.com"
date: "April 2026"
abstract: |
  This document is a clinical companion to a technical report on machine-learning prediction of antimicrobial resistance from bacterial genome assemblies. It presents the biological findings and their clinical consequences for readers whose expertise lies in microbiology, infectious disease, or clinical care rather than in machine learning. The statistical machinery of the main study is omitted here; the conclusions that matter for the clinic are not.
toc: false
---

## Why this question matters

Antibiotic resistance is the problem that, over the next few decades, may come to matter as much as cancer. The World Health Organization has estimated that by 2050 drug-resistant infections could cause more deaths per year than any single cancer. Modern medicine's routine operations — hip replacements, chemotherapy, caesarean sections — all rest on a quiet assumption that we can still kill bacteria when we need to. At the same time, a revolution has happened in the laboratory: sequencing the complete DNA of a bacterium, which cost hundreds of thousands of dollars twenty years ago, now costs tens of dollars and takes a few hours.

This raises a simple question with a surprisingly complicated answer. Can we read the DNA of a bacterium and predict whether a given antibiotic will work? The short answer is: mostly yes, with caveats that are more interesting than the prediction itself. The technical companion to this document reports a benchmark on 170 000 bacterial genomes across six important pathogen–drug combinations. Here I want to explain, without machine-learning jargon, what we found and why three of those findings matter more for the clinic than the prediction accuracies do.

## What "resistant" actually means

It is tempting to assume the word *resistant* describes a fixed fact about a bacterium: either the drug kills it or it does not. The truth is both more pragmatic and more fragile. When a clinical microbiology lab tests whether, say, an *Escherichia coli* isolate is resistant to ciprofloxacin, it measures the smallest concentration of the drug that stops the bacteria from growing — a number called the Minimum Inhibitory Concentration, or MIC. That number is a real physical measurement. But whether the lab then reports the isolate as *resistant* or *susceptible* depends on comparing the MIC against a threshold — a *breakpoint* — set by an expert committee. In the United States this committee is the Clinical and Laboratory Standards Institute (CLSI). In Europe it is the European Committee on Antimicrobial Susceptibility Testing (EUCAST). Both committees revise their breakpoints every year, sometimes substantially, as new evidence arrives about how drugs behave in patients.

The consequence is that the label *resistant* is partly a biological fact about the bacterium and partly a moving clinical definition. In the CABBAGE database we use in this study — a collection of roughly 170 000 bacterial genomes curated by the European Molecular Biology Laboratory's European Bioinformatics Institute — applying the 2025 EUCAST breakpoints to the raw historical MIC measurements changes the resistant/susceptible label for about 1.8 % of records, or 31 306 samples, compared to the label that was originally published alongside the genome. Nothing changed about the bacteria. What changed was the definition.

## What we did

We took those 170 000 bacterial genomes from CABBAGE and ran each one through a standard tool that scans DNA for a catalogue of roughly 400 known resistance genes and mutations. This gave us a large binary table: one row per isolate, one column per resistance determinant, with a 1 if the gene or mutation is present and a 0 if it is absent. We trained two standard statistical classifiers — ElasticNet, a regularised logistic regression, and LightGBM, a gradient-boosted decision-tree ensemble — to predict whether each isolate was resistant or susceptible for six pathogen–drug combinations of real clinical importance. These included methicillin resistance in *Staphylococcus aureus* (the bacterium known to clinicians as MRSA), vancomycin resistance in *Enterococcus faecium* (VRE), ciprofloxacin resistance in *Neisseria gonorrhoeae*, and three further combinations covering ampicillin and tetracycline resistance in *Salmonella enterica* and *Escherichia coli*. Crucially, we reserved one entire data source — 14 251 genomes from a cohort independent of the training data — and did not touch it until the final evaluation. That out-of-cohort test is the closest available proxy for what happens when a model trained in one laboratory is deployed on samples from another.

## Finding A: Where the bacterium was isolated matters enormously

The first and most clinically consequential result is not a property of the algorithm at all; it is a property of the data. The same species, tested against the same drug, shows resistance rates that differ by up to 91 percentage points across different data sources. *Escherichia coli* isolates from the United States Centers for Disease Control outbreak-tracking database are 94 % resistant to ampicillin. *Escherichia coli* isolates from the National Antimicrobial Resistance Monitoring System (NARMS), which samples the food supply and agricultural animals, are only 3.6 % resistant to the same drug.

This is not an artefact of sloppy data curation. It reflects ecology. The CDC collection is enriched for outbreak strains precisely because outbreak strains are, almost by definition, difficult to treat — that is usually what draws epidemiological attention to them in the first place. Food-animal isolates face a very different antibiotic selection pressure, dominated by veterinary use patterns that do not necessarily mirror clinical use in people. A useful predictive model has to know which cohort a sample came from, because the cohort label already encodes, before any genetics is considered, most of what an informed microbiologist would guess about the probability of resistance. Any published resistance rate that does not identify the collection it came from is close to uninterpretable.

## Finding B: Resistance is nearly monogenic for most well-studied combinations

The second finding surprised me, because I had expected the 402 features in our model to matter more than they did. For *Neisseria gonorrhoeae* tested against ciprofloxacin, a single point mutation in the gene *gyrA* — specifically the substitution S91F — accounts for roughly 42 percentage points of the prediction accuracy on its own. For ampicillin resistance in *Escherichia coli*, four genes carry almost all of the signal: *blaTEM-1*, *blaCMY-2*, and two variants of *blaCTX-M*. These are β-lactamase enzymes, proteins that the bacterium secretes to break the antibiotic apart before it can reach its target. For tetracycline resistance in *Salmonella*, two efflux-pump genes — *tet(A)* and *tet(B)* — do most of the work. They encode molecular pumps that push the drug back out of the bacterial cell faster than it can get in.

The practical implication is striking. A machine-learning model nominally built from 402 features is, in practice, a one- to four-gene classifier with a small residual. For these well-studied combinations, the microbiology community has already done the hard work over the past four decades of identifying the genes that matter. The clinical consequence is important: rapid point-of-care diagnostic tests that detect one, two, or four specific genes are scientifically justified for these pathogens. They are not a crude approximation of a sophisticated algorithm; they are nearly the algorithm.

## Finding C: The labels we try to predict change over time

The third finding recalls the distinction drawn earlier between measurement and definition. When we took the 2015 EUCAST breakpoints and compared the resistant/susceptible calls they would have assigned to the 2025 EUCAST calls for the same MIC measurements, roughly 1.8 % of calls flipped. Most of the flips came from breakpoint values moving by one or two doubling dilutions. The underlying biology did not change; the cutoff did.

This matters because it means any benchmark that trains a model on labels generated under old breakpoints and tests on labels generated under new breakpoints is mixing two different things: genuine biology, which is what we want to learn, and clinical-definition drift, which is noise from the model's point of view. A bacterium that was called *susceptible* in a 2018 publication may legitimately be called *resistant* in a 2025 publication without having acquired a single new resistance determinant. Any researcher comparing AMR prevalence across decades, and any AI system being validated against historical labels, has to account for this.

## What this means for machine learning

Put these three findings together and a picture emerges of where machine learning in bacterial AMR genomics actually adds value. For most of the six combinations in this study, a simple handwritten rule — *if* mecA *is present, call MRSA; otherwise call susceptible* — achieves an area under the receiver-operating-characteristic curve that is statistically indistinguishable from what our 402-feature ElasticNet model achieves. This is not a failure of machine learning. It is a demonstration that four decades of careful microbiology have already compressed the signal into a small catalogue of known determinants. Machine learning earns its keep in the minority of combinations where resistance is genuinely polygenic, such as ampicillin resistance in *Salmonella enterica*, where co-carriage of several genes on the same plasmid matters more than any single gene. For the rest, the ceiling is set by the quality of the gene catalogue, not by the sophistication of the classifier.

## Implications for the clinic and for public health

Several practical conclusions follow. First, rapid diagnostic tests built around small panels of canonical resistance genes are scientifically justified for the well-studied pathogen–drug combinations in this study. A clinician waiting on a blood culture does not need a neural network; a well-designed PCR assay targeting the right three or four genes will do almost as well. Second, any commercial or academic AI-driven resistance-prediction product must report its performance across multiple, ecologically distinct data sources, not merely within the cohort in which it was trained. A model that works well on outbreak isolates may be almost useless on routine surveillance samples, and vice versa. Third, published resistance prevalence figures systematically overestimate true community prevalence when they come from outbreak-enriched or publication-biased sources and underestimate it in agricultural surveillance. Policy decisions drawn from a single source are fragile. Fourth, when committees such as CLSI or EUCAST update breakpoints, the update should ideally be accompanied by a published quantification of how many retrospective classifications would flip under the new rules. This is relatively easy to compute and would help clinicians, epidemiologists, and researchers avoid false trend lines.

## What we did not do

This study has honest limits. Our pipeline cannot currently predict resistance in *Mycobacterium tuberculosis* because the gene-annotation tool we use does not capture the characteristic point mutations in *rpoB*, *katG*, *embB*, and *pncA* that drive tuberculosis resistance. This is a coverage gap in the bioinformatic tool we rely on, not a fundamental limit of genome-based prediction for TB. Closing it requires a specialised tuberculosis-aware annotator, and the community has several of those available. We also did not compute sequence-type information from multilocus sequence typing (MLST) or plasmid-replicon profiles, either of which would likely close the residual ten-percentage-point gap we see on methicillin resistance in *Staphylococcus aureus* and vancomycin resistance in *Enterococcus faecium*. That work is the natural next step.

## Closing

The punchline of this study is that bacterial antimicrobial resistance, for the well-studied pathogens that account for most clinical decisions, is a surprisingly parsimonious prediction problem. A small number of canonical genes, one or two well-characterised point mutations, and the clinical threshold of the year together explain the great majority of the signal. Machine learning plays an honest role, but it also reveals how little algorithmic complexity is required. The frontier of useful progress in this field is not bigger models. It is better features — full genome sequence, plasmid structure, MLST background — and clearer labels, meaning consistent and well-documented breakpoint annotation. Those are the unglamorous improvements that would move clinical AMR prediction forward.

## References {-}

[1] Jonathan Dickens et al. A comprehensive AMR genotype–phenotype database (CABBAGE). *bioRxiv*, 2025.11.12.688105, 2026.

[2] World Health Organization. WHO bacterial priority pathogens list 2024: bacterial pathogens of public health importance to guide research, development, and strategies to prevent and control antimicrobial resistance. WHO, Geneva, 2024.

[3] European Committee on Antimicrobial Susceptibility Testing. Breakpoint tables for interpretation of MICs and zone diameters, version 15.0. EUCAST, 2025.

[4] Clinical and Laboratory Standards Institute. Performance standards for antimicrobial susceptibility testing, 35th edition (M100). CLSI, 2025.

**Code and data.** <https://github.com/maher-coder/amr-benchmark>.
