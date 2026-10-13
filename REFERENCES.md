# References and data provenance

## Sequence data

All sequences retrieved from NCBI on 13 October 2026 by `data/fetch.sh`.

### Targets — desert locust, *Schistocerca gregaria*

| Accession | Gene | Role in the design | Length |
| --- | --- | --- | --- |
| [XM_049990443.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_049990443.1) | vitellogenin-1-like, transcript variant X1 | reproduction | 5,540 nt |
| [HQ634702.1](https://www.ncbi.nlm.nih.gov/nuccore/HQ634702.1) | juvenile hormone acid methyltransferase (JHAMT), complete cds | development | 927 nt |
| [XM_049985731.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_049985731.1) | V-type proton ATPase catalytic subunit A, variant X1 | digestion | 3,169 nt |
| [XM_049999838.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_049999838.1) | 5-hydroxytryptamine receptor 2A, variant X1 | swarming behaviour | 3,275 nt |

### Non-targets — RefSeq RNA (transcriptomes)

| Species | Assembly | Sequence screened |
| --- | --- | --- |
| *Apis mellifera* (honey bee) | GCF_003254395.2 Amel_HAv3.1 | 98,645,700 nt across 27,887 transcripts |
| *Danaus plexippus* (monarch butterfly) | GCF_018135715.1 MEX_DaPlex | 64,638,557 nt across 23,381 transcripts |
| *Homo sapiens* (human) | GCF_000001405.40 GRCh38.p14 | 718,571,139 nt across 186,185 transcripts |

## Literature

1. Mogren CL, Lundgren JG. In silico identification of off-target pesticidal dsRNA binding in honey bees (*Apis mellifera*). PeerJ. 2017;5:e4131. <https://doi.org/10.7717/peerj.4131>
2. Zhang H, Li HC, Miao XX. Feasibility, limitation and possible solutions of RNAi-based technology for insect pest control. Insect Science. 2013;20(1):15-30. <https://doi.org/10.1111/j.1744-7917.2012.01513.x>
3. Fire A, Xu S, Montgomery MK, Kostas SA, Driver SE, Mello CC. Potent and specific genetic interference by double-stranded RNA in *Caenorhabditis elegans*. Nature. 1998;391(6669):806-11. <https://doi.org/10.1038/35888>
4. Anstey ML, Rogers SM, Ott SR, Burrows M, Simpson SJ. Serotonin mediates behavioral gregarization underlying swarm formation in desert locusts. Science. 2009;323(5914):627-30. <https://doi.org/10.1126/science.1165939>
5. United States Environmental Protection Agency. EPA registers novel pesticide technology for potato crops. 22 December 2023. <https://www.epa.gov/pesticides/epa-registers-novel-pesticide-technology-potato-crops>

## Method

6. Aho AV, Corasick MJ. Efficient string matching: an aid to bibliographic search. Communications of the ACM. 1975;18(6):333-40. <https://doi.org/10.1145/360825.360855>

Implemented via the `pyahocorasick` package.

## Thresholds

A shared exact match of 21 nucleotides or more is the criterion used in regulatory
assessment of RNA-based pesticides, and matches reference 1's framing of the
question. The 16-nucleotide tier is reported alongside it as a control: in a
transcriptome of this size, matches of that length are expected by chance, and the
observed counts fall at or below that expectation, which is what a threshold
dominated by random matching looks like.

## Note on authorship

The analysis, choice of species, thresholds and interpretation are Abeer Israque's.
The code was written with AI assistance.
