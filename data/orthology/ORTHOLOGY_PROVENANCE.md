# Orthology input

## Input used by this analysis

The cross-species analysis uses the following pairwise OrthoFinder table:

`Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv`

- File size: 1,023,678 bytes
- SHA-256: `512758AF98B6C20295B47AD59E04858918748874CBBB2DA2DC34D8CA4D9A3BF7`
- *Acomys* annotation: Ensembl Rapid Release for *Acomys dimidiatus*, assembly GCA_907164435.1
- *Mus* annotation: Ensembl release 102 for *Mus musculus*, assembly GRCm38

The manuscript analysis removes Ensembl version suffixes, enumerates *Acomys-Mus* gene pairs, counts the number of partners for every gene, and retains only pairs in which each *Acomys* gene has one *Mus* partner and each *Mus* gene has one *Acomys* partner. The count matrices are then restricted to pairs represented in both species and filtered to genes with a count of at least 10 in at least 3 of the 18 samples.

The supplied table was generated with OrthoFinder (2.5.4) using primary protein transcripts generated from the Ensembl .pep annotations processed with Orthofinder's `primary_transcript.py`.

```
orthofinder.py -f primary_transcripts -t 32 -a 32 -M msa
```

## Acknowledgement

We thank Huayun Chen for originally generating and providing the *Acomys dimidiatus–Mus musculus* OrthoFinder orthology table used to construct the cross-species gene universe.
