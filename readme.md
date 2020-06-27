# Ontoclass
Ontology-based transporter substrate annotation for benchmark datasets

The tool takes UniProt identifiers of proteins as input. For each protein, it determines whether the protein has transporter-related GO MF annotations in the Swiss-Prot database. If a transporter-related GO MF is found, it looks for the ChEBI identifier of the transported substrates in the GO annotation. This ChEBI-ID and its ancestors in the ChEBI ontology are used to find the most specific substrate class according to a predetermined list of substrate classes and ChEBI-IDs. The tool outputs the final substrate class of each protein along with additional information. Details regarding the steps are available in:


<a id="1">[1]</a> 
Alballa, Munira, and Gregory Butler (2019).
Ontology-based transporter substrate annotation for benchmark datasets.
IEEE International Conference on Bioinformatics and Biomedicine (BIBM).
