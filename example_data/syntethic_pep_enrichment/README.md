# Sythethic peptide erichment analysis example

- mock differential analysis on peptide level
- pathway annotations on protein level, i.e. more than one peptide matches to a protein

## Stragegy to match

- create a unique peptide_protein identifier (peptides are unique themselves, but their 
association to a protein is highlighted)
- expand annotations for each protein to all it's peptides using the peptide_protein identifier
- perform per pathway enrichtment to all peptide_protein matching to that pathway in 
  the foreground and background.
