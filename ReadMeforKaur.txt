##### ReadMe 

# Potentially useful Modules
reconstructCDS() -> add/remove alternative exonic segments into/from reference CDS and set an ORF coordinate

testNMD() -> Given a query Transcript and its CDS, test for classical (default) and/or non-classical NMD features

classifyAltSegments() -> Extend indentifyAddedRemovedRegions function by adding in information on splicing types

resizeTranscripts() -> Resize a GRanges object containing exon coordinates of a transcript

matching()
ORF()

# Secondary Modules
testTXforStart() -> tests whether a query transcript contain the same annotated start codon as from a reference CDS

reconstructCDSstart() -> for a given transcript that do not contain an annotated start codon, attempt to predict an in-frame start codon

# Workflow/Pipeline Modules
functions are found in main.R and workflow_functions.R