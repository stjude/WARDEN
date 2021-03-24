<!-- dx-header -->
This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->




# Instructions :

1. For the FASTQ_FILES input filed:

    Select the field and identify all gzipped fastq files (*.fq.gz or *.fastq.gz) to be used.  There is a a limit of 20 samples.

    There should be at least 2 samples for each phenotype/condition tested

1. For the sample_list input field:
Identify sample sheet that adheres to the following format:

    The sample sheet is white space delimited (i.e. tab-delimited or space-delimited).  There should be a header line with at least columns.  The column names do not matter, but to be consistent the following can be used:
    SampleName	Phenotype	ReadFile1	ReadFile2
    Additional lines contain the above information for each sample.  For single-end reads use an '-'
    For any differentional expression analysis between groups, include a line starting with #COMPARISON=  or #COMPARISONS=   (either upper or lower case is fine).  Then provide a comma limited list with two phenotypes being compared, separated by a minus sign.  Additional columns will be ignored.
    Sample names and phenotypes should be alphanumeric with no spaces.  Dashes and underscores are ok.
    Example sample sheets:


```
#COMPARISONS=SJBALL101-SJPHALL005
SampleName	Phenotype	ReadFile1	ReadFile2
SJBALL101_1	SJBALL101	SJBALL101_10K_R1.fastq.gz	SJBALL101_10K_R2.fastq.gz
SJBALL101_2	SJBALL101	SJBALL101_10K_2_R1.fastq.gz	SJBALL101_10K_2_R2.fastq.gz
```

1. Select the 'WARDEN WORKFLOW GENERATOR" input box(Which brings you to the current menu)

    1. COMMON parameters

        1.  Enter an 'Output Folder'.  A folder must be specified here or from the 'Workflow Actions' pulldown menu from the workflow page.   The folder selected should be empty.  The folder does not need to currently exist.  (If selecting from the WORKFLOW menu, a folder can be created in the process)

        1.  Enter an analysisName.  This is required and the workflow created will be labelled with this name

        1.  Select 'sequencingStrandedness'.  This must not be an option starting with '***'. Determine the strandedness from your sequencing facility or source of sequence files.  Alternatively one can run a single sample through the pipeline with the three options and determine the option with the best alignment scores.

        1.  Select 'Genome'.  Choose any value not starting with '***'.  

        1.  allFASTQSInSampleSheet.  Not implemented.  For non-workflow mode, all FASTQ files provided should be in the sample sheet,  and all fastq files in the sample sheet should be provided.

    1.  LIMMA parameters.  Leave as is unless you are an advanced user. 
        1. calcNormFactorsMethod.  (See https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors).  

        1. filterCountType  rawCounts: Total number of counts for a gene across all samples must be greater than or equal to this number.  For CPM  the following is used (where filterCount is the cutoff):

        1. filterCount.  The cutoff used for filterCountType
        1. pValueAdjust.  Multiple testing adjustment for p-value.
           sSee http://svitsrv25.epfl.ch/R-doc/library/stats/html/p.adjust.html


