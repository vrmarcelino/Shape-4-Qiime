Shape-4-qiime
==============

Python scripts to adapt fastQ files from non-model sequencing designs (e.g. different amplicons or amplification strategies) to the Qiime pipeline.

Our PCR design consists of a two-step dual-indexing PCR protocol (similar to Nextera), which drastically reduces the costs of ordering multiple indexed primers. The following scripts are being used to demultiplex the samples, separate the amplicons, trimm the adaptors and etc...


NOTE (May 2016): We now use Heroen's perl pipeline to analyse the amplicon data. It is much faster as it spreads the job through multiple cores, and uses UPARSE to cluster OTUs.



