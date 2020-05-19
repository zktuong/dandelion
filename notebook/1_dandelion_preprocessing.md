# *dandelion* Notebook-1

![dandelion_logo](img/dandelion_logo.png)
## Foreword
***dandelion*** is written in `python==3.7.6` and it is primarily a single-cell BCR-seq analysis package. It makes use of some tools from the fantastic [*immcantation suite*](https://immcantation.readthedocs.io/) and the main idea is that it implements a workflow for the pre-processing and exploratory stages with integrated use of tools from *immcantation* for the BCR side of things and analysis tools from [*scanpy*](https://scanpy.readthedocs.io/) for the RNA-seq side of things. I hope to be able to introduce some new single-cell BCR-seq exploratory tools down the road through *dandelion*. 


## Pre-processing - Part 1
This notebook will cover the initial pre-processing of files after 10X's `cellranger vdj` immune profiling data analysis pipeline. The directory structure of a 10x output folder will typically look like this:
```console
(dandelion) mib113557i:Pan_Immune_BCR kt16$ tree Pan_T7918901
Pan_T7918901
├── all_contig.fasta
├── all_contig.fasta.fai
├── clonotypes.csv
├── concat_ref.fasta
├── concat_ref.fasta.fai
├── consensus.fasta
├── consensus.fasta.fai
├── consensus_annotations.csv
├── filtered_contig.fasta
├── filtered_contig_annotations.csv
└── metrics_summary.csv
```

At this stage, ***dandelion*** only needs the fasta files to start, particularly either *all_contig.fasta* or *filtered_contig.fasta*.

In this notebook, I'm running everything with the *all_contig.fasta* files to get a sense of how long it would take when the files are considerably larger. I'm using a standard laptop for the analysis here: entry level 2017 Macbook Pro with 2.3 GHz Intel Core i5 processor and 16 GB 2133 MHz LPDDR3 ram.

#### Before starting, a couple of environmental variables need to be set up to make it run smoothly:

In ***shell***, export the path to the database folders like as follows:
```bash
echo "export GERMLINE=/Users/kt16/Documents/Github/dandelion/database/germlines/" >> ~/.bash_profile
echo "export IGDATA=/Users/kt16/Documents/Github/dandelion/database/igblast/" >> ~/.bash_profile
echo "export BLASTDB=/Users/kt16/Documents/Github/dandelion/database/blast/" >> ~/.bash_profile
echo "export PATH=/Users/kt16/Documents/Github/dandelion/bin:$PATH" >> ~/.bash_profile
source ~/.bash_profile
```
The databases for igblast are basically setup using [changeo's instructions](https://changeo.readthedocs.io/en/stable/examples/igblast.html). The instruction for setting up blast database is simpler and will be covered later in this notebook.

Also check if the softwares can be found:
```console
(mypython3) mib113557i:~ kt16$ conda activate dandelion
(dandelion) mib113557i:~ kt16$ which makeblastdb
/Users/kt16/Documents/Github/dandelion/bin/makeblastdb
(dandelion) mib113557i:~ kt16$ which blastn
/Users/kt16/Documents/Github/dandelion/bin/blastn
(dandelion) mib113557i:~ kt16$ which igblastn
/Users/kt16/Documents/Github/dandelion/bin/igblastn
(dandelion) mib113557i:~ kt16$ which tigger-genotype.R
/Users/kt16/Documents/Github/dandelion/bin/tigger-genotype.R
```

If you don't have the softwares, download [blast+](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/) and [igblast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). For tigger-genotype, you can download it [here](https://bitbucket.org/kleinstein/immcantation/src/default/pipelines/). Just note that I made some minor modifications to this file, hence there's version specific to this package.


```python
# import modules
import os
os.chdir(os.path.expanduser('/Users/kt16/Documents/Github/dandelion'))
import dandelion as ddl
```

    /Users/kt16/miniconda3/envs/dandelion/lib/python3.7/site-packages/anndata/_core/anndata.py:21: FutureWarning: pandas.core.index is deprecated and will be removed in a future version.  The public classes are available in the top-level namespace.
      from pandas.core.index import RangeIndex



```python
# change directory to somewhere more workable
os.chdir(os.path.expanduser('/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/PIP/Pan_Immune_BCR/'))
# print current working directory
os.getcwd()
```




    '/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/PIP/Pan_Immune_BCR'



### Step 1:
#### Formatting the headers of the cellranger fasta file
This step immediately below is optional and is just a lazy way to make a dictionary from an external file using a utility function `utl.dict_from_table`.


```python
# prepare a dictionary from a meta data file.
sampledict = ddl.utl.dict_from_table('/Users/kt16/Documents/Clatworthy_scRNAseq/Ondrej/dandelion_files/meta/PIP_sampleInfo_kt16.txt', columns = ('SANGER SAMPLE ID', 'GEX_SAMPLE_ID')) # optional
```

I'm adding a sample prefix to the headers of each contig in the fasta files using the dictionary created above, via the function `pp.format_fasta(s)`. The prefix is basically just the folder name, so in this case it's `Pan_T7918901`. The function will also create subfolders where the new fasta file and all subsequent files will be located. The file structure should look something like this later on if the settings are left as default.
```console
(dandelion) mib113557i:Pan_Immune_BCR kt16$ tree Pan_T7918901
Pan_T7918901
├── all_contig.fasta
├── all_contig.fasta.fai
├── clonotypes.csv
├── concat_ref.fasta
├── concat_ref.fasta.fai
├── consensus.fasta
├── consensus.fasta.fai
├── consensus_annotations.csv
├── dandelion
│   └── data
│       ├── all_contig.fasta
│       ├── all_contig_igblast.tsv
│       ├── all_contig_igblast_db-pass.tsv
│       └── tmp
│           ├── all_contig_igblast.blastsummary.txt
│           ├── all_contig_igblast.fmt7
│           └── all_contig_igblast.xml
├── filtered_contig.fasta
├── filtered_contig_annotations.csv
└── metrics_summary.csv
```


```python
# the first option is a list of fasta files to format and the second option is the prefix to add to each file.
samples = ['Pan_T7918901', 'Pan_T7918902', 'Pan_T7918903', 'Pan_T7918904', 'Pan_T7918905', 'Pan_T7918906', 'Pan_T7918907', 'Pan_T7918908', 'Pan_T7918909', 'Pan_T7918910', 'Pan_T7918912', 'Pan_T7918913', 'Pan_T7918914']
ddl.pp.format_fastas([str(s)+'/all_contig.fasta' for s in samples], [sampledict[s] for s in samples])
```

    Formating fasta(s) : 100%|██████████| 13/13 [00:38<00:00,  2.92s/it]


The function above is just a wrapper for a for-loop:
```python
for s in samples:
    filePath = s+'/all_contig.fasta'
    ddl.pp.format_fasta(filePath, sampledict[s])
```

### Step 2:
#### Reannotate the V/D/J genes with *igblastn*.

`pp.reannotate_genes` uses [*changeo*](https://changeo.readthedocs.io/en/stable/examples/10x.html)'s scripts to call *igblastn* to reannotate the fasta files. Depending on the file format option, it will parse out as either an `airr` (default) or `changeo`-legacy TSV file. Importantly, with the recent update to changeo v1.0.0, all the column headers are now adhereing to the [*AIRR*](http://docs.airr-community.org/) standard (lowercase and some column name changes).


```python
# reannotate the vdj genes with igblastn and parses output to 'airr' (default) or 'changeo' tsv formats using changeo v1.0.0 scripts
ddl.pp.reannotate_genes(samples)
```

    Assigning genes : 100%|██████████| 13/13 [27:14<00:00, 125.76s/it]


Because the suffixes at the end of the files are different for the different output format, the output files shouldn't overwrite each other.

But anyway, either format should work for subsequent steps.


```python
# to write as changeo format
ddl.pp.reannotate_genes(samples, fileformat = 'changeo')
```

    Assigning genes : 100%|██████████| 13/13 [25:10<00:00, 116.19s/it]


### Step 3:
#### Assigning constant region calls

10x's \*annotation.csv file provides a *c_gene* column, but it is supposedly not accurate according to [hk6](https://twitter.com/hamish_king); It is also known that their v/d/j annotations are inaccurate too, hence above reannotation with *igblastn*. However, *igblastn* doesn't return the call for the constant genes.

Rather than simply relying on 10x's annotation, [hk6](https://twitter.com/hamish_king) reccomended using [*immcantation-presto*'s *MaskPrimers.py*](https://presto.readthedocs.io/en/version-0.5.3---license-change/tools/MaskPrimers.html) with his custom primer list and I tested that; worked well but it took ***20 min*** for the first file (~6k contigs). It also only calls the constant region for the heavy chains. 

So, I wrote a function, `pp.assign_isotype`, to use *blast* to annotate constant region calls for all contigs and retrieves the call, merging it with the tsv files.

This function will simply overwrite the output from previous steps and adds a *c_call* column at the end.

==========================

Before running, the there is a need to set up a database with IMGT constant gene fasta sequences using *makeblastdb*,
basically following the instructions from https://www.ncbi.nlm.nih.gov/books/NBK279688/.

The fasta files were downloaded from IMGT and only sequences corresponding to *CH1* region for each constant gene/allele were retained. The headers were trimmed to only keep the gene and allele information. Links to find the sequences can be found here : [***human***](http://www.imgt.org/genedb/GENElect?query=7.2+IGHC&species=Homo+sapiens) and [***mouse***](http://www.imgt.org/genedb/GENElect?query=7.2+IGHC&species=Mus).

I've written a utility function `utl.makeblastdb` to prep the fasta files/databases prior to running.
```python
ddl.utl.makeblastdb('/Users/kt16/Documents/Github/dandelion/database/blast/human/human_BCR_C.fasta')
```

We really only need to do it once and then the file path can be added as an environmental variable (like above). I've set it up that we only need to point to the blast folder and ***dandelion*** will have options to specify which organisms to use in specific functions.

```bash
echo "export BLASTDB=/Users/kt16/Documents/Github/dandelion/database/blast/" >> ~/.bash_profile
source ~/.bash_profile
```


```python
for s in samples:
    filePath = s+'/dandelion/data/all_contig.fasta'
    ddl.pp.assign_isotype(filePath)
```

    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 6356/6356 [01:37<00:00, 65.28it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 502/502 [00:00<00:00, 1171.12it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 3835/3835 [00:36<00:00, 103.97it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 8651/8651 [02:40<00:00, 53.96it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10213/10213 [03:53<00:00, 43.75it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 32052/32052 [40:02<00:00, 13.34it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10880/10880 [06:10<00:00, 29.39it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 3374/3374 [00:38<00:00, 87.88it/s] 
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 4656/4656 [01:46<00:00, 43.70it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 5352/5352 [01:57<00:00, 45.53it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10922/10922 [05:05<00:00, 35.76it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 16829/16829 [20:33<00:00, 13.65it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 18745/18745 [15:30<00:00, 20.14it/s]


Initial assessment by just comparing what isn't matching between this method and *MaskPrimers.py* was basically the light chains constant regions were not called by *MarkPrimers.py* (because lack of input IgK/L primer sequences). But as far as i could see, it caught most of the discrepancies that [hk6](https://twitter.com/hamish_king) was talking about (IGHA1/IGHA2 and IGHG2/IGHG4 mis-calling) and also reduced the processing time by 10x (***~20 min*** to ***~2 min***), so I'm happy with this.

This still takes a while when dealing with large files; the number of cpus to size of file isn't exactly linear, but I enabled parallelization as default because there were noticeable improvements in processing speeds with the smaller files. Maybe it will work better on a cluster with more cpus, rather than just a standard laptop. Other than 1 sample that took about ***~40min***, most ran within ***2-5min***.

Similarly, the function can be run with
```python
fileformat = 'changeo'
```


```python
for s in samples:
    filePath = s+'/dandelion/data/all_contig.fasta'
    ddl.pp.assign_isotype(filePath, fileformat = 'changeo')
```

    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 6356/6356 [02:49<00:00, 37.40it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 502/502 [00:00<00:00, 1032.69it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 3835/3835 [00:48<00:00, 78.56it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 8651/8651 [04:07<00:00, 34.97it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10213/10213 [06:48<00:00, 24.98it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 32052/32052 [47:39<00:00, 11.21it/s]  
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10880/10880 [06:04<00:00, 29.82it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 3374/3374 [00:42<00:00, 79.73it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 4656/4656 [01:24<00:00, 55.19it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 5352/5352 [01:32<00:00, 58.15it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 10922/10922 [03:57<00:00, 45.92it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 16829/16829 [14:41<00:00, 19.09it/s]
    Retrieving contant region calls, parallelizing with 4 cpus : 100%|██████████| 18745/18745 [12:15<00:00, 25.47it/s]


# Step 4 *(optional)*:
#### Reassigning heavy chain V gene alleles.

Last step for part one of pre-processing is to use *immcantation-tigger*'s method to reassign allelic calls with `pp.reassign_alleles`. 

This impact's on how ***dandelion*** picks contigs to go forward for finding clones so it's highly reccomended to run it. It's also important when condering to do mutational analysis. However, the main caveat is that this needs to be run on multiple samples from the same subject, allowing for more information to be used to confidently assign a genotype *v_call*. In this case, all the samples I was processing have been from different organs from a single patient. So while important, this step can be skipped if you don't have the samples, or the patience, to do this. Having said that, this step works pretty quickly. 

Unfortunately, I don't have any idea whether the same method will work on light chains. Currently it just runs everything, and assigns a genotyped heavy chain V call. The light chains V calls are then transferred back from the original calls from *igblastn*. It should be technically feasible to run it through with a light chain option but I will leave it for now.


```python
# this is also a for loop for multiple samples from the same subject
ddl.pp.reassign_alleles(samples, out_folder = 'A31', sample_dict = sampledict)
```

    Processing data file(s) : 100%|██████████| 13/13 [00:01<00:00,  8.44it/s]


    Concatenating objects
       Writing out concatenated object
          Reassigning alleles
       Reading genotyped object


       Returning light chain V calls: 100%|██████████| 36313/36313 [00:10<00:00, 3449.88it/s]


       Saving corrected genotyped object


    Writing out to individual folders : 100%|██████████| 13/13 [00:02<00:00,  4.59it/s]



```python
# doing the same thing for 'changeo' format
ddl.pp.reassign_alleles(samples, out_folder = 'A31', fileformat = 'changeo', sample_dict = sampledict)
```

    Processing data file(s) : 100%|██████████| 13/13 [00:01<00:00,  9.34it/s]


    Concatenating objects
       Writing out concatenated object
          Reassigning alleles
       Reading genotyped object


       Returning light chain V calls: 100%|██████████| 45213/45213 [00:11<00:00, 3980.27it/s]


       Saving corrected genotyped object


    Writing out to individual folders : 100%|██████████| 13/13 [00:01<00:00,  7.83it/s]



```python

```
