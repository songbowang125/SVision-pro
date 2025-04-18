## SVision-pro
SVision-pro detects germline, de novo and somatic structural variants (SV), including simple SV (SSV) and complex SV (CSV). 

Both long reads from HiFi or other technologies, and genome assemblies are acceptable.


SVision-pro implements genome-to-genome representation and comparative SV detection, genotyping and differentiating using a neural network based image instance segmentation framework. 

SVision-pro comprises two modules: a genome-to-genome representation module encoding genomic features from two samples to an image, from which a neural-network recognition module comparatively recognizes SVs as well as their inter-genome differences. In particular, SVision-pro formulates SV detection, genotyping and differentiating as a neural-network-based image instance segmentation task, facilitating the discovery of both de novo and somatic SV and CSV.

<div align=left><img width=80% height=80% src="https://github.com/songbowang125/SVision-pro/blob/main/src/pre_process/workflow.png"/></div> 

## Papers
You can cite both the 'pro' and 'orignial' versions:

SVision-pro: ["Nature Biotechnology (2024)"](https://www.nature.com/articles/s41587-024-02190-7)

SVision: ["Nature Methods (2022)"](https://www.nature.com/articles/s41592-022-01609-w)

SVision-popolation on reads: coming soon

SVision-popolation on assemblies: coming soon


## License

SVision-pro is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University.
For more information, please contact with Songbo Wang (songbowang125@163.com) or Kai Ye (kaiye@xjtu.edu.cn).

## Installation

### Operation systems
Linux-based systems required, including but are not limited to: 
* MacOS, Big Sur
* Ubuntu (including Windows Subsystem)
* CentOS

### Install from source

```
## Get the source code
git clone https://github.com/songbowang125/SVision-pro.git
cd SVision-pro

## Create a conda environment for SVision-pro
conda env create -f ./environment.yml 

## Install from source
conda activate svision-pro-env
python setup.py install

Note: Please make sure you have installed the same version of dependencies as in ./environment.yml. 
      We recommand that you create a new conda env and install by the command lines above. 
```

## Usage

The input of SVision-pro requires one target genome (aka case genome), N base genomes (aka control genomes, N>=0), the reference genome, the neural-network model (in ./src/pre_process), the sample name and the output path

There are four detection mode in SVision-pro: 

Command line format for 'germline' mode (N=0):
```commandline
## run SVision-pro
SVision-pro --target_path /path/to/target.bam  --genome_path /path/to/reference.fasta --model_path /path/to/model.pth --out_path /path/to/output/ --sample_name sample1 --detect_mode germline
```

Command line format for 'somatic' mode (N=1):
```commandline
## run SVision-pro
SVision-pro --target_path /path/to/tumor.bam  --base_path /path/to/normal.bam --genome_path /path/to/reference.fasta --model_path /path/to/model.pth --out_path /path/to/output/ --sample_name sample1 --detect_mode somatic

## extract somatic calls
python extract_op.py --input_vcf /path/to/svision_pro.vcf --extract somatic
```


Command line format for 'denovo' mode (N=2):
```commandline
## run SVision-pro
SVision-pro --target_path /path/to/child.bam  --base_path /path/to/father.bam /path/to/mother.bam --genome_path /path/to/reference.fasta --model_path /path/to/model.pth --out_path /path/to/output/ --sample_name sample1 --detect_mode denovo

## extract denovo calls
python extract_op.py --input_vcf /path/to/svision_pro.vcf --extract denovo
```

Command line format for 'genotype' mode (N=N):
```commandline
## run SVision-pro
SVision-pro --target_path /path/to/target.bam  --base_path /path/to/base1.bam /path/to/base2.bam /path/to/base3.bam... --genome_path /path/to/reference.fasta --model_path /path/to/model.pth --out_path /path/to/output/ --sample_name sample1 --detect_mode genotype --region chr1:1000-2000
```
Model (or image size) parameter:
```commandline
## For now, SVision-pro provides three different LiteUnet models (256/512/1024, in src/pre_process/model_liteunet_xxx_8_16_32_32_32.pth).
The three models correspond to three differenct image size settings (256, 512 and 1024), which determine the minimum detectable SV frequencies (0.04, 0.02 and 0.01). The description about this setting could be found in the Method section of our paper. 
   
In common, image size 256 (default) is capable for most conditions. 
So, you can directly use the 256 model (src/pre_process/model_liteunet_256_8_16_32_32_32.pth). 

While if you use other models/image sizes, please use this two parameters simultaneously:
SVision-pro ... ... --img_size 512 --model_path model_liteunet_512_8_16_32_32_32.pth, or
SVision-pro ... ... --img_size 1024 --model_path model_liteunet_1024_8_16_32_32_32.pth
```


Recommended parameters:

```commandline
## We have provided access BED files for human genomes (src/pre_process), which exclude centromere and heterochromatin regions.
--access_path   Absolute path to access BED file that contains accessable regions of the reference genome
```
Other optional parameters:
```commandline
Input parameters:
--process_num   Thread numbers
--preset        Sequence type, including hifi, error-prone and asm (for assembly-based calling)

Filter paramters:
--min_mapq      Minimum read map quality (defalut: 20)
--min_supp      Minimum supporting read number for considering a SV (defalut: 5)
--min_sv_size   Minimum SV size/length for considering a SV (defalut: 50)
--max_sv_size   Maximum SV size/length for considering a SV (defalut: 50000)
--interval_size The sliding genome region size for searching SVs (defalut: 10000000)
--max_coverage  Maximum read coverage of the searched genome regions (defalut: 500)
--skip_coverage_filter  SKip filtering genome regions by 'max_coverage'

Neural network parameters:
--img_size      The representation image sizes (default: 256, optional: 512, 1024)
--batch_size    The batch size for neural network prediction (default: 1)
--device        Using CPU or GPU for prediction (defalut: 'cpu')
--net_cpu_num   When using CPU, specific the CPU core number for each prediction sub-process (defalut: 2)
--gpu_id        When using GPU, specific the GPU id in your device (defalut: '0')
```

### Run demo data

The demo data located at './src/pre_process/', which contains a CSV deletion-inversion and is extracted from HiFi-sequenced AshkenazimTrio from Genome-In-A-Bottle.
* ``` demo.HG002.child.bam ```
* ``` demo.HG003.father.bam ```
* ``` demo.HG004.mother.bam ```

1. Download reference genome [GRCh38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)


2. Checking SVision-pro
  ``` commandline
SVision-pro -h
```

3. Run SVision-pro 'germline' mode
```commandline
SVision-pro --target_path ./src/pre_process/demo.HG002.child.bam  --genome_path /path/to/downloaded.GRCh38.fa --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.path --access_path ./src/pre_process/hg38.access.10M.bed --sample test --out_path ./test_out/ --detect_mode germline --process 8
```

4. Run SVision-pro 'denovo' mode

```commandline
SVision-pro --target_path ./src/pre_process/demo.HG002.child.bam  --base_path ./src/pre_process/demo.HG003.father.bam ./src/pre_process/demo.HG004.mother.bam --genome_path /path/to/downloaded.GRCh38.fa --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --access_path ./src/pre_process/hg38.access.10M.bed --sample test --out_path ./test_out/ --detect_mode denovo --process 8
```


5. Output files

* ``` *.vcf ``` The standard VCF output with genome-to-genome comparison info columns.


## Contact
If you have any questions, please feel free to contact: songbowang125@163.com
