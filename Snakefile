import os

def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f
            
SAMPLES = list(listdir_nohidden("./pairs"))

##############CONFIG################
BubblePath="/share/home/zliu/test/schicluster/BubbleCluster/"
nCPUS = 90
####################################


#############RULE_ALL###############
#this part decide how far you want to go
rule all:
    input:
        expand("absPairs/{sample}.abs",sample=SAMPLES),
        expand("contactMatrix/{sample}.conmat",sample=SAMPLES),
        expand("result.txt")
############END_rule_all############

rule pairs2abs:
    input:
        chrLength="{BubblePath}otherFiles/chr.len.hg19.tsv",
        pairs="pairs/{sample}",
    output: 
        abs="absPairs/{sample}.abs"
        
    shell:"""
        set +u
        source activate
        conda activate R
        set -u
        
        if [ ! -d absPairs  ];then
            mkdir absPairs
        fi

        Rscript {BubblePath}R_scripts/pairs2abs.R {input.chrLength} {input.pairs} {output.abs}

        set +u
        conda deactivate
        set -u
        """
        
rule abs2conMatrix:
    input:
        abs=rules.pairs2abs.output.abs,
    output: 
        conmat="contactMatrix/{sample}.conmat",
    params:
        resolution=1000000,
    shell:"""
        set +u
        source activate
        conda activate R
        set -u
        
        if [ ! -d contactMatrix  ];then
            mkdir contactMatrix
        fi
        
        Rscript {BubblePath}R_scripts/pairs2abs.R {params.resolution} {input.abs} {output.conmat}

        set +u
        conda deactivate
        set -u
        """

rule cluster:
    input:
        expand("contactMatrix/{sample}.conmat",sample=SAMPLES)
    output:
        ari="result.txt"
    threads:nCPUS
    shell:"""
        set +u 
        source activate
        conda activate schicluster
        set -u
        
        python {BubblePath}Python_scripts/bubbleCluster.py -b True -i contactMatrix/ -t $[{threads}-5] -l {BubblePath}otherFiles/chr.len.hg19.tsv

        set +u
        conda deactivate
        set -u
    """
    