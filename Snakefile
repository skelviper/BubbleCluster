import os

def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f
            
SAMPLES = list(listdir_nohidden("./pairs"))
##############CONFIG################
BubbleClusterPath = ""

#############RULE_ALL###############
#this part decide how far you want to go
rule all:
    input:
        expand("absPairs/{sample}.abs",sample=SAMPLES),
        expand("contactMatrix/{sample}.conmat",sample=SAMPLES),
        expand("ARIresult.txt")
############END_rule_all############

rule pairs2abs:
    input:
        chrLength="/share/home/zliu/test/schicluster/benchmark/scriptsAndFiles/chr.len.hg19.tsv",
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

        Rscript {BubbleClusterPath}}/R_scripts/pairs2abs.R {input.chrLength} {input.pairs} {output.abs}

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
        
        Rscript {BubbleClusterPath}}/R_scripts/abs2conMatrix.R {params.resolution} {input.abs} {output.conmat}

        set +u
        conda deactivate
        set -u
        """

rule cluster:
    input:
        expand("contactMatrix/{sample}.conmat",sample=SAMPLES)
    output:
        ari="ARIresult.txt"
    params:
        contactMatrixPath="contactMatrix/",
    threads: 90
    shell:"""
        set +u 
        source activate
        conda activate schicluster
        set -u
        
        python /share/home/zliu/test/schicluster/benchmark/scriptsAndFiles/bubbleCluster.py {params.contactMatrixPath}
        
        set +u
        conda deactivate
        set -u
    """
    
