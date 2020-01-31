#!usr/bin/en python3
#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV3.simg
#cp ~/git/MinION_HBV/Snakefile  /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR

#la ou sont mes datas
workdir : "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/"

(BARCODE,READ) = glob_wildcards('DATASET/{barcode}/{read}.fastq') 

rule merge:
    input:
        fastq_files = expand('DATASET/{barcode}/{read}.fastq',zip,barcode=BARCODE,read=READ)
    shell:
        "echo 'lolz'"    
