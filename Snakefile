import pandas as pd
import numpy as np
import os
jn = os.path.join

configfile: "metaAndconfig/config.yaml"
sample_mt = pd.read_csv(config["sample_mt"],
                                dtype=str, sep='\t').set_index("sequence_id", drop=False)
sample_mt.drop_duplicates(inplace=True)
n_samples = len(sample_mt)

IDS = "Dcoel"               
readnames = ["1P", "2P"]
kmer_ext = ["lower", "upper", "optimal"]
assembly = ["SOAP", "ABYSS", "DISCOVAR"]

rule all:
    input:
       #expand("assemblies/{id}_ABySS/{kgenie}KSize/{id}_abyss_{kgenie}K-stats.csv", id = IDS, kgenie = kmer_ext)
       expand("QUAST/{id}/report.tsv", id = IDS),
       "reports/trim/multiqc_report.html"

rule trimmomatic:
    input:
        f = "raw/{id}_2020_1.fastq.gz",
        r = "raw/{id}_2020_2.fastq.gz"
    output:
        fout = "trimmed/{id}_1P_trim.fastq",
        funp = "trimmed/{id}_1P_unpaired.fastq",
        rout = "trimmed/{id}_2P_trim.fastq",
        runp = "trimmed/{id}_2P_unpaired.fastq"
    threads: 24
    conda:
        "envs/trimmomatic.yml"
    shell:
        "trimmomatic PE -threads {threads} {input.f} {input.r} {output.fout} {output.funp} {output.rout} {output.runp} ILLUMINACLIP:adapterseq/Adapters_PE.fa:2:30:10: LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80"

rule fastqc_trim:
    input:
        fastq = "trimmed/{id}_{read}_trim.fastq"
    output:
        "reports/trim/{id}_{read}_trim_fastqc.html",
        "reports/trim/{id}_{read}_trim_fastqc.zip"
    conda:
        "envs/QCenv.yml"
    threads: 24
    shell:
        "fastqc -o reports/trim/ {input.fastq}"

rule multiqc_trim:
    input:
        expand("reports/trim/{id}_{read}_trim_fastqc.html", id = IDS, read = readnames)
    output:
        "reports/trim/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    shell:
        "multiqc reports/trim/ -o reports/trim/"

rule soap_config:
    input:
        fRead = "trimmed/{id}_1P_trim.fastq",
        rRead = "trimmed/{id}_2P_trim.fastq"
    output:
        soapconfig = "metaAndconfig/soapconfig/{id}_soapconfig"
    shell:
        "python scripts/writesoapconfig.py {output.soapconfig} {input.fRead} {input.rRead}"

rule KmerGenie_config:
    input:
        f="trimmed/{id}_1P_trim.fastq",
        r="trimmed/{id}_2P_trim.fastq"
    output:
        "{id}_KmerGenie_config.txt"
    shell:
        "ls -1 {input.f} {input.r} > {output}"

rule KmerGenie:
    input:
        infile = "{id}_KmerGenie_config.txt"
    output:
        "KmerGenie_reports/{id}/{id}_report.html",
        logfile = "KmerGenie_reports/{id}/{id}_kmergenie.log"
    params:
        outdir = "KmerGenie_reports/{id}/{id}"
    conda:
        "envs/kmergenie.yml"
    threads:24
    shell:
        "kmergenie {input.infile} -o {params.outdir} -s 6 -l 21 -k 121 --diploid | tee > {output.logfile}"

rule best_kmer:
    input:
        result_kmer = "KmerGenie_reports/{id}/{id}_kmergenie.log"
    output:
        lower = "KmerGenie_reports/{id}/{id}_lowerK.config",
        upper = "KmerGenie_reports/{id}/{id}_upperK.config",
        optimal = "KmerGenie_reports/{id}/{id}_optimalK.config"
    shell:
        "python scripts/kmerextractor.py {input.result_kmer} {output.lower} {output.upper} {output.optimal}"


rule soapdenovo:
    input:
        ksize = "KmerGenie_reports/{id}/{id}_{kgenie}K.config",
        config = "metaAndconfig/soapconfig/{id}_soapconfig"
    output:
        "assemblies/{id}_SOAPDENOVO/{kgenie}/{id}_soap_K{kgenie}.contig"
    params:
        outdir = "assemblies/{id}_SOAPDENOVO/{kgenie}/{id}_soap_K{kgenie}"
    conda:
        "envs/SoapDenovo.yml"
    threads: 24 
    shell:
        """
        file={input.ksize}
        kSIZE=$(cat "$file")
        echo $kSIZE
        if (($kSIZE <= 63)); 
        then
        SOAPdenovo-63mer all -p {threads} -s {input.config} -K $kSIZE -R -o {params.outdir}
        elif ((63 < $kSIZE <= 127));
        then
        SOAPdenovo-127mer all -p {threads} -s {input.config} -K $kSIZE -R -o {params.outdir}
        else
        echo "K-mer size has to be set between 1 and 127"
        fi
        """

rule abyss:
    input:
        f="trimmed/{id}_1P_trim.fastq",
        r="trimmed/{id}_2P_trim.fastq",
        kmer = "KmerGenie_reports/{id}/{id}_{kgenie}K.config"
    output:
        "assemblies/{id}_ABySS/{kgenie}KSize/{id}_abyss_K{kgenie}-contigs.fa",
        "assemblies/{id}_ABySS/{kgenie}KSize/{id}_abyss_K{kgenie}-stats.csv"
    params:
        pwd = os.getcwd(),
        name = "{id}_abyss_K{kgenie}",
        outdir ="assemblies/{id}_ABySS/{kgenie}KSize/",
        logfile = "assemblies/{id}_ABySS/{kgenie}KSize/{id}_abyss_K{kgenie}.log"
    threads: 24
    singularity:
        "docker://reslp/abyss:2.2.5"
    shell:
        """
        file={input.kmer}
        kSIZE=$(cat "$file")
        echo $kSIZE
        export TMPDIR={params.pwd}/{params.outdir}tmp
        abyss-pe -C {params.pwd}/{params.outdir} k=$kSIZE name={params.name} np={threads} in='{params.pwd}/{input.f} {params.pwd}/{input.r}' | tee {params.logfile}
        """

rule quast:
    input:
        abyss = expand("assemblies/{id}_ABySS/{kgenie}KSize/{id}_abyss_K{kgenie}-contigs.fa", id = IDS,  kgenie = kmer_ext),
        soapdenovo = expand("assemblies/{id}_SOAPDENOVO/{kgenie}/{id}_soap_K{kgenie}.contig", id = IDS,  kgenie = kmer_ext)
    output:
        "QUAST/{id}/report.tsv"
    params:
        outdir = "/QUAST/{id}",
    conda:
        "envs/quast.yml"
    shell:
        "quast.py -o {output} {input}"
