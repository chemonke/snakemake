import snakemake

rule count_words:
    input: "./input.txt"
    output: "./output.txt"
    shell: "wc -w {input} > {output}"