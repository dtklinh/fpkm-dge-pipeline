configfile: "config/config.yaml"

rule all:
    input:
        config["output"]["table"],
        config["output"]["plot"]

rule differential_expression:
    input:
        fpkm = config["input"]["fpkm"],
        meta = config["input"]["metadata"],
        script = "scripts/de.R"
    output:
        table = config["output"]["table"],
        plot = config["output"]["plot"]
    conda:
        "envs/environment.yaml"
    shell:
        """
        Rscript {input.script} \
            --input_fpkm {input.fpkm} \
            --metadata {input.meta} \
            --output_table {output.table} \
            --output_plot {output.plot} \
            --group_col {config[params][group_column]} \
            --control {config[params][control_level]} \
            --treatment {config[params][treatment_level]}
        """