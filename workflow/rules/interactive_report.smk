rule interactive_results_report:
    input:
        **interactive_results_report_input(),
    output:
        html="results/interactive_report.html",
    log:
        "logs/interactive_results_report.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        **interactive_results_report_params(),
    script:
        "../scripts/interactive_results_report.py"
