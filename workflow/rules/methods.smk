# Reference copy of the methods rules of snakemake-methods
# (https://github.com/niekwit/snakemake-methods). Copy into
# workflow/rules/methods.smk (same path as here) and include it after the
# other rules. See docs/adoption.md for the functions the workflow provides:
#
#   methods_envs()      conda environments (workflow/envs/{env}.yaml) used by
#                       the jobs of this configuration
#   methods_wrappers()  {name: (wrapper, [packages])} of the wrappers used,
#                       e.g. {"cutadapt": ("v9.18.0/bio/cutadapt/se", ["cutadapt"])}
#   methods_git()       {name: clone directory} of tools cloned at run time
#   methods_string_db() output files of the STRING queries ([] if none)
#   methods_targets()   all other targets (the methods describe the finished run)
#   env_packages(env)   package names of workflow/envs/{env}.yaml (see docs)
#   VERSION             workflow version
#
# Version files: results/software_versions/{name}.tsv (package, version,
# source), combined into methods/software_versions.tsv by the methods rule.
import snakemake.wrapper
from snakemake import __version__ as SNAKEMAKE_VERSION

# Methods template and database overrides of the workflow, as text (read
# while parsing). Text, not paths: source_path() gives a new runtime-cache path
# at every invocation, and a path in params would rerun the methods every time
with open(workflow.source_path("../methods/methods.md.j2"), encoding="utf-8") as fh:
    METHODS_TEMPLATE = fh.read()
with open(workflow.source_path("../methods/tools.yaml"), encoding="utf-8") as fh:
    METHODS_OVERRIDES = fh.read()

# Package names of the environment files (read while parsing)
ENV_PACKAGES = {env: env_packages(env) for env in methods_envs()}

# The installed packages of a conda environment (conda-meta records
# name-version-build.json) with these names, plus the Python and R
# interpreters if params.interpreters is 1 (the workflow's own environments,
# where its scripts run); params.source is the source column
CONDA_META_VERSIONS = (
    '(ls "$CONDA_PREFIX/conda-meta" | '
    "sed -nE 's/^(.+)-([^-]+)-([^-]+)\\.json$/\\1\\t\\2/p' | "
    "awk -F '\\t' -v OFS='\\t' -v packages='{params.packages}' "
    "-v source='{params.source}' -v interpreters='{params.interpreters}' "
    '\'BEGIN {{ n = split(packages, a, " "); for (i = 1; i <= n; i++) '
    "keep[a[i]] = 1 }} "
    '($1 in keep) || (interpreters && ($1 == "python" || $1 == "r-base")) '
    "{{ print $1, $2, source }}' | "
    "sort > {output}) 2> {log}"
)


# Version files of all sources used by this configuration
METHODS_VERSION_FILES = (
    [f"results/software_versions/{env}.tsv" for env in methods_envs()]
    + [f"results/software_versions/wrapper_{w}.tsv" for w in methods_wrappers()]
    + [f"results/software_versions/git_{name}.tsv" for name in methods_git()]
    + (["results/software_versions/string_db.tsv"] if methods_string_db() else [])
)


# Versions of the packages listed in an environment file of the workflow, as
# installed in its conda environment, without their dependencies
# -----------------------------------------------------
rule software_versions:
    output:
        temp("results/software_versions/{env}.tsv"),
    log:
        "logs/methods/software_versions_{env}.log",
    wildcard_constraints:
        env="|".join(methods_envs()),
    conda:
        "../envs/{env}.yaml"
    threads: 1
    resources:
        runtime=10,
        mem_mb=1000,
    params:
        packages=lambda w: " ".join(ENV_PACKAGES[w.env]),
        source="conda",
        interpreters=1,
    shell:
        CONDA_META_VERSIONS


# Versions of the packages of a Snakemake wrapper, from the wrapper's own
# conda environment: the conda directive is the wrapper's environment.yaml (and
# pin file), so this job runs in the environment the wrapper job ran in
# -----------------------------------------------------
for _name, (_wrapper, _packages) in methods_wrappers().items():

    rule:
        name:
            f"wrapper_versions_{_name}"
        output:
            temp(f"results/software_versions/wrapper_{_name}.tsv"),
        log:
            f"logs/methods/wrapper_versions_{_name}.log",
        conda:
            snakemake.wrapper.get_conda_env(
                _wrapper, prefix=workflow.workflow_settings.wrapper_prefix
            )
        threads: 1
        resources:
            runtime=10,
            mem_mb=1000,
        params:
            packages=" ".join(_packages),
            source="wrapper",
            # No workflow script runs in a wrapper's environment
            interpreters=0,
        shell:
            CONDA_META_VERSIONS


# Commit of a tool cloned from git at run time
# -----------------------------------------------------
rule git_versions:
    input:
        lambda w: methods_git()[w.name],
    output:
        temp("results/software_versions/git_{name}.tsv"),
    log:
        "logs/methods/git_versions_{name}.log",
    wildcard_constraints:
        name="|".join(methods_git()),
    conda:
        "../envs/methods.yaml"
    threads: 1
    resources:
        runtime=5,
        mem_mb=500,
    shell:
        'printf "%s\\t%s\\tgit\\n" {wildcards.name} '
        '"$(git -C {input} rev-parse HEAD)" > {output} 2> {log}'


# Version of STRING used by the queries (written by string_db.py when it
# queried STRING) and the date of the queries (the latest of the version
# files). Several versions are listed comma-separated (the methods then warn)
# -----------------------------------------------------
rule string_db_version:
    input:
        methods_string_db(),
    output:
        temp("results/software_versions/string_db.tsv"),
    log:
        "logs/methods/string_db_version.log",
    threads: 1
    resources:
        runtime=5,
        mem_mb=500,
    shell:
        "(printf 'string-db\\t%s\\tapi\\nstring-db-date\\t%s\\tapi\\n' "
        '"$(sort -u {input} | paste -sd , -)" '
        '"$(date -r "$(ls -t {input} | head -n 1)" +%F)" > {output}) 2> {log}'


# Methods section (Markdown), versions table and BibTeX of the parts of the
# workflow that were run. Runs after all other targets. Workflow-specific
# values for the template are added to params
# -----------------------------------------------------
rule methods:
    input:
        versions=METHODS_VERSION_FILES,
        # A function: the targets may be defined after this file is included.
        # Its result is flattened by one level (e.g. a list of expand() results)
        targets=lambda w: [
            path
            for target in methods_targets()
            for path in ([target] if isinstance(target, str) else target)
        ],
    output:
        methods="methods/methods.md",
        versions="methods/software_versions.tsv",
        bibtex="methods/references.bib",
    log:
        "logs/methods/methods.log",
    conda:
        "../envs/methods.yaml"
    threads: 1
    resources:
        runtime=10,
        mem_mb=1000,
    params:
        # The full config, so Snakemake's params-change tracking regenerates
        # methods.md when any config value changes, including ones not
        # otherwise read as their own param below (e.g. ngs_tracker.description)
        full_config=lambda w: config,
        # Functions: Snakemake would fill in "{...}" of a string param
        template_text=lambda w: METHODS_TEMPLATE,
        overrides_text=lambda w: METHODS_OVERRIDES,
        snakemake_version=SNAKEMAKE_VERSION,
        workflow_version=VERSION,
        samples=SAMPLES,
        cram=cram(),
        library_file=csv,
        cutadapt_extra=rules.cutadapt.params.extra,
        bowtie_extra=rules.count.params.extra,
        comparisons=[comparison_label(c) for c in COMPARISONS],
        comparisons_bagel2=[comparison_label(c) for c in COMPARISONS_BAGEL2],
        crisprcleanr_comparisons=[
            comparison_label(c) for c in crisprcleanr_comparisons()
        ],
        matrix_descriptions=matrix_descriptions(),
        mageck_control=mageck_control(),
        mageck_extra=extra_mageck_args(),
        string_db_data=STRING_DB_DATA,
        string_db_queries=bool(methods_string_db()),
    script:
        "../scripts/methods.py"
