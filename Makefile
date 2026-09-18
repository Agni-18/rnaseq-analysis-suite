.PHONY: envs refs bulk longread sc sc-validate validate real-example nfcore dag clean
envs:      ; mamba env create -f environments/env_bulk.yml || true
	mamba env create -f environments/env_singlecell.yml || true
	mamba env create -f environments/env_longread.yml || true
refs:        ; bash scripts/get_references.sh "chr20 chr21 chr22"
bulk:        ; snakemake -s arm1_bulk/Snakefile --use-conda -c4
longread:    ; snakemake -s arm3_longread/Snakefile --use-conda -c4
sc:          ; jupytext --to notebook arm2_singlecell/notebooks/01_scanpy_qc_to_annotation.py && jupyter lab
sc-validate: ; python arm2_singlecell/notebooks/validate_arm2.py
validate:    ; bash tests/run_all.sh
real-example: ; bash real_example/fetch_real_data.sh && Rscript real_example/real_arm1_de.R && python real_example/real_arm2_pbmc.py
nfcore:      ; bash nfcore_demo/run_test_profile.sh
dag:         ; snakemake -s arm1_bulk/Snakefile --dag | dot -Tpng > results/arm1_dag.png
clean:       ; rm -rf .snakemake results/*/tmp
