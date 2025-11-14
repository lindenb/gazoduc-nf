# Parabricks

Parabricks Pipeline

## Input / Ouput



| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA REFERENCE file. | `string` |  | True |
| `prefix` | prefix for the save files. A good prefix would be 'YYYY-MM-DD.project_name.'  <details><summary>Help</summary><small>Ouput file prefix</small></details>| `string` |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. <details><summary>Help</summary><small>output directory</small></details>| `string` |  |  |
| `help` |  | `boolean` |  |  |

## Mapping FASTQs



| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` | Path to samplesheet with the appropriate suffix (csv or tsv or json). <details><summary>Help</summary><small>path to fastq samplesheet.</small></details>| `string` |  |  |
| `mapper` | wich mapper to use. bwa,pb,parabricks,dragmap | `string` |  | True |
