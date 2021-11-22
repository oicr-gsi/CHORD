# CHORD

A WDL implementation of [UMCUGenetics CHORD](https://github.com/UMCUGenetics/CHORD).

# Usage

## Cromwell

Load cromwell:

`module load cromwell`

Then, to submit the workflow:

```
cromwell submit runChord.wdl --inputs /path/to/inputs.json --host http://cromwell-dev.hpc.oicr.on.ca:8000
```

## Inputs

**Required workflow parameters**

| Parameter   | Value | Description       |
| ----------- | ----- | ----------------- |
| `snvFile`   | File  | path to SNV VCF   |
| `indelFile` | File  | path to indel VCF |
| `svFile`    | File  | path to SV VCF    |
| `script`    | File  | path to R script  |

## Outputs

| Output        | Type | Description                                                        |
| ------------- | ---- | ------------------------------------------------------------------ |
| `contexts`    | File | merged_contexts.txt, txt file with extracted mutational signatires |
| `predictions` | File | chord_pred.txt, txt file with predictions                          |

## Debugging

To debug WDL execution, it is helpful to examine error logs located at

`/scratch2/groups/gsi/development/cromwell/cromwell-dev.hpc.oicr.on.ca/CHORD/{WORKFLOW-ID}/call-runCHORD/execution/`

Sometimes, there is an error with `ref.width` when running the script. This is indicative of an incorrect reference genome. Use `hg.19` for CHORD example data, and `hg.38` for OICR data.
