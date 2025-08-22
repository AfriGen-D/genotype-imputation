# CHIPIMPUTATION Workflow DAG

## Main Workflow Structure

```
┌─────────────────┐
│   Input VCFs    │
│  (samplesheet)  │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────────────────┐
│              PREPROCESS SUBWORKFLOW              │
├─────────────────────────────────────────────────┤
│                                                  │
│  CHECK_FILES ──► CHECK_CHROMOSOME                │
│       │                │                         │
│       ▼                ▼                         │
│  CHECK_GENOME_BUILD    │                         │
│       │                │                         │
│       ▼                ▼                         │
│  GET_CHROMOSOME ──► GENERATE_CHUNKS_VCF          │
│       │                │                         │
│       ▼                ▼                         │
│  SPLIT_TARGET_TO_CHUNK                           │
│       │                                          │
│       ├──► GET_REF_CHROMOSOMES                   │
│       │                                          │
│       ├──► GENERATE_CHUNK_MAP ──┐                │
│       │                         │                │
│       ├──► GENERATE_REF_MAP ────┤                │
│       │                         ▼                │
│       │                   CHECK_OVERLAP          │
│       │                         │                │
│       │               [Failed]  │  [Passed]      │
│       │                   ▼     │     │          │
│       │         MERGE_ADJACENT_CHUNKS │          │
│       │                   │           │          │
│       │                   └───────────┤          │
│       │                               ▼          │
│       ├────────────────────► CHECK_MISMATCH      │
│       │                               │          │
│       ▼                               ▼          │
│    QC_DUPL ──► SPLIT_MULTI_ALLELIC               │
│       │                │                         │
│       ▼                ▼                         │
│    FILTER_MIN_AC                                 │
│       │                                          │
│       ▼                                          │
│    QC_SITE_MISSINGNESS                           │
│       │                                          │
└───────┼─────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────┐
│               PHASE SUBWORKFLOW                  │
├─────────────────────────────────────────────────┤
│                                                  │
│    EAGLE_PHASING                                 │
│    (with chromosome-specific reference panels)   │
│                                                  │
└───────┬─────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────┐
│              IMPUTE SUBWORKFLOW                  │
├─────────────────────────────────────────────────┤
│                                                  │
│    IMPUTE_MINIMAC4                               │
│    (with chromosome-specific reference panels)   │
│                                                  │
└───────┬─────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────┐
│              REPORT SUBWORKFLOW                  │
├─────────────────────────────────────────────────┤
│                                                  │
│  ┌──────────────────────────┐                   │
│  │  Per-chunk processing:    │                   │
│  │  • FILTER_INFO_BY_TARGET  │                   │
│  │  • REPORT_WELL_IMPUTED    │                   │
│  │  • REPORT_SNP_DENSITY     │                   │
│  │  • PLOT_PERFORMANCE       │                   │
│  │  • REPORT_MAF_SCATTERPLOT │                   │
│  └──────────────┬───────────┘                   │
│                 │                                │
│                 ▼                                │
│  ┌──────────────────────────┐                   │
│  │  Aggregation:             │                   │
│  │  • COMBINE_IMPUTE         │                   │
│  │  • COMBINE_INFO           │                   │
│  │  • COMBINE_SUMMARY        │                   │
│  │  • GENERATE_REPORT        │                   │
│  └──────────────┬───────────┘                   │
│                 │                                │
└─────────────────┼────────────────────────────────┘
                  │
                  ▼
         ┌────────────────┐
         │  Final Output  │
         │   • VCFs       │
         │   • Reports    │
         │   • Plots      │
         └────────────────┘
```

## Key Features:

1. **PREPROCESS**: 
   - Validates input files and genome build
   - Chunks chromosomes for parallel processing
   - QC filtering (duplicates, multi-allelic, MAF, missingness)
   - Checks overlap with reference panel
   - Merges chunks with insufficient overlap

2. **PHASE**:
   - Eagle phasing with chromosome-matched reference panels
   - Fixed: Now correctly matches chr2 with chr2 reference (not chr1)

3. **IMPUTE**:
   - Minimac4 imputation with chromosome-matched reference panels
   - Fixed: Now correctly uses per-chromosome references

4. **REPORT**:
   - Per-chunk quality metrics
   - Combined summary reports
   - Performance visualization

## Recent Changes:
- ✅ Fixed chromosome-reference panel matching bug
- ✅ Removed unused SITES_ONLY process
- ✅ Both phasing and imputation now use correct chromosome-specific panels

