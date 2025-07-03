# Docker Directory Consolidation Log

## Date: 2024-07-03

## Summary
Consolidated `/docker/` directory into `/containers/` to eliminate duplication and confusion.

## Actions Taken

### 1. Analysis
- Found 4 tools in both locations: bcftools, impute5, crossmap, imputation_tools
- Most files were identical between locations
- Key differences found:
  - **bcftools**: `docker/bcftools/` had improved Alpine-based version with latest tools
  - **imputation_tools**: `docker/imputation_tools/` had additional variant files

### 2. Consolidation Steps
1. **bcftools**: Copied improved Alpine-based Dockerfile from `docker/bcftools/` → `containers/legacy/bcftools/`
2. **imputation_tools**: Copied variant files from `docker/imputation_tools/` → `containers/legacy/imputation_tools/`:
   - `Dockerfile.lite`
   - `Dockerfile.minimac4` 
   - `Dockerfile.fixed`
3. **impute5, crossmap**: Files were identical, no action needed

### 3. Files Preserved
- ✅ Improved bcftools Dockerfile (Alpine-based, latest versions)
- ✅ All imputation_tools variants
- ✅ All original legacy files maintained

### 4. Result
- **Before**: Duplicate containers in `/docker/` and `/containers/legacy/`
- **After**: Single consolidated location in `/containers/`
- **Benefits**: 
  - Eliminated confusion
  - Improved bcftools now in main location
  - All variant files preserved
  - Cleaner project structure

## Next Steps
- Remove `/docker/` directory (completed)
- Update any documentation references if needed
- Continue with systematic Alpine migration for other tools
