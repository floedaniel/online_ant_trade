# GABI AntMaps Data Integration Summary

## Overview
The script now integrates local occurrence data from the GABI AntMaps database with online sources (GBIF, BISON, ALA, iNat, iDigBio) to maximize data availability for species distribution modeling.

## Implementation Details

### 1. Data Loading (Lines 35-50)
**At script startup:**
- Loads `gabi_antmaps_data_clean.csv` once at the beginning
- File location: `./2_processed_data/gabi_antmaps_data_clean.csv`
- Includes error handling if file is missing
- Reports total number of records available

### 2. Species Matching (Lines 233-274)
**For each species in the loop:**

#### Matching Strategy:
The script searches GABI data using TWO criteria:
1. **By species name**: `valid_species_name == species_name`
2. **By GBIF key**: `acceptedKey == species_key` (if available)

This dual matching ensures maximum data retrieval even if names differ slightly.

#### Column Mapping:
GABI database columns are renamed to match the standard format:
```r
GABI Column          →  Standard Column
------------------------------------------
valid_species_name   →  species
dec_long             →  decimalLongitude
dec_lat              →  decimalLatitude
(constant)           →  source = "gabi_antmaps"
```

#### Data Validation:
- Filters out records with missing coordinates (`NA` values)
- Reports count before and after coordinate validation

### 3. Data Combination
**Merging process:**
- Online data sources downloaded first (GBIF, BISON, ALA, iNat, iDigBio)
- GABI records added using `bind_rows()`
- All source labels preserved in `source` column
- Duplicates handled by subsequent cleaning steps

### 4. Tracking & Reporting

#### Console Messages:
```
Raw occurrence records from online sources: 450
Searching for species in local GABI database...
Found 75 records in GABI database
GABI records with valid coordinates: 72
Total records after adding GABI data: 522
Raw occurrence records (total): 522 (online: 450, GABI: 72)
```

#### Summary File Enhancement:
New columns added to `summary_[species_key].xlsx`:

| Column | Description |
|--------|-------------|
| `raw_records` | Total records (online + GABI) |
| `raw_online` | Records from online sources only |
| `raw_gabi_antmaps` | Records from GABI AntMaps |
| `data_sources` | Detailed breakdown by source |

**Example `data_sources` value:**
```
gbif: 320; gabi_antmaps: 72; inat: 95; bison: 35
```

## Benefits

### 1. **Increased Data Coverage**
- Supplements online databases with curated local records
- Particularly valuable for:
  - Species with limited online presence
  - Historical records not digitized in GBIF
  - Museum specimens in GABI database

### 2. **Data Quality**
- GABI AntMaps is a curated, validated database
- Reduces reliance on potentially incomplete online sources
- Maintains source tracking for data provenance

### 3. **Transparency**
- Complete audit trail of data sources
- Can identify contribution of each source
- Facilitates citation and acknowledgment

### 4. **Fallback Protection**
- If online sources fail or have issues
- GABI data provides backup occurrence records
- Reduces species skipped due to insufficient data

## Data Flow Diagram

```
START
  ↓
Load GABI AntMaps data (once)
  ↓
For each species:
  ↓
Download online data (GBIF, BISON, etc.)
  ↓
Search GABI by species name/key
  ↓
Found matches? → YES → Validate coordinates
                  ↓
                Extract & rename columns
                  ↓
                Combine with online data
  ↓
Continue with cleaning, filtering, modeling...
  ↓
Save summary with source breakdown
```

## Example Output

### Species with GABI records:
```
Processing: Camponotus herculeanus
Raw occurrence records from online sources: 1,245
Found 89 records in GABI database
GABI records with valid coordinates: 87
Total records: 1,332 (online: 1,245, GABI: 87)
```

### Species without GABI records:
```
Processing: Rare species X
Raw occurrence records from online sources: 156
Searching for species in local GABI database...
No matching records found in GABI database
Total records: 156 (online: 156, GABI: 0)
```

## File Requirements

### GABI Data File Must Have:
- **Filename**: `gabi_antmaps_data_clean.csv`
- **Location**: `./2_processed_data/`
- **Required columns**:
  - `valid_species_name` (character)
  - `dec_long` (numeric)
  - `dec_lat` (numeric)
- **Optional column**:
  - `acceptedKey` (numeric) - enhances matching

### CSV Structure:
```csv
valid_species_name,dec_long,dec_lat,acceptedKey,...
Camponotus herculeanus,10.5,60.3,1234567,...
Formica rufa,8.2,59.1,2345678,...
```

## Error Handling

The script gracefully handles:
- Missing GABI file → Warning message, continues with online data only
- No matches for species → Message logged, continues without GABI data
- Missing coordinates → Records filtered out, count reported
- Column name variations → Uses flexible matching

## Backward Compatibility

✅ **Fully backward compatible:**
- If GABI file is missing → Script runs normally with online data only
- Existing functionality unchanged
- No breaking changes to downstream processes

## Performance Impact

- **Minimal**: GABI data loaded once at startup
- **Fast matching**: Simple filter operations per species
- **No network overhead**: Local file reading only

---

**Integration Date**: 2025-12-04
**Script Version**: Enhanced with GABI AntMaps integration
**Status**: Production ready
