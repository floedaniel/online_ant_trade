# GABI Species Name Matching Fix

## Problem Identified

**Issue**: Script was not finding GABI records even when they existed.

**Example**:
- Species in loop: `"Camponotus cruentatus (Latreille, 1802)"` (with author citation)
- Species in GABI: `"Camponotus cruentatus"` (without author citation)
- Result: **0 matches** (should be 853 records with coordinates)

## Root Cause

The original matching code was doing exact string comparison:
```r
filter(valid_species_name == species_name)
```

This failed because:
1. **Author citations**: Species names in the loop include author and year in parentheses
2. **Case sensitivity**: No case-insensitive matching
3. **Whitespace**: No trimming of leading/trailing spaces

## Solution Implemented

### 1. Remove Author Citations
```r
species_name_clean <- gsub("\\s*\\([^)]+\\)\\s*$", "", species_name)
species_name_clean <- trimws(species_name_clean)
```

**Examples:**
- `"Camponotus cruentatus (Latreille, 1802)"` → `"Camponotus cruentatus"`
- `"Formica rufa Linnaeus, 1761"` → `"Formica rufa Linnaeus, 1761"` (no parens, no change)
- `"Myrmica rubra (Linnaeus, 1758)"` → `"Myrmica rubra"`

### 2. Case-Insensitive Matching with Trimming
```r
filter(tolower(trimws(valid_species_name)) == tolower(species_name_clean))
```

**Benefits:**
- Handles case variations: `"Camponotus Cruentatus"` = `"camponotus cruentatus"`
- Removes whitespace: `" Camponotus cruentatus "` = `"Camponotus cruentatus"`
- More robust matching

### 3. Added Debug Message
```r
message("Matching species name: '", species_name_clean, "'")
```

Now you can see exactly what name is being searched for.

## Test Case: Camponotus cruentatus

**Before Fix:**
```
Searching for species in local GABI database...
No matching records found in GABI database
```

**After Fix:**
```
Searching for species in local GABI database...
Matching species name: 'Camponotus cruentatus'
Found 1018 records in GABI database
GABI records with valid coordinates: 853
Total records after adding GABI data: 1853
```

**Breakdown:**
- 1018 total GABI records
- 853 have valid coordinates (dec_lat and dec_long not empty)
- 165 records excluded (missing coordinates)

## Validation

### Manual Check:
```bash
# Count total records
grep -i "Camponotus cruentatus" gabi_antmaps_data_clean.csv | wc -l
# Result: 1018

# Count records with coordinates
grep -i "Camponotus cruentatus" gabi_antmaps_data_clean.csv |
  awk -F',' '{print $14","$15}' | grep -v "^,$" | wc -l
# Result: 853
```

✅ **Matches script output perfectly!**

## Code Changes

**File**: `updated_sdmtune_loop.R`
**Lines**: 233-247

**Key Changes:**
1. Added `gsub()` to remove author citations
2. Added `trimws()` to clean whitespace
3. Changed to `tolower()` comparison for case-insensitivity
4. Removed `acceptedKey` matching (not reliable in GABI data)
5. Added debug message showing cleaned name

## Impact

### Expected Improvements:
- ✅ **More data**: Previously missed GABI records now included
- ✅ **Better matching**: Handles name format variations
- ✅ **Transparency**: Debug message shows exact search term
- ✅ **Robust**: Works with various naming conventions

### Example Species That Will Benefit:
All species with author citations in parentheses:
- `"Camponotus cruentatus (Latreille, 1802)"`
- `"Formica aquilonia (Yarrow, 1955)"`
- `"Myrmica rubra (Linnaeus, 1758)"`
- `"Lasius niger (Linnaeus, 1758)"`

## Testing Recommendations

After running the script, check the console output for:

```
Matching species name: '[cleaned name]'
Found [N] records in GABI database
GABI records with valid coordinates: [M]
```

Where:
- `N` = total GABI records for species
- `M` = records with coordinates (M ≤ N)
- `M` should match your manual counts

## Summary Statistics

For production runs, the summary file will show:
```
raw_online: 1000
raw_gabi_antmaps: 853
raw_records: 1853
data_sources: gbif: 1000; gabi_antmaps: 853
```

---

**Fix Date**: 2025-12-04
**Status**: ✅ Implemented and tested
**Validation**: Confirmed with manual GABI data inspection
