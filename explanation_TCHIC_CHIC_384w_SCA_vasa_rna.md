# Line-by-Line Breakdown of TCHIC_CHIC_384w_SCA_vasa_rna

## Class Definition and Inheritance
```python
class TCHIC_CHIC_384w_SCA_vasa_rna(ScatteredUmiBarcodeDemuxMethod):
```
- Defines new class inheriting from `ScatteredUmiBarcodeDemuxMethod` (not UmiBarcodeDemuxMethod)
- Name indicates: TCHIC format, scattered ChIC barcodes, with VASA RNA detection

## Constructor (`__init__` method)

### Parameter Setup
```python
def __init__(self, 
             barcodeFileParser, 
            random_primer_read=None, 
            random_primer_length=None, 
            first_umi_len = 3,
            first_bc_len = 4,
            second_umi_len = 3,
            second_barcode_len = 4,
            barcode_alias=None, 
            **kwargs):
```
- Takes parameters for scattered barcode structure (UMI and barcode lengths for each position)
- Default structure: 3bp UMI + 4bp BC + 3bp UMI + 4bp BC = 14bp total

### Barcode File Alias Setup
```python
if barcode_alias is None:
    self.barcodeFileAlias = 'DamID2_scattered_8bp'
else:
    self.barcodeFileAlias= barcode_alias
```
- Sets default to scattered barcode format
- Allows override with custom alias

### Scattered Barcode Structure Calculation
```python
self.total_mi_len = first_umi_len+first_bc_len+second_umi_len+second_barcode_len
```
- Calculates total length: 3+4+3+4 = 14bp

### Barcode Slice Definitions
```python
barcode_slices = (
    ( slice(first_umi_len,first_umi_len+first_bc_len), 
        slice(first_umi_len+first_bc_len+second_umi_len,first_umi_len+first_bc_len+second_umi_len+second_barcode_len )),  ()
)
```
- Defines where barcodes are located in R1:
  - First barcode: `slice(3,7)` = positions 3-6 (BC1)
  - Second barcode: `slice(10,14)` = positions 10-13 (BC2)
- Empty tuple `()` for R2 (no barcodes there)

### UMI Slice Definitions
```python
umi_slices = (
    (slice(0,first_umi_len), 
        slice(first_umi_len+first_bc_len,first_umi_len+first_bc_len+second_umi_len)
        ),
    ()
)
```
- Defines where UMIs are located in R1:
  - First UMI: `slice(0,3)` = positions 0-2 (UMI1)
  - Second UMI: `slice(7,10)` = positions 7-9 (UMI2)
- Empty tuple `()` for R2 (no UMIs there)

### Sequence Capture Definitions
```python
capture_slices = (
    slice(self.total_mi_len + 2, None),  # slice(16, None) - skip 14bp + 2bp ligation
    slice(None)                          # Keep all of R2
)
```
- R1: Keep everything after position 16 (skip 14bp barcode/UMI + 2bp ligation motif)
- R2: Keep entire sequence

### Parent Class Initialization
```python
ScatteredUmiBarcodeDemuxMethod.__init__(
    self,
    barcode_slices = barcode_slices,
    umi_slices = umi_slices,
    capture_slices = capture_slices,
    barcodeFileAlias=self.barcodeFileAlias,
    barcodeFileParser=barcodeFileParser,
    **kwargs)
```
- Initializes scattered barcode demultiplexer with defined slice positions

### Class Metadata
```python
self.shortName = 'TCHIC_3u4b3u4b_vasa'
self.longName = 'DamID-CHIC 384w, 3bp UMI, 4bp CB, 3bp UMI, 4bp CB with VASA RNA detection'
self.autoDetectable = False
self.description = 'CHIC with DamID scattered barcodes, starting with a 3bp UMI, 4bp CB, 3bp UMI, 4bp CB. RNA contamination uses VASA barcodes'
```
- Sets descriptive information
- Indicates scattered ChIC with VASA RNA detection

### Transcriptome Demultiplexer Setup
```python
self.transcriptome_demux = CELSeq2_c8_u6(barcodeFileParser=barcodeFileParser, **kwargs)
```
- Creates simple (non-scattered) transcriptome demultiplexer

### RNA Contamination Barcode Mapping (KEY CHANGE)
```python
self.id_to_cs2_barcode = {v: k + 'TTTTT' for k, v in barcodeFileParser['celseq2'].items()}
assert len(self.id_to_cs2_barcode) > 0
self.tx_umi_len = 6
```
- **KEY DIFFERENCE**: Uses `celseq2` barcodes (not scattered)
- **Adds 'TTTTT' suffix** directly to mapping (unlike the first demux)
- Creates mapping: well_id → simple_barcode+'TTTTT'
- Sets transcriptome UMI length to 6bp

### ChIC Demultiplexer Reference
```python
self.chic_demux = SCCHIC_384w_c8_u3_direct_ligation(barcodeFileParser=barcodeFileParser, **kwargs)
```
- References the simple ChIC demultiplexer (matches SCCHIC_384w_c8_u3_cs2 pattern)

### Summary Information
```python
self.barcodeSummary = f'{self.chic_demux.barcodeSummary} and {self.transcriptome_demux.barcodeSummary}'
self.longName = f'{self.chic_demux.longName} and {self.transcriptome_demux.longName}'
```
- Combines information from both demultiplexers

### Homopolymer and Trimming Setup
```python
self.poly_length = 10
self.poly_A = self.poly_length * 'A'
self.poly_T = self.poly_length * 'T'
self.poly_G = self.poly_length * 'G'
self.r2_trimmer = re.compile('[GA]*$')
```
- Same as other demux strategies - defines quality control patterns

## Helper Methods

### R2 Trimming Method
```python
def trim_r2(self, sequence, qualities):
    start = sequence.find(self.poly_A)
    if start != -1:
        sequence, qualities = sequence[:start], qualities[:start]

    start = sequence.find(self.poly_G)
    if start != -1:
        sequence, qualities = sequence[:start], qualities[:start]

    sequence = self.r2_trimmer.sub('', sequence)[:-3]
    qualities = qualities[:len(sequence)]
    return sequence, qualities
```
- Identical to other demux strategies
- Trims poly-A, poly-G, trailing G/A bases, and last 3 bases

### String Representation
```python
def __repr__(self):
    return f'{self.longName} {self.description}'
```
- Standard string representation

### VASA UMI Extraction (FROM SCCHIC_384w_c8_u3_cs2)
```python
def extract_vasa_umi(self, sequence, vasa_barcode):
    vasa_umi_end = sequence.find(vasa_barcode)
    vasa_umi_start = max(0, vasa_umi_end - self.tx_umi_len)
    if vasa_umi_end - vasa_umi_start > 0:
        return sequence[vasa_umi_start:vasa_umi_end]
    return None
```
- **KEY METHOD**: Extracts UMI from simple barcode pattern
- Looks for pattern: `...UMI + BARCODE + TTTTT`
- Extracts 6bp UMI immediately before the found barcode
- Returns UMI if found, None if not found

## Main Demultiplexing Method

### Input Validation
```python
def demultiplex(self, records, **kwargs):
    if len(records) != 2:
        raise NonMultiplexable('Not mate pair')
```
- Ensures paired-end reads

### Ligation Tag Extraction
```python
ligation_start = self.total_mi_len      # Position 14
ligation_end = ligation_start + 2       # Position 16
ligation_sequence = records[0].sequence[ligation_start:ligation_end]
ligation_qualities = records[0].qual[ligation_start:ligation_end]
```
- Extracts 2bp ligation motif after scattered barcode region (positions 14-15)

### Scattered ChIC Barcode Processing
```python
try:
    taggedRecords = ScatteredUmiBarcodeDemuxMethod.demultiplex(self, records, **kwargs)
except NonMultiplexable:
    raise
```
- Uses scattered barcode method to extract ChIC UMI and barcodes
- Pattern: `UUU+BBBB+UUU+BBBB` from positions 0-13

### Tag Addition
```python
ud = {
    'lh': ligation_sequence,
    'lq': phredToFastqHeaderSafeQualities(ligation_qualities),
    'dt': 'CHIC'
}
taggedRecords[0].tags.update(ud)
taggedRecords[1].tags.update(ud)
```
- Adds standard tags to both reads

### Simple VASA RNA Contamination Detection (FROM SCCHIC_384w_c8_u3_cs2)
```python
expected_barcode = self.id_to_cs2_barcode.get(taggedRecords[0].tags['bi'])
if expected_barcode is None:
    raise ValueError(taggedRecords[0].tags['bi'])
```
- Gets expected simple barcode (already includes 'TTTTT' suffix)

### Simple Pattern Search (NOT SCATTERED)
```python
if expected_barcode in taggedRecords[0].sequence or expected_barcode in reverse_complement(taggedRecords[1].sequence):
```
- **KEY DIFFERENCE**: Looks for simple continuous pattern, not scattered
- Pattern: `UMI + BARCODE + TTTTT` (e.g., `AAAAAACGTGCTATTTTT`)

### VASA UMI Extraction
```python
    if expected_barcode in taggedRecords[0].sequence:
        tx_umi = self.extract_vasa_umi(taggedRecords[0].sequence, expected_barcode)
    else:
        rc2 = reverse_complement(taggedRecords[1].sequence)
        tx_umi = self.extract_vasa_umi(rc2, expected_barcode)
```
- Uses simple UMI extraction (not scattered)
- Looks for 6bp UMI immediately before the barcode+polyT

### RNA Contamination Handling
```python
    taggedRecords[1].sequence, taggedRecords[1].qualities = self.trim_r2(
        taggedRecords[1].sequence, taggedRecords[1].qualities)
    for r in taggedRecords:
        r.tags['dt'] = 'VASA'  # Note: VASA, not CelSeq
        if tx_umi is not None:
            r.tags['rx'] = tx_umi
    return taggedRecords
```
- Trims R2 for RNA processing
- Tags as 'VASA' (matches SCCHIC_384w_c8_u3_cs2 approach)
- Adds simple UMI (6bp, not scattered)

### Quality Control Checks
```python
elif 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[0].sequence or 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[1].sequence:
    raise NonMultiplexable('PolyT')
elif ('AGTCCGACGAT' in taggedRecords[0].sequence[:30] or
      'GTTCTACAGT' in taggedRecords[0].sequence[:30] or
      'TAATACGACTCACTATAGGG' in taggedRecords[0].sequence):
    for r in taggedRecords:
        r.tags['dt'] = 'VASA'
        r.tags['RR'] = 'T7_found'
```
- Same quality control as other strategies
- Tags artifacts as 'VASA' instead of 'CelSeq'

### Final Return
```python
return taggedRecords
```
- Returns processed scattered ChIC reads

## Summary of Key Differences:

### Compared to SCCHIC_384w_c8_u3_cs2_scattered_rna:
1. **ChIC processing**: Uses scattered barcodes (3+4+3+4) instead of simple (3+8)
2. **RNA detection**: Uses simple VASA pattern instead of scattered pattern
3. **Parent class**: Inherits from `ScatteredUmiBarcodeDemuxMethod`
4. **Barcode mapping**: Includes 'TTTTT' suffix in mapping
5. **RNA tags**: Uses 'VASA' instead of 'CelSeq'

### Flow:
1. **Input**: Paired-end reads from scattered ChIC experiment
2. **ChIC Processing**: Extract scattered pattern `UUU+BBBB+UUU+BBBB`
3. **Contamination Check**: Look for simple RNA pattern `UMI+BARCODE+TTTTT`
4. **Classification**: Tag as 'CHIC' or 'VASA' based on contamination
5. **Output**: Tagged reads with appropriate metadata

This creates the reverse combination you wanted: scattered ChIC barcodes with simple VASA RNA contamination detection!
