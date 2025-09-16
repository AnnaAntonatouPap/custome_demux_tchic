# Line-by-Line Breakdown of SCCHIC_384w_c8_u3_cs2_scattered_rna

## Class Definition and Inheritance
```python
class SCCHIC_384w_c8_u3_cs2_scattered_rna(UmiBarcodeDemuxMethod):
```
- Defines new class inheriting from `UmiBarcodeDemuxMethod`
- Name indicates: SCCHIC format, 384-well, 8bp cell barcode, 3bp UMI, CS2 chemistry, scattered RNA detection

## Constructor (`__init__` method)

### Basic Setup
```python
def __init__(self, barcodeFileParser, random_primer_read=None, random_primer_length=None, **kwargs):
    self.barcodeFileAlias = 'maya_384NLA'
```
- Constructor takes barcode file parser and optional parameters
- Sets barcode file alias to 'maya_384NLA' (the ChIC barcode set)

### Parent Class Initialization
```python
UmiBarcodeDemuxMethod.__init__(
    self,
    umiRead=0,           # UMI is in read 1 (R1)
    umiStart=0,          # UMI starts at position 0
    umiLength=3,         # UMI is 3 bases long
    barcodeRead=0,       # Barcode is in read 1 (R1)
    barcodeStart=3,      # Barcode starts at position 3 (after UMI)
    barcodeLength=8,     # Barcode is 8 bases long
    random_primer_read=random_primer_read,
    random_primer_length=random_primer_length,
    barcodeFileAlias=self.barcodeFileAlias,
    barcodeFileParser=barcodeFileParser,
    **kwargs)
```
- Calls parent constructor with ChIC barcode structure parameters
- Structure: `UUU BBBBBBBB` (3bp UMI + 8bp barcode at start of R1)

### Class Metadata
```python
self.description = '384 well format, mixed transcriptome and CHiC. scCHiC: 3bp umi followed by 8bp barcode and a single A. R2 has no random primer. Transcriptome uses scattered barcodes for contamination detection'
self.shortName = 'TCHIC_scattered_rna'
self.autoDetectable = False
```
- Sets descriptive text, short name, and disables auto-detection

### Transcriptome Demultiplexer Setup
```python
self.transcriptome_demux = T_new_SCA(barcodeFileParser=barcodeFileParser, **kwargs)
```
- Creates transcriptome demultiplexer using scattered barcode format

### RNA Contamination Barcode Mapping
```python
self.id_to_cs2_barcode = {v: k for k, v in barcodeFileParser['CS2_scattered_8bp'].items()}
assert len(self.id_to_cs2_barcode) > 0
self.tx_umi_len = 6
```
- Creates mapping from well IDs to scattered barcodes (8bp each)
- **Key change**: No 'TTTTT' suffix added here (done dynamically later)
- Sets transcriptome UMI length to 6bp (3+3 from scattered format)

### ChIC Demultiplexer Setup
```python
self.chic_demux = SCCHIC_384w_c8_u3_direct_ligation(barcodeFileParser=barcodeFileParser, **kwargs)
```
- Creates ChIC demultiplexer for direct ligation protocol

### Summary Information
```python
self.barcodeSummary = f'{self.chic_demux.barcodeSummary} and {self.transcriptome_demux.barcodeSummary}'
self.longName = f'{self.chic_demux.longName} and {self.transcriptome_demux.longName}'
```
- Combines summaries from both demultiplexers

### Homopolymer Definitions
```python
self.poly_length = 10
self.poly_A = self.poly_length * 'A'  # 'AAAAAAAAAA'
self.poly_T = self.poly_length * 'T'  # 'TTTTTTTTTT'
self.poly_G = self.poly_length * 'G'  # 'GGGGGGGGGG'
```
- Defines 10-base homopolymer stretches for quality control

### Trimming Pattern
```python
self.r2_trimmer = re.compile('[GA]*$')
```
- Regex to remove trailing G and A bases from R2 end

### Sequence Capture Definition
```python
self.sequenceCapture[0] = slice(self.barcodeLength + self.umiLength + 1, None)
```
- Defines what part of R1 to keep: everything after UMI+barcode+1bp
- `slice(12, None)` = skip first 12 bases (3 UMI + 8 barcode + 1 ligation)

## Helper Methods

### R2 Trimming Method
```python
def trim_r2(self, sequence, qualities):
    start = sequence.find(self.poly_A)
    if start != -1:
        sequence, qualities = sequence[:start], qualities[:start]
```
- Trims sequence at first poly-A stretch found

```python
    start = sequence.find(self.poly_G)
    if start != -1:
        sequence, qualities = sequence[:start], qualities[:start]
```
- Trims sequence at first poly-G stretch found

```python
    sequence = self.r2_trimmer.sub('', sequence)[:-3]
    qualities = qualities[:len(sequence)]
    return sequence, qualities
```
- Removes trailing G/A bases and removes last 3 bases
- Adjusts quality scores to match trimmed sequence length

### String Representation
```python
def __repr__(self):
    return f'{self.longName} {self.description}'
```
- Returns formatted string when object is printed

### Scattered RNA UMI Extraction
```python
def extract_scattered_celseq_umi_with_polyT(self, sequence, celseq2_barcode, polyT='TTTTT'):
```
- **Key method**: Extracts UMI from scattered barcode pattern

```python
    bc1 = celseq2_barcode[:4]   # First 4 bases of 8bp barcode
    bc2 = celseq2_barcode[4:]   # Last 4 bases of 8bp barcode
    block_len = 14 + len(polyT)  # 3+4+3+4+5 = 19 bases total
```
- Splits 8bp barcode into two 4bp parts
- Calculates expected pattern length

```python
    for i in range(len(sequence) - block_len + 1):
        block = sequence[i:i+block_len]
        if block[3:7] == bc1 and block[10:14] == bc2 and block[14:] == polyT:
```
- Slides through sequence looking for pattern
- Pattern: `UUU|BBBB|UUU|BBBB|TTTTT` (positions 3-7 and 10-14 must match barcodes)

```python
            umi1 = block[0:3]    # First UMI (positions 0-2)
            umi2 = block[7:10]   # Second UMI (positions 7-9)  
            return umi1 + umi2   # Return concatenated 6bp UMI
```
- Extracts both UMI parts and concatenates them
- Returns 6bp UMI if pattern found

```python
    return None
```
- Returns None if pattern not found anywhere

## Main Demultiplexing Method

### Input Validation
```python
def demultiplex(self, records, **kwargs):
    if len(records) != 2:
        raise NonMultiplexable('Not mate pair')
```
- Ensures exactly 2 reads (paired-end sequencing)

### Ligation Tag Extraction
```python
    ligation_start = self.barcodeLength + self.umiLength  # Position 11
    ligation_end = ligation_start + 2                     # Position 13
    ligation_sequence = records[0].sequence[ligation_start:ligation_end]
    ligation_qualities = records[0].qual[ligation_start:ligation_end]
```
- Extracts 2bp ligation motif after UMI+barcode
- Saves both sequence and quality scores

### ChIC Barcode Processing
```python
    try:
        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(self, records, **kwargs)
    except NonMultiplexable:
        raise
```
- Uses parent class method to extract ChIC UMI and barcode
- Re-raises exception if ChIC demultiplexing fails

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
- Adds tags to both reads:
  - `lh`: ligation sequence
  - `lq`: ligation qualities (encoded)
  - `dt`: data type = 'CHIC'

### RNA Contamination Detection
```python
    expected_barcode = self.id_to_cs2_barcode.get(taggedRecords[0].tags['bi'])
    if expected_barcode is None:
        raise ValueError(taggedRecords[0].tags['bi'])
```
- Gets expected RNA barcode for this ChIC barcode
- Raises error if mapping not found

### Scattered Barcode Search
```python
    tx_umi = self.extract_scattered_celseq_umi_with_polyT(taggedRecords[0].sequence, expected_barcode)
    if not tx_umi:
        rc2 = reverse_complement(taggedRecords[1].sequence)
        tx_umi = self.extract_scattered_celseq_umi_with_polyT(rc2, expected_barcode)
```
- First searches R1 for scattered RNA barcode pattern
- If not found, searches reverse complement of R2
- Looks for: `UUU+BBBB+UUU+BBBB+TTTTT` pattern

### RNA Contamination Handling
```python
    if tx_umi:
        taggedRecords[1].sequence, taggedRecords[1].qualities = self.trim_r2(
            taggedRecords[1].sequence, taggedRecords[1].qualities)
        for r in taggedRecords:
            r.tags['dt'] = 'CelSeq'  # Change data type
            r.tags['rx'] = tx_umi    # Add RNA UMI
        return taggedRecords
```
- If RNA contamination found:
  - Trims R2 to remove poly-A/G tails
  - Changes data type from 'CHIC' to 'CelSeq'
  - Adds RNA UMI as 'rx' tag
  - Returns contaminated reads

### Quality Control Checks
```python
    elif 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[0].sequence or 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[1].sequence:
        raise NonMultiplexable('PolyT')
```
- Rejects reads with 23+ consecutive T's (poly-T contamination)

```python
    elif ('AGTCCGACGAT' in taggedRecords[0].sequence[:30] or
          'GTTCTACAGT' in taggedRecords[0].sequence[:30] or
          'TAATACGACTCACTATAGGG' in taggedRecords[0].sequence):
```
- Checks first 30bp of R1 for sequencing artifacts:
  - `AGTCCGACGAT`: P5 adapter sequence
  - `GTTCTACAGT`: P7 adapter sequence  
  - `TAATACGACTCACTATAGGG`: T7 promoter sequence

```python
        for r in taggedRecords:
            r.tags['dt'] = 'CelSeq'
            r.tags['RR'] = 'T7_found'
```
- If artifacts found, marks as transcriptome reads with rejection reason

### Final Return
```python
    return taggedRecords
```
- Returns processed ChIC reads (if no contamination/artifacts found)

## Summary of Flow:
1. **Input**: Paired-end reads from ChIC experiment  
2. **ChIC Processing**: Extract 3bp UMI + 8bp barcode from R1
3. **Contamination Check**: Look for scattered RNA barcode pattern
4. **Classification**: Tag as 'CHIC' or 'CelSeq' based on contamination
5. **Output**: Tagged reads with appropriate metadata
