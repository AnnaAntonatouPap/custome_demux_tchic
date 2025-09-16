class SCCHIC_384w_c8_u3_cs2_scattered_rna(UmiBarcodeDemuxMethod):

    def __init__(self, barcodeFileParser, random_primer_read=None, random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)

        self.description = '384 well format, mixed transcriptome and CHiC. scCHiC: 3bp umi followed by 8bp barcode and a single A. R2 has no random primer. Transcriptome uses scattered barcodes for contamination detection'
        self.shortName = 'TCHIC_scattered_rna'

        self.autoDetectable = False

        # The demultiplexer used for the transcriptome reads:
        self.transcriptome_demux = CELSeq2_c8_u6(barcodeFileParser=barcodeFileParser, **kwargs)

        # CHANGED: Use scattered barcodes for RNA contamination detection instead of celseq2
        # Contains expected bleedthrough sequence from scattered barcodes
        self.id_to_cs2_barcode = {v: k + 'TTTTT' for k, v in barcodeFileParser['CS2_scattered_8bp'].items()}
        assert len(self.id_to_cs2_barcode) > 0
        self.tx_umi_len = 6

        # The demultiplexer used for the chic reads:
        self.chic_demux = SCCHIC_384w_c8_u3_direct_ligation(barcodeFileParser=barcodeFileParser, **kwargs)

        self.barcodeSummary = f'{self.chic_demux.barcodeSummary} and {self.transcriptome_demux.barcodeSummary}'
        self.longName = f'{self.chic_demux.longName} and {self.transcriptome_demux.longName}'

        self.poly_length = 10
        self.poly_A = self.poly_length * 'A'
        self.poly_T = self.poly_length * 'T'
        self.poly_G = self.poly_length * 'G'

        self.r2_trimmer = re.compile('[GA]*$')

        self.sequenceCapture[0] = slice(
            self.barcodeLength + self.umiLength + 1,
            None)  # dont capture the first base

    def trim_r2(self, sequence, qualities):
        start = sequence.find(self.poly_A)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        start = sequence.find(self.poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        # Trim any trailing A and G bases from the end and # Trim down 3 bases
        sequence = self.r2_trimmer.sub('', sequence)[:-3]
        qualities = qualities[:len(sequence)]

        return sequence, qualities

    def __repr__(self):
        return f'{self.longName} {self.description}'

    def extract_vasa_umi(self, sequence, vasa_barcode):
        # Extract vasa umi from sequence given the vasa barcode
        # Returns None when it cannot be extracted
        vasa_umi_end = sequence.find(vasa_barcode)
        vasa_umi_start = max(0, vasa_umi_end - self.tx_umi_len)
        if vasa_umi_end - vasa_umi_start > 0:
            return sequence[vasa_umi_start:vasa_umi_end]
        return None

    def extract_scattered_celseq_umi(self, sequence, celseq2_barcode):
        """
        Search for the scattered barcode (as a single 8bp string) anywhere in the sequence,
        and extract the 6bp UMI (3bp before bc1, 3bp before bc2).
        Returns the concatenated UMI, or None if not found.
        """
        # Remove the 'TTTTT' suffix that was added for contamination detection
        clean_barcode = celseq2_barcode.replace('TTTTT', '')
        bc1 = clean_barcode[:4]
        bc2 = clean_barcode[4:]
        block_len = 14  # 3umi + 4bc + 3umi + 4bc

        for i in range(len(sequence) - block_len + 1):
            block = sequence[i:i + block_len]
            if block[3:7] == bc1 and block[10:14] == bc2:
                umi1 = block[0:3]
                umi2 = block[7:10]
                return umi1 + umi2
        return None

    def demultiplex(self, records, **kwargs):
        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')

        # add first 2 bases as ligation tag:
        ligation_start = self.barcodeLength + self.umiLength
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]
        
        # Obtain the chic barcode and umi:
        try:
            taggedRecords = UmiBarcodeDemuxMethod.demultiplex(self, records, **kwargs)
        except NonMultiplexable:
            raise

        # Trim ligation motif:
        # add first 2 bases as ligation tag:
        ud = {
            'lh': ligation_sequence,
            'lq': phredToFastqHeaderSafeQualities(ligation_qualities),
            'dt': 'CHIC'
        }

        taggedRecords[0].tags.update(ud)
        taggedRecords[1].tags.update(ud)

        # CHANGED: Check for tx contamination using scattered barcodes
        expected_barcode = self.id_to_cs2_barcode.get(taggedRecords[0].tags['bi'])
        if expected_barcode is None:
            raise ValueError(taggedRecords[0].tags['bi'])

        # Try to detect scattered barcode contamination
        tx_umi = None
        
        # First check in R1
        tx_umi = self.extract_scattered_celseq_umi(taggedRecords[0].sequence, expected_barcode)
        
        # If not found in R1, check in reverse complement of R2
        if tx_umi is None:
            rc2 = reverse_complement(taggedRecords[1].sequence)
            tx_umi = self.extract_scattered_celseq_umi(rc2, expected_barcode)

        # If scattered barcode contamination detected
        if tx_umi is not None:
            # Trim read2 down
            taggedRecords[1].sequence, taggedRecords[1].qualities = self.trim_r2(
                taggedRecords[1].sequence, taggedRecords[1].qualities)
            for r in taggedRecords:
                r.tags['dt'] = 'CelSeq'  # Changed from 'VASA' to 'CelSeq'
                r.tags['rx'] = tx_umi
            return taggedRecords
        
        # Fallback: check for regular (non-scattered) contamination as well
        elif expected_barcode in taggedRecords[0].sequence or expected_barcode in reverse_complement(taggedRecords[1].sequence):
            # Regular tx contaminant:
            if expected_barcode in taggedRecords[0].sequence:
                tx_umi = self.extract_vasa_umi(taggedRecords[0].sequence, expected_barcode)
            else:
                rc2 = reverse_complement(taggedRecords[1].sequence)
                tx_umi = self.extract_vasa_umi(rc2, expected_barcode)

            # Trim read2 down
            taggedRecords[1].sequence, taggedRecords[1].qualities = self.trim_r2(
                taggedRecords[1].sequence, taggedRecords[1].qualities)
            for r in taggedRecords:
                r.tags['dt'] = 'CelSeq'
                if tx_umi is not None:
                    r.tags['rx'] = tx_umi
            return taggedRecords
        
        elif 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[0].sequence or 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[1].sequence:
            raise NonMultiplexable('PolyT')
        elif 'AGTCCGACGAT' in taggedRecords[0].sequence[:30] or 'GTTCTACAGT' in taggedRecords[0].sequence[:30] or 'TAATACGACTCACTATAGGG' in taggedRecords[0].sequence:
            for r in taggedRecords:
                r.tags['dt'] = 'CelSeq'  # Changed from 'VASA' to 'CelSeq'
                r.tags['RR'] = 'T7_found'

        return taggedRecords
