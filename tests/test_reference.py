"""
Unit tests for sacgf.genomics.reference.Reference

Created on Jun 20, 2012

@author: dlawrence
"""

from __future__ import print_function, absolute_import

import HTSeq
from inspect import getsourcefile
import os
from os.path import abspath
import six
import unittest

from pyreference import Reference, settings


class Test(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None  # Show all of diffs on error

        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        reference_dir = os.path.join(this_file_dir, "reference")

        genes_json = os.path.join(reference_dir, "hg19_chrY_300kb_genes.gtf.cdot.json.gz")
        genome_sequence_fasta = os.path.join(reference_dir, "hg19_chrY_300kb.fa")
        mature_mir_sequence_fasta = os.path.join(reference_dir, "mature_200ab_only.fa")

        self.reference = Reference(load_config_file=False,
                                   genome_accession='GRCh37',
                                   genes_json=genes_json,
                                   genome_sequence_fasta=genome_sequence_fasta,
                                   mature_mir_sequence_fasta=mature_mir_sequence_fasta)

    def tearDown(self):
        pass

    def test_regions(self):
        transcript_id = "NM_018390_2"
        transcript = self.reference.transcripts[transcript_id]

        for cds in transcript.get_features("CDS"):
            coding_regions = self.reference.get_region_names(cds.iv)
            self.assertEqual("coding", coding_regions)

        gene_regions = self.reference.get_regions_array(transcript.iv)
        self.assertTrue("intron" in gene_regions)
        self.assertTrue("coding" in gene_regions)
        self.assertTrue("5PUTR" in gene_regions)
        self.assertTrue("3PUTR" in gene_regions)

    def test_region_loading(self):
        # TODO: This call fails! Fix it!
        #        self.assertTrue(self.reference.regions is not None, "Fix me!")

        # Load AFTER call to transcripts works        
        self.reference.transcripts
        self.assertTrue(self.reference.regions is not None)

    def test_genes(self):
        test_cases = {  # Yes, I actually checked these in IGV
            "PLCXD1": {"transcript_id": "NM_018390_2",
                       "3p_utr": "CGGGACCCTT"},
            "PPP2R3B": {"transcript_id": "NM_013239",  # -'ve strand
                        "3p_utr": "CGCCGCCCGC"},
        }

        for test in six.itervalues(test_cases):
            transcript = self.reference.transcripts[test["transcript_id"]]

            utr = str(transcript.get_3putr_sequence()).upper()
            self.assertTrue(utr.startswith(test["3p_utr"]), "%s should start with %s" % (utr, test["3p_utr"]))

            m_rna = transcript.get_transcript_sequence()
            self.assertTrue(m_rna.find(test["3p_utr"]) > 1)

    def test_get_transcript_length(self):
        transcript_id = "NM_018390_2"
        transcript = self.reference.transcripts[transcript_id]

        seq = transcript.get_transcript_sequence()
        #        print "**** transcript_id = %s - sequence = %s" % (transcript_id, seq)
        sequence_length = len(seq)
        length = transcript.get_transcript_length()
        self.assertEqual(length, sequence_length)

    def test_mirna(self):
        miR = self.reference.get_mirna("hsa-miR-200a-3p")
        seed = miR.get_8mer_target()
        # print "seed = %s" % seed
        self.assertEqual(seed, "CAGTGTTA", "200a seed match")

        miR = self.reference.get_mirna("hsa-miR-200b-3p")
        seed = miR.get_8mer_target()
        # print "seed = %s" % seed
        self.assertEqual(seed, "CAGTATTA", "200b seed match")

    def test_promoter(self):
        tests = {"NM_018390_2": "CCGGGCAGCAGGGAAGATCT", "NM_013239": "CGCAGTGACGTGAACGCGGG"}

        for (tss_id, expected_promoter_sequence) in six.iteritems(tests):
            transcript = self.reference.transcripts[tss_id]
            promoter_sequence = str(transcript.get_promoter_sequence(10)).upper()
            # print "promoter sequence = %s" % promoter_sequence
            self.assertEqual(expected_promoter_sequence, promoter_sequence,
                             "%s strand promoter sequence" % transcript.iv.strand)

    def test_get_gene_names(self):
        intron = HTSeq.GenomicInterval("chrY", 144043, 144218, '+')
        gene_name = self.reference.get_gene_names(intron)
        self.assertEqual("PLCXD1", gene_name)

    def test_get_gene_region_names(self):
        intron = HTSeq.GenomicInterval("chrY", 144043, 144218, '+')
        region = self.reference.get_region_names(intron)
        self.assertEqual("intron", region)

    def test_gene_transcripts(self):
        plcxd1 = self.reference.get_gene("PLCXD1")
        expected_transcripts = ["NR_028057_2", "NM_018390_2"]  # chrY ones have _2

        plcxd1_transcript_ids = {pt.get_id() for pt in plcxd1.transcripts}
        message = "plcxd1 gene contains transcripts %s" % expected_transcripts
        for et_id in expected_transcripts:
            self.assertIn(et_id, plcxd1_transcript_ids, message)

    def test_has_chrom(self):
        self.assertTrue(self.reference.has_chr)

    def test_stranded_order(self):
        transcript = self.reference.transcripts["NM_018390_2"]  # + strand
        exons = transcript.get_features_in_stranded_order("exon")

        first_start = exons[0][settings.START]
        last_start = exons[-1][settings.END]
        self.assertGreater(last_start, first_start)  # genomic first comes first

        transcript = self.reference.transcripts["NM_013239"]  # - strand
        exons = transcript.get_features_in_stranded_order("exon")

        first_start = exons[0][settings.START]
        last_start = exons[-1][settings.END]
        self.assertGreater(first_start, last_start)  # - strand, genomic first comes last

    def test_get_features_positive_strand(self):
        """ We re-build features now from exons - test this matches GTF """
        gtf_cds = [
            # From GTF:
            # grep CDS.*NM_018390_2 tests/reference/hg19_chrY_300kb_genes.gtf | cut -d$'\t' -f 1,4,5,7
            ('chrY', 150855, 150981, '+'),
            ('chrY', 155400, 155536, '+'),
            ('chrY', 157315, 157443, '+'),
            ('chrY', 158166, 158321, '+'),
            ('chrY', 159702, 159885, '+'),
            ('chrY', 165764, 165999, '+'),
        ]
        expected_cds = []
        for (chrom, start, stop, strand) in gtf_cds:
            # Adjust start as GTF is 1-based
            expected_cds.append({"chrom": chrom, "start": start-1, "stop": stop, "strand": strand})

        transcript = self.reference.transcripts["NM_018390_2"]
        print(transcript._dict)

        cds_features = transcript.get_features_in_stranded_order("CDS")
        self.assertEqual(cds_features, expected_cds)

        expected_start_codon = [{"chrom": "chrY", 'start': 150854, 'stop': 150857, "strand": "+"}]
        start_codon = transcript.get_features_in_stranded_order("start_codon")
        self.assertEqual(start_codon, expected_start_codon)

        expected_stop_codon = [{"chrom": "chrY", 'start': 165999, 'stop': 166002, "strand": "+"}]
        stop_codon = transcript.get_features_in_stranded_order("stop_codon")
        self.assertEqual(stop_codon, expected_stop_codon)

    def test_get_features_negative_strand(self):
        """ We re-build features now from exons - test this matches GTF """
        gtf_cds = [
            # From GTF:
            # grep CDS.*NM_013239 tests/reference/hg19_chrY_300kb_genes.gtf | cut -d$'\t' -f 1,4,5,7
            ('chrY', 245105, 245252, '-'),
            ('chrY', 249339, 249445, '-'),
            ('chrY', 249513, 249631, '-'),
            ('chrY', 251500, 251675, '-'),
            ('chrY', 252042, 252131, '-'),
            ('chrY', 252618, 252666, '-'),
            ('chrY', 256251, 256407, '-'),
            ('chrY', 256909, 256995, '-'),
            ('chrY', 257436, 257510, '-'),
            ('chrY', 257969, 258071, '-'),
            ('chrY', 258325, 258428, '-'),
            ('chrY', 272140, 272325, '-'),
            ('chrY', 297103, 297426, '-'),
        ]
        # Reverse as NM_013239 is -'ve strand
        expected_cds = []
        for (chrom, start, stop, strand) in reversed(gtf_cds):
            # Adjust start as GTF is 1-based
            expected_cds.append({"chrom": chrom, "start": start-1, "stop": stop, "strand": strand})

        transcript = self.reference.transcripts["NM_013239"]
        print(transcript._dict)

        cds_features = transcript.get_features_in_stranded_order("CDS")
        self.assertEqual(cds_features, expected_cds)

        expected_start_codon = [{"chrom": "chrY", "start": 297423, "stop": 297426, "strand": "-"}]
        start_codon = transcript.get_features_in_stranded_order("start_codon")
        self.assertEqual(start_codon, expected_start_codon)

        expected_stop_codon = [{"chrom": "chrY", "start": 245101, "stop": 245104, "strand": "-"}]
        stop_codon = transcript.get_features_in_stranded_order("stop_codon")
        self.assertEqual(stop_codon, expected_stop_codon)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_name']
    unittest.main()
