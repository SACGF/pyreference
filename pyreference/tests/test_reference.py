from __future__ import print_function, absolute_import

'''
Unit tests for sacgf.genomics.reference.Reference

Created on Jun 20, 2012

@author: dlawrence
'''
import HTSeq
from inspect import getsourcefile
import os
from os.path import abspath
import unittest

from pyreference import Reference
from pyreference.utils.genomics_utils import last_base


class Test(unittest.TestCase):
    def setUp(self):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda:0)))
        reference_dir = os.path.join(this_file_dir, "reference")
        
        genes_json = os.path.join(reference_dir, "hg19_chrY_300kb_genes.gtf.json.gz")
        genome_sequence_fasta = os.path.join(reference_dir, "hg19_chrY_300kb.fa")
        self.reference = Reference(genes_json=genes_json,
                                   genome_sequence_fasta=genome_sequence_fasta)

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
        test_cases = { # Yes, I actually checked these in IGV
                        "PLCXD1":   {"transcript_id"  : "NM_018390_2",
                                     "3p_utr"          : "CGGGACCCTT" },
                        "PPP2R3B":  {"transcript_id"  : "NM_013239", # -'ve strand
                                     "3p_utr"          : "CGCCGCCCGC" },
        }

        for test in test_cases.itervalues():
            transcript = self.reference.transcripts[test["transcript_id"]]

            utr = str(transcript.get_3putr_sequence()).upper()
            self.assertTrue(utr.startswith(test["3p_utr"]), "%s should start with %s" % (utr, test["3p_utr"]))
            
            m_rna = transcript.get_transcript_sequence()
            self.assertTrue(m_rna.find(test["3p_utr"]) > 1)


    def test_mrna_position(self):
        transcript_id = "NR_028057_2"
        start = 142990
        transcript = self.reference.transcripts[transcript_id]

        # Verify that -1 is not on transcript
        try:
            pos = HTSeq.GenomicPosition("chrY", start - 1)
            mpos = transcript.get_transcript_position(pos)
            self.assertTrue(False, "Shouldn't be here!")
        except:
            self.assertTrue(True)

        for i in range(5):
            pos = HTSeq.GenomicPosition("chrY", start + i)
            mpos = transcript.get_transcript_position(pos)
            self.assertEquals(i, mpos)

        # 1st exon length 71 chrs
        exons = transcript.get_features("exon")
        offset = 0
        for e in exons:
            pos = e.iv.start_d_as_pos
            mpos = transcript.get_transcript_position(pos)
            self.assertEquals(offset, mpos, "Exon start position in transcript = sum of previous exons")
            offset += e.iv.length

    def test_mrna_position2(self):
        transcript_id = "NM_018390_2"
        transcript = self.reference.transcripts[transcript_id]

        # Start Codon
        start_codon = transcript.get_features("start_codon")[0]
        mpos = transcript.get_transcript_position(start_codon.iv.start_d_as_pos)
        fiveputr_len = len(transcript.get_5putr_sequence())
        self.assertEquals(mpos, fiveputr_len)

        mpos = transcript.get_transcript_position(last_base(transcript.iv))
        self.assertEquals(mpos, len(transcript.get_transcript_sequence())-1)


    def test_get_feature_extents(self):
        transcript_id = "NR_028057_2"
        transcript = self.reference.transcripts[transcript_id]

        exon1 = transcript.get_features("exon")[0].iv
        exon_intron_span = exon1.copy()
        exon_intron_span.length += 20

        (start, end) = transcript.get_first_and_last_genomic_position_on_transcript(exon_intron_span)
        self.assertTrue(exon1.contains(start))
        self.assertTrue(exon1.contains(end))

        self.assertEquals(start, exon1.start_as_pos)
        self.assertEquals(end, last_base(exon1))

        # Test smaller extent
        exon_sub_span = exon1.copy()
        exon_sub_span.length = 10
        (start, end) = transcript.get_first_and_last_genomic_position_on_transcript(exon_sub_span)
        self.assertTrue(exon_sub_span.contains(start))
        self.assertTrue(exon_sub_span.contains(end))
        
        self.assertEquals(start, exon_sub_span.start_as_pos)
        self.assertEquals(end, last_base(exon_sub_span))

        # Test intron...
        intron = HTSeq.GenomicInterval("chrY", 144043, 144218, '+')
        (start, end) = transcript.get_first_and_last_genomic_position_on_transcript(intron)
        self.assertEqual(start, None)
        self.assertEqual(end, None)

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
        #print "seed = %s" % seed
        self.assertEqual(seed, "CAGTGTTA", "200a seed match")

        miR = self.reference.get_mirna("hsa-miR-200b-3p")
        seed = miR.get_8mer_target()
        #print "seed = %s" % seed
        self.assertEqual(seed, "CAGTATTA", "200b seed match")

        
    def test_promoter(self):
        tests = {"NM_018390_2" : "CCGGGCAGCAGGGAAGATCT", "NM_013239" : "CGCAGTGACGTGAACGCGGG"}
        
        for (tss_id, expected_promoter_sequence) in tests.iteritems():
            transcript = self.reference.transcripts[tss_id]
            promoter_sequence = str(transcript.get_promoter_sequence(10)).upper()
            #print "promoter sequence = %s" % promoter_sequence
            self.assertEqual(expected_promoter_sequence, promoter_sequence, "%s strand promoter sequence" % transcript.iv.strand)

    def test_get_gene_names(self):
        t = self.reference.get_transcript_by_id("NR_028057_2")
        print("%s is_coding: %s" % (t, t.is_coding))
        
        intron = HTSeq.GenomicInterval("chrY", 144043, 144218, '+')
        gene_name = self.reference.get_gene_names(intron)
        self.assertEquals("PLCXD1", gene_name)

    def test_get_gene_region_names(self):
        intron = HTSeq.GenomicInterval("chrY", 144043, 144218, '+')
        region = self.reference.get_region_names(intron)
        self.assertEquals("intron", region)

    def test_gene_transcripts(self):
        plcxd1 = self.reference.get_gene("PLCXD1")
        expected_transcripts = ["NR_028057_2", "NM_018390_2"] # chrY ones have _2

        plcxd1_transcript_ids = {pt.get_id() for pt in plcxd1.transcripts}
        message = "plcxd1 gene contains transcripts %s" % expected_transcripts
        for et_id in expected_transcripts:
            self.assertIn(et_id, plcxd1_transcript_ids, message)

    def test_has_chrom(self):
        self.assertTrue(self.reference.has_chr)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_name']
    unittest.main()