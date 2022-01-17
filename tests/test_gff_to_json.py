
import os
from inspect import getsourcefile
import unittest
from bin.pyreference_gff_to_json import parser_factory


class Test(unittest.TestCase):
    this_file_dir = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    reference_dir = os.path.join(this_file_dir, "reference")
    ENSEMBL_GTF_FILENAME = os.path.join(reference_dir, "ensembl_test.GRCh38.104.gtf")
    REFSEQ_GFF3_FILENAME = os.path.join(reference_dir, "refseq_test.GRCh38.p13_genomic.109.20210514.gff")
    UCSC_GTF_FILENAME = os.path.join(reference_dir, "hg19_chrY_300kb_genes.gtf")

    def _test_exon_length(self, data, transcript_id, expected_length):
        transcript = data["transcripts_by_id"][transcript_id]
        length = sum([exon[1] - exon[0] for exon in transcript["exons"]])
        self.assertEquals(expected_length, length, "%s exons sum" % transcript_id)

    def test_ucsc_gtf(self):
        parser = parser_factory(gtf=self.UCSC_GTF_FILENAME)
        data = parser.get_data()
        self._test_exon_length(data, "NM_013239", 2426)

    def test_ensembl_gtf(self):
        parser = parser_factory(gtf=self.ENSEMBL_GTF_FILENAME)
        data = parser.get_data()
        self._test_exon_length(data, "ENST00000357654.9", 7088)

    def test_refseq_gff3(self):
        parser = parser_factory(gff3=self.REFSEQ_GFF3_FILENAME)
        data = parser.get_data()
        self._test_exon_length(data, "NM_007294.4", 7088)

    def test_exons_in_genomic_order(self):
        parser = parser_factory(gtf=self.ENSEMBL_GTF_FILENAME)
        data = parser.get_data()
        transcript = data["transcripts_by_id"]["ENST00000357654.9"]
        first_exon = transcript["exons"][0]
        last_exon = transcript["exons"][-1]
        self.assertGreater(last_exon[0], first_exon[0])

        parser = parser_factory(gff3=self.REFSEQ_GFF3_FILENAME)
        data = parser.get_data()
        transcript = data["transcripts_by_id"]["NM_007294.4"]
        first_exon = transcript["exons"][0]
        last_exon = transcript["exons"][-1]
        self.assertGreater(last_exon[0], first_exon[0])


