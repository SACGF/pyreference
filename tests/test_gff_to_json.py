
import os
import unittest
from bin.pyreference_gff_to_json import parser_factory


class Test(unittest.TestCase):
    base_dir = os.path.join(os.path.dirname(__file__), "..")
    ENSEMBL_GTF_FILENAME = os.path.join(base_dir, "pyreference", "tests", "reference", "ensembl_test.GRCh38.104.gtf")
    REFSEQ_GFF3_FILENAME = os.path.join(base_dir, "pyreference", "tests", "reference", "refseq_test.GRCh38.p13_genomic.109.20210514.gff")
    UCSC_GTF_FILENAME = os.path.join(base_dir, "pyreference", "tests", "reference", "hg19_chrY_300kb_genes.gtf")

    def _test_exon_length(self, data, transcript_id, expected_length):
        transcript = data["transcripts_by_id"][transcript_id]
        exons = transcript["features_by_type"]["exon"]
        length = sum([d["stop"] - d["start"] for d in exons])
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
