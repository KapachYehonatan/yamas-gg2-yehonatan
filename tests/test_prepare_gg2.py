import tempfile
import unittest
import zipfile
from pathlib import Path

from yamas.prepare_gg2 import prepare_reference


class PrepareGreengenes2Test(unittest.TestCase):
    def test_converts_qza_artifacts_to_sintax_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sequences = tmp / "sequences.qza"
            taxonomy = tmp / "taxonomy.qza"
            output = tmp / "gg2-sintax.fasta"

            with zipfile.ZipFile(sequences, "w") as artifact:
                artifact.writestr(
                    "uuid/data/dna-sequences.fasta",
                    ">seq1\nACGT\n>seq2\nTGCA\n",
                )
            with zipfile.ZipFile(taxonomy, "w") as artifact:
                artifact.writestr(
                    "uuid/data/taxonomy.tsv",
                    "Feature ID\tTaxon\n"
                    "seq1\td__Bacteria; p__Bacillota; c__Bacilli; o__; f__; g__; s__\n"
                    "seq2\td__Archaea; p__Thermoproteota; c__; o__; f__; g__; s__\n",
                )

            self.assertEqual(prepare_reference(sequences, taxonomy, output), 2)
            self.assertEqual(
                output.read_text(),
                ">seq1;tax=d:Bacteria,p:Bacillota,c:Bacilli;\nACGT\n"
                ">seq2;tax=d:Archaea,p:Thermoproteota;\nTGCA\n",
            )

    def test_removes_partial_output_when_taxonomy_is_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sequences = tmp / "sequences.qza"
            taxonomy = tmp / "taxonomy.qza"
            output = tmp / "gg2-sintax.fasta"

            with zipfile.ZipFile(sequences, "w") as artifact:
                artifact.writestr("uuid/data/dna-sequences.fasta", ">missing\nACGT\n")
            with zipfile.ZipFile(taxonomy, "w") as artifact:
                artifact.writestr(
                    "uuid/data/taxonomy.tsv",
                    "Feature ID\tTaxon\nseq1\td__Bacteria; p__; c__; o__; f__; g__; s__\n",
                )

            with self.assertRaisesRegex(ValueError, "no matching taxonomy"):
                prepare_reference(sequences, taxonomy, output)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
