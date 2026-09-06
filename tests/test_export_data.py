import tempfile
import unittest
from pathlib import Path

from yamas.export_data import _concatenate_files


class ConcatenateFilesTest(unittest.TestCase):
    def test_concatenates_many_files_in_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            inputs = []
            for index in range(1100):
                path = tmp / f"{index}.fasta"
                path.write_bytes(f">seq{index}\nACGT\n".encode())
                inputs.append(path)

            output = tmp / "all.fasta"
            _concatenate_files(inputs, output)

            data = output.read_bytes()
            self.assertTrue(data.startswith(b">seq0\nACGT\n"))
            self.assertTrue(data.endswith(b">seq1099\nACGT\n"))
            self.assertEqual(data.count(b">seq"), 1100)


if __name__ == "__main__":
    unittest.main()
