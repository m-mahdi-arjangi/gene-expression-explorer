"""Unit tests for the expression data loader."""

import sys
import tempfile
import unittest
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

from data_loader import load_expression_data


class LoadExpressionDataTests(unittest.TestCase):
    def write_csv(self, content: str) -> Path:
        temp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8")
        temp.write(content)
        temp.close()
        self.addCleanup(Path(temp.name).unlink, missing_ok=True)
        return Path(temp.name)

    def test_loads_valid_table(self):
        path = self.write_csv("gene_id,S1,S2\nGeneA,1,2\nGeneB,3,4\n")
        result = load_expression_data(path)
        self.assertEqual(result.shape, (2, 2))
        self.assertEqual(list(result.index), ["GeneA", "GeneB"])
        self.assertEqual(result.loc["GeneB", "S2"], 4)

    def test_rejects_missing_values(self):
        path = self.write_csv("gene_id,S1,S2\nGeneA,1,\n")
        with self.assertRaisesRegex(ValueError, "numeric and must not be missing"):
            load_expression_data(path)

    def test_rejects_duplicate_gene_ids(self):
        path = self.write_csv("gene_id,S1\nGeneA,1\nGeneA,2\n")
        with self.assertRaisesRegex(ValueError, "unique"):
            load_expression_data(path)

    def test_rejects_missing_gene_id_column(self):
        path = self.write_csv("S1,S2\n1,2\n")
        with self.assertRaisesRegex(ValueError, "gene_id"):
            load_expression_data(path)

    def test_rejects_infinite_values(self):
        path = self.write_csv("gene_id,S1\nGeneA,inf\n")
        with self.assertRaisesRegex(ValueError, "finite"):
            load_expression_data(path)


if __name__ == "__main__":
    unittest.main()
