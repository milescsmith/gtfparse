import unittest

from gtfparse import expand_attribute_strings


class TestRefseqGTF(unittest.TestCase):
    def setUp(self):
        self.attributes_in_quotes = [
            'gene_id "ENSG001"; tag "bogotron"; version "1";',
            'gene_id "ENSG002"; tag "wolfpuppy"; version "2";',
        ]
        self.attributes_without_quotes = [
            "gene_id ENSG001; tag bogotron; version 1;",
            "gene_id ENSG002; tag wolfpuppy; version 2",
        ]
        self.optional_attributes = [
            "gene_id ENSG001; sometimes-present bogotron;",
            "gene_id ENSG002;",
            "gene_id ENSG003; sometimes-present wolfpuppy;",
        ]

    def test_attributes_in_quotes(self):
        parsed_dict = expand_attribute_strings(self.attributes_in_quotes)
        self.assertEqual(list(parsed_dict), ["gene_id", "tag", "version"])
        self.assertEqual(parsed_dict["gene_id"], ["ENSG001", "ENSG002"])
        self.assertEqual(parsed_dict["tag"], ["bogotron", "wolfpuppy"])
        self.assertEqual(parsed_dict["version"], ["1", "2"])


    def test_attributes_without_quotes(self):
        parsed_dict = expand_attribute_strings(self.attributes_without_quotes)
        self.assertEqual(list(sorted(parsed_dict.keys())), ["gene_id", "tag", "version"])
        self.assertEqual(parsed_dict["gene_id"], ["ENSG001", "ENSG002"])
        self.assertEqual(parsed_dict["tag"], ["bogotron", "wolfpuppy"])
        self.assertEqual(parsed_dict["version"], ["1", "2"])


    def test_optional_attributes(self):
        parsed_dict = expand_attribute_strings(self.optional_attributes)
        self.assertEqual(list(parsed_dict), ["gene_id", "sometimes-present"])
        self.assertEqual(parsed_dict["gene_id"], ["ENSG001", "ENSG002", "ENSG003"])
        self.assertEqual(parsed_dict["sometimes-present"], ["bogotron", "", "wolfpuppy"])

if __name__ == "__main__":
    unittest.main()