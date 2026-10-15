"""Unit tests for validate_results.py parsing logic, using small synthetic
VCF content so tests run instantly without needing the real pipeline output."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from validate_results import count_variants, find_variant_near_position

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

def make_vcf(tmp_path, variant_lines):
    content = VCF_HEADER + "".join(variant_lines)
    f = tmp_path / "test.vcf"
    f.write_text(content)
    return str(f)

def test_count_variants_ignores_header():
    pass  # covered by the parametrized tests below via tmp_path

def test_count_variants_empty_vcf(tmp_path):
    path = make_vcf(tmp_path, [])
    assert count_variants(path) == 0

def test_count_variants_counts_only_data_lines(tmp_path):
    lines = [
        "chr1\t100\t.\tA\tG\t.\tPASS\t.\n",
        "chr1\t200\t.\tC\tT\t.\tPASS\t.\n",
        "chr1\t300\t.\tG\tA\t.\tPASS\t.\n",
    ]
    path = make_vcf(tmp_path, lines)
    assert count_variants(path) == 3

def test_find_variant_near_position_exact_match(tmp_path):
    lines = ["chr1\t2753768\t.\tC\tT\t.\tPASS\t.\n"]
    path = make_vcf(tmp_path, lines)
    result = find_variant_near_position(path, 2753768)
    assert result is not None
    assert "2753768" in result

def test_find_variant_near_position_within_window(tmp_path):
    lines = ["chr1\t2753770\t.\tC\tT\t.\tPASS\t.\n"]  # 2bp off, within default window of 5
    path = make_vcf(tmp_path, lines)
    result = find_variant_near_position(path, 2753768, window=5)
    assert result is not None

def test_find_variant_near_position_outside_window_returns_none(tmp_path):
    lines = ["chr1\t2753900\t.\tC\tT\t.\tPASS\t.\n"]  # far outside window
    path = make_vcf(tmp_path, lines)
    result = find_variant_near_position(path, 2753768, window=5)
    assert result is None

def test_find_variant_near_position_no_variants(tmp_path):
    path = make_vcf(tmp_path, [])
    result = find_variant_near_position(path, 2753768)
    assert result is None

def test_find_variant_ignores_malformed_lines(tmp_path):
    lines = [
        "chr1\tNOTANUMBER\t.\tC\tT\t.\tPASS\t.\n",  # malformed POS, should be skipped not crash
        "chr1\t2753768\t.\tC\tT\t.\tPASS\t.\n",
    ]
    path = make_vcf(tmp_path, lines)
    result = find_variant_near_position(path, 2753768)
    assert result is not None
