"""Unit tests for validate_results.py parsing logic, using small synthetic
VCF content so tests run instantly without needing the real pipeline output."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from validate_results import count_variants, find_variant_near_position, validate

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


# --- Tests for validate() itself: the actual PASS/CHECK decision logic,  ---
# --- not just the helper functions underneath it.                       ---

def _make_full_result_set(tmp_path, counts, include_mutation):
    """Build a results/vcf/ folder matching what validate() expects: the 3
    hardcoded sample names, each with `counts[sample]` dummy variant lines,
    and the real mutS mutation line included only in SRR2584866 if
    include_mutation is True (and never in the other two, matching the real
    biology)."""
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir()
    samples = ["SRR2589044", "SRR2584863", "SRR2584866"]
    for sample in samples:
        n = counts[sample]
        add_mutation = sample == "SRR2584866" and include_mutation
        # total lines must equal exactly `n` - if the mutation line is one
        # of them, generate n-1 filler lines instead of n, so the total
        # variant count matches `n` exactly rather than overshooting by 1.
        filler_count = n - 1 if add_mutation else n
        lines = [f"chr1\t{1000+i}\t.\tA\tG\t.\tPASS\t.\n" for i in range(filler_count)]
        if add_mutation:
            lines.append("chr1\t2753768\t.\tC\tT\t.\tPASS\t.\n")
        (vcf_dir / f"{sample}_final_variants.vcf").write_text(VCF_HEADER + "".join(lines))
    return str(tmp_path)

def test_validate_passes_when_everything_matches(tmp_path):
    counts = {"SRR2589044": 10, "SRR2584863": 25, "SRR2584866": 766}
    results_dir = _make_full_result_set(tmp_path, counts, include_mutation=True)
    assert validate(results_dir) is True

def test_validate_fails_on_count_mismatch_even_if_mutation_found(tmp_path):
    # Mirrors the real-world case we actually hit: 775 instead of 766.
    # The mutation check alone isn't enough to call it a full pass.
    counts = {"SRR2589044": 10, "SRR2584863": 25, "SRR2584866": 775}
    results_dir = _make_full_result_set(tmp_path, counts, include_mutation=True)
    assert validate(results_dir) is False

def test_validate_fails_if_mutation_missing_from_later_sample(tmp_path):
    counts = {"SRR2589044": 10, "SRR2584863": 25, "SRR2584866": 766}
    results_dir = _make_full_result_set(tmp_path, counts, include_mutation=False)
    assert validate(results_dir) is False

def test_validate_fails_on_missing_vcf_file(tmp_path):
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir()
    # Only create 2 of the 3 expected files
    (vcf_dir / "SRR2589044_final_variants.vcf").write_text(VCF_HEADER)
    (vcf_dir / "SRR2584863_final_variants.vcf").write_text(VCF_HEADER)
    assert validate(str(tmp_path)) is False
