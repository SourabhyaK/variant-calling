"""
Validates the variant-calling pipeline's output against two independent,
published sources of ground truth:

1. Total variant counts per sample (Data Carpentry "Wrangling Genomics"
   lesson, using this exact dataset) - a coarse sanity check.
2. Whether the specific, real causal mutation for the hypermutator
   phenotype (mutS Q606*, REL606 position 2,753,768; Favate et al.,
   PMC7523724, Table 1) is actually detected in the 50,000-generation
   sample's VCF, and absent from the earlier timepoints.

This second check is the real test of whether the pipeline works
correctly - matching a total count could happen by coincidence with a
flawed pipeline (wrong variants cancelling out to the right number);
finding the exact published mutation at the exact published position
cannot.
"""
import sys
import os

# Published benchmarks (Data Carpentry "Wrangling Genomics" lesson, same dataset)
EXPECTED_COUNTS = {
    "SRR2589044": 10,   # 5,000 generations
    "SRR2584863": 25,   # 15,000 generations
    "SRR2584866": 766,  # 50,000 generations
}

# The real, published causal mutation (Favate et al., PMC7523724, Table 1)
MUTS_POSITION = 2753768
MUTS_WINDOW = 5  # allow a small window for indel-shifted coordinate reporting

# Samples taken before generation 34,750 (when the mutS mutation arose)
# should NOT carry it; samples after should.
SAMPLE_GENERATIONS = {
    "SRR2589044": 5000,
    "SRR2584863": 15000,
    "SRR2584866": 50000,
}
MUTS_APPEARANCE_GENERATION = 34750


def count_variants(vcf_path):
    """Count non-header lines (actual variant calls) in a VCF file."""
    count = 0
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


def find_variant_near_position(vcf_path, position, window=MUTS_WINDOW):
    """Return the matching VCF line if a variant exists within `window` bp
    of `position`, else None. VCF POS is column 2 (1-indexed, tab-separated)."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            if abs(pos - position) <= window:
                return line.strip()
    return None


def validate(results_dir):
    all_passed = True
    print("=" * 70)
    print("VALIDATION 1: Variant counts vs. published benchmark")
    print("=" * 70)
    for sample, expected in EXPECTED_COUNTS.items():
        vcf_path = os.path.join(results_dir, "vcf", f"{sample}_final_variants.vcf")
        if not os.path.exists(vcf_path):
            print(f"  {sample}: MISSING ({vcf_path} not found)")
            all_passed = False
            continue
        actual = count_variants(vcf_path)
        status = "PASS" if actual == expected else "CHECK"
        if actual != expected:
            all_passed = False
        print(f"  {sample}: expected {expected}, got {actual}  [{status}]")

    print()
    print("=" * 70)
    print("VALIDATION 2: Does the pipeline detect the real causal mutation?")
    print(f"  (mutS Q606*, REL606 position {MUTS_POSITION}, arose at gen "
          f"{MUTS_APPEARANCE_GENERATION})")
    print("=" * 70)
    for sample, generation in SAMPLE_GENERATIONS.items():
        vcf_path = os.path.join(results_dir, "vcf", f"{sample}_final_variants.vcf")
        if not os.path.exists(vcf_path):
            print(f"  {sample}: MISSING")
            all_passed = False
            continue
        hit = find_variant_near_position(vcf_path, MUTS_POSITION)
        should_have_it = generation > MUTS_APPEARANCE_GENERATION
        has_it = hit is not None

        if should_have_it == has_it:
            status = "PASS"
        else:
            status = "CHECK"
            all_passed = False

        expectation = "should carry it" if should_have_it else "should NOT carry it"
        found = f"found: {hit}" if hit else "not found"
        print(f"  {sample} (gen {generation}, {expectation}): {found}  [{status}]")

    print()
    if all_passed:
        print("ALL VALIDATIONS PASSED - pipeline output matches published "
              "ground truth at both the aggregate and the single-mutation level.")
    else:
        print("Some validations did not match expectations - see CHECK lines "
              "above. A small count mismatch can be normal (aligner/caller "
              "version differences); a missing causal-mutation detection is "
              "more concerning and worth investigating.")
    return all_passed


if __name__ == "__main__":
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "results"
    success = validate(results_dir)
    sys.exit(0 if success else 1)
