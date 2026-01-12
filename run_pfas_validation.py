#!/usr/bin/env python3
"""
PFAScreen Validation - Using PFAS_feature_prioritization directly
Converts TSV data to PFAScreen format and runs PFAS_feature_prioritization
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys
import re
from typing import Dict

# Import the actual PFAScreen function
from PFAS_feature_prioritization import PFAS_feature_prioritization


def parse_formula(formula: str) -> Dict[str, int]:
    """Parse chemical formula to get element counts"""
    if pd.isna(formula) or not formula:
        return {}
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    elements = {}
    for element, count in matches:
        count = int(count) if count else 1
        elements[element] = elements.get(element, 0) + count
    return elements


def estimate_isotope_areas(precursor_intens: float, formula: str) -> tuple:
    """
    Estimate isotope peak areas from formula
    Returns (C12_area, C13_area, C14_area)
    """
    elements = parse_formula(formula)
    C = elements.get('C', 0)

    if C == 0:
        return precursor_intens, 0.0, 0.0

    # Natural abundance calculations
    # C13/C12 ratio ≈ 1.1% per carbon
    # C13 intensity = M0 * C * 0.011
    # C14 intensity ≈ M0 * C*(C-1)/2 * 0.011^2

    C12_area = precursor_intens
    C13_area = precursor_intens * C * 0.011145  # Exact natural abundance
    C14_area = precursor_intens * C * (C - 1) / 2 * 0.011145 ** 2

    return C12_area, C13_area, C14_area


def convert_tsv_to_pfascreen_format(df: pd.DataFrame) -> tuple:
    """
    Convert TSV format to PFAScreen's expected format

    Returns:
        (Df_FeatureData, Df_MS2RawData, idx_in_features)
    """
    print("Converting TSV data to PFAScreen format...")

    # Create Df_FeatureData (one row per compound)
    feature_data = []
    ms2_data = []
    idx_in_features = []

    for idx, row in df.iterrows():
        if idx % 10000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(df)} rows...")

        # Parse MS2 spectrum
        try:
            mzs = [float(x) for x in str(row['mzs']).split(',')]
            intensities = [float(x) for x in str(row['intensities']).split(',')]
        except:
            continue

        if len(mzs) == 0:
            continue

        # Get precursor info
        precursor_mz = row.get('precursor_mz', np.nan)
        formula = str(row.get('precursor_formula', ''))

        # Estimate isotope areas
        # Use max fragment intensity as proxy for precursor intensity
        precursor_intens = max(intensities) if intensities else 100.0
        C12_area, C13_area, C14_area = estimate_isotope_areas(precursor_intens, formula)

        # Create feature data entry
        feature_data.append({
            'mz': precursor_mz,
            'mz_area': C12_area,
            'mz+1_area': C13_area,
            'mz+2_area': C14_area,
            'rt': 60.0,  # Dummy RT (1 minute) since we don't have real RT
            'original_index': idx
        })

        # Create MS2 data entry
        ms2_data.append({
            'mz': precursor_mz,
            'rt': 60.0,
            'intens': precursor_intens,
            'ms2_spec_mz': np.array(mzs),
            'ms2_spec_intens': np.array(intensities),
            'idx_ms1': len(feature_data) - 1  # Maps to feature data index
        })

        # Track which MS1 features have MS2
        idx_in_features.append(len(feature_data) - 1)

    Df_FeatureData = pd.DataFrame(feature_data)
    Df_MS2RawData = pd.DataFrame(ms2_data)
    idx_in_features = np.array(idx_in_features)

    print(f"  Processed {len(df)}/{len(df)} rows...Done!")
    print(f"Created {len(Df_FeatureData)} features with {len(Df_MS2RawData)} MS2 spectra")

    return Df_FeatureData, Df_MS2RawData, idx_in_features


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Compute classification metrics"""
    tp = np.sum((y_true == True) & (y_pred == True))
    fp = np.sum((y_true == False) & (y_pred == True))
    tn = np.sum((y_true == False) & (y_pred == False))
    fn = np.sum((y_true == True) & (y_pred == False))

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    accuracy = (tp + tn) / (tp + fp + tn + fn) if (tp + fp + tn + fn) > 0 else 0.0

    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'accuracy': accuracy,
        'tp': tp,
        'fp': fp,
        'tn': tn,
        'fn': fn,
        'total': len(y_true)
    }


def main():
    parser = argparse.ArgumentParser(description='PFAScreen Validation using PFAS_feature_prioritization')
    parser.add_argument('--input', type=str,
                       default='~/Downloads/merged_massspec_nist20_nist_new_env_pfas_with_fold.tsv',
                       help='Input TSV file path')
    parser.add_argument('--fold', type=str, default='val',
                       help='Fold to evaluate (train/val/all)')
    parser.add_argument('--sample-size', type=int, default=None,
                       help='Sample size for testing (optional)')
    parser.add_argument('--output', type=str, default='pfas_results',
                       help='Output folder name')

    # PFAScreen parameters (matching PFAS_feature_prioritization defaults)
    parser.add_argument('--diffs', type=str, default='C2F4,C2F4,HF',
                       help='Comma-separated fragment differences (default: C2F4,C2F4,HF)')
    parser.add_argument('--number-of-fragments', type=int, default=1,
                       help='Minimum number of fragment differences (default: 1)')
    parser.add_argument('--mass-tolerance', type=float, default=0.002,
                       help='Mass tolerance for fragment matching in Da (default: 0.002)')
    parser.add_argument('--intensity-threshold', type=float, default=5,
                       help='Intensity threshold for MS2 filtering (default: 5)')
    parser.add_argument('--mdc-range', type=str, default='-0.5,0.5',
                       help='MD/C range comma-separated min,max (default: -0.5,0.5)')
    parser.add_argument('--mc-range', type=str, default='0,inf',
                       help='m/C range comma-separated min,max (default: 0,inf)')
    parser.add_argument('--diffs-kmd', type=str, default='CF2',
                       help='Comma-separated KMD differences (default: CF2)')
    parser.add_argument('--mass-tolerance-kmd', type=float, default=0.002,
                       help='Mass tolerance for KMD analysis in Da (default: 0.002)')
    parser.add_argument('--n-homologues', type=int, default=3,
                       help='Minimum homologues for KMD (default: 3)')
    parser.add_argument('--adducts', type=int, default=1,
                       help='Adducts for suspect screening (default: 1)')
    parser.add_argument('--md-range', type=str, default='-0.5,0.5',
                       help='MD range comma-separated min,max (default: -0.5,0.5)')

    args = parser.parse_args()

    # Create output folder and plots subfolder
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    plots_dir = os.path.join(args.output, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    # Load data
    print(f"Loading data from {args.input}...")
    input_path = os.path.expanduser(args.input)
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    print(f"Loaded {len(df)} spectra")

    # Filter by fold
    if args.fold != 'all':
        df = df[df['fold'] == args.fold]
        print(f"Filtered to {len(df)} spectra in {args.fold} fold")

    # Sample if requested
    if args.sample_size:
        df = df.sample(n=min(args.sample_size, len(df)), random_state=42)
        print(f"Sampled {len(df)} spectra for testing")

    # Keep ground truth labels
    ground_truth = df['is_PFAS'].values
    original_indices = df.index.values

    # Convert to PFAScreen format
    Df_FeatureData, Df_MS2RawData, idx_in_features = convert_tsv_to_pfascreen_format(df)

    # Parse parameters
    diffs = args.diffs.split(',')
    diffs_kmd = args.diffs_kmd.split(',')
    mdc_range = [float(x) if x != 'inf' else float('inf') for x in args.mdc_range.split(',')]
    mc_range = [float(x) if x != 'inf' else float('inf') for x in args.mc_range.split(',')]
    md_range = [float(x) if x != 'inf' else float('inf') for x in args.md_range.split(',')]

    print("\n" + "="*60)
    print("Running PFAS_feature_prioritization")
    print("="*60)
    print(f"Parameters:")
    print(f"  Fragment differences (diffs): {diffs}")
    print(f"  Number of fragments: {args.number_of_fragments}")
    print(f"  Mass tolerance FindPFAS: {args.mass_tolerance} Da")
    print(f"  Intensity threshold: {args.intensity_threshold}")
    print(f"  KMD differences: {diffs_kmd}")
    print(f"  Mass tolerance KMD: {args.mass_tolerance_kmd} Da")
    print(f"  Min homologues: {args.n_homologues}")
    print(f"  Adducts: {args.adducts}")
    print(f"  MD/C range: {mdc_range}")
    print(f"  m/C range: {mc_range}")
    print(f"  MD range: {md_range}")
    print()

    # Call actual PFAS_feature_prioritization function
    try:
        Df_FeatureData_Final, Df_FindPFAS = PFAS_feature_prioritization(
            Df_FeatureData,
            Df_MS2RawData,
            idx_in_features,
            Results_folder=args.output,
            diffs=diffs,
            number_of_fragments=args.number_of_fragments,
            mass_tolerance_FindPFAS=args.mass_tolerance,
            intensity_threshold=args.intensity_threshold,
            diffs_KMD=diffs_kmd,
            mass_tolerance_KMD=args.mass_tolerance_kmd,
            n_homologues=args.n_homologues,
            adducts=args.adducts,
            tol_suspect=args.mass_tolerance,
            save_MSMS_spectra=False,
            mC_range=mc_range,
            MDC_range=mdc_range,
            MD_range=md_range
        )

        print("\n" + "="*60)
        print("PFAS Detection Results")
        print("="*60)

        # Map results back to original data
        # Features are PFAS if they have:
        # 1. Passed filters (in Df_FeatureData_Final)
        # 2. Have PFAS indicators: n_diffs > 0 OR n_dias > 0 OR suspect hits OR HS detected

        predictions = []
        for i in range(len(Df_FeatureData)):
            if i in Df_FeatureData_Final.index:
                # Check PFAS indicators
                row = Df_FeatureData_Final.loc[i]
                n_diffs = row.get('n_diffs', 0)
                n_dias = row.get('n_dias', 0)
                has_suspect = row.get('compound_names', '') != ''
                has_hs = row.get('unique_homologues', False)

                # Classify as PFAS if any indicator is positive
                is_pfas = (n_diffs > 0) or (n_dias > 0) or has_suspect or has_hs
                predictions.append(is_pfas)
            else:
                # Filtered out = not PFAS
                predictions.append(False)

        predictions = np.array(predictions)

        # Compute metrics
        metrics = compute_metrics(ground_truth, predictions)

        print(f"\nOverall Metrics ({args.fold} fold):")
        print(f"  Precision: {metrics['precision']:.4f}")
        print(f"  Recall:    {metrics['recall']:.4f}")
        print(f"  F1 Score:  {metrics['f1']:.4f}")
        print(f"  Accuracy:  {metrics['accuracy']:.4f}")
        print(f"\nConfusion Matrix:")
        print(f"  TP: {metrics['tp']:6d}  FP: {metrics['fp']:6d}")
        print(f"  FN: {metrics['fn']:6d}  TN: {metrics['tn']:6d}")
        print(f"  Total: {metrics['total']}")

        # Save metrics
        metrics_file = os.path.join(args.output, 'validation_metrics.txt')
        with open(metrics_file, 'w') as f:
            f.write("PFAS Detection Validation Results\n")
            f.write("="*60 + "\n\n")
            f.write(f"Fold: {args.fold}\n")
            f.write(f"Total samples: {metrics['total']}\n\n")
            f.write(f"Precision: {metrics['precision']:.4f}\n")
            f.write(f"Recall:    {metrics['recall']:.4f}\n")
            f.write(f"F1 Score:  {metrics['f1']:.4f}\n")
            f.write(f"Accuracy:  {metrics['accuracy']:.4f}\n\n")
            f.write("Confusion Matrix:\n")
            f.write(f"  TP: {metrics['tp']:6d}  FP: {metrics['fp']:6d}\n")
            f.write(f"  FN: {metrics['fn']:6d}  TN: {metrics['tn']:6d}\n")

        print(f"\nResults saved to: {args.output}/")
        print(f"  - Results_{args.output}.xlsx (detailed feature table)")
        print(f"  - validation_metrics.txt (performance metrics)")

    except Exception as e:
        print(f"\nError running PFAS_feature_prioritization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
