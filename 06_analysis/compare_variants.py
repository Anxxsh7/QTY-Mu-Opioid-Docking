#!/usr/bin/env python3
"""
Step 6: Comparative Analysis - QTY vs WT Control
================================================
Compares binding results between QTY mutant and WT control.
Generates final ranking and assessment.

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import json
from pathlib import Path
from typing import List, Dict, Optional
import warnings

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def load_docking_results(results_dir: Path) -> Dict:
    """Load docking analysis results."""
    
    analysis_file = results_dir / 'docking_analysis.json'
    if analysis_file.exists():
        with open(analysis_file) as f:
            return json.load(f)
    return None


def load_md_results(results_dir: Path) -> Dict:
    """Load MD simulation results."""
    
    summary_file = results_dir / 'md_summary.json'
    if summary_file.exists():
        with open(summary_file) as f:
            return json.load(f)
    return None


def calculate_binding_score(docking_results: Dict, md_results: Dict) -> Dict:
    """Calculate composite binding score."""
    
    score = {
        'docking': {},
        'md': {},
        'composite': None,
        'assessment': None
    }
    
    # Docking metrics
    if docking_results:
        scores = docking_results.get('scores', {})
        score['docking'] = {
            'mean_affinity': scores.get('mean'),
            'best_affinity': scores.get('min'),
            'salt_bridge_fraction': docking_results.get('asp_salt_bridge', {}).get('fraction', 0)
        }
    
    # MD metrics
    if md_results:
        all_rmsd = []
        all_sb = []
        
        for result in md_results.get('results', []):
            analysis = result.get('analysis', {})
            
            rmsd_data = analysis.get('ligand_rmsd', {})
            if 'final_angstrom' in rmsd_data:
                all_rmsd.append(rmsd_data['final_angstrom'])
            
            sb_data = analysis.get('asp_salt_bridge', {})
            if 'intact_fraction' in sb_data:
                all_sb.append(sb_data['intact_fraction'])
        
        score['md'] = {
            'mean_final_rmsd': np.mean(all_rmsd) if all_rmsd else None,
            'rmsd_std': np.std(all_rmsd) if all_rmsd else None,
            'mean_salt_bridge_intact': np.mean(all_sb) if all_sb else None,
            'n_replicates': len(all_rmsd)
        }
    
    # Composite assessment
    issues = []
    good_signs = []
    
    # Check docking
    if score['docking'].get('salt_bridge_fraction', 0) < 0.5:
        issues.append("Salt bridge lost in >50% of docked poses")
    else:
        good_signs.append("Salt bridge preserved in docking")
    
    if score['docking'].get('mean_affinity'):
        if score['docking']['mean_affinity'] > -6.0:
            issues.append(f"Weak docking scores ({score['docking']['mean_affinity']:.1f} kcal/mol)")
        else:
            good_signs.append(f"Good docking scores ({score['docking']['mean_affinity']:.1f} kcal/mol)")
    
    # Check MD
    if score['md'].get('mean_final_rmsd'):
        if score['md']['mean_final_rmsd'] > 3.0:
            issues.append(f"Ligand drifts from pocket (RMSD={score['md']['mean_final_rmsd']:.1f}Å)")
        else:
            good_signs.append(f"Ligand stable in pocket (RMSD={score['md']['mean_final_rmsd']:.1f}Å)")
    
    if score['md'].get('mean_salt_bridge_intact'):
        if score['md']['mean_salt_bridge_intact'] < 0.5:
            issues.append("Salt bridge breaks during MD")
        else:
            good_signs.append(f"Salt bridge maintained ({score['md']['mean_salt_bridge_intact']*100:.0f}%)")
    
    # Overall assessment
    if len(issues) == 0:
        score['assessment'] = 'LIKELY_BINDER'
        score['confidence'] = 'HIGH'
    elif len(issues) <= 1 and len(good_signs) >= 2:
        score['assessment'] = 'POSSIBLE_BINDER'
        score['confidence'] = 'MEDIUM'
    else:
        score['assessment'] = 'UNLIKELY_BINDER'
        score['confidence'] = 'HIGH' if len(issues) >= 3 else 'MEDIUM'
    
    score['issues'] = issues
    score['good_signs'] = good_signs
    
    return score


def compare_variants(variants: Dict[str, Dict], output_dir: Path) -> Dict:
    """Compare multiple variants (QTY vs WT)."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    comparison = {
        'variants': {},
        'ranking': []
    }
    
    for name, dirs in variants.items():
        docking = load_docking_results(dirs.get('docking', Path()))
        md = load_md_results(dirs.get('md', Path()))
        
        score = calculate_binding_score(docking, md)
        score['name'] = name
        comparison['variants'][name] = score
    
    # Rank variants
    def rank_key(name):
        v = comparison['variants'][name]
        
        # Primary: assessment
        assessment_score = {
            'LIKELY_BINDER': 0,
            'POSSIBLE_BINDER': 1,
            'UNLIKELY_BINDER': 2
        }.get(v.get('assessment'), 3)
        
        # Secondary: docking score
        docking_score = v['docking'].get('mean_affinity', 0) or 0
        
        # Tertiary: MD stability
        md_rmsd = v['md'].get('mean_final_rmsd', 10) or 10
        
        return (assessment_score, docking_score, md_rmsd)
    
    comparison['ranking'] = sorted(
        comparison['variants'].keys(),
        key=rank_key
    )
    
    return comparison


def generate_report(comparison: Dict, output_dir: Path) -> str:
    """Generate final comparison report."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    report_lines = [
        "=" * 70,
        "QTY μOR - FENTANYL BINDING ASSESSMENT REPORT",
        "=" * 70,
        "",
        "SUMMARY",
        "-" * 70
    ]
    
    # Ranking
    report_lines.append("\nVariant Ranking (best to worst):")
    for i, name in enumerate(comparison['ranking'], 1):
        v = comparison['variants'][name]
        report_lines.append(f"  {i}. {name}: {v['assessment']} ({v['confidence']} confidence)")
    
    report_lines.append("\n" + "-" * 70)
    report_lines.append("DETAILED ANALYSIS")
    report_lines.append("-" * 70)
    
    for name in comparison['ranking']:
        v = comparison['variants'][name]
        
        report_lines.append(f"\n### {name} ###")
        report_lines.append(f"Assessment: {v['assessment']}")
        report_lines.append(f"Confidence: {v['confidence']}")
        
        if v['good_signs']:
            report_lines.append("\n✓ Positive indicators:")
            for sign in v['good_signs']:
                report_lines.append(f"    • {sign}")
        
        if v['issues']:
            report_lines.append("\n✗ Concerns:")
            for issue in v['issues']:
                report_lines.append(f"    • {issue}")
        
        # Docking details
        if v['docking']:
            report_lines.append("\n  Docking Results:")
            if v['docking'].get('mean_affinity'):
                report_lines.append(f"    Mean affinity: {v['docking']['mean_affinity']:.2f} kcal/mol")
            if v['docking'].get('best_affinity'):
                report_lines.append(f"    Best affinity: {v['docking']['best_affinity']:.2f} kcal/mol")
            if v['docking'].get('salt_bridge_fraction') is not None:
                report_lines.append(f"    Salt bridge preservation: {v['docking']['salt_bridge_fraction']*100:.1f}%")
        
        # MD details
        if v['md'] and v['md'].get('mean_final_rmsd'):
            report_lines.append("\n  MD Simulation Results:")
            report_lines.append(f"    Replicates: {v['md'].get('n_replicates', 'N/A')}")
            report_lines.append(f"    Mean final RMSD: {v['md']['mean_final_rmsd']:.2f} ± {v['md'].get('rmsd_std', 0):.2f} Å")
            if v['md'].get('mean_salt_bridge_intact') is not None:
                report_lines.append(f"    Salt bridge intact: {v['md']['mean_salt_bridge_intact']*100:.1f}%")
    
    # Conclusions
    report_lines.append("\n" + "=" * 70)
    report_lines.append("CONCLUSIONS")
    report_lines.append("=" * 70)
    
    qty_result = comparison['variants'].get('QTY_mutant', {})
    wt_result = comparison['variants'].get('WT_control', {})
    
    if qty_result and wt_result:
        qty_assessment = qty_result.get('assessment', 'UNKNOWN')
        wt_assessment = wt_result.get('assessment', 'UNKNOWN')
        
        if qty_assessment == 'LIKELY_BINDER' and wt_assessment == 'LIKELY_BINDER':
            report_lines.append("""
The QTY variant shows binding characteristics comparable to wild-type.
This suggests the QTY modifications do not significantly disrupt
the fentanyl binding pocket, supporting the hypothesis that the
water-soluble QTY μOR can still bind fentanyl.

RECOMMENDATION: Proceed with experimental validation.
""")
        elif qty_assessment == 'LIKELY_BINDER':
            report_lines.append("""
The QTY variant shows strong binding indicators.

RECOMMENDATION: Proceed with experimental validation.
""")
        elif qty_assessment == 'POSSIBLE_BINDER':
            report_lines.append("""
The QTY variant shows mixed binding indicators. Some concerns exist
but binding cannot be ruled out.

RECOMMENDATION: Consider additional computational analysis or
proceed cautiously with experimental validation.
""")
        else:
            report_lines.append("""
The QTY variant shows poor binding characteristics. The modifications
may have disrupted the fentanyl binding pocket.

RECOMMENDATION: Review mutation positions, particularly those near
the binding pocket. Consider alternative QTY designs that preserve
pocket residues.
""")
    
    report_lines.append("\n" + "=" * 70)
    
    report_text = '\n'.join(report_lines)
    
    # Save report
    report_file = output_dir / 'binding_assessment_report.txt'
    with open(report_file, 'w') as f:
        f.write(report_text)
    
    print(report_text)
    print(f"\nReport saved: {report_file}")
    
    return report_text


def create_comparison_plots(comparison: Dict, output_dir: Path):
    """Create comparison plots."""
    
    if not HAS_MATPLOTLIB:
        print("Matplotlib not available - skipping plots")
        return
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    variant_names = list(comparison['variants'].keys())
    
    # Prepare data
    docking_scores = []
    md_rmsd = []
    sb_docking = []
    sb_md = []
    
    for name in variant_names:
        v = comparison['variants'][name]
        docking_scores.append(v['docking'].get('mean_affinity') or 0)
        md_rmsd.append(v['md'].get('mean_final_rmsd') or 0)
        sb_docking.append((v['docking'].get('salt_bridge_fraction') or 0) * 100)
        sb_md.append((v['md'].get('mean_salt_bridge_intact') or 0) * 100)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    x = np.arange(len(variant_names))
    width = 0.6
    
    # Docking scores
    ax1 = axes[0, 0]
    bars1 = ax1.bar(x, docking_scores, width, color=['#2ecc71' if s < -7 else '#e74c3c' for s in docking_scores])
    ax1.set_ylabel('Docking Score (kcal/mol)')
    ax1.set_title('Docking Affinity')
    ax1.set_xticks(x)
    ax1.set_xticklabels(variant_names)
    ax1.axhline(y=-7, color='gray', linestyle='--', alpha=0.5, label='Good binding threshold')
    
    # MD RMSD
    ax2 = axes[0, 1]
    bars2 = ax2.bar(x, md_rmsd, width, color=['#2ecc71' if r < 3 else '#e74c3c' for r in md_rmsd])
    ax2.set_ylabel('Final Ligand RMSD (Å)')
    ax2.set_title('MD Stability')
    ax2.set_xticks(x)
    ax2.set_xticklabels(variant_names)
    ax2.axhline(y=3, color='gray', linestyle='--', alpha=0.5, label='Stability threshold')
    
    # Salt bridge in docking
    ax3 = axes[1, 0]
    bars3 = ax3.bar(x, sb_docking, width, color=['#2ecc71' if s > 50 else '#e74c3c' for s in sb_docking])
    ax3.set_ylabel('Salt Bridge Preserved (%)')
    ax3.set_title('Asp3.32 Interaction (Docking)')
    ax3.set_xticks(x)
    ax3.set_xticklabels(variant_names)
    ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax3.set_ylim(0, 100)
    
    # Salt bridge in MD
    ax4 = axes[1, 1]
    bars4 = ax4.bar(x, sb_md, width, color=['#2ecc71' if s > 50 else '#e74c3c' for s in sb_md])
    ax4.set_ylabel('Salt Bridge Intact (%)')
    ax4.set_title('Asp3.32 Interaction (MD)')
    ax4.set_xticks(x)
    ax4.set_xticklabels(variant_names)
    ax4.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax4.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'comparison_plots.png', dpi=150)
    plt.close()
    
    print(f"Comparison plots saved: {output_dir / 'comparison_plots.png'}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare QTY variant to WT control"
    )
    parser.add_argument('--qty-docking-dir',
                        help='QTY docking results directory')
    parser.add_argument('--qty-md-dir',
                        help='QTY MD results directory')
    parser.add_argument('--wt-docking-dir',
                        help='WT docking results directory')
    parser.add_argument('--wt-md-dir',
                        help='WT MD results directory')
    parser.add_argument('--output-dir', default='06_analysis/output',
                        help='Output directory')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 6: Comparative Analysis")
    print("=" * 60)
    
    # Set up variant paths
    variants = {}
    
    if args.qty_docking_dir or args.qty_md_dir:
        variants['QTY_mutant'] = {
            'docking': Path(args.qty_docking_dir) if args.qty_docking_dir else Path(),
            'md': Path(args.qty_md_dir) if args.qty_md_dir else Path()
        }
    
    if args.wt_docking_dir or args.wt_md_dir:
        variants['WT_control'] = {
            'docking': Path(args.wt_docking_dir) if args.wt_docking_dir else Path(),
            'md': Path(args.wt_md_dir) if args.wt_md_dir else Path()
        }
    
    if not variants:
        print("No variant directories provided. Using example paths...")
        variants = {
            'QTY_mutant': {
                'docking': Path('04_docking/output/qty'),
                'md': Path('05_md_simulation/output/qty')
            },
            'WT_control': {
                'docking': Path('04_docking/output/wt'),
                'md': Path('05_md_simulation/output/wt')
            }
        }
    
    # Run comparison
    comparison = compare_variants(variants, output_dir)
    
    # Save comparison data
    with open(output_dir / 'comparison_data.json', 'w') as f:
        json.dump(comparison, f, indent=2, default=str)
    
    # Generate report
    generate_report(comparison, output_dir)
    
    # Create plots
    create_comparison_plots(comparison, output_dir)


if __name__ == "__main__":
    main()
