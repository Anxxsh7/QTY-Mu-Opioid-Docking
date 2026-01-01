#!/usr/bin/env python3
"""
Analyze QTY mutations between WT and QTY-variant OPRM1 sequences.
Identifies all mutations and flags those near the binding pocket.
"""

# WT sequence (cleaned - remove spaces)
wt_seq = """
MDSSAAPTNA SNCTDALAYS SCSPAPSPGS WVNLSHLDGN LSDPCGPNRT
DLGGRDSLCP PTGSPSMITA ITIMALYSIV CVVGLFGNFL VMYVIVRYTK 
MKTATNIYIF NLALADALAT STLPFQSVNY LMGTWPFGTI LCKIVISIDY
YNMFTSIFTL CTMSVDRYIA VCHPVKALDF RTPRNAKIIN VCNWILSSAI
GLPVMFMATT KYRQGSIDCT LTFSHPTWYW ENLLKICVFI FAFIMPVLII
TVCYGLMILR LKSVRMLSGS KEKDRNLRRI TRMVLVVVAV FIVCWTPIHI
YVIIKALVTI PETTFQTVSW HFCIALGYTN SCLNPVLYAF LDENFKRCFR
EFCIPTSSNI EQQNSTRIRQ NTRDHPSTAN TVDRTNHQLE NLEAETAPLP
""".replace(" ", "").replace("\n", "")

# QTY variant sequence (cleaned)
qty_seq = """
MDSSAAPTNASNCTDALAYSSCSPAPSPGSWVNLSHLDGNLSDPCGPNRTDLGGRDSLCP
PTGSPSMITATTTMAQYSTTCTTGQYGNYQTMYVIVRYTKMKTATNTYTYNQAQADAQAT
STQPYQSTNYQMGTWPFGTILCKTTTSTDYYNMYTSTYTQCTMSTDRYIAVCHPVKALDF
RTPRNAKTTNTCNWTQSSATGQPTMYMATTKYRQGSIDCTLTFSHPTWYWENQQKTCTYT
YAYTMPTQTTTTCYGLMILRLKSVRMLSGSKEKDRNLRRTTRMTQTTTATYTTCWTPTHT
YTTTKALVTIPETTYQTTSWHYCTAQGYTNSCQNPTQYAFLDENFKRCFREFCIPTSSNI
EQQNSTRIRQNTRDHPSTANTVDRTNHQLENLEAETAPLP
""".replace(" ", "").replace("\n", "")

# QTY rules
QTY_RULES = {
    'L': 'Q',  # Leucine → Glutamine
    'I': 'T',  # Isoleucine → Threonine
    'V': 'T',  # Valine → Threonine
    'F': 'Y',  # Phenylalanine → Tyrosine
}

# Critical binding pocket residues (approximate positions in sequence)
# Based on μOR literature - Asp147 (Asp3.32) is the key residue
# Binding pocket residues from 8EF5 structure analysis
BINDING_POCKET_RESIDUES = {
    147: "Asp3.32 - CRITICAL salt bridge for opioids",
    148: "Tyr3.33 - pocket lining",
    149: "Met3.34 - pocket lining", 
    228: "Trp6.48 - toggle switch",
    229: "Ile6.49 - pocket lining",
    232: "His6.52 - pocket interaction",
    293: "Trp6.48 - toggle switch",
    297: "His6.52 - key contact",
    298: "Ile6.53 - pocket",
    301: "Tyr7.35 - pocket",
    302: "Trp7.35 - pocket",
}

# Transmembrane helix boundaries (approximate)
TM_REGIONS = {
    'TM1': (65, 95),
    'TM2': (100, 130),
    'TM3': (135, 170),
    'TM4': (185, 210),
    'TM5': (230, 265),
    'TM6': (275, 305),
    'TM7': (312, 345),
}

def get_tm_region(pos):
    """Get TM region for a position."""
    for tm, (start, end) in TM_REGIONS.items():
        if start <= pos <= end:
            return tm
    return "Loop/Terminal"

print("=" * 70)
print("QTY MUTATION ANALYSIS FOR OPRM1")
print("=" * 70)

print(f"\nWT sequence length: {len(wt_seq)}")
print(f"QTY sequence length: {len(qty_seq)}")

if len(wt_seq) != len(qty_seq):
    print(f"\n⚠️  WARNING: Sequence lengths differ!")
    print(f"    This may cause alignment issues.")

# Find all mutations
mutations = []
qty_mutations = []
non_qty_mutations = []

min_len = min(len(wt_seq), len(qty_seq))

for i in range(min_len):
    if wt_seq[i] != qty_seq[i]:
        pos = i + 1  # 1-based numbering
        wt_aa = wt_seq[i]
        qty_aa = qty_seq[i]
        
        is_qty = (wt_aa in QTY_RULES and QTY_RULES[wt_aa] == qty_aa)
        is_pocket = pos in BINDING_POCKET_RESIDUES
        tm_region = get_tm_region(pos)
        
        mut_info = {
            'position': pos,
            'wt': wt_aa,
            'qty': qty_aa,
            'string': f"{wt_aa}{pos}{qty_aa}",
            'is_qty_rule': is_qty,
            'is_pocket': is_pocket,
            'pocket_note': BINDING_POCKET_RESIDUES.get(pos, ""),
            'tm_region': tm_region
        }
        
        mutations.append(mut_info)
        if is_qty:
            qty_mutations.append(mut_info)
        else:
            non_qty_mutations.append(mut_info)

print(f"\n{'='*70}")
print(f"MUTATION SUMMARY")
print(f"{'='*70}")
print(f"Total mutations: {len(mutations)}")
print(f"QTY-rule mutations: {len(qty_mutations)}")
print(f"Non-QTY mutations: {len(non_qty_mutations)}")

# Group by TM region
print(f"\n{'='*70}")
print("MUTATIONS BY TRANSMEMBRANE REGION")
print(f"{'='*70}")

tm_counts = {}
for mut in mutations:
    tm = mut['tm_region']
    tm_counts[tm] = tm_counts.get(tm, 0) + 1

for tm in ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'Loop/Terminal']:
    if tm in tm_counts:
        print(f"  {tm}: {tm_counts[tm]} mutations")

# List all QTY mutations
print(f"\n{'='*70}")
print("ALL QTY-RULE MUTATIONS (L→Q, I→T, V→T, F→Y)")
print(f"{'='*70}")

for mut in qty_mutations:
    pocket_flag = " ⚠️  BINDING POCKET!" if mut['is_pocket'] else ""
    print(f"  {mut['string']:8s} ({mut['tm_region']}){pocket_flag}")
    if mut['pocket_note']:
        print(f"           → {mut['pocket_note']}")

# Check for non-QTY mutations (these might be errors or intentional)
if non_qty_mutations:
    print(f"\n{'='*70}")
    print("⚠️  NON-QTY MUTATIONS (verify these are correct)")
    print(f"{'='*70}")
    for mut in non_qty_mutations:
        print(f"  {mut['string']:8s} ({mut['tm_region']})")

# Binding pocket analysis
print(f"\n{'='*70}")
print("🔬 BINDING POCKET ANALYSIS")
print(f"{'='*70}")

pocket_mutations = [m for m in mutations if m['is_pocket']]
near_pocket = [m for m in mutations if 140 <= m['position'] <= 160 or 
               225 <= m['position'] <= 235 or 290 <= m['position'] <= 305]

if pocket_mutations:
    print("\n❌ CRITICAL: Mutations IN the binding pocket:")
    for mut in pocket_mutations:
        print(f"   {mut['string']} - {mut['pocket_note']}")
else:
    print("\n✓ No mutations directly in critical binding pocket residues")

# Check Asp147 specifically
asp147_mutated = any(m['position'] == 147 for m in mutations)
if asp147_mutated:
    print("\n❌ CRITICAL ERROR: Asp147 (Asp3.32) is mutated!")
    print("   This residue forms the essential salt bridge with fentanyl.")
    print("   Binding is very unlikely to be preserved!")
else:
    print("\n✓ Asp147 (Asp3.32) is preserved - essential for opioid binding")

# Near-pocket mutations
near_pocket_muts = [m for m in mutations if not m['is_pocket'] and
                   (140 <= m['position'] <= 160 or 
                    225 <= m['position'] <= 240 or 
                    290 <= m['position'] <= 310)]

if near_pocket_muts:
    print(f"\n⚠️  Mutations NEAR binding pocket (may affect binding):")
    for mut in near_pocket_muts:
        print(f"   {mut['string']} (position {mut['position']})")

# Save mutations to file
print(f"\n{'='*70}")
print("SAVING MUTATION FILES")
print(f"{'='*70}")

# Save mutation list
mutation_strings = [m['string'] for m in qty_mutations]
with open('qty_mutations.txt', 'w') as f:
    f.write("# QTY mutations for OPRM1\n")
    f.write("# Generated by analyze_qty_mutations.py\n")
    f.write(f"# Total: {len(mutation_strings)} mutations\n\n")
    for m in mutation_strings:
        f.write(f"{m}\n")
print(f"  Saved: qty_mutations.txt ({len(mutation_strings)} mutations)")

# Save full analysis as JSON
import json
analysis = {
    'wt_length': len(wt_seq),
    'qty_length': len(qty_seq),
    'total_mutations': len(mutations),
    'qty_rule_mutations': len(qty_mutations),
    'non_qty_mutations': len(non_qty_mutations),
    'pocket_mutations': len(pocket_mutations),
    'asp147_preserved': not asp147_mutated,
    'mutations': mutations,
    'mutation_list': mutation_strings
}

with open('qty_mutation_analysis.json', 'w') as f:
    json.dump(analysis, f, indent=2)
print(f"  Saved: qty_mutation_analysis.json")

# Print final assessment
print(f"\n{'='*70}")
print("📋 FINAL ASSESSMENT")
print(f"{'='*70}")

if asp147_mutated:
    print("""
❌ HIGH RISK: The critical Asp3.32 residue is mutated.
   This is the anchor point for opioid binding.
   Fentanyl binding is very unlikely to be preserved.
   
   RECOMMENDATION: Revise QTY design to preserve Asp147.
""")
elif pocket_mutations:
    print("""
⚠️  MODERATE RISK: Some binding pocket residues are mutated.
   This may affect binding affinity and pose.
   
   RECOMMENDATION: Proceed with docking but carefully monitor
   the binding pocket integrity during MD simulations.
""")
elif near_pocket_muts:
    print("""
⚠️  LOW-MODERATE RISK: Mutations near binding pocket.
   These may indirectly affect pocket geometry.
   
   RECOMMENDATION: Proceed with full pipeline.
   Pay attention to pocket RMSD during MD.
""")
else:
    print("""
✓ LOW RISK: No binding pocket mutations detected.
   QTY mutations are in transmembrane regions away from pocket.
   Fentanyl binding likely to be preserved.
   
   RECOMMENDATION: Proceed with full docking/MD pipeline.
""")

print(f"\nTo run the pipeline with these mutations:")
print(f"  python run_pipeline.py --mutations-file qty_mutations.txt")
print(f"\nOr use individual mutations:")
print(f"  python run_pipeline.py --mutations {' '.join(mutation_strings[:5])} ...")
