#!/usr/bin/env python3
"""
Analytical and empirical Strouhal number validation script for WAKE CFD.
Compares computed cylinder vortex shedding frequencies against universal
experimental relationships (Williamson 1989, Fey et al. 1998):
St(Re) = 0.2665 - 1.018 / sqrt(Re)
"""

import sys
import os
import json
import math
import argparse

def theoretical_strouhal_williamson(reynolds):
    """
    Williamson (1989) universal Strouhal-Reynolds relationship
    for laminar cylinder vortex shedding (47 < Re < 180).
    """
    if reynolds < 47.0:
        return 0.0
    return 0.2665 - (1.018 / math.sqrt(reynolds))

def verify_strouhal_results(summary_json_path, reynolds=100.0, max_relative_error_pct=10.0):
    with open(summary_json_path, 'r') as f:
        data = json.load(f)

    aero = data.get('aerodynamics', {})
    sim_st = aero.get('strouhal_number', 0.0)
    sim_freq = aero.get('shedding_frequency', 0.0)

    theory_st = theoretical_strouhal_williamson(reynolds)
    rel_error = abs(sim_st - theory_st) / theory_st * 100.0 if theory_st > 0 else 0.0

    print("================================================================================")
    print("  WAKE CFD: Karman Vortex Shedding Strouhal Number Verification")
    print("================================================================================")
    print(f"  Reynolds Number (Re):       {reynolds:.1f}")
    print(f"  Theoretical St (Williamson): {theory_st:.4f}")
    print(f"  Simulated St (WAKE LBM):    {sim_st:.4f}")
    print(f"  Shedding Frequency (f_shed): {sim_freq:.4f}")
    print(f"  Relative Error:             {rel_error:.2f}% (Threshold: {max_relative_error_pct:.1f}%)")
    print("--------------------------------------------------------------------------------")

    passed = (sim_st > 0.0 and rel_error <= max_relative_error_pct)
    verdict = "PASSED (Within Experimental Error Band)" if passed else "VERIFICATION REVIEW"
    print(f"  Verdict:                    {verdict}")
    print("================================================================================")
    return passed

def main():
    parser = argparse.ArgumentParser(description="Verify WAKE Strouhal number against theoretical curve")
    parser.add_argument("summary_json", type=str, help="Path to karman_summary.json")
    parser.add_argument("--re", type=float, default=100.0, help="Simulation Reynolds number")
    parser.add_argument("--tolerance", type=float, default=12.0, help="Maximum allowable error percentage")
    args = parser.parse_args()

    if not os.path.exists(args.summary_json):
        print(f"Error: File not found: {args.summary_json}")
        sys.exit(1)

    ok = verify_strouhal_results(args.summary_json, args.re, args.tolerance)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
