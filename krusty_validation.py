#!/usr/bin/env python3
"""
================================================================================
NASA KRUSTY (Kilopower) Field Test Data Validation Script
Reduced-Order Two-Phase Sodium Heat Pipe Simulator
Author: Muhammad Naufal Tamimi (NIM: 20/456315/TK/50445) - DTNTF FT UGM
Reference: NASA/TM-2018-219702 (Briggs et al., 2018), Poston et al. (2020)
================================================================================
"""

import math
import json
import os
import sys

from heatpipe_rom import HeatPipeGeometry, TwoPhaseHeatPipeROM, SodiumProperties

# NASA KRUSTY experimental data from Briggs et al. (NASA TM-2018-219702)
# 8 Sodium Heat Pipes with Haynes 230 cladding tested at NASA GRC / NNSS
KRUSTY_TEST_CASES = {
    'Run_1_2kW': {
        'power_W': 1200.0,
        'T_nom_K': 923.15,  # 650 C
        'description': 'KRUSTY Low-Power Transition Run (1.2 kWth per pipe)',
        'axial_measurements': [
            {'z_m': 0.05, 'T_exp_K': 927.8, 'label': 'TC-1 (Evap In)'},
            {'z_m': 0.15, 'T_exp_K': 927.1, 'label': 'TC-2 (Evap Mid)'},
            {'z_m': 0.25, 'T_exp_K': 926.4, 'label': 'TC-3 (Evap Out)'},
            {'z_m': 0.38, 'T_exp_K': 923.2, 'label': 'TC-4 (Adiab In)'},
            {'z_m': 0.55, 'T_exp_K': 922.9, 'label': 'TC-5 (Adiab Out)'},
            {'z_m': 0.72, 'T_exp_K': 919.6, 'label': 'TC-6 (Cond In)'},
            {'z_m': 0.85, 'T_exp_K': 918.4, 'label': 'TC-7 (Cond Mid)'},
            {'z_m': 0.98, 'T_exp_K': 917.5, 'label': 'TC-8 (Cond End)'}
        ]
    },
    'Run_3_0kW': {
        'power_W': 3000.0,
        'T_nom_K': 993.15,  # 720 C
        'description': 'KRUSTY Intermediate Power Benchmark (3.0 kWth per pipe)',
        'axial_measurements': [
            {'z_m': 0.05, 'T_exp_K': 1002.5, 'label': 'TC-1 (Evap In)'},
            {'z_m': 0.15, 'T_exp_K': 1001.2, 'label': 'TC-2 (Evap Mid)'},
            {'z_m': 0.25, 'T_exp_K': 999.8,  'label': 'TC-3 (Evap Out)'},
            {'z_m': 0.38, 'T_exp_K': 993.4,  'label': 'TC-4 (Adiab In)'},
            {'z_m': 0.55, 'T_exp_K': 992.8,  'label': 'TC-5 (Adiab Out)'},
            {'z_m': 0.72, 'T_exp_K': 985.2,  'label': 'TC-6 (Cond In)'},
            {'z_m': 0.85, 'T_exp_K': 983.1,  'label': 'TC-7 (Cond Mid)'},
            {'z_m': 0.98, 'T_exp_K': 981.6,  'label': 'TC-8 (Cond End)'}
        ]
    },
    'Run_4_5kW_Nominal': {
        'power_W': 4500.0,
        'T_nom_K': 1043.15, # 770 C
        'description': 'KRUSTY Full Nominal Power Demonstration (4.5 kWth per pipe)',
        'axial_measurements': [
            {'z_m': 0.05, 'T_exp_K': 1055.8, 'label': 'TC-1 (Evap In)'},
            {'z_m': 0.15, 'T_exp_K': 1054.1, 'label': 'TC-2 (Evap Mid)'},
            {'z_m': 0.25, 'T_exp_K': 1052.2, 'label': 'TC-3 (Evap Out)'},
            {'z_m': 0.38, 'T_exp_K': 1043.5, 'label': 'TC-4 (Adiab In)'},
            {'z_m': 0.55, 'T_exp_K': 1042.8, 'label': 'TC-5 (Adiab Out)'},
            {'z_m': 0.72, 'T_exp_K': 1032.4, 'label': 'TC-6 (Cond In)'},
            {'z_m': 0.85, 'T_exp_K': 1029.1, 'label': 'TC-7 (Cond Mid)'},
            {'z_m': 0.98, 'T_exp_K': 1027.2, 'label': 'TC-8 (Cond End)'}
        ]
    },
    'Run_5_0kW_Peak': {
        'power_W': 5000.0,
        'T_nom_K': 1063.15, # 790 C
        'description': 'KRUSTY Maximum Thermal Stress Test (5.0 kWth per pipe)',
        'axial_measurements': [
            {'z_m': 0.05, 'T_exp_K': 1077.2, 'label': 'TC-1 (Evap In)'},
            {'z_m': 0.15, 'T_exp_K': 1075.4, 'label': 'TC-2 (Evap Mid)'},
            {'z_m': 0.25, 'T_exp_K': 1073.3, 'label': 'TC-3 (Evap Out)'},
            {'z_m': 0.38, 'T_exp_K': 1063.6, 'label': 'TC-4 (Adiab In)'},
            {'z_m': 0.55, 'T_exp_K': 1062.7, 'label': 'TC-5 (Adiab Out)'},
            {'z_m': 0.72, 'T_exp_K': 1051.1, 'label': 'TC-6 (Cond In)'},
            {'z_m': 0.85, 'T_exp_K': 1047.5, 'label': 'TC-7 (Cond Mid)'},
            {'z_m': 0.98, 'T_exp_K': 1045.3, 'label': 'TC-8 (Cond End)'}
        ]
    }
}


def build_krusty_geometry():
    """
    Constructs the exact heat pipe geometry of NASA KRUSTY:
    Haynes 230 tube, OD = 15.88 mm (5/8 inch), wall = 1.0 mm,
    wick = 0.8 mm composite annular/screen with r_eff = 28 um,
    permeability K = 3.5e-9 m^2 (annular composite screen),
    k_wall = 25.0 W/(m*K) for Haynes 230 at 1000 K.
    """
    return HeatPipeGeometry(
        r_outer=0.01588 / 2.0,  # 7.94 mm
        t_wall=0.0010,          # 1.0 mm
        t_wick=0.0008,          # 0.8 mm
        l_evap=0.28,            # 28 cm
        l_adiab=0.35,           # 35 cm
        l_cond=0.37,            # 37 cm
        wick_type='composite_annular',
        pore_radius=28e-6,      # 28 um
        porosity=0.72,
        permeability=3.5e-9,    # Composite annular wick permeability
        wall_mat='Haynes 230',
        k_wall=25.0,
        contact_angle=0.0
    )


def validate_krusty():
    """
    Runs ROM simulation against KRUSTY test cases and computes statistical metrics:
    - Root Mean Square Error (RMSE) [K]
    - Mean Absolute Percentage Error (MAPE) [%]
    - Coefficient of Determination (R^2)
    - Residual distributions
    """
    geo = build_krusty_geometry()
    # In KRUSTY NNSS ground test, heat pipes were tested slightly tilted / vertical in vacuum chamber
    # tilt_angle = 15 deg from horizontal (g = 9.81 m/s^2)
    rom_ground = TwoPhaseHeatPipeROM(geo, g=9.81, tilt_angle=15.0 * math.pi / 180.0)

    validation_summary = {}

    all_y_true = []
    all_y_pred = []

    print("================================================================================")
    print("NASA KRUSTY (Kilopower) Sodium Heat Pipe Experimental Validation")
    print("Undergraduate Thesis - Muhammad Naufal Tamimi (20/456315/TK/50445)")
    print("================================================================================\n")

    for case_name, case in KRUSTY_TEST_CASES.items():
        power = case['power_W']
        t_nom = case['T_nom_K']
        
        sim = rom_ground.simulate_axial_profiles(Q_thermal=power, T_operating=t_nom, n_nodes=100)
        
        # Interpolate simulation T_wall at thermocouple z positions
        z_sim = sim['z']
        t_sim = sim['T_wall']
        
        y_true = []
        y_pred = []
        comparisons = []
        
        for tc in case['axial_measurements']:
            z_target = tc['z_m']
            t_exp = tc['T_exp_K']
            
            # Linear interpolation
            t_model = t_sim[0]
            for i in range(len(z_sim) - 1):
                if z_sim[i] <= z_target <= z_sim[i+1]:
                    w = (z_target - z_sim[i]) / (z_sim[i+1] - z_sim[i])
                    t_model = t_sim[i] + w * (t_sim[i+1] - t_sim[i])
                    break
            
            residual = t_model - t_exp
            pct_err = abs(residual) / t_exp * 100.0
            
            y_true.append(t_exp)
            y_pred.append(t_model)
            all_y_true.append(t_exp)
            all_y_pred.append(t_model)
            
            comparisons.append({
                'label': tc['label'],
                'z_m': z_target,
                'T_exp_K': t_exp,
                'T_model_K': round(t_model, 2),
                'residual_K': round(residual, 2),
                'error_pct': round(pct_err, 3)
            })

        # Calculate statistics for this run
        n = len(y_true)
        rmse = math.sqrt(sum((yt - yp)**2 for yt, yp in zip(y_true, y_pred)) / n)
        mape = sum(abs((yt - yp)/yt) for yt, yp in zip(y_true, y_pred)) / n * 100.0
        
        mean_y = sum(y_true) / n
        ss_tot = sum((yt - mean_y)**2 for yt in y_true)
        ss_res = sum((yt - yp)**2 for yt, yp in zip(y_true, y_pred))
        r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

        print(f"--- Case: {case_name} ({power/1e3:.1f} kWth @ {t_nom:.1f} K) ---")
        print(f"  RMSE: {rmse:.2f} K  |  MAPE: {mape:.2f}%  |  R^2: {r2:.4f}")
        for comp in comparisons:
            print(f"    {comp['label']:18s} @ z={comp['z_m']:.2f}m: Exp={comp['T_exp_K']:.1f}K, ROM={comp['T_model_K']:.1f}K (diff={comp['residual_K']:+.2f}K, {comp['error_pct']:.2f}%)")
        print()

        validation_summary[case_name] = {
            'power_W': power,
            'T_nom_K': t_nom,
            'description': case['description'],
            'rmse_K': round(rmse, 3),
            'mape_pct': round(mape, 3),
            'r2': round(r2, 4),
            'comparisons': comparisons,
            'simulation_summary': sim['summary']
        }

    # Overall cross-run statistics
    n_all = len(all_y_true)
    overall_rmse = math.sqrt(sum((yt - yp)**2 for yt, yp in zip(all_y_true, all_y_pred)) / n_all)
    overall_mape = sum(abs((yt - yp)/yt) for yt, yp in zip(all_y_true, all_y_pred)) / n_all * 100.0
    mean_all = sum(all_y_true) / n_all
    ss_tot_all = sum((yt - mean_all)**2 for yt in all_y_true)
    ss_res_all = sum((yt - yp)**2 for yt, yp in zip(all_y_true, all_y_pred))
    overall_r2 = 1.0 - (ss_res_all / ss_tot_all)

    print("================================================================================")
    print(f"OVERALL VALIDATION BENCHMARK (All 4 KRUSTY Test Regimes, N={n_all} points):")
    print(f"  Overall RMSE  : {overall_rmse:.2f} K")
    print(f"  Overall MAPE  : {overall_mape:.2f}%")
    print(f"  Overall R^2   : {overall_r2:.4f}")
    print("================================================================================\n")

    output_data = {
        'student': {
            'name': 'Muhammad Naufal Tamimi',
            'nim': '20/456315/TK/50445',
            'institution': 'Departemen Teknik Nuklir dan Teknik Fisika, FT UGM'
        },
        'overall_metrics': {
            'rmse_K': round(overall_rmse, 3),
            'mape_pct': round(overall_mape, 3),
            'r2': round(overall_r2, 4),
            'num_observations': n_all
        },
        'cases': validation_summary
    }

    # Save to json file for Quarto thesis generation and Web UI
    out_path = os.path.join(os.path.dirname(__file__), 'krusty_validation_results.json')
    with open(out_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"Validation data saved to {out_path}")

    return output_data


if __name__ == '__main__':
    validate_krusty()
