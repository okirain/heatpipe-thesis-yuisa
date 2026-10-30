#!/usr/bin/env python3
"""
================================================================================
Indonesian Space Reactor (ISR) Heat Pipe Core Prototype Simulation
Author: Muhammad Naufal Tamimi (NIM: 20/456315/TK/50445)
Department: Departemen Teknik Nuklir dan Teknik Fisika, Universitas Gadjah Mada
Reference: Khandaq, Harto, & Agung (2020), Progress in Nuclear Energy 118, 103109
================================================================================
"""

import math
import json
import os
import sys

# Ensure local module is importable
sys.path.insert(0, os.path.dirname(__file__))
from heatpipe_rom import HeatPipeGeometry, TwoPhaseHeatPipeROM, SodiumProperties, run_operational_limit_envelope

def simulate_isr_core(core_power_kW=500.0, T_core_operating_K=1050.0, heatpipe_od_mm=18.0):
    """
    Simulates the full 61-pin core of the Indonesian Space Reactor (ISR)
    operating in microgravity (g = 0.0 m/s^2) with Sodium Heat Pipes.
    
    heatpipe_od_mm: Outer diameter of the heat pipe inserted into the fuel pin hollow channel.
                    The fuel pin cladding has a 32.65 mm ID cylindrical hole (Khandaq et al. 2020).
                    OD = 18 mm allows 6.3 mm UN annular fuel thickness while providing 12.4+ kW capacity.
    """
    n_pins = 61
    q_avg_kW = core_power_kW / n_pins  # 8.197 kWth per pin at 500 kWth
    
    # Radial ring power peaking factors based on Khandaq et al. (2020) fast spectrum neutronic profile
    ring_specs = [
        {'ring': 1, 'n_pins': 1,  'radius_cm': 0.0,  'peaking_factor': 1.38, 'name': 'Central Fuel Pin'},
        {'ring': 2, 'n_pins': 6,  'radius_cm': 3.4,  'peaking_factor': 1.25, 'name': 'Inner Core Ring 2'},
        {'ring': 3, 'n_pins': 12, 'radius_cm': 6.8,  'peaking_factor': 1.08, 'name': 'Intermediate Ring 3'},
        {'ring': 4, 'n_pins': 18, 'radius_cm': 10.2, 'peaking_factor': 0.88, 'name': 'Outer Ring 4'},
        {'ring': 5, 'n_pins': 24, 'radius_cm': 13.6, 'peaking_factor': 0.76, 'name': 'Peripheral Ring 5 (Reflector Edge)'}
    ]

    r_outer = (heatpipe_od_mm / 2.0) * 1e-3
    t_wall = 0.0010   # 1.0 mm Mo-30Re
    t_wick = 0.0008   # 0.8 mm composite capillary wick

    # Prototype geometry
    geo_isr = HeatPipeGeometry(
        r_outer=r_outer,
        t_wall=t_wall,
        t_wick=t_wick,
        l_evap=0.25,          # 25 cm active core fuel height (Khandaq et al. 2020)
        l_adiab=0.25,         # 25 cm upper axial reflector & shield
        l_cond=0.50,          # 50 cm condenser interface to power conversion
        pore_radius=25e-6,    # 25 um
        porosity=0.70,
        permeability=3.5e-9,  # High-permeability composite screen/annular wick
        wall_mat='Mo-30Re',
        k_wall=65.0,
        contact_angle=0.0
    )

    # Microgravity space environment: g = 0.0 m/s^2
    rom_isr = TwoPhaseHeatPipeROM(geo_isr, g=0.0, tilt_angle=0.0)

    # Operational limits at T_core_operating_K
    limits = rom_isr.calculate_operational_limits(T_core_operating_K)

    ring_results = []
    all_pins = []

    pin_counter = 1
    for r in ring_specs:
        q_pin_kW = q_avg_kW * r['peaking_factor']
        q_pin_W = q_pin_kW * 1000.0
        
        sim_res = rom_isr.simulate_axial_profiles(
            Q_thermal=q_pin_W,
            T_operating=T_core_operating_K,
            n_nodes=60
        )
        
        # Overall thermal design margin
        q_limit_env = limits['q_envelope']
        margin_overall = (q_limit_env - q_pin_W) / q_limit_env * 100.0
        safety_factor = q_limit_env / q_pin_W if q_pin_W > 0 else 999.0

        r_summary = {
            'ring_id': r['ring'],
            'name': r['name'],
            'n_pins': r['n_pins'],
            'peaking_factor': r['peaking_factor'],
            'power_per_pin_kW': round(q_pin_kW, 3),
            'total_ring_power_kW': round(q_pin_kW * r['n_pins'], 2),
            'max_evaporator_T_K': round(max(sim_res['T_wall']), 1),
            'min_condenser_T_K': round(min(sim_res['T_wall']), 1),
            'delta_T_total_K': round(max(sim_res['T_wall']) - min(sim_res['T_wall']), 1),
            'vapor_pressure_drop_Pa': round(sim_res['summary']['delta_P_v_Pa'], 1),
            'liquid_pressure_drop_Pa': round(sim_res['summary']['delta_P_l_Pa'], 1),
            'total_pressure_drop_Pa': round(sim_res['summary']['total_delta_P_Pa'], 1),
            'capillary_margin_ratio_pct': round(sim_res['summary']['capillary_margin_ratio'] * 100.0, 1),
            'max_mach': round(sim_res['summary']['max_mach'], 4),
            'safety_factor': round(safety_factor, 2),
            'margin_overall_pct': round(margin_overall, 1),
            'status': 'OPERABLE_SAFE' if safety_factor >= 1.10 else 'WARNING_MARGIN'
        }
        ring_results.append(r_summary)

        for p_idx in range(r['n_pins']):
            all_pins.append({
                'pin_id': pin_counter,
                'ring_id': r['ring'],
                'power_kW': round(q_pin_kW, 3),
                'safety_factor': round(safety_factor, 2),
                'status': r_summary['status']
            })
            pin_counter += 1

    # Mass comparison: Original Na-78K Pump System vs Heat Pipe Passive Core
    mass_analysis = {
        'original_khandaq_core_mass_kg': 447.12,
        'pumped_auxiliary_mass_kg': 73.5,
        'heat_pipe_cooling_system_mass_kg': 36.7,
        'net_mass_saving_kg': 36.8,
        'new_heatpipe_reactor_mass_kg': 410.32,
        'mass_reduction_percent': 8.23,
        'reliability_advantages': [
            'Elimination of moving mechanical parts or EM pump coils prone to radiation degradation',
            'Zero acoustic/vibrational noise interference for sensitive spacecraft instrumentation',
            'Complete immunity to single-point coolant pump failure (graceful degradation if 1 pin fails)',
            'Natural passive start-up and shutdown driven solely by capillary surface tension forces',
            'Higher operating temperature capability (1050 K vs 950 K), yielding higher thermodynamic efficiency'
        ]
    }

    limit_envelope = run_operational_limit_envelope(geo_isr, g=0.0, t_start=800.0, t_end=1250.0, t_steps=46)

    isr_summary = {
        'core_power_kWth': core_power_kW,
        'core_operating_temp_K': T_core_operating_K,
        'heatpipe_od_mm': heatpipe_od_mm,
        'total_pins': n_pins,
        'average_pin_power_kW': round(q_avg_kW, 3),
        'peak_pin_power_kW': round(q_avg_kW * ring_specs[0]['peaking_factor'], 3),
        'limits_at_operating_temp': limits,
        'rings': ring_results,
        'mass_comparison': mass_analysis,
        'limit_envelope': limit_envelope
    }

    out_path = os.path.join(os.path.dirname(__file__), 'isr_core_simulation_results.json')
    with open(out_path, 'w') as f:
        json.dump(isr_summary, f, indent=2)
    print(f"ISR simulation data saved to {out_path}")

    return isr_summary


if __name__ == '__main__':
    print("================================================================================")
    print("Indonesian Space Reactor (ISR) 500 kWth Core Heat Pipe Simulation (OD = 18 mm)")
    print("Undergraduate Thesis - Muhammad Naufal Tamimi (20/456315/TK/50445)")
    print("================================================================================\n")
    res = simulate_isr_core(core_power_kW=500.0, T_core_operating_K=1050.0, heatpipe_od_mm=18.0)
    print(f"Total Core Power        : {res['core_power_kWth']} kWth")
    print(f"Total Fuel Pins         : {res['total_pins']}")
    print(f"Nominal Power Per Pin   : {res['average_pin_power_kW']} kWth")
    print(f"Peaked Central Pin Power: {res['peak_pin_power_kW']} kWth\n")
    print("Ring-by-Ring Thermal Hydraulic Results:")
    for r in res['rings']:
        print(f"  Ring {r['ring_id']} ({r['n_pins']} pins): P={r['power_per_pin_kW']:.2f} kW/pin, "
              f"T_evap={r['max_evaporator_T_K']} K, Delta_P={r['total_pressure_drop_Pa']} Pa, "
              f"Safety Factor={r['safety_factor']}x, Status={r['status']}")
    print(f"\nNet Mass Reduction: {res['mass_comparison']['net_mass_saving_kg']:.2f} kg ({res['mass_comparison']['mass_reduction_percent']:.2f}%)")
