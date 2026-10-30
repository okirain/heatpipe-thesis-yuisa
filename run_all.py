#!/usr/bin/env python3
"""
================================================================================
Master Execution Script: Two-Phase Sodium Heat Pipe Simulation & Thesis Analysis
Student: Muhammad Naufal Tamimi (NIM: 20/456315/TK/50445)
Departemen Teknik Nuklir dan Teknik Fisika, Universitas Gadjah Mada
================================================================================
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

import heatpipe_rom
import krusty_validation
import isr_core_prototype

def main():
    print("================================================================================")
    print("RUNNING MASTER HEAT PIPE NUMERICAL SIMULATION & THESIS DATA GENERATION")
    print("Author: Muhammad Naufal Tamimi (20/456315/TK/50445)")
    print("Departemen Teknik Nuklir dan Teknik Fisika, Universitas Gadjah Mada")
    print("================================================================================\n")
    
    print("[1/2] Executing NASA KRUSTY Field Test Validation...")
    val_data = krusty_validation.validate_krusty()
    print(">> Validation completed successfully. Overall R^2 =", val_data['overall_metrics']['r2'])
    
    print("\n[2/2] Executing Indonesian Space Reactor (ISR 500 kWth) Core Simulation...")
    isr_data = isr_core_prototype.simulate_isr_core(core_power_kW=500.0, T_core_operating_K=1050.0, heatpipe_od_mm=18.0)
    print(">> Core simulation completed successfully. 61 pins simulated.")
    print(">> Net reactor mass reduction =", isr_data['mass_comparison']['net_mass_saving_kg'], "kg")
    
    print("\nALL SIMULATION DATA READY FOR QUARTO THESIS AND WEB INTERFACE!")

if __name__ == '__main__':
    main()
