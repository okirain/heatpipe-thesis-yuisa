#!/usr/bin/env python3
"""
================================================================================
Reduced-Order Numerical Model (ROM) for Liquid Metal Heat Pipes in Space Reactors
Author: Muhammad Naufal Tamimi (NIM: 20/456315/TK/50445)
Department: Departemen Teknik Nuklir dan Teknik Fisika, Universitas Gadjah Mada
Application: Microgravity Thermal-Hydraulics for Space Nuclear Power
             (KRUSTY Benchmark & Indonesian Space Reactor / ISR Core Integration)
================================================================================
"""

import math
import json
import sys
import os

class SodiumProperties:
    """
    Thermophysical properties of liquid and vapor sodium as a function of
    temperature (K), based on Fink & Leibowitz (1995, ANL-RE-95-2), Faghri (1995),
    and nuclear thermal hydraulics data.
    """
    @staticmethod
    def p_sat(T):
        """Saturation vapor pressure [Pa] from Fink & Leibowitz (1995)"""
        # ln(P_MPa) = 11.9463 - 12633.73/T - 0.5 * ln(T)
        T_clamped = max(500.0, min(2400.0, float(T)))
        ln_p_mpa = 11.9463 - (12633.73 / T_clamped) - 0.5 * math.log(T_clamped)
        p_mpa = math.exp(ln_p_mpa)
        return p_mpa * 1.0e6  # Pa (1 MPa = 1e6 Pa)

    @staticmethod
    def t_sat(P):
        """Saturation temperature [K] given vapor pressure [Pa] via Newton-Raphson"""
        T = 1000.0
        P_target = max(10.0, float(P))
        for _ in range(25):
            p_calc = SodiumProperties.p_sat(T)
            # Derivative: d(ln P)/dT = 12633.73 / T^2 - 0.5 / T
            dp_dT = p_calc * (12633.73 / (T * T) - 0.5 / T)
            diff = p_calc - P_target
            if abs(diff) < 1.0 or abs(dp_dT) < 1e-12:
                break
            T = T - diff / dp_dT
            T = max(600.0, min(2000.0, T))
        return T

    @staticmethod
    def h_fg(T):
        """Latent heat of vaporization [J/kg]"""
        T_c = 2503.7  # Critical temperature of sodium [K]
        T_clamped = max(500.0, min(T_c - 1.0, float(T)))
        theta = max(1e-4, 1.0 - T_clamped / T_c)
        return 4.76e6 * (theta ** 0.38)

    @staticmethod
    def rho_l(T):
        """Liquid density [kg/m^3] from ANL-RE-95-2"""
        T_c = float(T)
        return 1011.8 - 0.22054 * T_c - 1.9226e-5 * (T_c**2) + 5.637e-9 * (T_c**3)

    @staticmethod
    def rho_v(T, P):
        """Vapor density [kg/m^3] from ideal gas law with sodium molecular weight"""
        R_na = 8314.46 / 22.98977  # 361.66 J/(kg*K)
        return float(P) / (R_na * float(T))

    @staticmethod
    def mu_l(T):
        """Liquid dynamic viscosity [Pa*s]"""
        T_c = float(T)
        return math.exp(-6.4406 - 0.3958 * math.log(T_c) + 556.835 / T_c)

    @staticmethod
    def mu_v(T):
        """Vapor dynamic viscosity [Pa*s]"""
        return 1.0e-7 * math.sqrt(float(T))

    @staticmethod
    def sigma(T):
        """Surface tension [N/m]"""
        T_c = 2503.7
        T_clamped = max(500.0, min(T_c - 1.0, float(T)))
        theta = max(1e-4, 1.0 - T_clamped / T_c)
        return 0.2405 * (theta ** 1.126)

    @staticmethod
    def k_l(T):
        """Liquid thermal conductivity [W/(m*K)]"""
        return 104.0 - 0.047 * float(T)

    @staticmethod
    def cp_l(T):
        """Liquid specific heat capacity [J/(kg*K)]"""
        T_c = float(T)
        return 1436.1 - 0.5802 * T_c + 4.621e-4 * (T_c**2)


class HeatPipeGeometry:
    """
    Geometry, wick, and wall material parameters for sodium space reactor heat pipe.
    """
    def __init__(self,
                 r_outer=0.008,        # Outer radius of envelope [m] (e.g., 16 mm OD)
                 t_wall=0.001,         # Cladding wall thickness [m]
                 t_wick=0.0008,        # Wick structure thickness [m]
                 l_evap=0.25,          # Evaporator length [m] (ISR active fuel height = 25 cm)
                 l_adiab=0.25,         # Adiabatic transport section length [m]
                 l_cond=0.50,          # Condenser length [m]
                 wick_type='screen',   # 'screen', 'sintered', 'groove'
                 pore_radius=25e-6,    # Effective capillary pore radius [m] (25 um)
                 porosity=0.70,        # Wick porosity epsilon
                 permeability=1.5e-10, # Wick permeability K [m^2]
                 wall_mat='Mo-30Re',   # Wall material
                 k_wall=65.0,          # Wall thermal conductivity [W/(m*K)]
                 contact_angle=0.0):   # Contact angle [rad] (0 for perfect wetting)
        self.r_outer = float(r_outer)
        self.t_wall = float(t_wall)
        self.r_inner = self.r_outer - self.t_wall
        self.t_wick = float(t_wick)
        self.r_vapor = self.r_inner - self.t_wick
        
        self.l_evap = float(l_evap)
        self.l_adiab = float(l_adiab)
        self.l_cond = float(l_cond)
        self.l_total = self.l_evap + self.l_adiab + self.l_cond
        self.l_eff = 0.5 * self.l_evap + self.l_adiab + 0.5 * self.l_cond
        
        self.wick_type = wick_type
        self.r_eff = float(pore_radius)
        self.porosity = float(porosity)
        self.permeability = float(permeability)
        self.wall_mat = wall_mat
        self.k_wall = float(k_wall)
        self.contact_angle = float(contact_angle)

        self.a_vapor = math.pi * (self.r_vapor ** 2)
        self.a_wick = math.pi * (self.r_inner**2 - self.r_vapor**2)
        self.d_vapor = 2.0 * self.r_vapor

    def effective_wick_k(self, k_liquid):
        """Effective thermal conductivity of liquid-saturated wick structure"""
        eps = self.porosity
        kw = self.k_wall
        kl = k_liquid
        if self.wick_type == 'screen':
            num = kl * ((kl + kw) - (1.0 - eps) * (kl - kw))
            den = (kl + kw) + (1.0 - eps) * (kl - kw)
            return num / den
        else:
            return (kw ** (1.0 - eps)) * (kl ** eps)


class TwoPhaseHeatPipeROM:
    """
    1D Reduced-Order Model for counter-current two-phase flow in microgravity/variable-g,
    solving coupled vapor core compressible/frictional dynamics, porous wick Darcy flow,
    and calculating the 5 fundamental heat transport operational limits.
    """
    def __init__(self, geometry: HeatPipeGeometry, g=0.0, tilt_angle=0.0):
        self.geo = geometry
        self.g = float(g)
        self.tilt_angle = float(tilt_angle)

    def calculate_operational_limits(self, T_sat):
        """
        Calculates all 5 operational limits [W] at a given saturation temperature T_sat [K]:
        1. Capillary Limit (Young-Laplace / Darcy / Poiseuille balance)
        2. Sonic Limit (Levy / Busse choked flow criterion)
        3. Entrainment Limit (Cotter Weber number shear criterion)
        4. Viscous Limit (Busse vapor pressure choked limit)
        5. Boiling Limit (Nucleation critical superheat)
        """
        geo = self.geo
        p_sat = SodiumProperties.p_sat(T_sat)
        h_fg = SodiumProperties.h_fg(T_sat)
        rho_l = SodiumProperties.rho_l(T_sat)
        rho_v = SodiumProperties.rho_v(T_sat, p_sat)
        mu_l = SodiumProperties.mu_l(T_sat)
        mu_v = SodiumProperties.mu_v(T_sat)
        sigma = SodiumProperties.sigma(T_sat)
        k_l = SodiumProperties.k_l(T_sat)
        k_eff = geo.effective_wick_k(k_l)
        gamma = 1.67  # Monatomic sodium vapor ratio of specific heats
        R_na = 361.66 # Gas constant J/(kg*K)

        # 1. Capillary Limit
        delta_p_c_max = (2.0 * sigma * math.cos(geo.contact_angle)) / geo.r_eff
        delta_p_g = rho_l * self.g * geo.l_total * math.sin(self.tilt_angle)
        f_l = mu_l / (rho_l * geo.a_wick * geo.permeability * h_fg)
        f_v = (8.0 * mu_v) / (math.pi * rho_v * (geo.r_vapor ** 4) * h_fg)

        available_capillary = delta_p_c_max - delta_p_g
        if available_capillary > 0:
            q_cap = available_capillary / ((f_l + f_v) * geo.l_eff)
        else:
            q_cap = 0.0

        # 2. Sonic Limit
        q_sonic = geo.a_vapor * rho_v * h_fg * math.sqrt((gamma * R_na * T_sat) / (2.0 * (gamma + 1.0)))

        # 3. Entrainment Limit
        r_hs = geo.r_eff * 0.5
        q_ent = geo.a_vapor * h_fg * math.sqrt((sigma * rho_v) / (2.0 * r_hs))

        # 4. Viscous Limit
        q_visc = (geo.r_vapor ** 2 * h_fg * rho_v * p_sat * geo.a_vapor) / (16.0 * mu_v * geo.l_eff)

        # 5. Boiling Limit
        r_n = 2.5e-6  # Cavity radius ~ 2.5 um
        delta_T_crit = (2.0 * sigma * T_sat) / (rho_v * h_fg * r_n)
        q_boil = (2.0 * math.pi * geo.l_evap * k_eff * delta_T_crit) / math.log(geo.r_inner / geo.r_vapor)

        # Governing envelope
        limits = [
            ('Capillary', q_cap),
            ('Sonic', q_sonic),
            ('Entrainment', q_ent),
            ('Viscous', q_visc),
            ('Boiling', q_boil)
        ]
        gov_limit, q_max = min(limits, key=lambda x: x[1])

        return {
            'T_sat': T_sat,
            'P_sat': p_sat,
            'q_cap': q_cap,
            'q_sonic': q_sonic,
            'q_ent': q_ent,
            'q_visc': q_visc,
            'q_boil': q_boil,
            'q_envelope': q_max,
            'governing_limit': gov_limit,
            'delta_p_c_max': delta_p_c_max,
            'delta_p_g': delta_p_g
        }

    def simulate_axial_profiles(self, Q_thermal, T_operating, n_nodes=60):
        """
        Discretizes heat pipe axially into n_nodes and computes profiles for:
        m_dot(z), u_vapor(z), Mach(z), P_v(z), P_l(z), T_v(z), T_wall(z), Capillary margin.
        """
        geo = self.geo
        dz = geo.l_total / (n_nodes - 1)
        
        p_ref = SodiumProperties.p_sat(T_operating)
        h_fg = SodiumProperties.h_fg(T_operating)
        rho_l = SodiumProperties.rho_l(T_operating)
        mu_l = SodiumProperties.mu_l(T_operating)
        mu_v = SodiumProperties.mu_v(T_operating)
        sigma = SodiumProperties.sigma(T_operating)
        k_l = SodiumProperties.k_l(T_operating)
        k_eff = geo.effective_wick_k(k_l)
        R_na = 361.66
        gamma = 1.67

        m_dot_max = float(Q_thermal) / h_fg

        z_vals = [i * dz for i in range(n_nodes)]
        m_dot_vals = []
        u_v_vals = []
        ma_vals = []
        p_v_vals = [0.0] * n_nodes
        p_l_vals = [0.0] * n_nodes
        t_v_vals = [0.0] * n_nodes
        t_wall_vals = [0.0] * n_nodes
        section_types = []

        # 1. Mass flow and velocity
        rho_v_ref = SodiumProperties.rho_v(T_operating, p_ref)
        c_sound_ref = math.sqrt(gamma * R_na * T_operating)

        for z in z_vals:
            if z <= geo.l_evap:
                frac = z / geo.l_evap if geo.l_evap > 0 else 1.0
                m_dot = m_dot_max * frac
                st = 'evaporator'
            elif z <= (geo.l_evap + geo.l_adiab):
                m_dot = m_dot_max
                st = 'adiabatic'
            else:
                z_cond = z - (geo.l_evap + geo.l_adiab)
                frac = max(0.0, 1.0 - z_cond / geo.l_cond)
                m_dot = m_dot_max * frac
                st = 'condenser'
            
            m_dot_vals.append(m_dot)
            section_types.append(st)
            u_v = m_dot / (rho_v_ref * geo.a_vapor)
            u_v_vals.append(u_v)
            ma_vals.append(u_v / c_sound_ref)

        # 2. Vapor pressure drop integration
        p_v_vals[0] = p_ref
        for i in range(1, n_nodes):
            m_mid = 0.5 * (m_dot_vals[i-1] + m_dot_vals[i])
            T_loc = max(600.0, SodiumProperties.t_sat(p_v_vals[i-1]))
            rho_v_loc = SodiumProperties.rho_v(T_loc, p_v_vals[i-1])
            u_v_loc = m_mid / (rho_v_loc * geo.a_vapor) if geo.a_vapor > 0 else 0.0
            
            re_v = max(1.0, (rho_v_loc * abs(u_v_loc) * geo.d_vapor) / mu_v)
            if re_v < 2300:
                f_v = 64.0 / re_v
            else:
                f_v = 0.3164 * (re_v ** -0.25)

            dp_fric = (f_v * rho_v_loc * (u_v_loc**2) / (2.0 * geo.d_vapor)) * dz
            
            if section_types[i] == 'evaporator':
                dp_inert = rho_v_loc * (u_v_vals[i]**2 - u_v_vals[i-1]**2)
                dp_step = dp_fric + dp_inert
            elif section_types[i] == 'adiabatic':
                dp_step = dp_fric
            else:
                # Condenser pressure recovery factor ~ 0.65
                dp_inert = 0.65 * rho_v_loc * (u_v_vals[i]**2 - u_v_vals[i-1]**2)
                dp_step = dp_fric + dp_inert
            
            p_v_vals[i] = max(100.0, p_v_vals[i-1] - dp_step)

        # 3. Liquid pressure drop integration (Darcy flow backwards from condenser tip)
        p_l_vals[-1] = p_v_vals[-1]
        for i in range(n_nodes - 2, -1, -1):
            m_mid = 0.5 * (m_dot_vals[i] + m_dot_vals[i+1])
            dp_fric_l = (mu_l * m_mid / (rho_l * geo.a_wick * geo.permeability)) * dz
            dp_grav_l = rho_l * self.g * math.sin(self.tilt_angle) * dz
            p_l_vals[i] = p_l_vals[i+1] - (dp_fric_l + dp_grav_l)

        # 4. Temperature profiles & capillary margins
        r_wick_cond = math.log(geo.r_inner / geo.r_vapor) / (2.0 * math.pi * k_eff)
        r_wall_cond = math.log(geo.r_outer / geo.r_inner) / (2.0 * math.pi * geo.k_wall)
        r_radial_e = (r_wick_cond + r_wall_cond) / geo.l_evap
        r_radial_c = (r_wick_cond + r_wall_cond) / geo.l_cond

        delta_p_c_avail = (2.0 * sigma * math.cos(geo.contact_angle)) / geo.r_eff
        capillary_margins = []

        for i in range(n_nodes):
            t_v = SodiumProperties.t_sat(p_v_vals[i])
            t_v_vals[i] = t_v
            
            if section_types[i] == 'evaporator':
                t_wall_vals[i] = t_v + Q_thermal * r_radial_e
            elif section_types[i] == 'condenser':
                t_wall_vals[i] = t_v - Q_thermal * r_radial_c
            else:
                t_wall_vals[i] = t_v

            delta_p_req = p_v_vals[i] - p_l_vals[i]
            capillary_margins.append(delta_p_c_avail - delta_p_req)

        total_delta_p_v = p_v_vals[0] - min(p_v_vals)
        total_delta_p_l = p_l_vals[-1] - p_l_vals[0]
        total_p_drop = total_delta_p_v + total_delta_p_l
        min_capillary_margin = min(capillary_margins)
        capillary_margin_ratio = min_capillary_margin / delta_p_c_avail if delta_p_c_avail > 0 else 0.0
        max_mach = max(ma_vals)

        return {
            'z': z_vals,
            'section': section_types,
            'm_dot': m_dot_vals,
            'u_vapor': u_v_vals,
            'mach': ma_vals,
            'P_vapor': p_v_vals,
            'P_liquid': p_l_vals,
            'T_vapor': t_v_vals,
            'T_wall': t_wall_vals,
            'capillary_margin': capillary_margins,
            'summary': {
                'Q_thermal_W': Q_thermal,
                'T_operating_K': T_operating,
                'delta_P_v_Pa': total_delta_p_v,
                'delta_P_l_Pa': total_delta_p_l,
                'total_delta_P_Pa': total_p_drop,
                'delta_P_c_max_Pa': delta_p_c_avail,
                'min_capillary_margin_Pa': min_capillary_margin,
                'capillary_margin_ratio': capillary_margin_ratio,
                'max_mach': max_mach,
                'is_operable': min_capillary_margin > 0.0 and max_mach < 1.0,
                'effective_thermal_resistance_K_W': (max(t_wall_vals) - min(t_wall_vals)) / Q_thermal if Q_thermal > 0 else 0.0
            }
        }


def run_operational_limit_envelope(geo, g=0.0, t_start=750.0, t_end=1250.0, t_steps=51):
    """Generates operating limit curves across temperature range [K]"""
    rom = TwoPhaseHeatPipeROM(geo, g=g)
    results = []
    dt = (t_end - t_start) / (t_steps - 1)
    for i in range(t_steps):
        T = t_start + i * dt
        res = rom.calculate_operational_limits(T)
        results.append(res)
    return results


if __name__ == '__main__':
    print("================================================================================")
    print("Liquid Metal Sodium Heat Pipe Two-Phase Flow Reduced-Order Numerical Simulator")
    print("Author: Muhammad Naufal Tamimi (NIM: 20/456315/TK/50445) - DTNTF FT UGM")
    print("================================================================================")
    
    # Test ISR prototype single heat pipe
    # For ISR 500 kWth with 61 pins:
    # Outer envelope OD = 16 mm (r_outer = 8 mm), inner vapor diameter = 12.4 mm (r_vapor = 6.2 mm)
    geo_isr = HeatPipeGeometry(
        r_outer=0.008,        # 16 mm OD
        t_wall=0.001,         # 1 mm Mo-30Re
        t_wick=0.0008,        # 0.8 mm porous screen
        l_evap=0.25,          # 25 cm active core fuel height (Khandaq et al. 2020)
        l_adiab=0.25,         # 25 cm upper shield/reflector section
        l_cond=0.50,          # 50 cm heat exchanger / Stirling head
        pore_radius=25e-6,
        permeability=1.5e-10,
        wall_mat='Mo-30Re',
        k_wall=65.0
    )
    
    rom_isr = TwoPhaseHeatPipeROM(geo_isr, g=0.0) # Microgravity g=0
    limits_1050K = rom_isr.calculate_operational_limits(1050.0)
    print(f"\n--- ISR Heat Pipe Operational Limits at Nominal Sodium Temperature (1050 K / 777 C) ---")
    print(f"Capillary Limit  : {limits_1050K['q_cap']/1e3:.2f} kWth")
    print(f"Sonic Limit      : {limits_1050K['q_sonic']/1e3:.2f} kWth")
    print(f"Entrainment Limit: {limits_1050K['q_ent']/1e3:.2f} kWth")
    print(f"Viscous Limit    : {limits_1050K['q_visc']/1e3:.2f} kWth")
    print(f"Boiling Limit    : {limits_1050K['q_boil']/1e3:.2f} kWth")
    print(f"Envelope Limit   : {limits_1050K['q_envelope']/1e3:.2f} kWth (Governed by {limits_1050K['governing_limit']})")
    
    sim = rom_isr.simulate_axial_profiles(Q_thermal=8200.0, T_operating=1050.0, n_nodes=50) # 8.2 kW nominal per pin
    print(f"\n--- Simulation at Nominal Pin Power Q = 8.2 kWth (Microgravity g=0) ---")
    print(f"Vapor Pressure Drop  : {sim['summary']['delta_P_v_Pa']:.1f} Pa")
    print(f"Liquid Pressure Drop : {sim['summary']['delta_P_l_Pa']:.1f} Pa")
    print(f"Total Pressure Drop  : {sim['summary']['total_delta_P_Pa']:.1f} Pa")
    print(f"Max Capillary Head   : {sim['summary']['delta_P_c_max_Pa']:.1f} Pa")
    print(f"Capillary Margin     : {sim['summary']['min_capillary_margin_Pa']:.1f} Pa ({sim['summary']['capillary_margin_ratio']*100:.1f}%)")
    print(f"Max Vapor Mach No.   : {sim['summary']['max_mach']:.4f}")
    print(f"Effective Resistance : {sim['summary']['effective_thermal_resistance_K_W']:.5f} K/W")
    print(f"Operability Status   : {'PASSED (SAFE)' if sim['summary']['is_operable'] else 'FAILED'}")
