#!/usr/bin/env python3
"""
================================================================================
    Script para reproducir la Figura 8 de Winske & Leroy (1984)
    Lee datos del código Fortran hyb1d_rrkfft.f90
    
    DETECTA AUTOMÁTICAMENTE el formato del archivo energy.d12
    
    Genera 3 figuras:
    1. Figura 8 (4 paneles, escala log en Wf)
    2. Figura 8 (4 paneles, escala lineal)
    3. Energías globales del sistema (1 panel)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import struct

IHPARM = 18


def read_header_fortran(filename, nsp_guess=2):
    """Lee el header del archivo binario Fortran."""
    
    with open(filename, 'rb') as f:
        # Leer record marker inicial
        rec_marker = struct.unpack('i', f.read(4))[0]
        
        # Determinar nsp basado en el tamaño del header
        # header tiene IHPARM + 8*nsp floats
        # Probar diferentes combinaciones
        for nsp in [2, 1]:
            nheadr = IHPARM + 8 * nsp
            for float_size in [4, 8]:  # float32 o float64
                expected_size = nheadr * float_size
                if rec_marker == expected_size:
                    print(f"Header detectado: nsp={nsp}, float{float_size*8}")
                    f.seek(4)  # Volver al inicio de datos
                    
                    dtype = np.float32 if float_size == 4 else np.float64
                    header = np.frombuffer(f.read(nheadr * float_size), dtype=dtype)
                    
                    # Verificar marker final
                    rec_end = struct.unpack('i', f.read(4))[0]
                    
                    params = parse_header(header, nsp)
                    params['float_size'] = float_size
                    params['dtype'] = dtype
                    header_bytes = 4 + nheadr * float_size + 4
                    
                    return params, header_bytes
        
        raise ValueError(f"No se pudo determinar formato del header. Record marker = {rec_marker}")


def parse_header(header, nsp):
    """Parsea el header en un diccionario."""
    params = {
        'dt': float(header[0]),
        'dx': float(header[1]),
        'nx': int(header[2]),
        'itmax': int(header[3]),
        'lfld': int(header[4]),
        'nsp': int(header[5]),
        'np': int(header[6]),
        'npg': int(header[7]),
        'betaen': float(header[8]),
        'wpiwci': float(header[9]),
        'bf0': [float(header[10]), float(header[11]), float(header[12])],
        'ifield': int(header[13]),
        'iparticles': int(header[14]),
        'ienergy': int(header[15]),
        'ifilter': int(header[16]),
        'itrestart': int(header[17]),
    }
    
    for is_ in range(nsp):
        indh = IHPARM + 8 * is_
        params[f'rdn_{is_+1}'] = float(header[indh])
        params[f'vdr_{is_+1}'] = float(header[indh + 1])
        params[f'betain_{is_+1}'] = float(header[indh + 2])
        params[f'anis_{is_+1}'] = float(header[indh + 3])
        params[f'qi_{is_+1}'] = float(header[indh + 4])
        params[f'ai_{is_+1}'] = float(header[indh + 5])
        params[f'nions_{is_+1}'] = int(header[indh + 6])
        params[f'gg_{is_+1}'] = float(header[indh + 7])
    
    params['xlen'] = params['dx'] * params['nx']
    
    return params


def detect_energy_format(filename, header_bytes, nsp, float_size):
    """
    Detecta el formato del registro de energía.
    
    Formatos posibles:
    1. Original: tim(2), tpal(nsp), tper(nsp), umx(nsp), enrgys(6)
       = 2 + nsp + nsp + nsp + 6 = 2 + 3*nsp + 6
    
    2. Con escalares: tim(2), tpal(nsp), tper(nsp), wkpal_total, wkper_total, umx(nsp), enrgys(6)
       = 2 + nsp + nsp + 1 + 1 + nsp + 6 = 4 + 3*nsp + 6
    
    3. Con arrays: tim(2), tpal(nsp), tper(nsp), wkpal(nsp), wkper(nsp), umx(nsp), enrgys(6)
       = 2 + 5*nsp + 6
    """
    with open(filename, 'rb') as f:
        f.seek(header_bytes)
        rec_marker = struct.unpack('i', f.read(4))[0]
    
    record_bytes = rec_marker
    record_floats = record_bytes // float_size
    
    # Determinar formato
    formats = {
        'original': 2 + 3*nsp + 6,           # tim, tpal, tper, umx, enrgys
        'scalar_wk': 2 + 3*nsp + 2 + 6,      # + wkpal_total, wkper_total
        'array_wk': 2 + 5*nsp + 6,           # + wkpal(nsp), wkper(nsp)
    }
    
    detected_format = None
    for fmt_name, expected_floats in formats.items():
        if record_floats == expected_floats:
            detected_format = fmt_name
            break
    
    if detected_format is None:
        print(f"ADVERTENCIA: Formato no reconocido. record_floats={record_floats}")
        print(f"  Esperados: {formats}")
        # Intentar adivinar
        if record_floats > formats['array_wk']:
            detected_format = 'array_wk'
        elif record_floats > formats['scalar_wk']:
            detected_format = 'scalar_wk'
        else:
            detected_format = 'original'
    
    print(f"Formato de energía detectado: '{detected_format}' ({record_floats} floats = {record_bytes} bytes)")
    
    return detected_format, record_bytes


def read_energy_file_fortran(filename):
    """Lee el archivo de energía con detección automática de formato."""
    
    # Leer header
    params, header_bytes = read_header_fortran(filename)
    nsp = params['nsp']
    float_size = params['float_size']
    dtype = params['dtype']
    
    # Detectar formato de registros de energía
    energy_format, record_bytes = detect_energy_format(filename, header_bytes, nsp, float_size)
    record_total = 4 + record_bytes + 4  # con markers
    
    # Calcular número de registros
    file_size = Path(filename).stat().st_size
    data_bytes = file_size - header_bytes
    n_records = data_bytes // record_total
    
    print(f"\nLeyendo archivo: {filename}")
    print(f"  nsp = {nsp}, dtype = {dtype}")
    print(f"  Header: {header_bytes} bytes")
    print(f"  Registros: {n_records}")
    print(f"  Formato: {energy_format}")
    
    # Inicializar arrays
    time = np.zeros(n_records)
    step = np.zeros(n_records, dtype=int)
    tpal = np.zeros((n_records, nsp))
    tper = np.zeros((n_records, nsp))
    umx = np.zeros((n_records, nsp))
    enrgys = np.zeros((n_records, 6))
    
    # Arrays adicionales según formato
    if energy_format == 'scalar_wk':
        wkpal_total = np.zeros(n_records)
        wkper_total = np.zeros(n_records)
        wkpal = None
        wkper = None
    elif energy_format == 'array_wk':
        wkpal = np.zeros((n_records, nsp))
        wkper = np.zeros((n_records, nsp))
        wkpal_total = None
        wkper_total = None
    else:
        wkpal = None
        wkper = None
        wkpal_total = None
        wkper_total = None
    
    records_read = 0
    with open(filename, 'rb') as f:
        f.seek(header_bytes)
        
        for i in range(n_records):
            # Leer record marker
            marker_data = f.read(4)
            if len(marker_data) < 4:
                break
            rec_start = struct.unpack('i', marker_data)[0]
            
            if rec_start != record_bytes:
                print(f"ADVERTENCIA registro {i}: marker={rec_start}, esperado={record_bytes}")
                break
            
            # Leer datos
            data_raw = np.frombuffer(f.read(record_bytes), dtype=dtype)
            
            # Parsear según formato
            idx = 0
            step[i] = int(data_raw[idx])
            time[i] = data_raw[idx + 1]
            idx += 2
            
            tpal[i, :] = data_raw[idx:idx + nsp]
            idx += nsp
            
            tper[i, :] = data_raw[idx:idx + nsp]
            idx += nsp
            
            if energy_format == 'scalar_wk':
                wkpal_total[i] = data_raw[idx]
                wkper_total[i] = data_raw[idx + 1]
                idx += 2
            elif energy_format == 'array_wk':
                wkpal[i, :] = data_raw[idx:idx + nsp]
                idx += nsp
                wkper[i, :] = data_raw[idx:idx + nsp]
                idx += nsp
            
            umx[i, :] = data_raw[idx:idx + nsp]
            idx += nsp
            
            enrgys[i, :] = data_raw[idx:idx + 6]
            
            # Leer marker final
            rec_end = struct.unpack('i', f.read(4))[0]
            
            records_read += 1
    
    print(f"  Registros leídos: {records_read}")
    
    # Truncar al número real de registros
    time = time[:records_read]
    step = step[:records_read]
    tpal = tpal[:records_read]
    tper = tper[:records_read]
    umx = umx[:records_read]
    enrgys = enrgys[:records_read]
    
    data = {
        'time': time,
        'step': step,
        'tpal': tpal,
        'tper': tper,
        'umx': umx,
        'enrgys': enrgys,
        'format': energy_format,
    }
    
    if energy_format == 'scalar_wk':
        data['wkpal_total'] = wkpal_total[:records_read]
        data['wkper_total'] = wkper_total[:records_read]
    elif energy_format == 'array_wk':
        data['wkpal'] = wkpal[:records_read]
        data['wkper'] = wkper[:records_read]
    
    return params, data


def plot_figure8_fortran(params, data, output_file="figure8_winske_fortran.png"):
    """Genera la Figura 8 de Winske & Leroy (1984)."""
    
    t = data['time']
    enrgys = data['enrgys']
    nsp = params['nsp']
    energy_format = data['format']
    
    # =========================================================================
    # EXTRAER ENERGÍAS SEGÚN FORMATO
    # =========================================================================
    
    if energy_format == 'array_wk':
        # Tenemos energías por especie
        wkpal = data['wkpal']
        wkper = data['wkper']
        
        if nsp >= 2:
            Wm_par = wkpal[:, 0]
            Wm_perp = 2.0 * wkper[:, 0]  # Corregir factor 2
            Wb_par = wkpal[:, 1]
            Wb_perp = 2.0 * wkper[:, 1]
        else:
            Wm_par = wkpal[:, 0]
            Wm_perp = 2.0 * wkper[:, 0]
            Wb_par = np.zeros_like(Wm_par)
            Wb_perp = np.zeros_like(Wm_par)
            
    elif energy_format == 'scalar_wk':
        # Solo tenemos totales, estimamos por especie usando enrgys
        wkpal_total = data['wkpal_total']
        wkper_total = data['wkper_total']
        
        # enrgys[2] = wkpalenrgy (total paralela)
        # enrgys[3] = wkperenrgy (total perpendicular)
        # Estos ya son los totales correctos
        
        # Para nsp=2, necesitamos estimar la división
        # Usamos la aproximación de que inicialmente Wb >> Wm para el beam
        if nsp >= 2:
            # Usar parámetros para estimar división inicial
            rdn_m = params.get('rdn_1', 0.95)
            rdn_b = params.get('rdn_2', 0.05)
            vdr_b = params.get('vdr_2', 0)
            beta_m = params.get('betain_1', 1.0)
            beta_b = params.get('betain_2', 1.0)
            
            # Energía térmica inicial proporcional a rdn * beta
            Wth_m = rdn_m * beta_m
            Wth_b = rdn_b * beta_b
            
            # Energía de drift del beam: (1/2) * rdn_b * vdr^2
            Wdrift_b = 0.5 * rdn_b * vdr_b**2
            
            # Fracción inicial del beam en energía paralela
            frac_b_par = (Wth_b/2 + Wdrift_b) / (Wth_m/2 + Wth_b/2 + Wdrift_b) if (Wth_m + Wth_b + 2*Wdrift_b) > 0 else 0.5
            frac_b_perp = Wth_b / (Wth_m + Wth_b) if (Wth_m + Wth_b) > 0 else 0.5
            
            # Esta es una aproximación gruesa - los valores evolucionarán
            Wb_par = frac_b_par * wkpal_total
            Wm_par = (1 - frac_b_par) * wkpal_total
            Wb_perp = frac_b_perp * 2.0 * wkper_total  # Factor 2 para corrección
            Wm_perp = (1 - frac_b_perp) * 2.0 * wkper_total
            
            print(f"\nADVERTENCIA: Formato 'scalar_wk' - estimando energías por especie")
            print(f"  frac_beam_par inicial ≈ {frac_b_par:.3f}")
            print(f"  frac_beam_perp inicial ≈ {frac_b_perp:.3f}")
        else:
            Wm_par = wkpal_total
            Wm_perp = 2.0 * wkper_total
            Wb_par = np.zeros_like(Wm_par)
            Wb_perp = np.zeros_like(Wm_par)
            
    else:  # 'original'
        # No tenemos energías cinéticas separadas, usar enrgys
        # enrgys = [bfldenrgy, efldenrgy, wkpalenrgy, wkperenrgy, wktotenrgy, totenrgy]
        
        wkpal_total = enrgys[:, 2]  # wkpalenrgy
        wkper_total = enrgys[:, 3]  # wkperenrgy
        
        if nsp >= 2:
            # Misma estimación que arriba
            rdn_m = params.get('rdn_1', 0.95)
            rdn_b = params.get('rdn_2', 0.05)
            vdr_b = params.get('vdr_2', 0)
            
            frac_b = rdn_b / (rdn_m + rdn_b)
            
            # Para el drift, toda la energía de drift está en el beam
            Wdrift = 0.5 * rdn_b * vdr_b**2
            
            Wb_par = frac_b * wkpal_total  # Aproximación
            Wm_par = (1 - frac_b) * wkpal_total
            Wb_perp = frac_b * 2.0 * wkper_total
            Wm_perp = (1 - frac_b) * 2.0 * wkper_total
            
            print(f"\nADVERTENCIA: Formato 'original' - no hay energías por especie")
            print(f"  Usando totales de enrgys[2:3]")
        else:
            Wm_par = wkpal_total
            Wm_perp = 2.0 * wkper_total
            Wb_par = np.zeros_like(Wm_par)
            Wb_perp = np.zeros_like(Wm_par)
    
    # Energías totales
    Wm = Wm_par + Wm_perp
    Wb = Wb_par + Wb_perp
    Wf = enrgys[:, 0]  # bfldenrgy
    
    # =========================================================================
    # NORMALIZACIÓN
    # =========================================================================
    xlen = params['xlen']
    W0_mag = xlen / 2.0
    W0_kin = Wm[0] + Wb[0]
    W0 = W0_kin if W0_kin > 0 else 1.0
    
    Wm_par_norm = Wm_par / W0
    Wm_perp_norm = Wm_perp / W0
    Wb_par_norm = Wb_par / W0 if nsp >= 2 else np.zeros_like(Wm_par_norm)
    Wb_perp_norm = Wb_perp / W0 if nsp >= 2 else np.zeros_like(Wm_par_norm)
    Wm_norm = Wm / W0
    Wb_norm = Wb / W0 if nsp >= 2 else np.zeros_like(Wm_norm)
    Wf_norm_mag = Wf / W0_mag
    
    # =========================================================================
    # DIAGNÓSTICO
    # =========================================================================
    print("\n" + "="*70)
    print("DIAGNÓSTICO DE ENERGÍAS - Código FORTRAN")
    print("="*70)
    print(f"Formato detectado: {energy_format}")
    print(f"Parámetros: xlen={params['xlen']:.1f}, dt={params['dt']:.4f}, nsp={nsp}")
    
    if nsp >= 2:
        print(f"\nEspecie 1 (Main): rdn={params.get('rdn_1', 'N/A')}, vdr={params.get('vdr_1', 'N/A')}")
        print(f"Especie 2 (Beam): rdn={params.get('rdn_2', 'N/A')}, vdr={params.get('vdr_2', 'N/A')}")
    
    print(f"\nRango temporal: t ∈ [{t[0]:.1f}, {t[-1]:.1f}]")
    print(f"W₀_kin = {W0_kin:.2e}, W₀_mag = {W0_mag:.1f}")
    
    print(f"\nEnergías iniciales (normalizadas):")
    print(f"  W_m∥/W₀ = {Wm_par_norm[0]:.4f}, W_m⊥/W₀ = {Wm_perp_norm[0]:.4f}")
    if nsp >= 2:
        print(f"  W_b∥/W₀ = {Wb_par_norm[0]:.4f}, W_b⊥/W₀ = {Wb_perp_norm[0]:.4f}")
    print(f"  W_f/W₀_mag = {Wf_norm_mag[0]:.2e}")
    
    print(f"\nEnergías finales:")
    print(f"  W_m∥/W₀ = {Wm_par_norm[-1]:.4f}, W_m⊥/W₀ = {Wm_perp_norm[-1]:.4f}")
    if nsp >= 2:
        print(f"  W_b∥/W₀ = {Wb_par_norm[-1]:.4f}, W_b⊥/W₀ = {Wb_perp_norm[-1]:.4f}")
    print(f"  W_f/W₀_mag max = {np.max(Wf_norm_mag):.2e}")
    print("="*70 + "\n")
    
    # =========================================================================
    # FIGURA 1: 4 paneles (escala log en Wf)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    colors = {'par': 'blue', 'perp': 'red', 'Wm': 'green', 'Wb': 'orange', 'Wf': 'purple'}
    
    # Panel (a): Beam
    ax = axes[0, 0]
    if nsp >= 2 and np.any(Wb_par_norm > 0):
        ax.plot(t, Wb_par_norm, color=colors['par'], lw=2, label=r'$W_{b\parallel}$')
        ax.plot(t, Wb_perp_norm, color=colors['perp'], lw=2, label=r'$W_{b\perp}$')
        ax.legend(loc='best')
    else:
        ax.text(0.5, 0.5, 'N/A (nsp=1)', ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel(r'$\Omega_i t$')
    ax.set_ylabel(r'$W / W_0$')
    ax.set_title('(a) Beam Energies')
    ax.grid(True, alpha=0.3)
    
    # Panel (b): Main ions
    ax = axes[0, 1]
    ax.plot(t, Wm_par_norm, color=colors['par'], lw=2, label=r'$W_{m\parallel}$')
    ax.plot(t, Wm_perp_norm, color=colors['perp'], lw=2, label=r'$W_{m\perp}$')
    ax.set_xlabel(r'$\Omega_i t$')
    ax.set_ylabel(r'$W / W_0$')
    ax.set_title('(b) Main Ion Energies')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel (c): Totales
    ax = axes[1, 0]
    if nsp >= 2:
        ax.plot(t, Wb_norm, color=colors['Wb'], lw=2, label=r'$W_b$')
    ax.plot(t, Wm_norm, color=colors['Wm'], lw=2, label=r'$W_m$')
    Wtot_norm = (Wm + Wb + Wf) / W0
    ax.plot(t, Wtot_norm, 'k--', lw=1.5, label=r'$W_{tot}$')
    ax.set_xlabel(r'$\Omega_i t$')
    ax.set_ylabel(r'$W / W_0$')
    ax.set_title('(c) Total Species Energies')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel (d): Fluctuaciones magnéticas
    ax = axes[1, 1]
    ax.semilogy(t, Wf_norm_mag, color=colors['Wf'], lw=2, label=r'$W_f$')
    ax.set_xlabel(r'$\Omega_i t$')
    ax.set_ylabel(r'$W_f / W_{0,mag}$')
    ax.set_title('(d) Magnetic Fluctuation Energy')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'Energy Histories - Winske & Leroy (1984)\nFortran data (format: {energy_format})', 
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"Figura guardada: {output_file}")
    
    # =========================================================================
    # FIGURA 2: Versión lineal
    # =========================================================================
    output_linear = output_file.replace('.png', '_linear.png')
    axes[1, 1].set_yscale('linear')
    plt.savefig(output_linear, dpi=200, bbox_inches='tight')
    print(f"Figura guardada: {output_linear}")
    
    return fig


def plot_global_energies(params, data, output_file="figure_global_energies_fortran.png"):
    """
    Genera una figura con las energías globales del sistema.
    
    Contenido de enrgys según el código Fortran:
        enrgys[:,0] = bfldenrgy   → Energía magnética de fluctuaciones (δB²/8π)
        enrgys[:,1] = efldenrgy   → Energía eléctrica (E²/8π)
        enrgys[:,2] = wkpalenrgy  → Energía cinética paralela total
        enrgys[:,3] = wkperenrgy  → Energía cinética perpendicular total
        enrgys[:,4] = wktotenrgy  → Energía cinética total
        enrgys[:,5] = totenrgy    → Energía total del sistema
    """
    t = data['time']
    enrgys = data['enrgys']
    
    # Extraer cada componente
    W_B = enrgys[:, 0]      # Energía magnética de fluctuaciones
    W_E = enrgys[:, 1]      # Energía eléctrica
    K_par = enrgys[:, 2]    # Energía cinética paralela total
    K_perp = enrgys[:, 3]   # Energía cinética perpendicular total
    K_tot = enrgys[:, 4]    # Energía cinética total
    E_tot = enrgys[:, 5]    # Energía total del sistema
    
    # =========================================================================
    # FIGURA 3: Energías globales (1 panel)
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(t, W_B, lw=2, label=r'$W_B$ (magnética)', color='blue')
    ax.plot(t, W_E, lw=2, label=r'$W_E$ (eléctrica)', color='cyan')
    ax.plot(t, K_par, lw=2, label=r'$K_\parallel$ (cin. paralela)', color='red')
    ax.plot(t, K_perp, lw=2, label=r'$K_\perp$ (cin. perpendicular)', color='orange')
    ax.plot(t, K_tot, lw=2, label=r'$K_{tot}$ (cin. total)', color='green', linestyle='--')
    ax.plot(t, E_tot, lw=2, label=r'$E_{tot}$ (total sistema)', color='black', linestyle='-.')
    
    ax.set_xlabel(r'$\Omega_i t$', fontsize=14)
    ax.set_ylabel('Energía', fontsize=14)
    ax.set_title('Energías Globales del Sistema (Fortran)', fontsize=16, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"Figura guardada: {output_file}")
    
    # =========================================================================
    # DIAGNÓSTICO DE ENERGÍAS GLOBALES
    # =========================================================================
    print("\n" + "="*50)
    print("DIAGNÓSTICO DE ENERGÍAS GLOBALES (Fortran)")
    print("="*50)
    print("Energías iniciales:")
    print(f"  W_B(0)    = {W_B[0]:.6e}")
    print(f"  W_E(0)    = {W_E[0]:.6e}")
    print(f"  K_∥(0)    = {K_par[0]:.6e}")
    print(f"  K_⊥(0)    = {K_perp[0]:.6e}")
    print(f"  K_tot(0)  = {K_tot[0]:.6e}")
    print(f"  E_tot(0)  = {E_tot[0]:.6e}")
    print("\nEnergías finales:")
    print(f"  W_B(f)    = {W_B[-1]:.6e}")
    print(f"  W_E(f)    = {W_E[-1]:.6e}")
    print(f"  K_∥(f)    = {K_par[-1]:.6e}")
    print(f"  K_⊥(f)    = {K_perp[-1]:.6e}")
    print(f"  K_tot(f)  = {K_tot[-1]:.6e}")
    print(f"  E_tot(f)  = {E_tot[-1]:.6e}")
    print("\nConservación de energía:")
    dE = (E_tot[-1] - E_tot[0]) / E_tot[0] * 100
    print(f"  ΔE/E₀ = {dE:.4f} %")
    print("="*50 + "\n")
    
    return fig


def main():
    import sys
    
    energy_file = sys.argv[1] if len(sys.argv) > 1 else "energy.d12"
    
    if not Path(energy_file).exists():
        print(f"ERROR: No se encontró '{energy_file}'")
        print("Uso: python3 plotenergy_beam.py [energy.d12]")
        return
    
    params, data = read_energy_file_fortran(energy_file)
    
    # Generar las 3 figuras
    plot_figure8_fortran(params, data)
    plot_global_energies(params, data)
    
    plt.show()
    
    print("\n¡Figuras generadas!")


if __name__ == "__main__":
    main()
