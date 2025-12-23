#!/usr/bin/env julia
#=
================================================================================
    Script para reproducir la Figura 8 de Winske & Leroy (1984)
    "Diffuse Ions Produced by Electromagnetic Ion Beam Instabilities"
    J. Geophys. Res., Vol. 89, No. A5, pp. 2673-2688
    
    Figura 8: Historias temporales de energía para modo resonante
    
    Panel (a): W_b∥ y W_b⊥  - Energías cinéticas del beam (paralela y perpendicular)
    Panel (b): W_m∥ y W_m⊥  - Energías cinéticas de los iones principales
    Panel (c): W_b y W_m    - Energías cinéticas totales por especie
    Panel (d): W_f          - Energía de fluctuaciones magnéticas
    
    Todas las energías normalizadas por W_0 = B_0² L / 8π
================================================================================
=#





using Printf # Para escritura científica con cifras significativas y notación científica, y uso de latex para labels
using CairoMakie # Librería de ploteo. Tiene mejores detalles que Plots, y usa comandos similares.

const nsp = 2 #número de partículas
const ihparm = 18 # bits en binario del header de la simulación. Basado en el código Fortran original. Se eliminará en una futura versión
const nheadr = ihparm + 8 * nsp # dimensiones del header

"""
    read_header(filename) -> Dict

Lee el header del archivo binario y retorna un diccionario con los parámetros.
"""

function read_header(filename::String)
    header = zeros(Float64, nheadr)
    open(filename, "r") do f
        read!(f, header)
    end
    
    params = Dict{String, Any}(
        "dt" => header[1],
        "dx" => header[2],
        "nx" => Int(header[3]),
        "itmax" => Int(header[4]),
        "lfld" => Int(header[5]),
        "nsp" => Int(header[6]),
        "np" => Int(header[7]),
        "npg" => Int(header[8]),
        "betaen" => header[9],
        "wpiwci" => header[10],
        "bf0" => [header[11], header[12], header[13]],
        "ifield" => Int(header[14]),
        "iparticles" => Int(header[15]),
        "ienergy" => Int(header[16]),
        "ifilter" => Int(header[17]),
        "itrestart" => Int(header[18])
    )
    
    for is in 1:nsp
        indh = ihparm + 8 * (is - 1)
        params["rdn_$is"] = header[indh + 1]
        params["vdr_$is"] = header[indh + 2]
        params["betain_$is"] = header[indh + 3]
        params["anis_$is"] = header[indh + 4]
        params["qi_$is"] = header[indh + 5]
        params["ai_$is"] = header[indh + 6]
        params["nions_$is"] = Int(header[indh + 7])
        params["gg_$is"] = header[indh + 8]
    end
    
    params["xlen"] = params["dx"] * params["nx"]
    
    return params
end

function read_energy_file(filename::String)
    params = read_header(filename)
    nsp_file = params["nsp"]
    
    record_size = 2 + nsp_file + nsp_file + nsp_file + 6 + nsp_file + nsp_file
    
    filesize_bytes = filesize(filename)
    header_bytes = nheadr * 8
    data_bytes = filesize_bytes - header_bytes
    n_records = div(data_bytes, record_size * 8)
    
    println("Leyendo archivo: $filename")
    println("  Registros: $n_records")
    
    time = zeros(Float64, n_records)
    step = zeros(Int, n_records)
    tpal = zeros(Float64, n_records, nsp_file)
    tper = zeros(Float64, n_records, nsp_file)
    umx = zeros(Float64, n_records, nsp_file)
    enrgys = zeros(Float64, n_records, 6)
    wkpal_sp = zeros(Float64, n_records, nsp_file)
    wkper_sp = zeros(Float64, n_records, nsp_file)
    
    open(filename, "r") do f
        seek(f, header_bytes)
        
        for i in 1:n_records
            tim = zeros(Float64, 2)
            read!(f, tim)
            step[i] = Int(tim[1])
            time[i] = tim[2]
            
            tpal_rec = zeros(Float64, nsp_file)
            read!(f, tpal_rec)
            tpal[i, :] = tpal_rec
            
            tper_rec = zeros(Float64, nsp_file)
            read!(f, tper_rec)
            tper[i, :] = tper_rec
            
            umx_rec = zeros(Float64, nsp_file)
            read!(f, umx_rec)
            umx[i, :] = umx_rec
            
            enrgys_rec = zeros(Float64, 6)
            read!(f, enrgys_rec)
            enrgys[i, :] = enrgys_rec
            
            wkpal_rec = zeros(Float64, nsp_file)
            read!(f, wkpal_rec)
            wkpal_sp[i, :] = wkpal_rec
            
            wkper_rec = zeros(Float64, nsp_file)
            read!(f, wkper_rec)
            wkper_sp[i, :] = wkper_rec
        end
    end
    
    data = Dict(
        "time" => time,
        "step" => step,
        "tpal" => tpal,
        "tper" => tper,
        "umx" => umx,
        "enrgys" => enrgys,
        "wkpal_sp" => wkpal_sp,
        "wkper_sp" => wkper_sp
    )
    
    return params, data
end

"""
    plot_figure8_corrected(params, data; output_file="figure8_winske.png")

Versión CORREGIDA con:
1. Normalización correcta: W₀ = energía cinética total inicial
2. Factor 2 aplicado a wkper_sp (corrige el bug del código original de Viñas)
"""
function plot_figure8_corrected(params, data; output_file="figure8_winske_corrected.png")
    t = data["time"]
    wkpal_sp = data["wkpal_sp"]
    wkper_sp = data["wkper_sp"]
    enrgys = data["enrgys"]
    
    # =========================================================================
    # CORRECCIÓN 1: Factor 2 en energía perpendicular
    # El código de simulación calcula wkper_sp = 0.5 * (suma de vy² + vz²) * fac2
    # La energía perpendicular real es el doble
    # =========================================================================
    Wm_par = wkpal_sp[:, 1]           # W_m∥ (correcta)
    Wm_perp = 2.0 .* wkper_sp[:, 1]   # W_m⊥ (CORREGIDA: factor 2)
    Wb_par = wkpal_sp[:, 2]           # W_b∥ (correcta)
    Wb_perp = 2.0 .* wkper_sp[:, 2]   # W_b⊥ (CORREGIDA: factor 2)
    
    # Energías totales por especie (ya no necesitan factor 2 adicional)
    Wm = Wm_par .+ Wm_perp   # W_m total
    Wb = Wb_par .+ Wb_perp   # W_b total
    
    # Energía de fluctuaciones magnéticas
    Wf = enrgys[:, 1]
    
    # =========================================================================
    # CORRECCIÓN 2: Normalización apropiada
    # Opción A: W₀ = B₀²L/8π ≈ xlen/2 (en unidades normalizadas)
    # Opción B: W₀ = energía cinética total inicial (más útil para comparar)
    # =========================================================================
    xlen = params["xlen"]
    
    # Opción A: Energía magnética del campo de fondo
    W0_mag = xlen / 2.0
    
    # Opción B: Energía cinética total inicial (recomendada)
    W0_kin = Wm[1] + Wb[1]
    
    # Usamos W0_kin para que las energías cinéticas sean O(1)
    W0 = W0_kin
    
    # También calculamos Wf normalizada por W0_mag para comparar con Winske
    Wf_norm_mag = Wf ./ W0_mag
    
    # Normalizar todas las energías
    Wm_par_norm = Wm_par ./ W0
    Wm_perp_norm = Wm_perp ./ W0
    Wb_par_norm = Wb_par ./ W0
    Wb_perp_norm = Wb_perp ./ W0
    Wm_norm = Wm ./ W0
    Wb_norm = Wb ./ W0
    Wf_norm = Wf ./ W0
    
    # =========================================================================
    # DIAGNÓSTICO
    # =========================================================================
    println("\n" * "="^70)
    println("DIAGNÓSTICO DE ENERGÍAS (Figura 8 - Winske & Leroy 1984)")
    println("="^70)
    println("Parámetros de simulación:")
    println("  xlen = $(params["xlen"]), dx = $(params["dx"]), nx = $(params["nx"])")
    println("  dt = $(params["dt"])")
    
    println("\nParámetros por especie:")
    for is in 1:nsp
        name = is == 1 ? "Main ions" : "Beam"
        println("  Especie $is ($name):")
        println("    rdn = $(params["rdn_$is"]), vdr = $(params["vdr_$is"]), β = $(params["betain_$is"])")
    end
    
    # Parámetros de Winske
    f = params["rdn_2"]  # fracción del beam
    V = params["vdr_2"]  # drift del beam en V_A
    F = f * V^2 / 2      # parámetro F de Winske
    println("\nParámetros de Winske:")
    println("  f (beam fraction) = $f")
    println("  V (beam drift/V_A) = $V")
    println("  F = fV²/2 = $F")
    println("  Criterio no-resonante (f/2)^(2/3)*V = $((f/2)^(2/3) * V)")
    
    println("\nEnergías de referencia:")
    println("  W₀_mag = B₀²L/8π ≈ xlen/2 = $W0_mag")
    println("  W₀_kin = Wm(0) + Wb(0) = $W0_kin")
    println("  W₀ usado = $W0")
    
    println("\nEnergías iniciales (normalizadas por W₀_kin):")
    println("  W_b∥/W₀ = $(Wb_par_norm[1])")
    println("  W_b⊥/W₀ = $(Wb_perp_norm[1])")
    println("  W_m∥/W₀ = $(Wm_par_norm[1])")
    println("  W_m⊥/W₀ = $(Wm_perp_norm[1])")
    println("  W_f/W₀_mag = $(Wf_norm_mag[1])")
    
    println("\nEnergías finales (normalizadas):")
    println("  W_b∥/W₀ = $(Wb_par_norm[end])")
    println("  W_b⊥/W₀ = $(Wb_perp_norm[end])")
    println("  W_m∥/W₀ = $(Wm_par_norm[end])")
    println("  W_m⊥/W₀ = $(Wm_perp_norm[end])")
    println("  W_f/W₀_mag = $(Wf_norm_mag[end]) (máx = $(maximum(Wf_norm_mag)))")
    
    # Predicción teórica de δB/B₀ máximo
    dB_B0_pred = sqrt(f) * V
    println("\nPredicción teórica (Winske eq. 16):")
    println("  (δB/B₀)_max ≈ √f × V = $(dB_B0_pred)")
    println("  W_f/W₀_mag max ≈ (δB/B₀)² = $(dB_B0_pred^2)")
    
    # Conservación de energía
    Etot_init = Wm[1] + Wb[1] + Wf[1]
    Etot_final = Wm[end] + Wb[end] + Wf[end]
    println("\nConservación de energía:")
    println("  E_tot(0) = $Etot_init")
    println("  E_tot(final) = $Etot_final")
    println("  ΔE/E = $((Etot_final - Etot_init)/Etot_init * 100) %")
    println("="^70 * "\n")
    
    # =========================================================================
    # CREAR FIGURA
    # =========================================================================
    fig = Figure(size = (1000, 900), fontsize = 14)
    
    color_par = :blue
    color_perp = :red
    color_total_m = :green
    color_total_b = :orange
    color_Wf = :purple
    
    # Panel (a): Beam Energies
    ax_a = Axis(fig[1, 1],
        xlabel = "Ωᵢt",
        ylabel = "W / W₀",
        title = "(a) Beam Energies"
    )
    lines!(ax_a, t, Wb_par_norm, color = color_par, linewidth = 2, label = "W_b∥")
    lines!(ax_a, t, Wb_perp_norm, color = color_perp, linewidth = 2, label = "W_b⊥")
    axislegend(ax_a, position = :rt)
    
    # Panel (b): Main Ion Energies
    ax_b = Axis(fig[1, 2],
        xlabel = "Ωᵢt",
        ylabel = "W / W₀",
        title = "(b) Main Ion Energies"
    )
    lines!(ax_b, t, Wm_par_norm, color = color_par, linewidth = 2, label = "W_m∥")
    lines!(ax_b, t, Wm_perp_norm, color = color_perp, linewidth = 2, label = "W_m⊥")
    axislegend(ax_b, position = :rb)
    
    # Panel (c): Total Species Energies
    ax_c = Axis(fig[2, 1],
        xlabel = "Ωᵢt",
        ylabel = "W / W₀",
        title = "(c) Total Species Energies"
    )
    lines!(ax_c, t, Wb_norm, color = color_total_b, linewidth = 2, label = "W_b")
    lines!(ax_c, t, Wm_norm, color = color_total_m, linewidth = 2, label = "W_m")
    # Añadir energía total
    Wtot_norm = (Wm .+ Wb .+ Wf) ./ W0
    lines!(ax_c, t, Wtot_norm, color = :black, linewidth = 1.5, linestyle = :dash, label = "W_tot")
    axislegend(ax_c, position = :rt)
    
    # Panel (d): Magnetic Fluctuation Energy (normalizado por W0_mag)
    ax_d = Axis(fig[2, 2],
        xlabel = "Ωᵢt",
        ylabel = "W_f / W₀_mag",
        title = "(d) Magnetic Fluctuation Energy"#,
#        yscale = log10
    )
    lines!(ax_d, t, Wf_norm_mag, color = color_Wf, linewidth = 2, label = "W_f")
    # Añadir predicción teórica
#    hlines!(ax_d, [dB_B0_pred^2], color = :gray, linewidth = 1, linestyle = :dash, 
 #           label = "fV² (pred.)")
    axislegend(ax_d, position = :rb)
    
    # Título
    Label(fig[0, :], 
        text = "Energy Histories - Winske & Leroy (1984) Figure 8 (Corrected)",
        fontsize = 18, font = :bold
    )
    
    save(output_file, fig, px_per_unit = 2)
    println("Figura guardada en: $output_file")
    
    return fig
end

"""
Versión con escala logarítmica en todos los paneles
"""
function plot_figure8_log(params, data; output_file="figure8_winske_log_corrected.png")
    t = data["time"]
    wkpal_sp = data["wkpal_sp"]
    wkper_sp = data["wkper_sp"]
    enrgys = data["enrgys"]
    
    # Correcciones
    Wm_par = wkpal_sp[:, 1]
    Wm_perp = 2.0 .* wkper_sp[:, 1]
    Wb_par = wkpal_sp[:, 2]
    Wb_perp = 2.0 .* wkper_sp[:, 2]
    
    Wm = Wm_par .+ Wm_perp
    Wb = Wb_par .+ Wb_perp
    Wf = enrgys[:, 1]
    
    xlen = params["xlen"]
    W0_mag = xlen / 2.0
    W0_kin = Wm[1] + Wb[1]
    W0 = W0_kin
    
    Wm_par_norm = Wm_par ./ W0
    Wm_perp_norm = Wm_perp ./ W0
    Wb_par_norm = Wb_par ./ W0
    Wb_perp_norm = Wb_perp ./ W0
    Wm_norm = Wm ./ W0
    Wb_norm = Wb ./ W0
    Wf_norm_mag = Wf ./ W0_mag
    
    fig = Figure(size = (1000, 900), fontsize = 14)
    
    color_par = :blue
    color_perp = :red
    color_total_m = :green
    color_total_b = :orange
    color_Wf = :purple
    
    ax_a = Axis(fig[1, 1], xlabel = "Ωᵢt", ylabel = "W / W₀", 
                title = "(a) Beam Energies", yscale = log10)
    lines!(ax_a, t, Wb_par_norm, color = color_par, linewidth = 2, label = "W_b∥")
    lines!(ax_a, t, Wb_perp_norm, color = color_perp, linewidth = 2, label = "W_b⊥")
    axislegend(ax_a, position = :rt)
    
    ax_b = Axis(fig[1, 2], xlabel = "Ωᵢt", ylabel = "W / W₀", 
                title = "(b) Main Ion Energies", yscale = log10)
    lines!(ax_b, t, Wm_par_norm, color = color_par, linewidth = 2, label = "W_m∥")
    lines!(ax_b, t, Wm_perp_norm, color = color_perp, linewidth = 2, label = "W_m⊥")
    axislegend(ax_b, position = :rb)
    
    ax_c = Axis(fig[2, 1], xlabel = "Ωᵢt", ylabel = "W / W₀", 
                title = "(c) Total Species Energies", yscale = log10)
    lines!(ax_c, t, Wb_norm, color = color_total_b, linewidth = 2, label = "W_b")
    lines!(ax_c, t, Wm_norm, color = color_total_m, linewidth = 2, label = "W_m")
    axislegend(ax_c, position = :rt)
    
    ax_d = Axis(fig[2, 2], xlabel = "Ωᵢt", ylabel = "W_f / W₀_mag", 
                title = "(d) Magnetic Fluctuation Energy", yscale = log10)
    lines!(ax_d, t, Wf_norm_mag, color = color_Wf, linewidth = 2, label = "W_f")
    axislegend(ax_d, position = :rb)
    
    Label(fig[0, :], 
        text = "Energy Histories - Winske & Leroy (1984) Figure 8 (Log Scale)",
        fontsize = 18, font = :bold
    )
    
    save(output_file, fig, px_per_unit = 2)
    println("Figura guardada en: $output_file")
    
    return fig
end

#=
Solo las energías totales
=#

function total_energy_plots(params, data; output_file="figure8_winske_tot.png")
    t = data["time"]
    enrgys = data["enrgys"]

    wbtot = enrgys[:,1]
    wetot = enrgys[:,2]
    kpal = enrgys[:,3]
    kper = enrgys[:,4]
    wktot = enrgys[:,5]
    etot = enrgys[:,6]

    fig = Figure(size = (1000, 900), fontsize = 14)

    ax_a = Axis(fig[1, 1], xlabel = "Ωᵢt", ylabel = "E", 
                title = "Total system Energies")
    lines!(ax_a, t, enrgys[:,1], linewidth = 2, label = "Wb")
    lines!(ax_a, t, enrgys[:,2], linewidth = 2, label = "We")
    lines!(ax_a, t, enrgys[:,3], linewidth = 2, label = "K_∥")
    lines!(ax_a, t, enrgys[:,4], linewidth = 2, label = "K_⊥")
    lines!(ax_a, t, enrgys[:,5], linewidth = 2, label = "K_t")
    lines!(ax_a, t, enrgys[:,6], linewidth = 2, label = "E")
 #   lines!(ax_a, t, Wb_perp_norm, color = color_perp, linewidth = 2, label = "W_b⊥")
    axislegend(ax_a, position = :rt)

    save(output_file, fig, px_per_unit = 2)
    println("Figura guardada en: $output_file")

    return fig
end

# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================
function main()
    energy_file = "energy.d12"
    
    if !isfile(energy_file)
        println("ERROR: No se encontró el archivo '$energy_file'")
        return
    end
    
    params, data = read_energy_file(energy_file)
    
    plot_figure8_corrected(params, data, output_file = "figure8_winske_corrected.png")
    plot_figure8_log(params, data, output_file = "figure8_winske_log_corrected.png")
    total_energy_plots(params, data; output_file="figure8_winske_tot.png")
    
    println("\n¡Figuras generadas exitosamente!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 1
        energy_file = ARGS[1]
        if isfile(energy_file)
            params, data = read_energy_file(energy_file)
            plot_figure8_corrected(params, data)
            plot_figure8_log(params, data)
            total_energy_plots(params, data)
        else
            println("ERROR: No se encontró el archivo '$energy_file'")
        end
    else
        main()
    end
end



