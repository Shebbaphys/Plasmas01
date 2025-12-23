using Printf, LinearAlgebra, TOML
using Random, ProgressBars
using FFTW
using UnPack: @unpack

#*********************************
# MÓDULO DE DIMENSIONES (constantes globales)
#*********************************
const lpx = 9
const nx = 2^lpx
const nxp1 = nx + 1
const nxp2 = nx + 2
const nxh = div(nx, 2)

const nsp = 2
const npg1 = 512
const npg2 = 512
const npg = npg1 + npg2
const np = nx * npg

const ihparm = 18
const nheadr = ihparm + 8 * nsp

#*********************************
# MÓDULO DE DATOS COMPARTIDOS
#*********************************
mutable struct SimState
    # --- Parámetros de Configuración ---
    home::String
    xlen::Float64
    dx::Float64
    dt::Float64
    tmax::Float64
    betaen::Float64
    wpiwci::Float64
    game::Float64
    dus::Float64
    bnamp::Float64
    slx::Float64
    eta::Float64
    bx::Float64
    
    # --- Control ---
    itmax::Int
    lfld::Int
    ifield::Int
    iparticles::Int
    ienergy::Int
    ifilter::Int
    itrestart::Int
    it1::Int
    it2::Int
    it3::Int
    it4::Int
    
    # --- Arrays por Especie (tamaño nsp) ---
    bf0::Vector{Float64}
    nions::Vector{Int}
    rdn::Vector{Float64}
    betain::Vector{Float64}
    vdr::Vector{Float64}
    anis::Vector{Float64}
    qi::Vector{Float64}
    ai::Vector{Float64}
    qm::Vector{Float64}
    gg::Vector{Float64}
    vthpal::Vector{Float64}
    vthper::Vector{Float64}
    tpal::Vector{Float64}
    tper::Vector{Float64}
    umx::Vector{Float64}
    
    # --- Campos (tamaño nxp2 para ghosts) ---
    by::Vector{Float64}
    bz::Vector{Float64}
    ex::Vector{Float64}
    ey::Vector{Float64}
    ez::Vector{Float64}
    xi::Vector{Float64}
    
    # --- Momentos (Grilla con ghosts) ---
    denh::Vector{Float64}
    denp1::Vector{Float64}
    uxh::Vector{Float64}
    uyh::Vector{Float64}
    uzh::Vector{Float64}
    
    # --- Partículas (tamaño np) ---
    xp0::Vector{Float64}
    xp1::Vector{Float64}
    vxp::Vector{Float64}
    vyp::Vector{Float64}
    vzp::Vector{Float64}
    qmr::Vector{Float64}
    
    # Índices y Shapes de Partículas
    imx::Vector{Int}
    iox::Vector{Int}
    ipx::Vector{Int}
    smm::Vector{Float64}
    smo::Vector{Float64}
    smp::Vector{Float64}
    
    # --- Momentos Acumuladores (nx x nsp) ---
    dnsp0::Matrix{Float64}
    dnsp1::Matrix{Float64}
    dnsh::Matrix{Float64}
    uxsh::Matrix{Float64}
    uysh::Matrix{Float64}
    uzsh::Matrix{Float64}
    
    # --- Vectores k y atenuación ---
    ek::Vector{Float64}
    attenu::Vector{Float64}
    
    # --- Workspaces (Pre-asignados) ---
    run::Array{Float64, 3}
    fnp1::Matrix{Float64}
    un::Matrix{Float64}
    g1::Matrix{Float64}
    g2::Matrix{Float64}
    vpp::Matrix{Float64}
    
    grdnx::Vector{Float64}
    cury::Vector{Float64}
    curz::Vector{Float64}
    wrk1::Vector{Float64}
    wrk2::Vector{Float64}
    wrk3::Vector{Float64}
    
    # rrkfld workspaces
    gy1::Vector{Float64}
    gz1::Vector{Float64}
    gy2::Vector{Float64}
    gz2::Vector{Float64}
    byd::Vector{Float64}
    bzd::Vector{Float64}
    
    # FFT Workspaces (usar Any para compatibilidad)
    fft_buffer::Vector{ComplexF64}
    fft_work::Vector{ComplexF64}
    fft_plan::Any
    ifft_plan::Any

    # Arrays de trabajo para filter
    g::Vector{Float64}
    gy::Vector{Float64}
    gz::Vector{Float64}
    
    # Salidas
    enrgys::Vector{Float64}
    header::Vector{Float64}
    wkpal_sp::Vector{Float64}
    wkper_sp::Vector{Float64}
end


function SimState()
    # FFT - crear planes
    fft_buffer = zeros(ComplexF64, nx)
    fft_work = zeros(ComplexF64, nx)
    fft_plan = plan_fft(fft_buffer; flags=FFTW.MEASURE)
    ifft_plan = plan_ifft(fft_buffer; flags=FFTW.MEASURE)

    SimState(
        # Parámetros escalares (13 campos)
        "./",       # home
        256.0,      # xlen
        0.5,        # dx
        0.1,        # dt
        100.0,      # tmax
        1.0,        # betaen
        1e4,        # wpiwci
        1.0,        # game
        10.0,       # dus
        1e-6,       # bnamp
        256.0,      # slx
        0.0,        # eta
        1.0,        # bx
        
        # Control (11 campos)
        0,          # itmax
        1,          # lfld
        8,        # ifield
        8,       # iparticles
        10,         # ienergy
        2,          # ifilter
        0,          # itrestart
        1,          # it1
        1,          # it2
        1,          # it3
        0,          # it4
        
        # Arrays por especie (15 campos)
        zeros(3),           # bf0
        zeros(Int, nsp),    # nions
        zeros(nsp),         # rdn
        zeros(nsp),         # betain
        zeros(nsp),         # vdr
        zeros(nsp),         # anis
        zeros(nsp),         # qi
        zeros(nsp),         # ai
        zeros(nsp),         # qm
        zeros(nsp),         # gg
        zeros(nsp),         # vthpal
        zeros(nsp),         # vthper
        zeros(nsp),         # tpal
        zeros(nsp),         # tper
        zeros(nsp),         # umx
        
        # Campos (6 campos)
        zeros(nxp2),        # by
        zeros(nxp2),        # bz
        zeros(nxp2),        # ex
        zeros(nxp2),        # ey
        zeros(nxp2),        # ez
        zeros(nxp2),        # xi
        
        # Momentos grilla (5 campos)
        zeros(nxp2),        # denh
        zeros(nxp2),        # denp1
        zeros(nxp2),        # uxh
        zeros(nxp2),        # uyh
        zeros(nxp2),        # uzh
        
        # Partículas (6 campos)
        zeros(np),          # xp0
        zeros(np),          # xp1
        zeros(np),          # vxp
        zeros(np),          # vyp
        zeros(np),          # vzp
        zeros(np),          # qmr
        
        # Índices y shapes (6 campos)
        zeros(Int, np),     # imx
        zeros(Int, np),     # iox
        zeros(Int, np),     # ipx
        zeros(np),          # smm
        zeros(np),          # smo
        zeros(np),          # smp
        
        # Momentos por especie (6 campos)
        zeros(nx, nsp),     # dnsp0
        zeros(nx, nsp),     # dnsp1
        zeros(nx, nsp),     # dnsh
        zeros(nx, nsp),     # uxsh
        zeros(nx, nsp),     # uysh
        zeros(nx, nsp),     # uzsh
        
        # Vectores k (2 campos)
        zeros(nx),          # ek
        ones(nx),           # attenu
        
        # Workspaces pushv (12 campos)
        zeros(nx, 3, nsp),  # run
        zeros(nx, 3),       # fnp1
        zeros(nx, 3),       # un
        zeros(np, 3),       # g1
        zeros(np, 3),       # g2
        zeros(np, 3),       # vpp
        zeros(nx),          # grdnx
        zeros(nx),          # cury
        zeros(nx),          # curz
        zeros(nx),          # wrk1
        zeros(nx),          # wrk2
        zeros(nx),          # wrk3
        
        # rrkfld workspaces (6 campos)
        zeros(nx),          # gy1
        zeros(nx),          # gz1
        zeros(nx),          # gy2
        zeros(nx),          # gz2
        zeros(nx),          # byd
        zeros(nx),          # bzd
        
        # FFT (4 campos)
        fft_buffer,
        fft_work,
        fft_plan,
        ifft_plan,
        
        # filter workspaces (3 campos)
        zeros(nxp2),        # g
        zeros(nxp2),        # gy
        zeros(nxp2),        # gz
        
        # Salidas (2 campos)
        zeros(6),           # enrgys #original tenía 6
        zeros(nheadr),      # header
        zeros(nsp),         # wkpal_sp
        zeros(nsp),         # wkper_sp
    )
end


function generate_header!(s::SimState)
    @unpack dt, dx, itmax, lfld, betaen, wpiwci, bf0 = s
    @unpack ifield, iparticles, ienergy, ifilter, itrestart = s
    @unpack rdn, vdr, betain, anis, qi, ai, nions, gg, header = s
    
    header[1] = dt
    header[2] = dx
    header[3] = Float64(nx)
    header[4] = Float64(itmax)
    header[5] = Float64(lfld)
    header[6] = Float64(nsp)
    header[7] = Float64(np)
    header[8] = Float64(npg)
    header[9] = betaen
    header[10] = wpiwci
    header[11] = bf0[1]
    header[12] = bf0[2]
    header[13] = bf0[3]
    header[14] = Float64(ifield)
    header[15] = Float64(iparticles)
    header[16] = Float64(ienergy)
    header[17] = Float64(ifilter)
    header[18] = Float64(itrestart)

    @inbounds for is in 1:nsp
        indh = ihparm + 8 * (is - 1)
        header[indh+1] = rdn[is]
        header[indh+2] = vdr[is]
        header[indh+3] = betain[is]
        header[indh+4] = anis[is]
        header[indh+5] = qi[is]
        header[indh+6] = ai[is]
        header[indh+7] = Float64(nions[is])
        header[indh+8] = gg[is]
    end
end


#*******************************************
# INICIALIZACIÓN
#*******************************************
function init!(s::SimState)
    config_path = "config.toml"
    params = TOML.parsefile(config_path)
    
    sim_params = params["simulation"]
    ion_params = params["ions"]
    
    # Asignar directamente a la estructura
    s.home = sim_params["home"]
    s.xlen = sim_params["xlen"]
    s.dt = sim_params["dt"]
    s.tmax = sim_params["tmax"]
    s.betaen = sim_params["betaen"]
    s.wpiwci = sim_params["wpiwci"]
    s.game = sim_params["game"]
    s.dus = sim_params["dus"]
    s.bnamp = sim_params["bnamp"]
    s.bf0 .= sim_params["bf0"]
    s.itrestart = sim_params["itrestart"]
    s.ifield = sim_params["ifield"]
    s.iparticles = sim_params["iparticles"]
    s.ienergy = sim_params["ienergy"]
    s.ifilter = sim_params["ifilter"]
    s.eta = sim_params["eta"]
    
    if ion_params["nisp"] != nsp
        error("Discrepancia de nsp en config.toml")
    end
    
    # Manejo de valores escalares vs arrays
    qi_val = ion_params["qi"]
    ai_val = ion_params["ai"]
    betain_val = ion_params["betain"]
    vdr_val = ion_params["vdr"]
    rdn_val = ion_params["rdn"]
    anis_val = ion_params["anis"]
    
    if nsp == 1
        s.qi[1] = isa(qi_val, AbstractArray) ? qi_val[1] : qi_val
        s.ai[1] = isa(ai_val, AbstractArray) ? ai_val[1] : ai_val
        s.betain[1] = isa(betain_val, AbstractArray) ? betain_val[1] : betain_val
        s.vdr[1] = isa(vdr_val, AbstractArray) ? vdr_val[1] : vdr_val
        s.rdn[1] = isa(rdn_val, AbstractArray) ? rdn_val[1] : rdn_val
        s.anis[1] = isa(anis_val, AbstractArray) ? anis_val[1] : anis_val
    else
        s.qi .= qi_val
        s.ai .= ai_val
        s.betain .= betain_val
        s.vdr .= vdr_val
        s.rdn .= rdn_val
        s.anis .= anis_val
    end

    println("Parámetros cargados desde '$config_path'")
    
    s.lfld = 1
    s.dx = s.xlen / Float64(nx)
    s.slx = Float64(nx) * s.dx
    s.itmax = Int(floor(s.tmax / s.dt)) + 1

    println("*** tmax & itmax *** $(s.tmax)  $(s.itmax)")

    # Crear superpartículas
    betan = 0.0
    @inbounds for is in 1:nsp
        s.nions[is] = div(np, nsp)
        s.gg[is] = (s.rdn[is] / s.qi[is]) * Float64(nx) / s.nions[is]
        betan += s.betain[is]
    end
    betan += s.betaen

    # Razón carga-masa
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 = n2 + s.nions[is]
        s.qm[is] = s.qi[is] / s.ai[is]
        qm_is = s.qm[is]
        @simd for n in n1:n2
            s.qmr[n] = qm_is
        end
    end

    # Velocidades térmicas
    @inbounds for is in 1:nsp
        s.vthpal[is] = sqrt(s.betain[is] / s.ai[is])
        s.vthper[is] = s.vthpal[is] * sqrt(s.anis[is])
    end

    println("Inicializando posiciones...")
    println("Inicializando velocidades...")

    # Vector de ondas
    freqs = FFTW.fftfreq(nx, 1.0/s.dx)
    @inbounds for i in 1:nx
        s.ek[i] = 2.0 * π * freqs[i]
    end
    s.ek[nxh+1] = 0.0
    
    fill!(s.attenu, 1.0)
    
    @inbounds @simd for ix in 1:nxp2
        s.xi[ix] = Float64(ix - 2) * s.dx
    end
    
    s.bx = s.bf0[1]

    println("Inicialización completada.\n")
end


#*******************************************
# INITP
#*******************************************
function initp!(s::SimState)
    @unpack dx, nions, xp0, xp1, anis, ai, betain, vthpal, vthper, vdr = s
    @unpack dnsp1, vxp, vyp, vzp, game, gg, denp1 = s
    @unpack imx, iox, ipx, smm, smo, smp = s
    
    xlen_local = s.xlen
    xl = Float64(nx) * dx

    # Inicializa posiciones
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        for n in n1:n2
            xp0[n] = xl * rand()
            xp1[n] = xp0[n]
        end
    end

    # Shape functions
    shapefun!(s, xp0)

    # Inicializa velocidades
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        
        vthpal[is] = sqrt(betain[is] / ai[is])
        vthper[is] = vthpal[is] * sqrt(anis[is])
        
        vdrpal = vdr[is]
        
        for n in n1:n2
            vxp[n] = vthpal[is] * randn() + vdrpal
            vyp[n] = vthper[is] * randn()
            vzp[n] = vthper[is] * randn()
        end
    end

    # Densidad inicial
    fill!(dnsp1, 0.0)

    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        for n in n1:n2
            jm = mod1(imx[n] - 1, nx)
            jo = mod1(iox[n] - 1, nx)
            jp = mod1(ipx[n] - 1, nx)
            
            dnsp1[jm, is] += smm[n]
            dnsp1[jo, is] += smo[n]
            dnsp1[jp, is] += smp[n]
        end
    end

    @inbounds for is in 1:nsp
        @simd for ix in 1:nx
            dnsp1[ix, is] *= gg[is]
        end
    end

    gam = game - 1.0
    @inbounds for ix in 1:nx
        sumd = 0.0
        for is in 1:nsp
            sumd += dnsp1[ix, is] * ai[is]
        end
        denp1[ix+1] = sumd
    end
    
    denp1[1] = denp1[nxp1]
    denp1[nxp2] = denp1[2]
end


#*******************************************
# INITF
#*******************************************
function initf!(s::SimState)
    @unpack ex, ey, ez, bf0, by, bz, bnamp = s

    fill!(ex, 0.0)
    fill!(ey, 0.0)
    fill!(ez, 0.0)

    s.bx = bf0[1]
    @inbounds for ix in 2:nxp1
        by[ix] = bnamp * (rand() - 0.5) * 2.0 + bf0[2]
        bz[ix] = bnamp * (rand() - 0.5) * 2.0 + bf0[3]
    end

    sumy = sum(by[2:nxp1] .- bf0[2]) / Float64(nx)
    sumz = sum(bz[2:nxp1] .- bf0[3]) / Float64(nx)

    @inbounds @simd for ix in 2:nxp1
        by[ix] -= sumy
        bz[ix] -= sumz
    end

    # BC corregidas
    by[1] = by[nxp1]
    by[nxp2] = by[2]
    bz[1] = bz[nxp1]
    bz[nxp2] = bz[2]
end


#*******************************************
# SHAPEFUN
#*******************************************
function shapefun!(s::SimState, xp_array::Vector{Float64})
    @unpack dx, xlen, imx, iox, ipx, smm, smo, smp = s

    dxi = 1.0 / dx
    
    @inbounds for n in 1:np
        x_periodic = mod(xp_array[n], xlen)
       
        rpos = x_periodic * dxi
        ix_nearest = round(Int, rpos)
        dlx = rpos - Float64(ix_nearest)
        base_idx = ix_nearest + 2
  
        iox_n = base_idx
        imx_n = base_idx - 1
        ipx_n = base_idx + 1
        
        # Condiciones periódicas
        if iox_n > nxp1; iox_n -= nx; end
        if iox_n < 2; iox_n += nx; end
        if imx_n < 2; imx_n += nx; end
        if ipx_n > nxp1; ipx_n -= nx; end
        
        iox[n] = iox_n
        imx[n] = imx_n
        ipx[n] = ipx_n
        
        # Spline cuadrático
        dlx2 = dlx * dlx
        smm[n] = 0.5 * dlx2 - 0.5 * dlx + 0.125
        smo[n] = -dlx2 + 0.750
        smp[n] = 0.5 * dlx2 + 0.5 * dlx + 0.125
    end
end


#*******************************************
# PUSHX
#*******************************************
function pushx!(s::SimState)
    @unpack dt, xlen, xp0, xp1, vxp = s

    xmax = xlen
    
    @inbounds @simd for n in 1:np
        xnew = xp0[n] + dt * vxp[n]
        if xnew >= xmax
            xnew -= xmax
        elseif xnew < 0.0
            xnew += xmax
        end
        xp1[n] = xnew
    end
end


#*******************************************
# REPX
#*******************************************
function repx!(s::SimState)
    @unpack xp0, xp1 = s

    @inbounds @simd for n in 1:np
        xp0[n] = xp1[n]
    end
end


#*******************************************
# MOMENT
#*******************************************
function moment!(s::SimState, it::Int)
    @unpack nions, qi, gg, dx, xp1, ifilter = s
    @unpack dnsp0, dnsp1, dnsh, uxsh, uysh, uzsh = s
    @unpack imx, iox, ipx, smm, smo, smp, vxp, vyp, vzp = s
    @unpack denh, denp1, uxh, uyh, uzh = s

    # Guardar y reinicializar
    @inbounds for is in 1:nsp
        @simd for ix in 1:nx
            dnsp0[ix, is] = dnsp1[ix, is]
            dnsp1[ix, is] = 0.0
            uxsh[ix, is] = 0.0
            uysh[ix, is] = 0.0
            uzsh[ix, is] = 0.0
        end
    end
    
    # Primera acumulación
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        
        for n in n1:n2
            jm = mod1(imx[n] - 1, nx)
            jo = mod1(iox[n] - 1, nx)
            jp = mod1(ipx[n] - 1, nx)
            
            uxsh[jm, is] += smm[n] * vxp[n]
            uxsh[jo, is] += smo[n] * vxp[n]
            uxsh[jp, is] += smp[n] * vxp[n]
            
            uysh[jm, is] += smm[n] * vyp[n]
            uysh[jo, is] += smo[n] * vyp[n]
            uysh[jp, is] += smp[n] * vyp[n]
            
            uzsh[jm, is] += smm[n] * vzp[n]
            uzsh[jo, is] += smo[n] * vzp[n]
            uzsh[jp, is] += smp[n] * vzp[n]
        end
    end
    
    shapefun!(s, xp1)
    
    # Segunda acumulación
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        
        for n in n1:n2
            jm = mod1(imx[n] - 1, nx)
            jo = mod1(iox[n] - 1, nx)
            jp = mod1(ipx[n] - 1, nx)
            
            dnsp1[jm, is] += smm[n]
            dnsp1[jo, is] += smo[n]
            dnsp1[jp, is] += smp[n]
            
            uxsh[jm, is] += smm[n] * vxp[n]
            uxsh[jo, is] += smo[n] * vxp[n]
            uxsh[jp, is] += smp[n] * vxp[n]
            
            uysh[jm, is] += smm[n] * vyp[n]
            uysh[jo, is] += smo[n] * vyp[n]
            uysh[jp, is] += smp[n] * vyp[n]
            
            uzsh[jm, is] += smm[n] * vzp[n]
            uzsh[jo, is] += smo[n] * vzp[n]
            uzsh[jp, is] += smp[n] * vzp[n]
        end
    end
    
    # Normalizar
    @inbounds for is in 1:nsp
        factor = gg[is]
        factor_v = factor * 0.5
        @simd for ix in 1:nx
            dnsp1[ix, is] *= factor
            uxsh[ix, is] *= factor_v
            uysh[ix, is] *= factor_v
            uzsh[ix, is] *= factor_v
            dnsh[ix, is] = (dnsp1[ix, is] + dnsp0[ix, is]) * 0.5
        end
    end
    
    # Momentos totales
    @inbounds for ix in 1:nx
        sumdh = 0.0
        sumd1 = 0.0
        sumvx = 0.0
        sumvy = 0.0
        sumvz = 0.0
        
        for is in 1:nsp
            q = qi[is]
            sumdh += dnsh[ix, is] * q
            sumd1 += dnsp1[ix, is] * q
            sumvx += uxsh[ix, is] * q
            sumvy += uysh[ix, is] * q
            sumvz += uzsh[ix, is] * q
        end
        
        denh[ix+1] = sumdh
        denp1[ix+1] = sumd1
        
        if sumdh > 0.0
            inv_den = 1.0 / sumdh
            uxh[ix+1] = sumvx * inv_den
            uyh[ix+1] = sumvy * inv_den
            uzh[ix+1] = sumvz * inv_den
        else
            uxh[ix+1] = 0.0
            uyh[ix+1] = 0.0
            uzh[ix+1] = 0.0
        end
    end

    # Filtros
    filter!(s, denh)
    filter!(s, denp1)
    filter!(s, uxh)
    filter!(s, uyh)
    filter!(s, uzh)
    
    # Condiciones de borde
    denh[1] = denh[nxp1]
    denh[nxp2] = denh[2]
    denp1[1] = denp1[nxp1]
    denp1[nxp2] = denp1[2]
    uxh[1] = uxh[nxp1]
    uxh[nxp2] = uxh[2]
    uyh[1] = uyh[nxp1]
    uyh[nxp2] = uyh[2]
    uzh[1] = uzh[nxp1]
    uzh[nxp2] = uzh[2]
end


#*******************************************
# CURL2 y GRAD2 con FFT
#*******************************************
@inline function curl2!(s::SimState, wrk2::Vector{Float64}, wrk3::Vector{Float64}, 
                        cury::Vector{Float64}, curz::Vector{Float64})
    @unpack fft_buffer, fft_work, fft_plan, ifft_plan, ek = s

    # cury = ∂Bz/∂x
    @inbounds @simd for i in 1:nx
        fft_buffer[i] = ComplexF64(wrk3[i], 0.0)
    end
    mul!(fft_work, fft_plan, fft_buffer)
    @inbounds @simd for i in 1:nx
        fft_work[i] *= im * ek[i]
    end
    mul!(fft_buffer, ifft_plan, fft_work)
    @inbounds @simd for i in 1:nx
        cury[i] = -real(fft_buffer[i])
    end
    
    # curz = -∂By/∂x
    @inbounds @simd for i in 1:nx
        fft_buffer[i] = ComplexF64(wrk2[i], 0.0)
    end
    mul!(fft_work, fft_plan, fft_buffer)
    @inbounds @simd for i in 1:nx
        fft_work[i] *= im * ek[i]
    end
    mul!(fft_buffer, ifft_plan, fft_work)
    @inbounds @simd for i in 1:nx
        curz[i] = real(fft_buffer[i])
    end
end

@inline function grad2!(s::SimState, x::Vector{Float64}, rx::Vector{Float64})
    @unpack fft_buffer, fft_work, fft_plan, ifft_plan, ek = s

    @inbounds @simd for i in 1:nx
        fft_buffer[i] = ComplexF64(x[i], 0.0)
    end
    mul!(fft_work, fft_plan, fft_buffer)
    @inbounds @simd for i in 1:nx
        fft_work[i] *= im * ek[i]
    end
    mul!(fft_buffer, ifft_plan, fft_work)
    @inbounds @simd for i in 1:nx
        rx[i] = real(fft_buffer[i])
    end
end


#*******************************************
# PUSHV
#*******************************************
function pushv!(s::SimState)
    @unpack dt, nions, qi, qmr, eta, betaen, bx, gg = s
    @unpack denp1, by, bz, ex, ey, ez = s
    @unpack vxp, vyp, vzp, uxh, uyh, uzh = s
    @unpack imx, iox, ipx, smm, smo, smp = s
    @unpack run, fnp1, un, g1, g2, vpp = s
    @unpack wrk1, wrk2, wrk3, grdnx, cury, curz = s

    fill!(run, 0.0)

    # Gradiente de densidad
    @inbounds @simd for ix in 1:nx
        wrk1[ix] = denp1[ix+1]
    end
    grad2!(s, wrk1, grdnx)
    
    # Curl de B
    @inbounds @simd for ix in 1:nx
        wrk2[ix] = by[ix+1]
        wrk3[ix] = bz[ix+1]
    end
    curl2!(s, wrk2, wrk3, cury, curz)
    
    # Términos conocidos fnp1
    @inbounds for ix in 1:nx
        dval = denp1[ix+1]
        gdnnp1 = (dval > 0.0) ? 1.0 / dval : 0.0
        
        fnp1[ix, 1] = (-betaen * grdnx[ix] * 0.5 + 
                        cury[ix] * bz[ix+1] - curz[ix] * by[ix+1]) * gdnnp1
        fnp1[ix, 2] = (curz[ix] * bx) * gdnnp1 + eta * cury[ix]
        fnp1[ix, 3] = (-cury[ix] * bx) * gdnnp1 + eta * curz[ix]
    end
    
    # RRK Paso 1: Predictor
    @inbounds for n in 1:np
        jm = mod1(imx[n] - 1, nx)
        jo = mod1(iox[n] - 1, nx)
        jp = mod1(ipx[n] - 1, nx)
        
        smm_n = smm[n]
        smo_n = smo[n]
        smp_n = smp[n]
        qmr_n = qmr[n]
        
        for k in 1:3
            g2[n, k] = (fnp1[jm, k] * smm_n + fnp1[jo, k] * smo_n + fnp1[jp, k] * smp_n) * qmr_n
        end
    end
    
    # g1 con fuerza de Lorentz
    @inbounds for n in 1:np
        imx_n = imx[n]
        iox_n = iox[n]
        ipx_n = ipx[n]
        
        smm_n = smm[n]
        smo_n = smo[n]
        smp_n = smp[n]
        qmr_n = qmr[n]
        
        vxp_n = vxp[n]
        vyp_n = vyp[n]
        vzp_n = vzp[n]
        
        # Componente X
        term_im = (vyp_n - uyh[imx_n]) * bz[imx_n] - (vzp_n - uzh[imx_n]) * by[imx_n]
        term_io = (vyp_n - uyh[iox_n]) * bz[iox_n] - (vzp_n - uzh[iox_n]) * by[iox_n]
        term_ip = (vyp_n - uyh[ipx_n]) * bz[ipx_n] - (vzp_n - uzh[ipx_n]) * by[ipx_n]
        g1[n, 1] = g2[n, 1] + qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
        
        # Componente Y
        term_im = (vzp_n - uzh[imx_n]) * bx - (vxp_n - uxh[imx_n]) * bz[imx_n]
        term_io = (vzp_n - uzh[iox_n]) * bx - (vxp_n - uxh[iox_n]) * bz[iox_n]
        term_ip = (vzp_n - uzh[ipx_n]) * bx - (vxp_n - uxh[ipx_n]) * bz[ipx_n]
        g1[n, 2] = g2[n, 2] + qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
        
        # Componente Z
        term_im = (vxp_n - uxh[imx_n]) * by[imx_n] - (vyp_n - uyh[imx_n]) * bx
        term_io = (vxp_n - uxh[iox_n]) * by[iox_n] - (vyp_n - uyh[iox_n]) * bx
        term_ip = (vxp_n - uxh[ipx_n]) * by[ipx_n] - (vyp_n - uyh[ipx_n]) * bx
        g1[n, 3] = g2[n, 3] + qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
    end
    
    # Velocidad predicha
    dt_half = 0.5 * dt
    @inbounds for n in 1:np
        vpp[n, 1] = vxp[n] + dt_half * g1[n, 1]
        vpp[n, 2] = vyp[n] + dt_half * g1[n, 2]
        vpp[n, 3] = vzp[n] + dt_half * g1[n, 3]
    end
    
    # Acumular corriente
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 += nions[is]
        q = qi[is]
        
        for n in n1:n2
            jm = mod1(imx[n] - 1, nx)
            jo = mod1(iox[n] - 1, nx)
            jp = mod1(ipx[n] - 1, nx)
            
            smm_n = smm[n]
            smo_n = smo[n]
            smp_n = smp[n]
            
            for k in 1:3
                val = q * vpp[n, k]
                run[jm, k, is] += val * smm_n
                run[jo, k, is] += val * smo_n
                run[jp, k, is] += val * smp_n
            end
        end
    end
    
    # Velocidad promedio
    @inbounds for ix in 1:nx
        dval = denp1[ix+1]
        gdnnp1 = (dval > 0.0) ? 1.0 / dval : 0.0
        
        rvsumx = 0.0
        rvsumy = 0.0
        rvsumz = 0.0
        
        for is in 1:nsp
            gg_is = gg[is]
            rvsumx += run[ix, 1, is] * gg_is
            rvsumy += run[ix, 2, is] * gg_is
            rvsumz += run[ix, 3, is] * gg_is
        end
        
        un[ix, 1] = rvsumx * gdnnp1
        un[ix, 2] = rvsumy * gdnnp1
        un[ix, 3] = rvsumz * gdnnp1
    end
    
    # Campo eléctrico
    @inbounds for ix in 1:nx
        ex[ix+1] = -un[ix, 2] * bz[ix+1] + un[ix, 3] * by[ix+1] + fnp1[ix, 1]
        ey[ix+1] = -un[ix, 3] * bx + un[ix, 1] * bz[ix+1] + fnp1[ix, 2]
        ez[ix+1] = -un[ix, 1] * by[ix+1] + un[ix, 2] * bx + fnp1[ix, 3]
    end
    
    # RRK Paso 2: Corrector
    @inbounds for n in 1:np
        jm = mod1(imx[n] - 1, nx)
        jo = mod1(iox[n] - 1, nx)
        jp = mod1(ipx[n] - 1, nx)
        
        imx_n = imx[n]
        iox_n = iox[n]
        ipx_n = ipx[n]
        
        smm_n = smm[n]
        smo_n = smo[n]
        smp_n = smp[n]
        qmr_n = qmr[n]
        
        vpp1 = vpp[n, 1]
        vpp2 = vpp[n, 2]
        vpp3 = vpp[n, 3]
        
        # Componente X
        term_im = (vpp2 - un[jm, 2]) * bz[imx_n] - (vpp3 - un[jm, 3]) * by[imx_n]
        term_io = (vpp2 - un[jo, 2]) * bz[iox_n] - (vpp3 - un[jo, 3]) * by[iox_n]
        term_ip = (vpp2 - un[jp, 2]) * bz[ipx_n] - (vpp3 - un[jp, 3]) * by[ipx_n]
        g2[n, 1] += qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
        
        # Componente Y
        term_im = (vpp3 - un[jm, 3]) * bx - (vpp1 - un[jm, 1]) * bz[imx_n]
        term_io = (vpp3 - un[jo, 3]) * bx - (vpp1 - un[jo, 1]) * bz[iox_n]
        term_ip = (vpp3 - un[jp, 3]) * bx - (vpp1 - un[jp, 1]) * bz[ipx_n]
        g2[n, 2] += qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
        
        # Componente Z
        term_im = (vpp1 - un[jm, 1]) * by[imx_n] - (vpp2 - un[jm, 2]) * bx
        term_io = (vpp1 - un[jo, 1]) * by[iox_n] - (vpp2 - un[jo, 2]) * bx
        term_ip = (vpp1 - un[jp, 1]) * by[ipx_n] - (vpp2 - un[jp, 2]) * bx
        g2[n, 3] += qmr_n * (term_im * smm_n + term_io * smo_n + term_ip * smp_n)
    end
    
    # Productos escalares
    g11 = 0.0
    g22 = 0.0
    g12 = 0.0
    
    @inbounds for k in 1:3
        for n in 1:np
            val_g1 = g1[n, k]
            val_g2 = 2.0 * val_g1 - g2[n, k]
            g2[n, k] = val_g2
            
            g11 += val_g1 * val_g1
            g22 += val_g2 * val_g2
            g12 += val_g1 * val_g2
        end
    end
    
    # Actualizar velocidades
    if g22 > 0.0
        factor = dt / g22
        @inbounds for n in 1:np
            vxp[n] += factor * (2.0 * g12 * g1[n, 1] - g11 * g2[n, 1])
            vyp[n] += factor * (2.0 * g12 * g1[n, 2] - g11 * g2[n, 2])
            vzp[n] += factor * (2.0 * g12 * g1[n, 3] - g11 * g2[n, 3])
        end
    end
end


#*******************************************
# RRKSUB y RRKFLD
#*******************************************
function rrksub!(s::SimState, bxd::Float64, byd::Vector{Float64}, 
                 bzd::Vector{Float64}, cury_out::Vector{Float64}, curz_out::Vector{Float64})
    @unpack wrk2, wrk3, eta, denh, uxh, uyh, uzh = s
    
    # Calcular curl(B)
    @inbounds @simd for ix in 1:nx
        wrk2[ix] = byd[ix]
        wrk3[ix] = bzd[ix]
    end
    curl2!(s, wrk2, wrk3, cury_out, curz_out)
    
    # Calcular -E
    @inbounds for ix in 1:nx
        if denh[ix+1] > 0.0
            inv_denh = 1.0 / denh[ix+1]
            wrk2[ix] = uzh[ix+1] * bxd - uxh[ix+1] * bzd[ix] -
                       curz_out[ix] * bxd * inv_denh - eta * cury_out[ix]
            wrk3[ix] = uxh[ix+1] * byd[ix] - uyh[ix+1] * bxd +
                       cury_out[ix] * bxd * inv_denh - eta * curz_out[ix]
        else
            wrk2[ix] = 0.0
            wrk3[ix] = 0.0
        end
    end
    
    # Calcular -curl(E) = ∂B/∂t
    curl2!(s, wrk2, wrk3, cury_out, curz_out)
end


function rrkfld!(s::SimState, it::Int)
    @unpack lfld, dt, bx, by, bz = s
    @unpack gy1, gz1, gy2, gz2, byd, bzd = s
    
    if lfld <= 0
        return
    end
    
    ddt = dt / Float64(lfld)
    bxd = bx
    
    for l in 1:lfld
        @inbounds @simd for ix in 1:nx
            byd[ix] = by[ix+1]
            bzd[ix] = bz[ix+1]
        end
        
        rrksub!(s, bxd, byd, bzd, gy1, gz1)
        
        @inbounds @simd for ix in 1:nx
            byd[ix] = by[ix+1] + ddt * 0.5 * gy1[ix]
            bzd[ix] = bz[ix+1] + ddt * 0.5 * gz1[ix]
        end
        
        rrksub!(s, bxd, byd, bzd, gy2, gz2)
        
        g11 = 0.0
        g12 = 0.0
        g22 = 0.0
        
        @inbounds for ix in 1:nx
            gy1_ix = gy1[ix]
            gz1_ix = gz1[ix]
            gy2_ix = 2.0 * gy1_ix - gy2[ix]
            gz2_ix = 2.0 * gz1_ix - gz2[ix]
            
            gy2[ix] = gy2_ix
            gz2[ix] = gz2_ix
            
            g11 += gy1_ix * gy1_ix + gz1_ix * gz1_ix
            g22 += gy2_ix * gy2_ix + gz2_ix * gz2_ix
            g12 += gy1_ix * gy2_ix + gz1_ix * gz2_ix
        end
        
        if g22 > 0.0
            gg22 = ddt / g22
            @inbounds @simd for ix in 1:nx
                by[ix+1] += gg22 * (2.0 * g12 * gy1[ix] - g11 * gy2[ix])
                bz[ix+1] += gg22 * (2.0 * g12 * gz1[ix] - g11 * gz2[ix])
            end
        end
    end
    
    by[1] = by[nxp1]
    bz[1] = bz[nxp1]
    by[nxp2] = by[2]
    bz[nxp2] = bz[2]
end


#*******************************************
# FILTROS
#*******************************************
function filter!(s::SimState, a::Vector{Float64})
    @unpack g = s
    
    @inbounds for i in 2:nxp1
        g[i] = (2.0 * a[i] + a[i+1] + a[i-1]) * 0.25
    end
    
    g[1] = g[nxp1]
    g[nxp2] = g[2]
    
    @inbounds for i in 2:nxp1
        a[i] = (6.0 * g[i] - g[i+1] - g[i-1]) * 0.25
    end
    
    a[1] = a[nxp1]
    a[nxp2] = a[2]
end


function filter2!(s::SimState)
    @unpack gy, gz, by, bz = s
    
    @inbounds for i in 2:nxp1
        gy[i] = (2.0 * by[i] + by[i+1] + by[i-1]) * 0.25
        gz[i] = (2.0 * bz[i] + bz[i+1] + bz[i-1]) * 0.25
    end
    
    gy[1] = gy[nxp1]
    gy[nxp2] = gy[2]
    gz[1] = gz[nxp1]
    gz[nxp2] = gz[2]
    
    @inbounds for i in 2:nxp1
        by[i] = (6.0 * gy[i] - gy[i+1] - gy[i-1]) * 0.25
        bz[i] = (6.0 * gz[i] - gz[i+1] - gz[i-1]) * 0.25
    end
    
    by[1] = by[nxp1]
    by[nxp2] = by[2]
    bz[1] = bz[nxp1]
    bz[nxp2] = bz[2]
end


#*******************************************
# ENERGÍA
#*******************************************
function energy!(s::SimState, it::Int, t::Float64)
    @unpack nions, gg, ai, bf0, ienergy = s
    @unpack vxp, vyp, vzp = s
    @unpack tpal, tper, umx, enrgys = s
    @unpack by, bz, ex, ey, ez = s
    @unpack wkpal_sp, wkper_sp = s
    
    wksum = 0.0
    wksum1 = 0.0
    wksum2 = 0.0
    
    n2 = 0
    @inbounds for is in 1:nsp
        n1 = n2 + 1
        n2 = n2 + nions[is]
        
        fac1 = 1.0 / Float64(nions[is])
        fac2 = 0.5 * ai[is] * gg[is]
        
        umx_is = 0.0
        umy_is = 0.0
        umz_is = 0.0
        
        @simd for n in n1:n2
            umx_is += vxp[n]
            umy_is += vyp[n]
            umz_is += vzp[n]
        end
        
        umx_is *= fac1
        umy_is *= fac1
        umz_is *= fac1
        umx[is] = umx_is
        
        wkpal_is = 0.0
        wkper_is = 0.0
        tpal_is = 0.0
        tpery_is = 0.0
        tperz_is = 0.0
#        vkenrgy_is = 0.0
        
        @simd for n in n1:n2
            vx = vxp[n]
            vy = vyp[n]
            vz = vzp[n]
            
            wkpal_is += vx * vx
            wkper_is += vy * vy + vz * vz
            
            dvx = vx - umx_is
            dvy = vy - umy_is
            dvz = vz - umz_is
            
            tpal_is += dvx * dvx
            tpery_is += dvy * dvy
            tperz_is += dvz * dvz
        end
        
        wkpal_is *= fac2
        wkper_is *= fac2
        tpal_is *= fac2
        tpery_is *= fac2
        tperz_is *= fac2
        
        tpal[is] = tpal_is
        tper[is] = 0.5 * (tpery_is + tperz_is)
        wkpal_sp[is] = wkpal_is
        wkper_sp[is] = 0.5 * wkper_is   

        wksum1 += wkpal_is
        wksum2 += 0.5 * wkper_is
        wksum += wkpal_is + wkper_is
    end
    
    # Energía de campos
    bfldy = 0.0
    bfldz = 0.0
    efldx = 0.0
    efldy = 0.0
    efldz = 0.0
    
    bf0_2 = bf0[2]
    bf0_3 = bf0[3]
    
    @inbounds @simd for ix in 2:nxp1
        db_y = by[ix] - bf0_2
        db_z = bz[ix] - bf0_3
        bfldy += db_y * db_y
        bfldz += db_z * db_z
        efldx += ex[ix] * ex[ix]
        efldy += ey[ix] * ey[ix]
        efldz += ez[ix] * ez[ix]
    end
    
    bfldenrgy = 0.5 * (bfldy + bfldz)
    efldenrgy = 0.5 * (efldx + efldy + efldz)
    
    totenrgy = bfldenrgy + efldenrgy + wksum
    
    enrgys[1] = bfldenrgy
    enrgys[2] = efldenrgy
#    enrgys[3] = wkper_sp
#    enrgys[4] = wkpal_sp
    enrgys[3] = wksum1
    enrgys[4] = wksum2
    enrgys[5] = wksum
    enrgys[6] = totenrgy
    
    if mod(it, ienergy) == 0
        println(@sprintf("t=%.2e: Emag=%.6e Eelec=%.6e Ekin_par=%.6e Ekin_perp=%.6e Etot=%.6e", 
                        t, bfldenrgy, efldenrgy, wksum1, wksum2, totenrgy))
        @inbounds for is in 1:nsp
            println(@sprintf("  Especie %d: T_par=%.6e T_perp=%.6e <vx>=%.6e", 
                           is, tpal[is], tper[is], umx[is]))
        end
    end
end


#*******************************************
# SALIDA
#*******************************************
function outdat!(s::SimState, it::Int, t::Float64, 
                 fields_file::IO, particles_file::IO, energy_file::IO)
    @unpack it1, it2, it3, ifield, iparticles, ienergy = s
    @unpack by, bz, ex, ey, ez, uxsh, uysh, uzsh, dnsh = s
    @unpack tpal, tper, xp1, vxp, vyp, vzp, umx, enrgys = s
    @unpack home, header, tmax = s
    @unpack wkpal_sp, wkper_sp = s

    tim = Float64[Float64(it), Float64(t)]

    # GUARDAR CAMPOS
    if ifield > 0 && mod(it, ifield) == it1
        println(" *** Writing Field Data *** step = $it")
        
        write(fields_file, tim)
        write(fields_file, by[2:nxp1])
        write(fields_file, bz[2:nxp1])
        write(fields_file, ex[1:nxp2])
        write(fields_file, ey[1:nxp2])
        write(fields_file, ez[1:nxp2])
        write(fields_file, uxsh)
        write(fields_file, uysh)
        write(fields_file, uzsh)
        write(fields_file, dnsh)
        write(fields_file, tpal)
        write(fields_file, tper)
        
        flush(fields_file)
    end

    # GUARDAR PARTÍCULAS
    if iparticles > 0 && mod(it, iparticles) == it2
        println(" *** Writing Particle Data *** step = $it")
        
        write(particles_file, tim)
        write(particles_file, xp1)
        write(particles_file, vxp)
        write(particles_file, vyp)
        write(particles_file, vzp)
        
        flush(particles_file)
    end
    
    # GUARDAR ENERGÍAS
    if ienergy > 0 && mod(it, ienergy) == it3
        write(energy_file, tim)
        write(energy_file, tpal)
        write(energy_file, tper)
        write(energy_file, umx)
        write(energy_file, enrgys)
        write(energy_file, wkpal_sp)
        write(energy_file, wkper_sp)        
        flush(energy_file)
    end
    
    # ARCHIVO DE REINICIO
#    if t >= tmax
#        println(" *** Writing Restart File ***")
#        restart_path = joinpath(home, "restart.d13")
#        
#        open(restart_path, "w") do restart_file
#            write(restart_file, header)
#            write(restart_file, tim)
#            write(restart_file, by[2:nxp1])
#            write(restart_file, bz[2:nxp1])
#            write(restart_file, ex)
#            write(restart_file, ey)
#            write(restart_file, ez)
#            write(restart_file, uxsh)
#            write(restart_file, uysh)
#            write(restart_file, uzsh)
#            write(restart_file, dnsh)
#            write(restart_file, tpal)
#            write(restart_file, tper)
#            write(restart_file, xp1)
#            write(restart_file, vxp)
#            write(restart_file, vyp)
#            write(restart_file, vzp)
#            flush(restart_file)
#       end
#        println(" *** Restart file saved successfully ***")
 #   end

    if mod(it, ifield) == 0 && it > 0
        println(@sprintf("Output at step %d, t=%.3e", it, t))
    end
end


#function restart!(t_ref)
#    println("Restart not implemented.")
#end


#=============================================================================
# PROGRAMA PRINCIPAL
=============================================================================#
function hyb1d_rrkfft()
    t_begin = Base.time()
    it = 0
    t = 0.0
    
    println("╔════════════════════════════════════════════════════════════════╗")
    println("║     Hybrid PIC 1D - RRK Method (Rational Runge-Kutta)          ║")
    println("║     Created by A.F. Viñas (NASA-GSFC, 1993)                    ║")
    println("║     Julia version by Sebastián Pons (UdeC, 2025)               ║")
    println("╚════════════════════════════════════════════════════════════════╝")
    println()
    
    s = SimState()

    init!(s)
    initp!(s)
    initf!(s)


    # Abrir archivos binarios
    fields_file = open(joinpath(s.home, "fields.d10"), "w")
    particles_file = open(joinpath(s.home, "partcls.d11"), "w")
    energy_file = open(joinpath(s.home, "energy.d12"), "w")

    # Escribir headers
    generate_header!(s)
    write(fields_file, s.header)
    write(particles_file, s.header)
    write(energy_file, s.header)
    
    println("System dimensions:")
    println("  nx = $nx, np = $np, nsp = $nsp")
    println("  itmax = $(s.itmax), dt = $(s.dt)\n")

    s.it1 = (s.ifield == 1) ? 0 : 1
    s.it2 = (s.iparticles == 1) ? 0 : 1
    s.it3 = (s.ienergy == 1) ? 0 : 1
    s.it4 = 0
    
    energy!(s, it, t)
    outdat!(s, it, t, fields_file, particles_file, energy_file)
    
    pushx!(s)
    moment!(s, it)
    rrkfld!(s, it)
    
    if s.ifilter > 0 && mod(it, s.ifilter) == s.it4
        filter2!(s)
    end
    
    println("Starting main loop...")
    println("─" ^ 64)
    
    # Cachear valores locales para el loop
    dt_local = s.dt
    itmax_local = s.itmax
    
    for it in ProgressBar(1:itmax_local)
        t = it * dt_local
        
        if mod(it, s.ienergy) == s.it3
            energy!(s, it, t)
        end
        
        if mod(it, s.ifield) == s.it1
            outdat!(s, it, t, fields_file, particles_file, energy_file)
        end
        
        pushv!(s)
        repx!(s)
        pushx!(s)
        moment!(s, it)
        rrkfld!(s, it)
        
        if s.ifilter > 0 && mod(it, s.ifilter) == s.it4
            filter2!(s)
        end
        
        if mod(it, 10) == 0
            progress = it / itmax_local * 100
            println(@sprintf("  Step %6d / %6d  (%.1f%%)  t = %.3e", 
                           it, itmax_local, progress, t))
        end
    end
    
    it = itmax_local + 1
    t = it * dt_local

    if mod(it, s.ienergy) == s.it3
        energy!(s, it, t)
    end
    
    if mod(it, s.ifield) == s.it1
        outdat!(s, it, t, fields_file, particles_file, energy_file)
    end
    
    close(fields_file)
    close(particles_file)
    close(energy_file)

    println("─" ^ 64)
    println(" *** Simulation Completed ***")
    
    t_end = Base.time()
    tcpu = t_end - t_begin
    println(@sprintf(" CPU Time: %.2f seconds", tcpu))
    println()
end

#*********************************
# PUNTO DE ENTRADA
#*********************************
if abspath(PROGRAM_FILE) == @__FILE__
    hyb1d_rrkfft()
end
