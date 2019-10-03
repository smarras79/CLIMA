"""
    Diagnostics

Accumulate mean fields and covariance statistics on the computational grid

"""

# TODO: add functional ways to iterate over the mesh for these statistics
# TODO: add a way to reduce communication: compute averages locally only

module Diagnostics

using JLD2
using MPI
using StaticArrays

using ..Mesh.Topologies
using ..Mesh.Grids
using ..MPIStateArrays
using ..Atmos
using ..VariableTemplates
using ..PlanetParameters

export gather_diagnostics

function gather_diagnostics(mpicomm, dg, Q, current_time_string, κ, out_dir)
    mpirank = MPI.Comm_rank(mpicomm)
    nranks = MPI.Comm_size(mpicomm)

    bl = dg.balancelaw
    grid = dg.grid
    topology = grid.topology
    N = polynomialorder(grid)
    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    nrealelem = length(topology.realelems)
    nvertelem = topology.stacksize
    nhorzelem = div(nrealelem, nvertelem)

    if Array ∈ typeof(Q).parameters
        localQ = Q.realdata
        localaux = dg.auxstate.realdata
        localvgeo = grid.vgeo
    else
        localQ = Array(Q.realdata)
        localaux = Array(dg.auxstate.realdata)
        localvgeo = Array(grid.vgeo)
    end
    FT = eltype(localQ)

    nstate = num_state(bl, FT)
    nauxstate = num_aux(bl, FT)

    nstate = 6
    nthermo = 7

    thermoQ = zeros(Nq * Nq * Nqk, nthermo, nrealelem)
    Zvals = zeros(Nqk, nvertelem)

    for e in 1:nrealelem  
        for i in 1:Nq*Nq*Nqk
            z = localvgeo[i,grid.x3id,e]

            rho_node  = localQ[i,1,e]
            u_node    = localQ[i,2,e] / rho_node
            v_node    = localQ[i,3,e] / rho_node
            w_node    = localQ[i,4,e] / rho_node
            etot_node = localQ[i,5,e] / rho_node
            qt_node   = localQ[i,6,e] / rho_node

            e_int = etot_node - 1//2 * (u_node^2 + v_node^2 + w_node^2) - grav * z

            ts = PhaseEquil(e_int, qt_node, rho_node)
            Phpart = PhasePartition(ts)
            thermoQ[i,1,e] = Phpart.liq
            thermoQ[i,2,e] = Phpart.ice
            thermoQ[i,3,e] = qt_node-Phpart.liq-Phpart.ice
            thermoQ[i,4,e] = ts.T
            thermoQ[i,5,e] = liquid_ice_pottemp(ts)
            thermoQ[i,6,e] = dry_pottemp(ts)
            thermoQ[i,7,e] = virtual_pottemp(ts)
        end
    end

    fluctT = zeros(Nq * Nq * Nqk, nthermo, nrealelem)    
    fluctQ = zeros(Nq * Nq * Nqk, nstate, nrealelem)
    varQ   = zeros(Nq * Nq * Nqk, nstate, nrealelem)

    rho_localtot = sum(localQ[:,1,:]) / (size(localQ, 1) * size(localQ, 3))
    U_localtot   = sum(localQ[:,2,:])
    V_localtot   = sum(localQ[:,3,:])
    W_localtot   = sum(localQ[:,4,:]) / (size(localQ, 1) * size(localQ, 3))
    e_localtot   = sum(localQ[:,5,:])
    qt_localtot  = sum(localQ[:,6,:])

    qliq_localtot  = sum(thermoQ[:,1,:])
    qice_localtot  = sum(thermoQ[:,2,:])
    qvap_localtot  = sum(thermoQ[:,3,:])
    T_localtot     = sum(thermoQ[:,4,:])
    theta_localtot = sum(thermoQ[:,5,:])

    rho_tot   = MPI.Reduce(rho_localtot, +, 0, mpicomm)
    U_tot     = MPI.Reduce(U_localtot, +, 0, mpicomm)
    V_tot     = MPI.Reduce(V_localtot, +, 0, mpicomm)
    W_tot     = MPI.Reduce(W_localtot, +, 0, mpicomm)
    e_tot     = MPI.Reduce(e_localtot, +, 0, mpicomm)
    qt_tot    = MPI.Reduce(qt_localtot, +, 0, mpicomm)
    qliq_tot  = MPI.Reduce(qliq_localtot, +, 0, mpicomm)
    qice_tot  = MPI.Reduce(qice_localtot, +, 0, mpicomm)
    qvap_tot  = MPI.Reduce(qvap_localtot, +, 0, mpicomm)
    T_tot     = MPI.Reduce(T_localtot, +, 0, mpicomm)
    theta_tot = MPI.Reduce(theta_localtot, +, 0, mpicomm)

    if mpirank == 0
        # TODO: use the integral of the domain instead and delete this hack
        tot_num_elems = MPI.Allreduce(size(localQ, 3), +, mpicomm)
        tot_points = dofs_per_element(grid) * tot_num_elems

        rho_avg   = rho_tot / nranks
        U_avg     = (U_tot / tot_points) / rho_avg
        V_avg     = (V_tot / tot_points) / rho_avg
        W_avg     = (W_tot / (nranks * rho_avg)) 
        e_avg     = (e_tot / tot_points) / rho_avg
        qt_avg    = (qt_tot / tot_points) / rho_avg
        qliq_avg  = qliq_tot / tot_points
        qice_avg  = qice_tot / tot_points
        qvap_avg  = qvap_tot / tot_points
        T_avg     = T_tot / tot_points
        theta_avg = theta_tot / tot_points

        @debug "ρ average = $(rho_avg)"
        @debug "U average = $(U_avg)"
        @debug "V average = $(V_avg)"
        @debug "W average = $(W_avg)"
        @debug "e average = $(e_avg)"
        @debug "qt average = $(qt_avg)"
        @debug "qliq average = $(qliq_avg)"
        @debug "qice average = $(qice_avg)"
        @debug "qvap average = $(qvap_avg)"
        @debug "T average = $(T_avg)"
        @debug "theta average = $(theta_avg)"
    end

    AVG = SVector(rho_avg, U_avg, V_avg, W_avg, e_avg, qt_avg)
    AVG_T = SVector(qliq_avg, qice_avg, qvap_avg, T_avg, theta_avg)

    # Horizontal averages
    horz_avgs = zeros(Nqk, nvertelem, 10)
    LWP_local = 0
    for eh in 1:nhorzelem
        for ev in 1:nvertelem
            e = ev + (eh - 1) * nvertelem

            for k in 1:Nqk
                for j in 1:Nq
                    for i in 1:Nq
                        ijk = i + Nq * ((j-1) + Nq * (k-1)) 
                        horz_avgs[k,ev,1] += localQ[ijk,1,e] #density average
                        horz_avgs[k,ev,2] += localQ[ijk,2,e] #U average 
                        horz_avgs[k,ev,3] += localQ[ijk,3,e] #V average
                        horz_avgs[k,ev,4] += localQ[ijk,4,e] #W average
                        horz_avgs[k,ev,5] += thermoQ[ijk,5,e] # theta l  average
                        horz_avgs[k,ev,6] += thermoQ[ijk,1,e] #qliq average
                        horz_avgs[k,ev,7] += thermoQ[ijk,3,e] #qvap average 
                        horz_avgs[k,ev,8] += thermoQ[ijk,6,e] #dry theta average
                        horz_avgs[k,ev,9] += localQ[ijk,6,e] #qt average
                        horz_avgs[k,ev,10] += thermoQ[ijk,7,e] #virtual theta average

                        # liquid water path
                        if ev == floor(nvertelem/2) && k==floor(Nqk/2)
                            LWP_local += (localaux[ijk,1,e] + localaux[ijk,2,e])
                                         / κ / (Nq * Nq * nhorzelem)
                        end
                    end
                end
            end
        end
    end

    horz_avgs_tot = zeros(Nqk,nvertelem,10)
    for s in 1:10
        for ev in 1:nvertelem
            for k in 1:Nqk
                horz_avgs[k,ev,s] = horz_avgs[k,ev,s] /  (Nq * Nq * nhorzelem)
                horz_avgs_tot[k,ev,s] = MPI.Reduce(horz_avgs[k,ev,s], +, 0, mpicomm)
                if mpirank == 0
                    horz_avgs_tot[k,ev,s] = horz_avgs_tot[k,ev,s] / (nranks)
                end
            end
        end
    end

    if mpirank == 0
        LWP=MPI.Reduce(LWP_local, +, 0, mpicomm) / (nranks) 
    end

    S=zeros(Nqk,nvertelem,18)
    for eh in 1:nhorzelem
        for ev in 1:nvertelem
            e = ev + (eh - 1) * nvertelem

            for k in 1:Nqk
                for j in 1:Nq
                    for i in 1:Nq
                        ijk = i + Nq * ((j-1) + Nq * (k-1)) 
                        S[k,ev,1] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (thermoQ[ijk,6,e] - horz_avgs_tot[k,ev,8]) #fluctQ[ijk,4,e] * fluctT[ijk,5,e]
                        S[k,ev,2] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (thermoQ[ijk,3,e] - horz_avgs_tot[k,ev,7]) #fluctQ[ijk,4,e] * fluctT[ijk,3,e]
                        S[k,ev,3] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (localQ[ijk,2,e] / localQ[ijk,1,e] - horz_avgs_tot[k,ev,2] / horz_avgs_tot[k,ev,1])
                        #fluctQ[ijk,4,e] * fluctQ[ijk,2,e]
                        S[k,ev,4] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (localQ[ijk,3,e] / localQ[ijk,1,e] - horz_avgs_tot[k,ev,3] / horz_avgs_tot[k,ev,1])  
                        #fluctQ[ijk,4,e] * fluctQ[ijk,3,e]
                        S[k,ev,5] += (localQ[ijk,4,e]/localQ[ijk,1,e]-horz_avgs_tot[k,ev,4]/horz_avgs[k,ev,1])^2  #fluctQ[ijk,4,e] * fluctQ[ijk,4,e]
                        S[k,ev,6] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (localQ[ijk,1,e] - horz_avgs_tot[k,ev,1])
                        #fluctQ[ijk,4,e] * fluctQ[ijk,1,e]
                        S[k,ev,7] = horz_avgs[k,ev,6]
                        S[k,ev,8] += (localQ[ijk,4,e]/localQ[ijk,1,e]-horz_avgs_tot[k,ev,4]/horz_avgs_tot[k,ev,1]) * (thermoQ[ijk,1,e]-horz_avgs_tot[k,ev,6])  #fluctQ[ijk,4,e] * fluctT[ijk,1,e]
                        S[k,ev,9] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1])^3
                        #fluctQ[ijk,4,e] * fluctQ[ijk,4,e] * fluctQ[ijk,4,e]
                        S[k,ev,10] += (localQ[ijk,2,e]/localQ[ijk,1,e]-horz_avgs_tot[k,ev,2]/horz_avgs[k,ev,1])^2  #fluctQ[ijk,2,e] * fluctQ[ijk,2,e]
                        S[k,ev,11] += (localQ[ijk,3,e]/localQ[ijk,1,e]-horz_avgs_tot[k,ev,3]/horz_avgs[k,ev,1])^2  #fluctQ[ijk,3,e] * fluctQ[ijk,3,e]
                        Zvals[k,ev] = localvgeo[ijk,grid.x3id,e]
                        S[k,ev,12] = (horz_avgs_tot[k,ev,8])
                        S[k,ev,13] = (horz_avgs_tot[k,ev,9]/horz_avgs[k,ev,1])
                        S[k,ev,14] += S[k,ev,15]*(localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1])
                        S[k,ev,15] = (horz_avgs_tot[k,ev,5])
                        S[k,ev,16] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (thermoQ[ijk,7,e] - horz_avgs_tot[k,ev,10])
                        S[k,ev,17] += 0.5 * (S[k,ev,5] + S[k,ev,10] + S[k,ev,11])
                        S[k,ev,18] += (localQ[ijk,4,e]/localQ[ijk,1,e] - horz_avgs_tot[k,ev,4] / horz_avgs_tot[k,ev,1]) * (thermoQ[ijk,5,e] - horz_avgs_tot[k,ev,5])
                    end
                end
            end
        end
    end

    # See Outputs below for what S[k,ev,:] are respectively.

    S_avg = zeros(Nqk,nvertelem,18)
    for s in 1:18
        for ev in 1:nvertelem
            for k in 1:Nqk
                S[k,ev,s] = S[k,ev,s] /  (Nq * Nq * nhorzelem)
                S_avg[k,ev,s] = MPI.Reduce(S[k,ev,s], +, 0, mpicomm)

                if mpirank == 0
                    S_avg[k,ev,s] = S_avg[k,ev,s] / (nranks)
                end
            end
        end
    end

    OutputWtheta = zeros(nvertelem * Nqk)
    OutputWQVAP = zeros(nvertelem * Nqk)
    OutputWU = zeros(nvertelem * Nqk)
    OutputWV = zeros(nvertelem * Nqk)
    OutputWW = zeros(nvertelem * Nqk)
    OutputWRHO = zeros(nvertelem * Nqk)
    OutputQLIQ = zeros(nvertelem * Nqk)
    OutputWQLIQ = zeros(nvertelem * Nqk)
    OutputWWW = zeros(nvertelem * Nqk )
    OutputUU = zeros(nvertelem * Nqk )
    OutputVV = zeros(nvertelem * Nqk )
    OutputZ = zeros(nvertelem * Nqk)
    OutputU = zeros(nvertelem * Nqk)
    OutputV = zeros(nvertelem * Nqk)
    Outputtheta = zeros(nvertelem * Nqk)
    Outputqt = zeros(nvertelem * Nqk)
    OutputWqt = zeros(nvertelem * Nqk)
    Outputthetaliq = zeros(nvertelem * Nqk)
    OutputWthetav = zeros(nvertelem * Nqk)
    OutputTKE = zeros(nvertelem * Nqk)
    OutputWthetal = zeros(nvertelem * Nqk)

    for ev in 1:nvertelem
        for k in 1:Nqk
            i=k + Nqk * (ev - 1)
            OutputWtheta[i] = S_avg[k,ev,1] # <w'theta'>
            OutputWQVAP[i] = S_avg[k,ev,2] # <w'qvap'>
            OutputWU[i] = S_avg[k,ev,3] # <w'u'>
            OutputWV[i] = S_avg[k,ev,4] # <w'v'>
            OutputWW[i] = S_avg[k,ev,5] # <w'w'>
            OutputWRHO[i] = S_avg[k,ev,6] #<w'rho'>
            OutputQLIQ[i] = S_avg[k,ev,7] # qliq 
            OutputWQLIQ[i] = S_avg[k,ev,8] # <w'qliq'>
            OutputWWW[i] = S_avg[k,ev,9]  #<w'w'w'>
            OutputUU[i] = S_avg[k,ev,10] #<u'u'>
            OutputVV[i] = S_avg[k,ev,11] #<v'v'>
            OutputZ[i] = Zvals[k,ev] # Height
            OutputU[i] =  horz_avgs_tot[k,ev,2] / horz_avgs_tot[k,ev,1] # <u> 
            OutputV[i] =  horz_avgs_tot[k,ev,3] / horz_avgs_tot[k,ev,1] # <v>
            Outputtheta[i] = S_avg[k,ev,12] # <theta> 
            Outputqt[i] = S_avg[k,ev,13]  # <qt>
            OutputWqt[i] = S_avg[k,ev,14] #<w'qt'> 
            Outputthetaliq[i] = S_avg[k,ev,15] # <thetaliq>
            OutputWthetav[i] = S_avg[k,ev,16] # <w'thetav'>
            OutputTKE[i] = S_avg[k,ev,17] # <TKE>
            OutputWthetal[i] = S_avg[k,ev,18] # <w'thetal'>
        end
    end

if mpirank == 0

  io = open(diagnostics_fileout, "a")
     current_time_str = string(current_time_string, "\n")
     write(io, current_time_str)
     writedlm(io, [OutputWtheta OutputWQVAP OutputWU OutputWV OutputWW OutputWRHO OutputQLIQ OutputWQLIQ OutputWWW OutputUU OutputVV OutputU OutputV Outputtheta Outputqt OutputWqt Outputthetaliq OutputWthetav OutputTKE OutputWthetal OutputZ])
  close(io)
  io = open(LWP_fileout, "a")
     current_time_str = string(current_time_string, "\n")
     write(io, current_time_str)
     writedlm(io, LWP)
  close(io)
end
end # function gather_diagnostics

end # module Diagnostics
