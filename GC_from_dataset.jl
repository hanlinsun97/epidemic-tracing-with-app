
#####################  GC_from_dataset  #############################################
# 
# This code can be redistributed and/or modified
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#  
# This program is distributed ny the authors in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  
# If you use this code please cite 
#
#  Ginestra Bianconi, Hanlin Sun, Giacomo Rapisardi, and Alex Arenas.
#  "Message-passing approach to epidemic tracing and mitigation with apps."
#  Physical Review Research 3, no. 1 (2021): L012014.
#
# (c) 
#     Hanlin Sun (hanlin.sun@qmul.ac.uk)
#     Ginestra Bianconi (g.bianconi@qmul.ac.uk) 
#     
####################################################################################

# This program implements the (averaged) message passing algorithm and Monte Carlo simulation on real datasets
# (Friendship networks from the music streaming site Deezer in the countries of Romania, Hungary and Croatia)


using DelimitedFiles
using Plots

function network_from_data(name)

    # Read data from dataset.

    # Parameter: name (Chossiing from "Hungary", "Romania", "Croatia")
    # Return: 
    #   N: total number of node
    #   k1: degree sequence
    #   knn1:  knn1 returns a matrix of neighbor information. knn1[i, j] returns the ID of the j-th neighbor of node i.
    #   sigmaN and sigmaT: random initialization of messages.

    if name == "Hungary"
        file = "REF9_HU_edges.txt"
    elseif name == "Romania"
        file = "REF9_RO_edges.txt"
    elseif name == "Croatia"
        file = "REF9_HR_edges.txt"
    end
    data = readdlm(file, Int64) .+ 1   # Make the node ID starts from 1
    N = maximum(data)
    L = maximum(size(data)) # Number of links
    k1 = zeros(Int64, N)
    knn1 = Array{Int64}(undef, N, N)
    sigmaT = Array{Float64}(undef, N, N)
    sigmaN = Array{Float64}(undef, N, N)
    for n = 1: L
        line = data[n, :]
        i = line[1]
        j = line[2]
        k1[i] += 1
        k1[j] += 1
        knn1[i, k1[i]] = j
        knn1[j, k1[j]] = i
        sigmaT[i, j] = rand()
        sigmaT[j ,i] = rand()
        sigmaN[i, j] = rand()
        sigmaN[j, i] = rand()
    end
    return N, k1, knn1,sigmaT,sigmaN
end

function messagePassing(N::Int64, kc::Int64, k1::Vector{Int64}, 
                        knn1::Matrix{Int64}, rho::Float64, sigmaT::Matrix{Float64}, 
                        sigmaN::Matrix{Float64}, err_tol::Float64)
    
    #  Implementation of the averaged message passing algorithm (Eq.7-9 in the paper)
    #  err_tol: Error control of message passing algorithm.

    print("Message Passing Algorithm Starts......\n")

    nsum1 = 0
    T = rho * ones(Float64, N)      # Store the probability of adopting the app.
    p_list = zeros(Float64, 101)    # Store the probability p of retaining a link.
    R_list = zeros(Float64, 101)    # Store the fraction of the node in the giant component.


    for i = 1 : N
        if k1[i] > kc
            T[i] = 1
        end
    end

    for counter = 0 : 100
        p = 1 - (counter / 100)
        print("p = ", p, "\n")
        p_list[counter+1] = p
        nsumold1 = 1

        while (abs(nsum1 - nsumold1) > err_tol)
            nsumold1 = nsum1
            nsum1 = 0
            for i = 1 : N
                for n = 1 : k1[i]
                    j = knn1[i, n]
                    aus1 = 1.0
                    aus2 = 1.0
                    for np = 1 : k1[i]
                        if (np != n)
                            aus1 *= (1 - sigmaT[knn1[i, np], i] - sigmaN[knn1[i, np], i])
                            aus2 *= (1 - sigmaN[knn1[i, np], i])
                        end
                    end

                    nsum1 -= (sigmaT[i, j] + sigmaN[i, j])

                    sigmaT[i, j] = p * T[i] * (1 - aus2)
                    sigmaN[i, j] = p * (1 - T[i]) * (1 - aus1)

                    nsum1 += (sigmaT[i, j] + sigmaN[i, j])
                end
            end
        end

        # After convering, calculate the fraction of node in the giant component

        R = 0  #  Fraction of node in the giant component
        for i = 1 : N
            if k1[i] > 0
                aus1 = 1
                for n = 1 : k1[i]
                    aus1 *= (1 - sigmaT[knn1[i, n], i] - sigmaN[knn1[i, n], i])
                end
                R += (1 - aus1)
            end
        end

        R_list[counter+1] =  R / N
    end
 
    return p_list, R_list
end

function Recurrence(i::Int64, cluster_size::Int64, clusterID::Int64, clusterID_map::Vector{Int64},
                    k1::Vector{Int64}, knn1::Matrix{Int64}, y::Matrix{Int64})

    #  Recurrence is a subrutine for calulating the giant component in link percolation
    #  A unique cluter ID is asigned to each clusters.
    #  clusterID_map[i] gives the cluster ID that node i is contained.

    cluster_size += 1
    clusterID_map[i] = clusterID
    for n = 1: k1[i]
        j = knn1[i, n]
        if y[i, j] == 1 && clusterID_map[j] == 0
            aus_cluster_size = Recurrence(j, cluster_size, clusterID, clusterID_map, k1, knn1, y)
            cluster_size = aus_cluster_size
        end
    end
    return cluster_size
end

function recursive_method(N::Int64, Nrunmax::Int64, kc::Int64, k1::Vector{Int64}, knn1::Matrix{Int64}, rho::Float64)

    print("Monte Carlo simulation starts......\n")

    T = rho * ones(Float64, N)           # A node has adopted the app or not
    x = Array{Float64}(undef, N, N)         # A link is retained or damaged
    y = Array{Int64}(undef, N, N)
    GC_list = zeros(Float64, 101)
    p_list = zeros(Float64, 101)
    for i = 1 : N
        if k1[i] > kc
            T[i] = 1
        end
    end

    for nrun = 1: Nrunmax
        print("Iteration:", nrun, "\n")
        for i = 1: N
            for n = 1: k1[i]
                j = knn1[i, n]
                x[i, j] = rand()
                x[j, i] = x[i, j]  # p_ij = p_ji.
            end
        end

        for counter = 1:100

            p = 1 - (counter / 100)
            p_list[counter+1] = p

            for i = 1: N
                for n = 1: k1[i]
                    j = knn1[i, n]
                    y[i, j] = 0
                    y[j, i] = 0
                    if x[i, j] < p && T[i] * T[j] != 1
                        y[i ,j] = 1
                        y[j, i] = 1
                    end
                end
            end

            clusterID = 0
            max_ClusterSize = 0
            max_ClusterID = 0
            clusterID_map = zeros(Int64, N)

            for n = 1: N

                if clusterID_map[n] == 0
                    cluster_size = 0
                    clusterID += 1

                    cluster_size = Recurrence(n, cluster_size, clusterID, clusterID_map, k1, knn1, y)

                    if cluster_size > max_ClusterSize

                        max_ClusterSize = cluster_size
                        max_ClusterID = clusterID
                    end
                end


            end

            Nc1 = max_ClusterID
            GC = 0

            for i = 1: N
                if clusterID_map[i] == Nc1
                    GC += 1
                end

                if T[i] == 1 && clusterID_map[i] != Nc1
                    ausi = 0
                    for n = 1: k1[i]
                        j = knn1[i, n]
                        if clusterID_map[j] == Nc1 && T[j] == 1 && x[i, j] < p
                            ausi = 1
                        end
                    end
                    GC += ausi
                end
            end
            GC_list[counter+1] += GC
        end
    end
    return p_list, GC_list
end


function main()
 
    Nrunmax = 5  # Number of iterations of Monte Carlo simulatiion
    kc = 10
    rho = 0.
    err_tol = 0.01  #  Error control in message passing
    N, k1, knn1, sigmaT, sigmaN = network_from_data("Hungary")
    p_list1, R_list = messagePassing(N, kc, k1, knn1, rho, sigmaT, sigmaN, err_tol)
    p_list2, GC_list = recursive_method(N, Nrunmax, kc, k1, knn1, rho)

    plot(p_list1, R_list)
    scatter!(p_list2, GC_list / (N * Nrunmax))

end

@time main()