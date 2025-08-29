using Distributed

# Define the number of parallel threads.
n_thr = 3
@info "How many parallel threads are required? (Default is 3)"
try
	global n_thr = parse(Int64,readline())
catch eee
	@info "Default number of threads used: 3."
end

@info "Loading optimizer..."

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

# If the directory "data" does not exists, then create it.
if !isdir("data/")
	mkdir("data/")
end

@everywhere using FFTW, LinearAlgebra, DelimitedFiles

#################################################
@everywhere function run_SALO_par_GMO(id::String, L::Matrix{Float64}, Mh::Matrix{Float64}, Dh::Matrix{Float64}, gen::Vector{Int64}, loa::Vector{Int64}, Xs::Matrix{Float64}, τ::Float64, ls::Vector{Int64}, ks::Vector{Int64}, b::Float64=0., μ::Float64=0.5e-4, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(round(nn/2,digits=0))
	N = NN-1
	Lr = L[gen,gen] - L[gen,loa]*pinv(L[loa,loa])*L[loa,gen]
	Sigmagl = I(n) + L[gen,loa]*pinv(L[loa,loa])*pinv(L[loa,loa])*L[loa,gen]
	Sgl1 = inv(Sigmagl)
	# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[:,2:end] - Xs[:,1:end-1])/τ
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end
	args = Array{Tuple{String,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Float64,Float64,Float64},1}()

	for l in ls
		for k in ks
			push!(args,(id,Lr,L,Mh,Dh,Sgl1,gen,loa,x,Dx,xt,Dxt,l,k,b,μ,bp))
		end
	end

	pmap(Lmax_SALO_par_GMO,args)
end

# likelihood maximum using generator measurements only (GMO)
@everywhere function Lmax_SALO_par_GMO(tups::Tuple{String,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Float64,Float64,Float64})
	id,Lr,L,Mh,Dh,Sigmagl,gen,loa,x,Dx,xt,Dxt,l,kl,b,μ,bp = tups

	@info "===================================================================================="
	@info "Computing "*id*" with SALO_GMO optimization: ℓ = $l, k = $kl."
	@info "===================================================================================="
				m=(Mh)#.^(-1)
				d=(Dh)
				nn,N = size(x)
				n = Int(nn/2)
				Gammal = Vector{Float64}
				if(l <= length(gen))
					Gammal = (I(length(gen))[:,l])
				else
					Gammal = -L[gen,loa]*pinv(L[loa,loa])*(I(length(loa))[:,maximum([1,l-length(gen)])])
				end

				Sgl1 = Sigmagl

				# Definition of the needed parameters.
				xtk = xt[:,kl+1]		# THE FIRST COLUMN IS f=0, WHICH WE DON'T WANT TO TREAT.
				Dxtk = Dxt[:,kl+1]		# THE FIRST COLUMN IS f=0, WHICH WE DON'T WANT TO TREAT.

				Fk = real.(xtk*xtk')
				DFk = real.(Dxtk*Dxtk')
				DFk2 = real.(Dxtk*xtk')
				C = Sigmagl*Gammal*Gammal'*Sigmagl

				coef1 = Gammal'*Sigmagl*Gammal
				coef2 = 2.0/(sqrt(N))*sqrt((tr(Fk[n+1:end,n+1:end]*d*C*d )  + 2.0*tr(Fk[1:n,n+1:end]*d*C*Lr) + 2.0*tr(DFk[n+1:end,1:n]*d*C*m) + tr(Fk[1:n,1:n]*Lr*C*Lr) +2.0*tr(DFk2[n+1:end,1:n]*Lr*C*m) + tr(DFk[n+1:end,n+1:end]*m*C*m)))


				
				if((coef1)>0.0 && abs(coef1)>1.0e-10)
					gamma_hat = coef2/coef1
					writedlm("data/"*id*"GMO_$(l).$(kl)_obj.csv",-(gamma_hat^2*0.5*coef1 - gamma_hat*coef2),',')

				else
					writedlm("data/"*id*"GMO_$(l).$(kl)_obj.csv",0.0,',')
				end
end
############################################################################

@info "Loaded."