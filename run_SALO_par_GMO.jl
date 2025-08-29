include("Methods_GMO.jl")
include("SALO_par_GMO.jl")
lo=Int64[]
fr=Int64[]

sp=0

for gamma in [10.0]# amplitude of the forcing
    for f in [1.0]
        kf = collect(Int64(f*t2)-1:1:Int64(f*t2)+1)

        ii = 55 # location of the forcing
        println("DATA generation")
        RES4 = simulation2(L,M,D,Ngen,Nloa,ii,gamma,f,0.0,DeltaT,t1,t2+0*20)[:,1:40:end-1]
        writedlm("IEEE57_data_$(ii)_$(f)_$(gamma)",RES4,',')
        println("DATA generated")


        MH1, DH1, AA = Pre_Inference()

        run_SALO_par_GMO("$(ii)_$(f)_$(gamma)_IEEE57_", L, MH1, DH1, collect(1:Ngen), collect((Ngen+1):(Ngen+Nloa)), RES4[:,1:1:end], 40.0*DeltaT, collect(1:57),kf)

    end
end
