using LinearAlgebra, Statistics, PyPlot, Random, DelimitedFiles


function Pre_Inference()

    Xs = simulation3(L,M,D,Ngen,Nloa,0,0.0,0,phi,DeltaT,t1,900)
    x = Xs[:,1:end-1]
    nn,N=size(x)
    nn=Int64(nn/2)
	Dx = (Xs[:,2:end] - Xs[:,1:end-1])/(40.0*(DeltaT))
    S1 = (Dx*x')/N
    S0 = (x*x')/N

    hatA = S1*inv(S0)
    A1 = hatA[nn+1:end,1:nn]
    Lr=KronReductionMAT(L,Ngen,Nloa)
    hatM = zeros(nn,nn)
    for i in 1:nn
        hatM[i,i] = -Lr[i,i]/A1[i,i]
    end
    A2 = hatA[nn+1:end,nn+1:end]
    hatD = zeros(nn,nn)
    for i in 1:nn
        hatD[i,i] = -A2[i,i]*hatM[i,i]
    end

    return hatM, hatD, hatA

end

function KronReductionMAT(L,ngen,nloa)
    LR=L[1:ngen,1:ngen]-L[1:ngen,(ngen+1):(ngen+nloa)]*inv(L[(ngen+1):(ngen+nloa),(ngen+1):(ngen+nloa)])*L[(ngen+1):(ngen+nloa),1:ngen]
    return LR
end

function KronReductionVec(V,L,ngen,nloa)
    VR=V[1:ngen]-L[1:ngen,(ngen+1):(ngen+nloa)]*inv(L[(ngen+1):ngen+nloa,(ngen+1):ngen+nloa])*V[(ngen+1):(ngen+nloa)]
    return VR
end

function AR(L,M,D,ngen,nloa)
    Ar=zeros(2*ngen,2*ngen)
    Z=zeros(ngen,ngen)
    Id=1.0*Matrix(I,ngen,ngen) 
    LR=KronReductionMAT(L,ngen,nloa)
    B1=-inv(M)*LR
    B2=-1.0*inv(M)*D
    Ar[1:ngen,1:ngen]=Z
    Ar[1:ngen,(ngen+1):2*ngen]=Id
    Ar[(ngen+1):2*ngen,1:ngen]=B1
    Ar[(ngen+1):2*ngen,(ngen+1):2*ngen]=B2
    return Ar
end


phi=0 #phase
DeltaT=0.0005#1#0.02 #time step
t1=0 #initial time
t2=400#200 #final time

Ngen=7#3 #number of generators
Nloa=50#10 #number of loads
Random.seed!(109800)

L = 1.0*readdlm("adjIEEE57",',')


M = 2.5*Matrix(I(7))
D = 1.0*Matrix(I(7))

Random.seed!(1090000)

M[3,3]=4.0
M[5,5]=1.5
D[3,3]=1.6
D[5,5]=1.2



        
function simulation2(L,M,D,Ngen,Nloa,l,gamma,f,phi,DeltaT,t1,t2)
   Random.seed!(l*10*1000980)

    A=AR(L,M,D,Ngen,Nloa)
    Y=[zeros(2*Ngen) zeros(2*Ngen)]
    Mat=I+DeltaT*A

    Ydisc1 = zeros(2*Ngen,1)#Y[:,end]
    data = zeros(2*Ngen,Int64(t2/DeltaT)+1)
    data[:,1]=Ydisc1
    noise= zeros(Ngen,Int64(t2/DeltaT)+1)
    jsp=1
    for j in 1-19999:(Int64(t2/DeltaT))

        Ydisc1=Mat*Ydisc1 #[Id+\delta t A]\bar{Y_t}

        Ft=DeltaT*gamma*cos(2*pi*(f*j*DeltaT+phi))*I[1:(Ngen+Nloa), l]

        V2=vcat(randn(Ngen)*1.0*sqrt(DeltaT),1.0*randn(Nloa)*1.0*sqrt(DeltaT))

       if(j<0)
            Ydisc1=Ydisc1+vcat(zeros(Ngen),inv(M)*KronReductionVec(1.0*V2,L,Ngen,Nloa))
    else
            Ydisc1=Ydisc1+vcat(zeros(Ngen),inv(M)*KronReductionVec(Ft+1.0*V2,L,Ngen,Nloa))
    end


        if(mod(j,1)==0 && j>0)
        data[:,jsp+0] = Ydisc1
        jsp+=1
        end

    end
    return data
end

      
function simulation3(L,M,D,Ngen,Nloa,l,gamma,f,phi,DeltaT,t1,t2)
    Random.seed!(10760000)
 
     A=AR(L,M,D,Ngen,Nloa)
     Y=[zeros(2*Ngen) zeros(2*Ngen)]
     Mat=I+DeltaT*A

     Ydisc1 = zeros(2*Ngen,1)#Y[:,end]
     data = []#zeros(2*Ngen,Int64(t2/DeltaT)+1)
     data = push!(data, Ydisc1)
     noise= zeros(Ngen,Int64(t2/DeltaT)+1)
     jsp=1
     for j in 1-999:(Int64(t2/DeltaT))

         Ydisc1=Mat*Ydisc1 #[Id+\delta t A]\bar{Y_t}
         Ft=DeltaT*gamma*cos(2*pi*(f*(0+j)*DeltaT+phi))*I[1:Ngen+Nloa, l]

         V2=vcat(randn(Ngen)*1.0*sqrt(DeltaT),1.0.*randn(Nloa)*1.0*sqrt(DeltaT))

        Ydisc1=Ydisc1+vcat(zeros(Ngen),1.0*inv(M)*KronReductionVec(V2,L,Ngen,Nloa))

         if(mod(j,40)==0 && j>0)
            data = push!(data, Ydisc1)
            jsp+=1
         end

     end
     datar = zeros(2*size(M,1),size(data,1))
     for i in 1:size(data,1)
        datar[:,i] = data[i]
     end
     return datar
 end