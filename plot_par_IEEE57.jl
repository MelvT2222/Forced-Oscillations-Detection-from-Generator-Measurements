gamma=10.0
f=1.0
s=55
kf = collect(Int64(f*t2)-1:1:Int64(f*t2)+1) # pick frequencies around the forcing

#figure(1919)
for i in 1:57#[12]#1:5#7#7:42#[1,13]
    fr=Float64[]
    for k in kf
        a=readdlm("data/$(s)_$(f)_$(gamma)_IEEE57_GMO_$(i).$(k)_obj.csv",',')

        push!(fr,a[1])
    end

    if(i==s) 
        plot(kf./(t2) ,fr, "go")
    else
        plot(kf./(t2) ,fr,"grey", alpha=0.3,label="$i")
    end
    xlabel(L"f[Hz]",fontsize=16)
    ylabel(L"-L",fontsize=16)
end
tight_layout()
savefig("Figure_FO_GMO_$(s).pdf")
