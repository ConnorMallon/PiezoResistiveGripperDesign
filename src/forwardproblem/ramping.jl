function ramp(opd,nls,f_in,k∞)
    # d0 = solve(nls,opd(1e7,k∞))
    # solve!(d0,nls,opd(2e7,k∞))
    # solve!(d0,nls,opd(5e7,k∞))
    # solve!(d0,nls,opd(6e7,k∞))
    # solve!(d0,nls,opd(7.3e7,k∞))
    # solve!(d0,nls,opd(7.9e7,k∞))
    # solve!(d0,nls,opd(8.4e7,k∞))
    # solve!(d0,nls,opd(8.7e7,k∞))
    # solve!(d0,nls,opd(9e7,k∞))
    # solve!(d0,nls,opd(1e8,k∞))
    # solve!(d0,nls,opd(1.5e8,k∞))

    # f = 1.5e8
    # while f < f_in
    #     solve!(d0,nls,opd(f,k∞))
    #     f+=0.25e8
    #     if f<=f_in
    #     elseif f>f_in 
    #         println("f=f_in")
    #         f=f_in
    #     end
    #     @show f 
    # end
    # # solve!(d0,nls,opd(f,2e11))
    # # solve!(d0,nls,opd(f,5e11))
    # @show f 
    # d0

    # println("ramp / 1000")
    d0 = solve(nls,opd(f_in,k∞/100))
    # println("ramp / 750")
    # d0 = solve(nls,opd(f_in,k∞/750))
    # println("ramp / 500")
    # solve!(d0,nls,opd(f_in,k∞/500))
    # println("ramp / 100")
    #solve!(d0,nls,opd(f_in,k∞/100))
    println("ramp / 10")
    solve!(d0,nls,opd(f_in,k∞/10))
    println("ramp / 1")

    @show k∞
    solve!(d0,nls,opd(f_in,k∞))
    println("ramp done")
    d0



end
 