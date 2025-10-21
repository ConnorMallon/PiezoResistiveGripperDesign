function electrode_locations(hsnf)
    # Electrode positions
     points = []
    # #Region 1: x = 0:hs:0.085, y = 0.005
    # for x in 0:hs:0.085
    #     push!(points, (x, 0.0045))
    # end

    if hsnf == 1 
        push!(points, (0.085/2, 0.0045))
    elseif hsnf ==2 
        push!(points, (0.085/2, 0.0045))

        push!(points, (0.0, 0.005/2))
    elseif hsnf == 3
        push!(points, (0.085/2, 0.0045))

        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))
    elseif hsnf == 4

        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))

    elseif hsnf == 5 
        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))

    elseif hsnf == 6 
        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))

        push!(points, (0.13, 0.07))


    elseif hsnf ==7 
        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))

        push!(points, (0.13, 0.07))
        push!(points, (0.15, 0.065))

    elseif hsnf == 8 
        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))

        push!(points, (0.13, 0.07))
        push!(points, (0.15, 0.065))

        push!(points, (0.15, 0.07))

    elseif hsnf == 9

        push!(points, (0.085/2, 0.0045))
        
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))

        push!(points, (0.13, 0.07))
        push!(points, (0.15, 0.065))

        push!(points, (0.15, 0.07))

        push!(points, (0.121, 0.07))

    elseif hsnf > 9 

        #push!(points, (0.085/2, 0.0045))
        push!(points, (0.0, 0.005/2))
        push!(points, (0.085, 0.0045/2))
        
        x = hsnf - 8
        b = 0.085/(x+1)
        xpoints = b:b:(x*b) 
        for xp in xpoints
            push!(points, (xp, 0.0045))
        end

        push!(points, (0.12, 0.065))
        push!(points, (0.14, 0.07))
        push!(points, (0.13, 0.07))
        push!(points, (0.15, 0.065))
        push!(points, (0.15, 0.07))
        push!(points, (0.121, 0.07))


    end


         
    # # turn it into a DataFrame
    df = DataFrame(x = first.(points), y = last.(points))
    #CSV.write("points.csv", df)
    
    obs_points = [Point(x, y) for (x,y) in points]
end
